#!/usr/bin/env Rscript
# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: Aug 20, 2021
# Updated: Jan 29, 2022

# Estimtate the motif break score
options(stringsAsFactors = FALSE)

# Libraries
library(ggrepel)
library(magrittr)
library(data.table)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(tidyverse)
library(MotifDb) # Database

merge_active <- TRUE
if (merge_active) {
  disease_group <- c("active", "post")
  condition_map <- c("active" = "Hospitalized", "post" = "Convalescent")
} else {
  disease_group <- c("Severe", "Mild", "Post")
  condition_map <- c("Severe" = "Seve.", "Mild" = "Mild", "Post" = "Conv.")
}

celltype_order <- c("cMono", "CD4T", "CD8T", "NK", "B", "ncMono")

# Workspace
pjdir <- "~/Documents/projects/wp_covid19_mhh50"
wkdir <- str_glue("{pjdir}/outputs/scATAC-seq/motifbreakr")

asoc_snp_file <- str_glue("{pjdir}/outputs/scATAC-seq/summary/pseudo-bulk_per-celltype_asoc-snps_annotated.csv")
asoc_snp_tab <- fread(asoc_snp_file, verbose = F, showProgress = F)


# Motifs annotated by scATAC-seq results (chromVAR, CIS-BP 2.0)
load_motifs <- TRUE
if (load_motifs) {
  tar_motif_lst <- NULL
  for (cond in disease_group) {
    ocr_snp_file <- str_glue("{pjdir}/inputs/scATAC-seq/combined_{cond}.txt")
    if (is.null(tar_motif_lst)) {
      tar_motif_lst <- fread(ocr_snp_file, verbose = F, showProgress = F) %>%
        dplyr::mutate(condition = cond)
    } else {
      tar_motif_lst <- fread(ocr_snp_file, verbose = F, showProgress = F) %>%
        dplyr::mutate(condition = cond)
      rbind(tar_motif_lst)
    }
  }
  tar_motif_lst %<>%
    dplyr::filter(SNP_in_motif != "no_motif") %$%
    SNP_in_motif %>%
    str_c(collapse = ";") %>%
    str_split(pattern = ";", simplify = T) %>%
    str_remove(pattern = "_.*$") %>%
    unique() %>%
    discard(function(x) x == "")
}

# Enriched motifs per TF
atac_motif_lst <- c(
  "SPIC", "SPIB", "POU2F1", "POU2F3", "POU5F1B", "POU2F2", "POU3F4", "ARID3A",
  "NR2F1", "BCL11A", "BCL11B", "STAT2", "IRF4", "SPI1", "RUNX1", "RUNX2", "SP2",
  "KLF5", "KLF4", "ZNF148", "SP1", "SP3", "WT1", "FOS", "FOSL2", "JUN", "JUNB",
  "FOSL1", "JUND", "CEBPD", "CEBPA"
)

# Estimate motif break score
# NOTE:
#   1. jaspar2018 has the largest overlap with scATAC-seq motifs
mb_file <- str_glue("{wkdir}/motifbreakR.csv")
if (file.exists(mb_file)) {
  mb_est <- data.table::fread(mb_file, verbose = F, showProgress = F) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
} else {
  mb_est <- NULL
  asoc_snp_rsid <- asoc_snp_tab %$% variantID %>% unique()

  # Fetch ASoC SNP postion information
  asoc_snp_gr <- snps.from.rsid(
    rsid = asoc_snp_rsid,
    dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,
    search.genome = BSgenome.Hsapiens.UCSC.hg38
  )
  tryCatch(
    {
      # 6 hours
      tm_start <- proc.time()
      pwm_lst <- subset(MotifDb,
        subset = organism == "Hsapiens" &
          dataSource == "jaspar2018" &
          geneSymbol %in% c(tar_motif_lst, atac_motif_lst)
      )

      # Estimate disruption motif break effects.
      mb_est <- motifbreakR(
        snpList = asoc_snp_gr, filterp = T, pwmList = pwm_lst,
        threshold = 1e-4, method = "log",
        bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
        BPPARAM = BiocParallel::SnowParam(workers = 6)
      )
      print(proc.time() - tm_start)
    },
    error = function(e) {
      print(e)
      BiocParallel::bpstop()
    }
  )

  # Save motif break score into disk
  if (!is.null(mb_est)) {
    mb_est %>%
      as.data.table() %>%
      fwrite(mb_file, sep = ",", row.names = F, quote = F)
  }
}

# Calculate disruption score by quantifying the difference between scoreAlt() and scoreRef()
# NOTE:
#   1. Alleles were aligned base on risk alleles of COVID19_HGI_A2_ALL-leave_23andme_20210607.txt.gz
#   2. Only SNPs with predicted "strong" effect were used for the correlation.
cor_mthd <- "spearman"
cor_est_file <- str_glue("{wkdir}/motifbreakR_mbs_vs_ai_{cor_mthd}.csv")
if (file.exists(cor_est_file)) {
  cor_est_tab <- fread(cor_est_file, showProgress = F, verbose = F)
} else {
  # COVID19_HGI_B2_ALL_leave_23andme_20210607_risk_allele.csv
  gwas_ra_file <- str_glue("{pjdir}/outputs/covid19_gwas/COVID19_HGI_B2_ALL-leave_23andme_20210607_risk_allele.csv")
  gwas_ra_tab <- fread(gwas_ra_file, verbose = F, showProgress = F, tmpdir = dirname(gwas_ra_file))
  min_comp_pairs <- 5
  cor_est_tab <- NULL
  for (cur_ctype in c("cMono", "CD4T", "CD8T", "B", "NK", "ncMono")) {
    cor_est_tab <- mb_est %>%
      as.data.table() %>%
      dplyr::filter(effect %in% c("strong")) %>%
      dplyr::group_by(geneSymbol) %>%
      dplyr::summarise(cor = {
        mbscore <- cur_data() %>%
          as.data.table() %>%
          dplyr::select(SNP_id, REF, ALT, alleleDiff)

        cor_res <- NULL
        for (cur_cond in disease_group) {
          rc2ad_tab <- asoc_snp_tab %>%
            dplyr::filter(celltype %in% cur_ctype & condition %in% cur_cond) %>%
            dplyr::select(variantID, refAllele, altAllele, refCount, altCount) %>%
            dplyr::inner_join(mbscore, by = c("variantID" = "SNP_id")) %>%
            dplyr::inner_join(gwas_ra_tab, by = c("variantID" = "RSID")) %>%
            dplyr::mutate(
              rc_log2or = case_when(
                altAllele == RISK_ALLELE & refAllele == NON_RISK_ALLELE ~ log2(altCount / refCount),
                refAllele == RISK_ALLELE & altAllele == NON_RISK_ALLELE ~ log2(refCount / altCount),
                T ~ NA_real_
              ),
              alleleDiff = case_when(
                ALT == RISK_ALLELE & REF == NON_RISK_ALLELE ~ alleleDiff,
                REF == RISK_ALLELE & ALT == NON_RISK_ALLELE ~ -alleleDiff,
                T ~ NA_real_
              )
            ) %>%
            dplyr::select(variantID, rc_log2or, alleleDiff) %>%
            dplyr::distinct()

          tryCatch(
            {
              ad <- rc2ad_tab %$% alleleDiff
              rc <- rc2ad_tab %$% rc_log2or

              ht <- list(p.value = 1, estimate = 0)
              if (length(ad) == length(rc) && length(ad) >= min_comp_pairs) {
                ht <- cor.test(ad, rc, method = cor_mthd)
              }
              cor_res <- data.frame(p = ht$p.value, r = ht$estimate, cc = cur_cond, ct = cur_ctype) %>%
                rbind(cor_res)
            },
            error = function(e) {
              genesym <- cur_group() %$% geneSymbol
              cat(e$message, ", ", genesym, ", ", cur_cond, "\n", sep = "")
            }
          )

          print(paste(cur_cond, cur_ctype, cur_group(), sep = ", "))
        }
        cor_res
      }) %>%
      as.data.table() %>%
      dplyr::rename("p_value" = "cor.p", "coef" = "cor.r", "condition" = "cor.cc", "celltype" = "cor.ct") %>%
      dplyr::mutate(
        fdr = p.adjust(p_value, method = "fdr"),
        signif = factor(if_else(fdr < 0.05, "FDR<0.05", ""), levels = c("FDR<0.05", "")),
        condition = factor(condition, levels = disease_group)
      ) %>%
      rbind(cor_est_tab)
  }
  # Save the estimation results to the disk
  cor_est_tab %>% fwrite(cor_est_file, sep = ",", row.names = F, quote = F)
}


cor_est_for_plot <- cor_est_tab %>%
  dplyr::filter(geneSymbol %in% (dplyr::filter(., p_value < 0.01) %$% geneSymbol)) %>%
  dplyr::mutate(
    geneSymbol = factor(
      geneSymbol,
      levels = dplyr::filter(., condition == "active" & celltype == "cMono") %>%
        dplyr::arrange(coef) %$% geneSymbol
    ),
    celltype = factor(celltype, levels = celltype_order),
    condition = factor(condition, levels = rev(disease_group))
  )

g_hmap <- cor_est_for_plot %>%
  ggplot(aes(x = geneSymbol, y = condition)) +
  geom_tile(aes(fill = coef), color = "white") +
  geom_point(
    data = (dplyr::filter(cor_est_for_plot, signif == "FDR<0.05")),
    mapping = aes(shape = signif)
  ) +
  scale_shape_manual(name = NULL, values = c(4)) +
  scale_fill_gradient2(
    name = "Rho", # quote("Polerized" ~ -Log[10] ~ "(p-value)"),
    low = "darkblue", mid = "gray97", high = "darkred"
  ) +
  scale_y_discrete(labels = condition_map) +
  facet_grid(celltype ~ .) +
  labs(y = NULL, x = NULL) +
  theme_classic() +
  theme(
    axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14)
  )

g_hmap_save_to <- file.path(wkdir, "motifbreakR_mbs_vs_ai.pdf")
ggsave(g_hmap_save_to, plot = g_hmap, width = 14, height = 6)


# Estimate the motif break score at position of rs7255545
# rs7255545 <- snps.from.rsid(
#   rsid = "rs7255545",
#   dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,
#   search.genome = BSgenome.Hsapiens.UCSC.hg38
# )

sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 11 (bullseye)
#
# Matrix products: default
# BLAS:   /home/umcg-zzhang/tools/R/lib/R/lib/libRblas.so
# LAPACK: /home/umcg-zzhang/tools/R/lib/R/lib/libRlapack.so
#
# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
#  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats4    grid      stats     graphics  grDevices utils     datasets
# [8] methods   base
#
# other attached packages:
#  [1] forcats_0.5.1
#  [2] stringr_1.4.0
#  [3] dplyr_1.0.7
#  [4] purrr_0.3.4
#  [5] readr_2.1.1
#  [6] tidyr_1.1.4
#  [7] tibble_3.1.6
#  [8] tidyverse_1.3.1
#  [9] SNPlocs.Hsapiens.dbSNP144.GRCh38_0.99.20
# [10] BSgenome.Hsapiens.UCSC.hg38_1.4.3
# [11] BSgenome_1.60.0
# [12] rtracklayer_1.52.1
# [13] motifbreakR_2.6.1
# [14] MotifDb_1.34.0
# [15] Biostrings_2.62.0
# [16] XVector_0.34.0
# [17] GenomicRanges_1.46.1
# [18] GenomeInfoDb_1.30.0
# [19] IRanges_2.28.0
# [20] S4Vectors_0.32.3
# [21] BiocGenerics_0.40.0
# [22] data.table_1.14.2
# [23] magrittr_2.0.1
# [24] ggrepel_0.9.1
# [25] ggplot2_3.3.5
# [26] devtools_2.4.3
# [27] usethis_2.1.5
#
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1                backports_1.4.1
#   [3] Hmisc_4.6-0                 BiocFileCache_2.0.0
#   [5] lazyeval_0.2.2              splines_4.1.2
#   [7] BiocParallel_1.28.3         digest_0.6.29
#   [9] ensembldb_2.16.4            htmltools_0.5.2
#  [11] fansi_1.0.2                 checkmate_2.0.0
#  [13] memoise_2.0.1               cluster_2.1.2
#  [15] tzdb_0.2.0                  remotes_2.4.2
#  [17] modelr_0.1.8                matrixStats_0.61.0
#  [19] prettyunits_1.1.1           jpeg_0.1-9
#  [21] colorspace_2.0-2            rvest_1.0.2
#  [23] blob_1.2.2                  rappdirs_0.3.3
#  [25] haven_2.4.3                 xfun_0.29
#  [27] jsonlite_1.7.3              callr_3.7.0
#  [29] crayon_1.4.2                RCurl_1.98-1.5
#  [31] TFMPvalue_0.0.8             survival_3.2-13
#  [33] VariantAnnotation_1.38.0    glue_1.6.1
#  [35] gtable_0.3.0                zlibbioc_1.40.0
#  [37] DelayedArray_0.20.0         pkgbuild_1.3.1
#  [39] scales_1.1.1                DBI_1.1.2
#  [41] Rcpp_1.0.8                  progress_1.2.2
#  [43] htmlTable_2.4.0             foreign_0.8-82
#  [45] bit_4.0.4                   Formula_1.2-4
#  [47] htmlwidgets_1.5.4           httr_1.4.2
#  [49] RColorBrewer_1.1-2          ellipsis_0.3.2
#  [51] farver_2.1.0                pkgconfig_2.0.3
#  [53] XML_3.99-0.8                Gviz_1.36.2
#  [55] nnet_7.3-17                 dbplyr_2.1.1
#  [57] utf8_1.2.2                  labeling_0.4.2
#  [59] tidyselect_1.1.1            rlang_0.4.12
#  [61] AnnotationDbi_1.56.2        cellranger_1.1.0
#  [63] munsell_0.5.0               tools_4.1.2
#  [65] cachem_1.0.6                cli_3.1.1
#  [67] generics_0.1.1              RSQLite_2.2.9
#  [69] ade4_1.7-18                 broom_0.7.11
#  [71] fastmap_1.1.0               yaml_2.2.1
#  [73] processx_3.5.2              knitr_1.37
#  [75] bit64_4.0.5                 fs_1.5.2
#  [77] KEGGREST_1.34.0             AnnotationFilter_1.16.0
#  [79] splitstackshape_1.4.8       xml2_1.3.3
#  [81] biomaRt_2.48.3              brio_1.1.3
#  [83] compiler_4.1.2              rstudioapi_0.13
#  [85] filelock_1.0.2              curl_4.3.2
#  [87] png_0.1-7                   testthat_3.1.2
#  [89] reprex_2.0.1                stringi_1.7.6
#  [91] ps_1.6.0                    GenomicFeatures_1.44.2
#  [93] desc_1.4.0                  lattice_0.20-45
#  [95] ProtGenerics_1.24.0         Matrix_1.4-0
#  [97] vctrs_0.3.8                 pillar_1.6.4
#  [99] lifecycle_1.0.1             bitops_1.0-7
# [101] R6_2.5.1                    BiocIO_1.2.0
# [103] latticeExtra_0.6-29         gridExtra_2.3
# [105] motifStack_1.36.1           sessioninfo_1.2.2
# [107] dichromat_2.0-0             MASS_7.3-55
# [109] assertthat_0.2.1            pkgload_1.2.4
# [111] SummarizedExperiment_1.24.0 rprojroot_2.0.2
# [113] rjson_0.2.21                withr_2.4.3
# [115] GenomicAlignments_1.28.0    Rsamtools_2.8.0
# [117] GenomeInfoDbData_1.2.7      parallel_4.1.2
# [119] hms_1.1.1                   rpart_4.1-15
# [121] MatrixGenerics_1.6.0        biovizBase_1.40.0
# [123] lubridate_1.8.0             Biobase_2.54.0
# [125] base64enc_0.1-3             restfulr_0.0.13
