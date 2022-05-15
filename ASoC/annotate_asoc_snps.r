#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 27, 2021
# Updated: Jan 28, 2022

library(biomaRt)
library(data.table)
library(magrittr)
library(progress)
library(tidyverse)

# Working dir, input and output
pjdir <- "~/Documents/projects/wp_covid19_mhh50"
wkdir <- file.path(pjdir, "outputs", "scATAC-seq", "summary")

#
## Assign chromatin state based on cell type
#
# The epigenetic features to be included
promoters <- c("1_TssA", "2_PromU", "3_PromD1", "4_PromD2")
enhancers <- c(
  "9_TxReg", "10_TxEnh5", "11_TxEnh3", "12_TxEnhW", "13_EnhA1", "14_EnhA2",
  "15_EnhAF", "16_EnhW1", "17_EnhW2", "18_EnhAc", "19_DNase"
)
chosen_states <- c(promoters, enhancers)

condition_map <- c("active" = "active", "post" = "convalescent")

#
## Assign a SNP to genes if the SNP is located in the promoter region (-2kb, +1kb) of the gene.
#
## Compile a list of promoter regions
## 2K bp upstream to the TSS and 1K bp downstream to TSS.
## Here Ensembl general features (release 104) were retrieved. Reference genome is GRCh38
if (F) listEnsembl()
#         biomart                version
# 1         genes      Ensembl Genes 104
# 2 mouse_strains      Mouse strains 104
# 3          snps  Ensembl Variation 104
# 4    regulation Ensembl Regulation 104
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Here we use GENCODE basic and Ensembl Canonical transcripts
# FIXME: Please note, only protein coding genes on autosomals were included.
tss_tab <- getBM(
  attributes = c(
    "chromosome_name", "transcript_start", "transcript_end", "strand",
    "transcription_start_site", "ensembl_gene_id", "external_gene_name"
  ),
  filters = c(
    "chromosome_name", "biotype", "transcript_gencode_basic",
    "transcript_is_canonical"
  ),
  values = list(as.character(1:22), c("protein_coding"), T, T),
  mart = ensembl
)

tss_tab <- tss_tab %>%
  arrange(chromosome_name, transcription_start_site) %>%
  mutate(
    promoter_start = mapply(
      function(pos, strand) ifelse(strand == "1", pos - 2000, pos - 1000),
      pos = transcription_start_site, strand = strand
    ),
    promoter_end = mapply(
      function(pos, strand) ifelse(strand == "1", pos + 1000, pos + 2000),
      pos = transcription_start_site, strand = strand
    ),
    strand = case_when(strand == "-1" ~ "-", strand == "1" ~ "+")
  ) %>%
  dplyr::filter(external_gene_name != "")


###############################################################################
# Assign ASoC SNPs to DE gene,
###############################################################################
# NOTE: If the TSS of the DE gene is located in the 25kbp window of the SNP.
degene_file <- file.path(pjdir, "inputs", "scRNA-seq", "DEG", "DEtable.2cond.txt")
degene_tab <- fread(degene_file) %>%
  dplyr::filter(!cell %in% c("cMono", "CD163.cMono")) %>%
  dplyr::rename("external_gene_name" = "gene") %>%
  dplyr::left_join(tss_tab, by = "external_gene_name") %>%
  dplyr::mutate(cell = case_when(
    cell == "CD8.T" ~ "CD8T",
    cell == "CD4.T" ~ "CD4T",
    cell == "merged.cMonos" ~ "cMono",
    T ~ cell
  )) %>%
  dplyr::filter(!is.na(chromosome_name)) %>%
  dplyr::arrange(chromosome_name, transcript_start, transcript_end) %>%
  dplyr::relocate(ensembl_gene_id)


###############################################################################
# Assign an ASoC SNP to genes
###############################################################################
# NOTE: If the TSS of a eQTL gene is located within the 250kbp of the SNP.
# NOTE: Compile a list of blood eQTL SNP genes according whold blood samples
# by GTEx v8 We kept the variant-gene pair if it's
# 1. The variant is a SNP
# 2. The MAF of the SNP is greater equal 0.1
# 3. The gene is in the tss_tab (check tss_tab for more information)
egene_file <- file.path(
  pjdir, "inputs", "GTEx", "Whole_Blood.v8.EUR.signif_pairs.txt.gz"
)
egene_tab <- fread(egene_file) %>%
  tidyr::separate(variant_id,
    into = c("contig", "position", "refAllele", "altAllele", "build")
  ) %>%
  dplyr::filter(nchar(refAllele) == 1 & nchar(altAllele) == 1 & maf >= 0.05) %>%
  dplyr::mutate(
    phenotype_id = mapply(function(x) str_split(x, "\\.", simplify = T)[1, 1], phenotype_id),
    position = as.integer(position),
    contig = as.integer(str_remove(contig, "chr"))
  ) %>%
  dplyr::inner_join(tss_tab, by = c("phenotype_id" = "ensembl_gene_id")) %>%
  dplyr::arrange(contig, position) %>%
  dplyr::rename("ensembl_gene_id" = phenotype_id) %>%
  dplyr::select(c(
    ensembl_gene_id, contig, position, refAllele, altAllele, external_gene_name
  )) %>%
  dplyr::filter(external_gene_name != "")

###############################################################################
# Assign chromatin_state to ASoC SNPs
###############################################################################
# NOTE: Using bash/intersect_asocsnps_and_epimarks.sh to make heter-snps_epiannot.csv
epimark_file <- file.path(wkdir, "heter-snps_epiannot.csv")
if (!file.exists(epimark_file)) {
  stop("Use bash/intersect_asocsnps_and_epimarks.sh to make heter-snps_epiannot.csv")
} else {
  epimark_tab <- fread(epimark_file, verbose = F, showProgress = F) %>%
    group_by(contig, position, variantID) %>%
    summarise(
      chrom_state = paste(chrom_start, chrom_stop, chrom_state, celltype,
        sep = "|", collapse = ";"
      )
    ) %>%
    distinct()
}

# Chromatin state evidence.
celltype_trans <- c(
  "Primary_mononuclear_cells_from_peripheral_blood" = "Plasmablast pDC",
  # "Primary_T_cells_from_peripheral_blood"="CD4T CD8T",
  # "Primary_T_cells_effector_memory_enriched_from_peripheral_blood"="CD4T CD8T",
  "Primary_T_regulatory_cells_from_peripheral_blood" = "CD4T",
  "Primary_T_helper_cells_from_peripheral_blood" = "CD4T",
  "Primary_T_helper_naive_cells_from_peripheral_blood" = "CD4T",
  "Primary_T_helper_memory_cells_from_peripheral_blood_1" = "CD4T",
  "Primary_T_helper_memory_cells_from_peripheral_blood_2" = "CD4T",
  "Primary_T_helper_naive_cells_from_peripheral_blood" = "CD4T",
  "Primary_T_CD8+_memory_cells_from_peripheral_blood" = "CD8T",
  "Primary_T_CD8+_naive_cells_from_peripheral_blood" = "CD8T",
  "Primary_monocytes_from_peripheral_blood" = "cMono ncMono",
  "Monocytes-CD14+_RO01746_Primary_Cells" = "cMono ncMono",
  "Primary_B_cells_from_peripheral_blood" = "B",
  "Primary_Natural_Killer_cells_from_peripheral_blood" = "NK"
)

# Assign annotations, including gene by promoter, gene by eQTL, gene by
# differential expression, and chromatin state including promoters & enhancers
wd_size <- 25000
heter_snp_file <- file.path(wkdir, "pseudo-bulk_per-celltype_heter-snps.csv")
heter_snp_tab <- fread(heter_snp_file) %>%
  dplyr::arrange(contig, position) %>%
  dplyr::select(contig, position, variantID, celltype, condition) %>%
  dplyr::left_join(egene_tab, by = c("contig", "position")) %>%
  dplyr::filter(variantID != "NA") %>%
  dplyr::group_by(contig, position, variantID, celltype, condition) %>%
  dplyr::summarise(by_eqtlgene = paste0(external_gene_name, collapse = "|"), n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(epimark_tab, by = c("contig", "position", "variantID")) %>%
  dplyr::mutate(
    by_eqtlgene = case_when(by_eqtlgene == "NA" ~ "", T ~ by_eqtlgene),
    up_bound = if_else(position - wd_size > 0, position - wd_size, 1),
    dw_bound = position + wd_size
  )

###############################################################################
# Assign
###############################################################################
# NOTE:
# 1. In the apply function, the variable should not have the same name in the
# column names of the data frame, otherwise, it's overrided by the column name.
heter_snp_annot_file <- str_glue("{wkdir}/pseudo-bulk_per-celltype_heter-snps_annotated.csv")
if (file.exists(heter_snp_annot_file)) {
  heter_snp_annot_tab <- fread(heter_snp_annot_file, verbose = F, showProgress = F)
} else {
  # Make a progess bar
  pbar <- progress_bar$new(total = nrow(heter_snp_tab), format = "[:bar] :current/:total :percent :eta")

  heter_snp_annot_tab <- heter_snp_tab %>%
    # head(2500) %>%
    dplyr::mutate(ANN = base::apply(., 1, function(x, .tss_tab, .degene_tab) {
      if (!pbar$finished) pbar$tick() # make sure the progress bar does not break the "loop"

      .contig <- as.integer(x["contig"])
      .position <- as.integer(x["position"])
      .up_bound <- as.integer(x["up_bound"])
      .dw_bound <- as.integer(x["dw_bound"])
      .celltype <- x["celltype"]
      .chrom_state <- x["chrom_state"]
      .condition <- condition_map[x["condition"]]

      # Assign gene by promoter to the SNP
      by_promoter <- ""
      .by_promoter <- .tss_tab %>% # Gene assigned by promoter
        dplyr::filter(chromosome_name == .contig &
          .position >= promoter_start &
          .position <= promoter_end &
          external_gene_name != "")
      if (nrow(.by_promoter) > 0) {
        by_promoter <- .by_promoter %$%
          external_gene_name %>%
          unique() %>%
          paste(collapse = "|")
      }

      # Assign DEG to the SNP
      by_degene <- ""
      .by_degene <- .degene_tab %>% # Gene assigned by DE
        dplyr::filter(chromosome_name == .contig &
          .up_bound <= transcription_start_site &
          transcription_start_site <= .dw_bound &
          cell %in% .celltype &
          grepl(.condition, comp, fixed = T) &
          external_gene_name != "")
      if (nrow(.by_degene) > 0) {
        by_degene <- .by_degene %>%
          dplyr::select(external_gene_name, comp, direction) %>%
          dplyr::mutate(degene = paste(external_gene_name, ":", direction, "@", comp, sep = "")) %$%
          degene %>%
          unique() %>%
          paste(collapse = "|")
      }

      # Assign annotation (promoter/enhancers) to the SNP
      annot <- ""
      annot_type <- ""
      .regs <- .chrom_state %>%
        str_split(pattern = ";", simplify = T) %>%
        str_split(pattern = "\\|", simplify = T) %>%
        as.data.table()
      if (ncol(.regs) == 4) {
        .regs <- .regs %>%
          rlang::set_names(c("cs_start", "cs_stop", "cs_annot", "cs_celltype")) %>%
          dplyr::filter(cs_annot %in% chosen_states) %>%
          dplyr::mutate(
            cs_start = as.integer(cs_start),
            cs_stop = as.integer(cs_stop),
            cs_celltype_ident = celltype_trans[cs_celltype]
          ) %>%
          dplyr::filter(grepl(.celltype, cs_celltype_ident, fixed = T)) %>%
          dplyr::select(cs_annot) %>%
          dplyr::mutate(cs_annot_type = case_when(
            cs_annot %in% promoters ~ "promoter",
            cs_annot %in% enhancers ~ "enhancer"
          ))

        annot <- .regs %$% cs_annot %>%
          unique() %>%
          paste(collapse = "|")
        annot_type <- .regs %$% cs_annot_type %>%
          unique() %>%
          paste(collapse = "|")
      }

      paste(by_promoter, by_degene, annot, annot_type, sep = ";")
    },
    .tss_tab = as.data.table(tss_tab),
    .degene_tab = as.data.table(degene_tab)
    )) %>%
    separate(ANN,
      into = c("by_promoter", "by_degene", "chrom_annot", "chrom_annot_type"),
      sep = ";"
    ) %>%
    rename(chrom_annot_evi = chrom_state) %>%
    relocate(chrom_annot_evi, .after = chrom_annot_type) %T>%
    fwrite(heter_snp_annot_file, row.names = F, sep = ",", quote = F)
}

asoc_annot_file <- file.path(
  wkdir, "pseudo-bulk_per-celltype_asoc-snps_annotated.csv"
)
if (!file.exists(asoc_annot_file)) {
  asoc_snp_file <- file.path(wkdir, "pseudo-bulk_per-celltype_asoc-snps.csv")
  tmp <- fread(asoc_snp_file) %>%
    left_join(heter_snp_annot_tab,
      by = c("contig", "position", "variantID", "condition", "celltype")
    ) %T>%
    fwrite(asoc_annot_file, row.names = F, sep = ",", quote = F)
}


sessionInfo()
# R version 4.0.0 (2020-04-24)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
#
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so
#
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#  [1] biomaRt_2.46.3      progress_1.2.2      ggVennDiagram_1.1.4
#  [4] data.table_1.14.0   forcats_0.5.1       stringr_1.4.0
#  [7] dplyr_1.0.5         purrr_0.3.4         readr_1.4.0
# [10] tidyr_1.1.3         tibble_3.1.0        tidyverse_1.3.1
# [13] magrittr_2.0.1      ggupset_0.3.0       ggrepel_0.9.1
# [16] ggplot2_3.3.3       devtools_2.4.2      usethis_2.0.1
#
# loaded via a namespace (and not attached):
#  [1] nlme_3.1-147         fs_1.5.0             sf_1.0-1
#  [4] bit64_4.0.5          lubridate_1.7.10     httr_1.4.2
#  [7] rprojroot_2.0.2      tools_4.0.0          backports_1.2.1
# [10] utf8_1.2.1           R6_2.5.0             KernSmooth_2.23-16
# [13] BiocGenerics_0.36.0  DBI_1.1.1            mgcv_1.8-31
# [16] colorspace_2.0-0     withr_2.4.1          tidyselect_1.1.0
# [19] prettyunits_1.1.1    processx_3.5.2       curl_4.3
# [22] bit_4.0.4            compiler_4.0.0       Biobase_2.50.0
# [25] cli_2.5.0            rvest_1.0.0          xml2_1.3.2
# [28] desc_1.3.0           labeling_0.4.2       scales_1.1.1
# [31] classInt_0.4-3       askpass_1.1          callr_3.7.0
# [34] proxy_0.4-26         rappdirs_0.3.3       digest_0.6.27
# [37] R.utils_2.10.1       pkgconfig_2.0.3      sessioninfo_1.1.1
# [40] dbplyr_2.1.1         fastmap_1.1.0        rlang_0.4.10
# [43] readxl_1.3.1         RSQLite_2.2.4        rstudioapi_0.13
# [46] farver_2.1.0         generics_0.1.0       jsonlite_1.7.2
# [49] R.oo_1.24.0          Matrix_1.3-4         S4Vectors_0.28.1
# [52] Rcpp_1.0.6           munsell_0.5.0        fansi_0.4.2
# [55] R.methodsS3_1.8.1    lifecycle_1.0.0      stringi_1.6.2
# [58] BiocFileCache_1.14.0 pkgbuild_1.2.0       blob_1.2.1
# [61] grid_4.0.0           parallel_4.0.0       crayon_1.4.1
# [64] lattice_0.20-41      haven_2.4.1          splines_4.0.0
# [67] hms_1.0.0            ps_1.6.0             pillar_1.6.1
# [70] stats4_4.0.0         pkgload_1.2.1        XML_3.99-0.6
# [73] reprex_2.0.0         glue_1.4.2           remotes_2.4.0
# [76] modelr_0.1.8         vctrs_0.3.8          testthat_3.0.2
# [79] cellranger_1.1.0     openssl_1.4.3        gtable_0.3.0
# [82] assertthat_0.2.1     cachem_1.0.4         broom_0.7.6
# [85] e1071_1.7-7          class_7.3-16         RVenn_1.1.0
# [88] IRanges_2.24.1       AnnotationDbi_1.52.0 memoise_2.0.0
# [91] units_0.7-2          ellipsis_0.3.2
