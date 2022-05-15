#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 27, 2021
# Updated: Feb 09, 2022

# Enrichment analysis
#   1. The backgroud is all heterozygous sites, including X (TODO: add NR. of
#      SNPs) SNPs, of the cohort

options(stringsAsFactors = F)

# Databases
library(BSgenome)
library(org.Hs.eg.db)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

# Packages
library(VariantAnnotation)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(tidyverse)
library(magrittr)
library(snpStats)
library(ggrepel)
library(ggupset)

#' Format hypothesis test results into a data.frame
#'
as.data.frame.htest <- function(htest, test = "fisher.test") {
  return(data.frame(
    p.value = as.vector(htest$p.value),
    ci.left = as.vector(htest$conf.int[1]),
    ci.right = as.vector(htest$conf.int[2]),
    odds.ratio = as.vector(htest$estimate),
    test = test
  ))
}


#' Assign GWAS information to the ASoC SNPs
#'
#' @param data The input SNP annotation results to be merged with GWAS
#' @param snp_tab The GWAS summary statistics.
#' @param pval_col The colnames of p-values in the GWAS summary statistics.
#' @return A `data.frame` including ASoC SNPs.
#' @examples
assign_gwas_snp <- function(data, snp_tab, pval_col = "p_value", alpha = 5e-5,
                            asoc_rsid_col = "variantID", rsid_col = NULL,
                            snplocs = NULL, chrom_col = NULL, pos_col = NULL,
                            pick_flank = 0) {
  snp_tab %<>%
    (function(d) {
      cat(base::paste0("I: ", nrow(d), " items in total.\n"))
      d
    }) %>%
    dplyr::filter(!!rlang::sym(pval_col) <= alpha) %>%
    (function(d) {
      cat(base::paste0("I: ", nrow(d), " remain at ", alpha, "\n"))
      d
    })

  if (is.null(rsid_col)) {
    if (any(is.null(snplocs), is.null(chrom_col), is.null(pos_col))) {
      stop("When no rsid is available, snplocs, chrom_col, and pos_col should be given.")
    }

    rsid_col <- "rsid"
    snp_tab %<>%
      dplyr::mutate(rsid = mapply(function(chrom, pos, add_chr = F) {
        if (add_chr) chrom <- base::paste0("chr", chrom)
        reg <- base::paste0(chrom, ":", pos, "-", pos)
        snp <- snpsByOverlaps(snplocs, reg)$RefSNP_id
        ifelse(is.character(snp), snp, "NA")
      }, chrom = !!rlang::sym(chrom_col), pos = !!rlang::sym(pos_col), SIMPLIFY = T))
  }

  names(asoc_rsid_col) <- rsid_col
  snp_tab %<>%
    dplyr::filter(!!rlang::sym(rsid_col) != "NA") %>%
    (function(d) {
      cat(base::paste0("I: ", nrow(d), " items have rsid.\n"))
      d
    }) %>%
    dplyr::inner_join(data, by = asoc_rsid_col, suffix = c("_gwas", ""))

  return(snp_tab)
}


#' Estimate the over-representation of SNPs in the compiled GWAS SNP clumps
#'
#' @param query_snps A `data.frame` or `data.table` including rsid and its chromosome.
#' @param snp_bin A set of SNPs tha compiled in to clumps.
#' @param pair_r2 The minimum R-sequre score to assigne a query SNP to clump. Only work when the query SNPs are not in the clump SNPs.
#' @param pair_dist The maximum physical distance in base-pair to exclude SNPs in the clumps. Only work when the query SNPs are not in the clump SNPs.
#' @param use_bin_snps_proxy Use the proxy SNPs of non-top SNP in the clump to assign the query SNP to the clump. If `TRUE` the query_snps must have a chrom column.
#' @param snp_locs The version of genomic coordinations.
#'
#' @return `data.frame`
#' @examples
enrich_snps <- function(query_snps, snp_bin, pair_r2 = 0.6, pair_dist = 250000,
                        n_bk_snps = 5e6, vcffile = NULL, snp_locs = NULL,
                        use_bin_snps_proxy = F) {
  # Add a column including flag that the query SNPs belongs the clump.
  vcfhandle <- VcfFile(vcffile)
  snp_bin %<>%
    dplyr::mutate(flag = apply(., 1, function(x, .query_snps) {
      .query_snp_ids <- .query_snps %$% rsid
      bq_snps <- NULL
      bq_count <- 0

      idx_snp_id <- x["SNP"]
      # message("Working on SNP clump indexed by ", idx_snp_id)

      # Top SNP, it's perfect!
      if (idx_snp_id %in% .query_snp_ids) {
        bq_snps <- idx_snp_id
        bq_count <- 1
        .query_snp_ids <- .query_snp_ids[.query_snp_ids != idx_snp_id]
        message("Query SNP is the index SNP of the clump: ", bq_snps)
      }

      # Not top SNP, but belongs to the compiled SNP clump, it's great!
      ld_snps <- x["SP2"] %>%
        str_remove_all("\\([0-9]\\)") %>%
        str_split(",") %>%
        unlist()
      overlaps <- .query_snp_ids %in% ld_snps
      if (any(overlaps)) {
        # print(.query_snp_ids[.query_snp_ids %in% ld_snps])
        bq_snps %<>%
          stringr::str_c(.query_snp_ids[overlaps], sep = ";", collapse = ";") %>%
          stringr::str_remove(pattern = "^;")
        bq_count %<>% sum(overlaps)
        .query_snp_ids <- .query_snp_ids[!overlaps]
        message("Query SNP(s) belong to non index SNPs of the clump: ", bq_snps)
      }

      # Not top SNP, not clump SNP, but in LD with top SNP, it's still great
      cur_chrom <- as.integer(x["CHR"])
      if (pair_r2 != 0) {
        query_snps_cur_chr <- .query_snps %>% dplyr::filter(chrom == cur_chrom) %$% rsid
        snpsreg <- snpsById(snp_locs, ids = c(query_snps_cur_chr, idx_snp_id), ifnotfound = "drop")
        if (length(snpsreg) > 0) {
          snpsvcf <- readVcf(vcfhandle, "GRCh37", param = ScanVcfParam(which = snpsreg))
          snpsmat <- genotypeToSnpMatrix(snpsvcf)
          snpsgen <- snpsmat$genotypes
          if (idx_snp_id %in% snpsmat$map$snp.names) {
            query_snps_cur_chr <- query_snps_cur_chr[query_snps_cur_chr %in% snpsmat$map$snp.names]
            snpsldm <- ld(snpsgen[, idx_snp_id], snpsgen[, query_snps_cur_chr],
              stats = "R.squared",
              depth = length(query_snps_cur_chr)
            )

            which_query_snp <- snpsldm > pair_r2
            if (any(which_query_snp)) {
              bq_snps %<>%
                stringr::str_c(colnames(snpsldm)[which_query_snp], sep = ";", collapse = ";") %>%
                stringr::str_remove(pattern = "^;")
              bq_count %<>% sum(sum(which_query_snp))
              .query_snp_ids <- .query_snp_ids[!which_query_snp]
              message("Query SNP(s) in LD with the index SNP of the clump")
            }
          }
        }
      }

      # Not top SNP, not clump SNP, not in LD with top SNP, but in LD with any clump SNP, it's fair.
      if (use_bin_snps_proxy) {
        query_snps_cur_chr <- .query_snps %>% dplyr::filter(chrom == cur_chrom) %$% rsid
        snpsreg <- snpsById(snp_locs, ids = c(query_snps_cur_chr, ld_snps), ifnotfound = "drop")
        if (length(snpsreg) > 0) {
          snpsvcf <- readVcf(vcfhandle, "GRCh37", param = ScanVcfParam(which = snpsreg))
          snpsmat <- genotypeToSnpMatrix(snpsvcf)
          snpsgen <- snpsmat$genotypes
          if (any(ld_snps %in% snpsmat$map$snp.names)) {
            query_snps_cur_chr <- query_snps_cur_chr[query_snps_cur_chr %in% snpsmat$map$snp.names]
            clump_snps_cur_chr <- ld_snps[ld_snps %in% snpsmat$map$snp.names]
            for (.rsid in clump_snps_cur_chr) {
              snpsldm <- ld(snpsgen[, .rsid], snpsgen[, query_snps_cur_chr],
                stats = "R.squared", depth = length(query_snps_cur_chr)
              )

              which_query_snp <- snpsldm > pair_r2
              if (any(which_query_snp)) {
                bq_snps %<>%
                  stringr::str_c(colnames(snpsldm)[which_query_snp], sep = ";", collapse = ";") %>%
                  stringr::str_remove(pattern = "^;")
                bq_count %<>% sum(sum(which_query_snp))
                .query_snp_ids <- .query_snp_ids[!which_query_snp]

                message("Query SNP(s) in LD with non index SNPs of the clump")
              }
            }
          }
        }
      }

      record <- stringr::str_c(bq_count, bq_snps, sep = "|")
      record
    }, .query_snps = query_snps)) %>%
    tidyr::separate(flag, into = c("bq_count", "bq_snps"), sep = "\\|") %>%
    dplyr::mutate(bq_count = as.integer(bq_count))

  n_binned_query_snps <- snp_bin %$% bq_count %>% sum()
  n_binned_snps <- snp_bin %$% SP2 %>%
    paste(collapse = ";") %>%
    str_split(",") %>%
    unlist() %>%
    unique() %>%
    length()
  n_query_snps <- query_snps %>% nrow()

  if (n_binned_query_snps != 0) {
    # (n_binned_query_snps / n_bk_snps) / ((n_binned_snps / n_bk_snps) * (n_query_snps / n_bk_snps))
    ratio <- (n_binned_query_snps / n_query_snps) / (n_binned_snps / n_bk_snps)
    ppt <- prop.test(n_binned_query_snps, n_query_snps, n_binned_snps / n_bk_snps)
    p_value <- ppt$p.value
    ci_left <- ppt$conf.int[1]
    ci_right <- ppt$conf.int[2]
  } else {
    ratio <- 0
    p_value <- 1
    ci_left <- 0
    ci_right <- 1
  }

  return(list(
    ratio = ratio,
    p_value = p_value,
    estimate = n_binned_snps / n_bk_snps,
    ci_left = ci_left,
    ci_right = ci_right,
    data = list(
      n_binned_query_snps = n_binned_query_snps,
      n_binned_snps = n_binned_snps,
      n_query_snps = n_query_snps
    ),
    snp_clumps = snp_bin
  ))
}



###############################################################################
# Important meta information
###############################################################################
promoter_annot <- c("1_TssA", "2_PromU", "3_PromD1", "4_PromD2")
enhancer_annot <- c(
  "9_TxReg", "10_TxEnh5", "11_TxEnh3", "12_TxEnhW",
  "13_EnhA1", "14_EnhA2", "15_EnhAF", "16_EnhW1",
  "17_EnhW2", "18_EnhAc", "19_DNase"
)

merge_patient <- TRUE
if (merge_patient) {
  disease_group <- c("active", "post")
  condition_map <- c(
    # "Active" = "Active", "Post" = "Convalescent",
    # "active" = "Active", "post" = "Convalescent"
    "Active" = "Hospitalized", "Post" = "Convalescent",
    "active" = "Hospitalized", "post" = "Convalescent"
  ) # Convalescent
} else {
  disease_group <- c("Severe", "Mild", "Post") # COVID19 conditions
  condition_map <- c("Severe" = "Seve.", "Mild" = "Mild", "Post" = "Conv.")
}

celltypes <- c("cMono", "ncMono", "CD4T", "CD8T", "NK", "B")

# Working dir, input and output
pjdir <- "~/Documents/projects/wp_covid19_mhh50"
wkdir <- str_glue("{pjdir}/outputs/scATAC-seq/enrichment")

# All heterozygous sites.
heter_snp_file <- str_glue("{pjdir}/outputs/scATAC-seq/summary/pseudo-bulk_per-celltype_heter-snps_annotated.csv")
heter_snp_tab <- fread(heter_snp_file, verbose = F, showProgress = F)

# ASoC SNPs
asoc_snp_file <- str_glue("{pjdir}/outputs/scATAC-seq/summary/pseudo-bulk_per-celltype_asoc-snps_annotated.csv")
asoc_snp_tab <- fread(asoc_snp_file, verbose = F, showProgress = F) %>%
  dplyr::filter(refCount > 2 & altCount > 2)

# FDR threshold
p_val_th <- 0.05



##############################################################################
#  Ratio of condition specificity
##############################################################################
# We grouped ASoC SNPs into cell type specific and shared SNPs (check annotate_asoc_snps.r).
enrich_for_cond_spec <- TRUE
if (enrich_for_cond_spec) {
  asoc_snp_annot_count_tab <- asoc_snp_tab %>%
    dplyr::filter(p_value_adj < p_val_th) %>%
    dplyr::select(variantID, condition, celltype, by_eqtlgene, by_degene, chrom_annot_type) %>%
    dplyr::group_by(variantID, condition) %>%
    dplyr::summarise(ann = {
      cur_data() %>%
        apply(1, function(x) {
          annot <- list(
            "Promoter" = 0, "Enhancer" = 0, "Prom_and_Enh" = 0, "Non_prom_or_Enh" = 0,
            "DE_gene" = 0, "eQTL_gene" = 0, "DE_and_eQTL_gene" = 0
          )

          ann_type <- x["chrom_annot_type"]
          is_prom <- grepl("promoter", ann_type)
          is_enh <- grepl("enhancer", ann_type)
          if (is_prom && is_enh) {
            annot[["Prom_and_Enh"]] <- 1
          } else if (is_prom) {
            annot[["Promoter"]] <- 1
          } else if (is_enh) {
            annot[["Enhancer"]] <- 1
          } else {
            annot[["Non_prom_or_Enh"]] <- 1
          }

          de_gene <- ""
          if (!x["by_degene"] %in% c("", NA)) {
            de_gene <- x["by_degene"] %>%
              str_split(pattern = "[|@:]", simplify = T) %>%
              as.data.frame() %>%
              set_colnames(c("gene", "dire", "comp")) %$%
              gene
          }

          eq_gene <- ""
          if (!x["by_eqtlgene"] %in% c("", NA)) {
            eq_gene <- x["by_eqtlgene"] %>%
              str_split(pattern = "\\|") %>%
              unlist()
          }

          has_de_gene <- ifelse(sum(de_gene != ""), T, F)
          has_eq_gene <- ifelse(sum(eq_gene != ""), T, F)
          share_de_eq <- any(intersect(de_gene, eq_gene) != "")
          if (share_de_eq) {
            annot[["DE_and_eQTL_gene"]] <- 1
          } else if (has_de_gene || has_eq_gene) {
            if (has_de_gene) annot[["DE_gene"]] <- 1
            if (has_eq_gene) annot[["eQTL_gene"]] <- 1
          }

          annot
        }) %>%
        do.call(rbind.data.frame, .) %>%
        dplyr::summarise_all(max)
    }) %>%
    as.data.table() %>%
    tidyr::pivot_longer(
      cols = starts_with("ann."), names_to = "Annotation", values_to = "Check"
    )

  anncounts_for_upset <- asoc_snp_annot_count_tab %>%
    dplyr::filter(Check == 1) %>%
    dplyr::mutate(
      Annotation = str_remove_all(Annotation, "ann."),
      condition = if_else(condition == "post", "Convalescent", "Hospitalized")
    ) %>%
    dplyr::group_by(variantID, condition) %>%
    dplyr::summarise(Annotation = list(Annotation))

  # FIXME: A special case affected by cell type. It's assigned to only Enhancer
  #        in cMono, but assigned to both enhancer and promoter in CD4T

  #  variantID  condition count
  #  <chr>      <chr>     <chr>
  #  rs5995473  active    Enhancer;Prom_and_Enh;DE_gene;eQTL_gene
  asoc_snp_annot_count_tab %>%
    dplyr::group_by(variantID, condition) %>%
    dplyr::mutate(Annotation = str_remove_all(Annotation, "ann.")) %>%
    dplyr::summarise(annot = {
      cur_data() %>%
        dplyr::filter(Check == 1) %$%
        Annotation %>%
        paste0(collapse = ";")
    }) %>%
    dplyr::filter(
      base::grepl("eQTL_gene", annot) & base::grepl("Enhancer", annot) &
        base::grepl("DE_gene", annot) & base::grepl("Prom_and_Enh", annot)
    )

  asoc_snp_annot_count_tab %>%
    dplyr::mutate(Annotation = str_remove_all(Annotation, "ann.")) %>%
    tidyr::pivot_wider(names_from = "Annotation", values_from = "Check") %>%
    dplyr::mutate(
      Prom_DE_EQ = if_else((Promoter == 1 | Prom_and_Enh == 1) & DE_and_eQTL_gene == 1, 1, 0)
    ) %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(across(-variantID, sum)) %>%
    dplyr::mutate(
      nr_asoc_snp = c(292, 86),
      Prom_DE_EQ_pct = Prom_DE_EQ / nr_asoc_snp,
      Promoter_pct = (Promoter + Prom_and_Enh) / nr_asoc_snp,
      Enhancer_pct = (Enhancer + Prom_and_Enh) / nr_asoc_snp,
      eQTL_gene_pct = (eQTL_gene + DE_and_eQTL_gene) / nr_asoc_snp,
      DE_gene_pct = (DE_gene + DE_and_eQTL_gene) / nr_asoc_snp
    ) %>%
    as.data.frame()


  g <- anncounts_for_upset %>%
    ggplot(aes(x = Annotation)) +
    geom_bar() +
    geom_text(aes(label = after_stat(count)), size = 4.5, stat = "count", vjust = -1) +
    scale_x_upset(order = "degree") +
    scale_y_continuous(breaks = NULL, lim = c(NA, 75), name = "") +
    facet_wrap(~condition, nrow = 2) +
    theme_classic() +
    theme(
      plot.margin = margin(0.5, 0.5, 0.5, 3.5, "cm"),
      text = element_text(size = 16),
      axis.text.y = element_text(size = 12)
    )

  save_fig_to <- str_glue("{wkdir}/regulation_annotation.pdf")
  ggsave(save_fig_to, plot = g, width = 7, height = 5.5)


  #############################################################################
  # Check promoter enrichment
  #############################################################################
  bk_tab <- heter_snp_tab %>%
    dplyr::select(condition, variantID, chrom_annot_type) %>%
    dplyr::group_by(condition, variantID) %>%
    dplyr::summarise(chrom_annot_type = {
      nr_ann <- cur_data() %>% base::nrow()
      ifelse(nr_ann > 1, paste(chrom_annot_type, collapse = "|"), chrom_annot_type)
    }) %>%
    dplyr::distinct() %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(
      nr_enh_snp = sum(grepl("enhancer", chrom_annot_type)),
      nr_prom_snp = sum(grepl("promoter", chrom_annot_type)),
      nr_obs_snp = cur_data() %$% variantID %>% unique() %>% length(),
      type = "Bkg"
    )

  ob_tab <- asoc_snp_tab %>%
    dplyr::filter(p_value_adj < p_val_th) %>%
    dplyr::select(condition, variantID, chrom_annot_type) %>%
    dplyr::distinct() %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(
      nr_enh_snp = sum(grepl("enhancer", chrom_annot_type)),
      nr_prom_snp = sum(grepl("promoter", chrom_annot_type)),
      nr_obs_snp = cur_data() %$% variantID %>% unique() %>% length(),
      type = "Obs"
    )

  rbind(ob_tab, bk_tab) %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(ppt = {
      ob_count <- cur_data() %>% dplyr::filter(type == "Obs")
      n <- ob_count["nr_obs_snp"] %>% unlist()
      x <- ob_count["nr_enh_snp"] %>% unlist()

      bk_count <- cur_data() %>% dplyr::filter(type == "Bkg")
      N <- bk_count["nr_obs_snp"] %>% unlist()
      X <- bk_count["nr_enh_snp"] %>% unlist()

      ppt_enh <- prop.test(x, n, p = X / N) %>%
        as.data.frame() %>%
        dplyr::mutate(
          observed = x / n, x = x, n = n, expected = X / N, X = X, N = N,
          epig = "enhancer"
        )

      x <- ob_count["nr_prom_snp"] %>% unlist()
      X <- bk_count["nr_prom_snp"] %>% unlist()
      ppt_prom <- prop.test(x, n, p = X / N) %>%
        as.data.frame() %>%
        dplyr::mutate(
          observed = x / n, x = x, n = n, expected = X / N, X = X, N = N,
          epig = "promoter"
        )

      rbind(ppt_enh, ppt_prom)
    }) %>%
    as.data.table() %>%
    dplyr::mutate(p_value_adj = p.adjust(ppt.p.value, method = "holm"))
}


###############################################################################
# Enrichment of ASoC SNPs assigned to DE/eQTL genes
###############################################################################
enrich_for_de_genes <- TRUE
if (enrich_for_de_genes) {
  nr_heter_snp <- heter_snp_tab %>%
    dplyr::select(-c(chrom_annot_evi)) %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(
      n_snps = n(),
      n_snps_degene = sum(by_degene != ""),
      n_snps_eqtlgene = sum(by_eqtlgene != "")
    )

  nr_asoc_snp <- asoc_snp_tab %>%
    dplyr::select(-c(chrom_annot_evi)) %>%
    dplyr::group_by(condition, celltype) %>%
    dplyr::filter(p_value_adj < p_val_th) %>%
    dplyr::summarise(
      n_snps = n(),
      n_snps_degene = sum(by_degene != ""),
      n_snps_eqtlgene = sum(by_eqtlgene != "")
    )

  ph_dtfm <- data.frame(
    condition = c("post", "active", "post"),
    celltype = c("B", "ncMono", "ncMono"),
    p.value.degene = rep(1, 3),
    ci.left.degene = rep(0, 3),
    ci.right.degene = rep(0, 3),
    odds.ratio.degene = rep(0, 3),
    observed.degene = rep(0, 3),
    n.degene = rep(0, 3),
    x.degene = rep(0, 3),
    p.value.eqtlgene = rep(1, 3),
    ci.left.eqtlgene = rep(0, 3),
    ci.right.eqtlgene = rep(0, 3),
    odds.ratio.eqtlgene = rep(0, 3),
    observed.eqtlgene = rep(0, 3),
    n.eqtlgene = rep(0, 3),
    x.eqtlgene = rep(0, 3),
    expected = rep(0, 3),
    or.degene = rep(0, 3),
    log2.fc.degene = rep(0, 3),
    or.eqtlgene = rep(0, 3),
    log2.fc.eqtlgene = rep(0, 3),
    degene_label = rep("<NA>", 3),
    eqtlgene_label = rep("<NA>", 3)
  )


  # proportion test for Nr. of SNPs belongs to a DE gene by scRNA-seq
  ppt <- full_join(nr_asoc_snp, nr_heter_snp,
    by = "condition", suffix = c(".asoc", ".heter")
  ) %>%
    dplyr::group_by(condition, celltype) %>%
    dplyr::summarise(ht = {
      X <- cur_data() %$% n_snps_degene.heter
      N <- cur_data() %$% n_snps.heter

      x <- cur_data() %$% n_snps_degene.asoc
      n <- cur_data() %$% n_snps.asoc
      ppt_degene <- prop.test(x, n, X / N) %>%
        as.data.frame() %>%
        dplyr::mutate(observed = x / n, n = n, x = x) %>%
        dplyr::rename_with(.fn = function(n) paste0(n, ".degene"))

      x <- cur_data() %$% n_snps_eqtlgene.asoc
      n <- cur_data() %$% n_snps.asoc
      ppt_eqtlgene <- prop.test(x, n, X / N) %>%
        as.data.frame() %>%
        dplyr::mutate(observed = x / n, n = n, x = x) %>%
        dplyr::rename_with(.fn = function(n) paste0(n, ".eqtlgene"))

      cbind(ppt_degene, ppt_eqtlgene) %>%
        dplyr::select(-c(starts_with("test"))) %>%
        dplyr::mutate(
          expected = X / N,
          or.degene = observed.degene / expected,
          log2.fc.degene = log2(or.degene),
          or.eqtlgene = observed.eqtlgene / expected,
          log2.fc.eqtlgene = log2(or.eqtlgene)
        )
    }) %>%
    as.data.table() %>%
    dplyr::rename_with(.fn = function(e) str_remove_all(e, "^ht.")) %>%
    dplyr::mutate(
      p.value.degene = p.adjust(p.value.degene, method = "fdr"),
      p.value.eqtlgene = p.adjust(p.value.eqtlgene, method = "fdr"),
      degene_label = paste0(
        if_else(p.value.degene < 0.01,
          format(p.value.degene, trim = TRUE, digits = 2),
          format(round(p.value.degene, 2), nsmall = 2, scientific = FALSE)
        ),
        ", ", x.degene, "/", n.degene
      ),
      eqtlgene_label = paste0(
        if_else(p.value.eqtlgene < 0.01,
          format(p.value.eqtlgene, trim = TRUE, digits = 2),
          format(round(p.value.eqtlgene, 2), nsmall = 2, scientific = FALSE)
        ),
        ", ", x.eqtlgene, "/", n.eqtlgene
      ),
      log2.fc.degene = if_else(
        is.infinite(log2.fc.degene), sign(log2.fc.degene) * 0.5, log2.fc.degene
      ),
      log2.fc.eqtlgene = if_else(
        is.infinite(log2.fc.eqtlgene), sign(log2.fc.eqtlgene) * 0.5, log2.fc.eqtlgene
      )
    ) %>%
    rbind(ph_dtfm) %>%
    dplyr::mutate(celltype = factor(celltype, celltypes))

  # Enrichement for DE genes from scRNA-seq of the current study.
  g_ppt <- ppt %>%
    ggplot(aes(x = celltype, fill = condition)) +
    geom_col(aes(y = or.degene), position = position_dodge(0.9)) +
    geom_text(aes(y = 0, label = degene_label),
      angle = 90, size = 8, hjust = 0,
      position = position_dodge(0.9)
    ) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    scale_fill_discrete(name = "Condition", labels = condition_map) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 15),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
  save_fig_to <- str_glue("{wkdir}/enrich_asoc_snp_for_degene.pdf")
  ggsave(save_fig_to, plot = g_ppt, width = 9, height = 6)

  # Enrichment for eQTL genes from eQTLGen
  g_ppt <- ppt %>%
    ggplot(aes(x = celltype, fill = condition)) +
    geom_col(aes(y = or.eqtlgene), position = position_dodge(0.9)) +
    geom_text(aes(y = 0, label = eqtlgene_label),
      angle = 90, size = 8, hjust = 0,
      position = position_dodge(0.9)
    ) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    scale_fill_discrete(name = "Condition", labels = condition_map) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 15),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
  save_fig_to <- str_glue("{wkdir}/enrich_asoc_snp_for_eqtlgene.pdf")
  ggsave(save_fig_to, plot = g_ppt, width = 9, height = 6)
}



###############################################################################
#  Enrichment for COVID19 GWAS SNPs
###############################################################################
enrich_for_covid19_gwas <- TRUE
if (enrich_for_covid19_gwas) {
  # COVID19_HGI_A2_ALL_leave_23andme_20210607_5e-2.txt: Very severe respiratory confirmed covid vs. population
  # COVID19_HGI_B1_ALL_leave_23andme_20210607_5e-2.txt: Hospitalized covid vs.not hospitalized covid
  # COVID19_HGI_B2_ALL_leave_23andme_20210607_5e-2.txt: Hospitalized covid vs. population
  # COVID19_HGI_C2_ALL_leave_23andme_20210607_5e-2.txt: Covid vs. population
  gwas_ver <- "B2" # Different COVID-19 summary statistics.
  gwas_snp_file <- str_glue("{pjdir}/inputs/covid19_gwas/covid19hg.org/COVID19_HGI_{gwas_ver}_ALL_leave_23andme_20210607.txt.gz") # Covid vs. population)
  gwas_snp_tab <- fread(gwas_snp_file, tmpdir = dirname(as.character(gwas_snp_file))) %>% dplyr::rename("CHR" = "#CHR")

  asoc_vs_covid19_tab <- asoc_snp_tab %>%
    dplyr::filter(!celltype %in% c("Plasmablast") & p_value < 0.05) %>%
    dplyr::select(c(contig, position, variantID, p_value, celltype, condition)) %>%
    dplyr::inner_join(dplyr::select(gwas_snp_tab, rsid, all_inv_var_meta_p), by = c("variantID" = "rsid")) %>%
    dplyr::rename("p_value_gwas" = "all_inv_var_meta_p", "p_value_asoc" = "p_value") %>%
    dplyr::mutate(condition = factor(condition, levels = disease_group), celltype = factor(celltype, levels = celltypes)) %>%
    dplyr::group_by(condition, celltype) %>%
    dplyr::arrange(p_value_gwas) %>%
    dplyr::mutate(
      observed = -log10(p_value_gwas),
      expected = -log10(ppoints(n())),
      ci_upper = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = seq(n()), shape2 = rev(seq(n())))),
      ci_lower = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = seq(n()), shape2 = rev(seq(n()))))
    ) %>%
    dplyr::ungroup()

  asoc_snp_af <- gwas_snp_tab %>%
    dplyr::select(rsid, all_meta_AF) %>%
    dplyr::filter(rsid %in% asoc_vs_covid19_tab$variantID) %>%
    dplyr::rename("variantID" = rsid)

  set.seed(31415926)
  asoc_vs_covid19_tab <- gwas_snp_tab %>%
    dplyr::select(rsid, all_meta_AF, all_inv_var_meta_p) %>%
    dplyr::rename("rand_gwas_snp" = "rsid", "rand_gwas_snp_p_value" = "all_inv_var_meta_p") %>%
    dplyr::filter((!rand_gwas_snp %in% asoc_snp_af$variantID) & all_meta_AF %in% unique(asoc_snp_af$all_meta_AF)) %>%
    dplyr::group_by(all_meta_AF) %>%
    dplyr::slice_sample(n = 1) %>%
    dplyr::mutate(condition = "random") %>%
    dplyr::ungroup() %>%
    dplyr::right_join(asoc_snp_af, by = "all_meta_AF", suffix = c("", ".x")) %>%
    dplyr::right_join(asoc_vs_covid19_tab, by = "variantID", suffix = c("", ".x")) %>%
    dplyr::select(-c(variantID, all_meta_AF, p_value_gwas:ci_lower), -ends_with(".x")) %>%
    dplyr::rename("variantID" = "rand_gwas_snp", "p_value_gwas" = "rand_gwas_snp_p_value") %>%
    dplyr::group_by(condition, celltype) %>%
    dplyr::arrange(p_value_gwas) %>%
    dplyr::mutate(
      observed = -log10(p_value_gwas),
      expected = -log10(ppoints(n())),
      ci_upper = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = seq(n()), shape2 = rev(seq(n())))),
      ci_lower = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = seq(n()), shape2 = rev(seq(n()))))
    ) %>%
    dplyr::select(
      contig, position, variantID, p_value_asoc, celltype, condition, p_value_gwas,
      observed:ci_lower
    ) %>%
    rbind(asoc_vs_covid19_tab) %>%
    dplyr::mutate(condition = factor(condition, levels = c(disease_group, "random")))


  # Are the distributions different?
  asoc_vs_covid19_tab %>%
    dplyr::group_by(condition, celltype) %>%
    dplyr::summarise(ks = {
      obs <- cur_data() %$% observed
      exp <- cur_data() %$% expected
      rdm <- dplyr::filter(., condition == "random") %$% observed
      kt_rdm <- ks.test(x = obs, y = rdm, "pbeta")
      kt_exp <- ks.test(x = obs, y = exp, "pbeta")
      data.frame(
        p_val_exp = kt_exp$p.value,
        stat_exp = kt_exp$statistic,
        p_val_rdm = kt_rdm$p.value,
        stat_rdm = kt_rdm$statistic
      )
    }) %>%
    as.data.table() %>%
    dplyr::rename_with(starts_with("ks."), .fn = function(x) str_remove(x, "ks."))


  # cMono
  g <- asoc_vs_covid19_tab %>%
    dplyr::filter(celltype == "cMono") %>%
    ggplot(aes(x = expected, y = observed)) +
    geom_point(aes(color = condition), size = 1) +
    geom_line(aes(color = condition), size = 0.05) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    # geom_ribbon(aes(ymax=ci_upper, ymin=ci_lower), fill="grey50", alpha=0.5) +
    geom_text_repel(aes(label = variantID), . %>% filter(observed >= -log10(1e-5)),
      min.segment.length = 0, nudge_x = 0.2, size = 5
    ) +
    labs(
      x = quote("Expected" ~ -Log[10] ~ "(p-value)"),
      y = quote("Observed" ~ -Log[10] ~ "(p-value) from COVID-19 GWAS")
    ) +
    scale_color_discrete(name = "Condition", labels = condition_map) +
    theme_classic()
  save_to <- str_glue("{wkdir}/qqplot.asoc_vs_COVID19_HGI_{gwas_ver}.p_value_gwas.cMono.pdf")
  ggsave(save_to, plot = g, width = 5.5, height = 5)

  # All celltype
  g <- asoc_vs_covid19_tab %>%
    ggplot(aes(x = expected, y = observed)) +
    geom_point(aes(color = condition), size = 1) +
    geom_line(aes(color = condition), size = 0.05) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    # geom_ribbon(aes(ymax=ci_upper, ymin=ci_lower), fill="grey50", alpha=0.5) +
    geom_text_repel(aes(label = variantID), . %>% filter(observed >= -log10(1e-5)),
      min.segment.length = 0, nudge_x = 0.2, size = 5
    ) +
    facet_wrap(~celltype, nrow = 2) +
    labs(
      x = quote("Expected" ~ -Log[10] ~ "(p-value)"),
      y = quote("Observed" ~ -Log[10] ~ "(p-value) from COVID-19 GWAS")
    ) +
    scale_color_discrete(name = "Condition", labels = condition_map) +
    theme_classic() +
    theme(
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 14)
    )
  save_to <- str_glue("{wkdir}/qqplot.asoc_vs_COVID19_HGI_{gwas_ver}.p_value_gwas.pdf")
  ggsave(save_to, plot = g, width = 10, height = 6)
}



###############################################################################
#  Enrichment of OCR
###############################################################################
enrich_for_ocr <- TRUE
# TODO: The heter_snps_in_peaks_uniq.txt should be updated.
if (enrich_for_ocr) {
  ocr_snp_tab <- NULL
  for (cond in disease_group) {
    ocr_snp_file <- str_glue("{pjdir}/inputs/scATAC-seq/combined_{cond}.txt")
    if (is.null(ocr_snp_tab)) {
      ocr_snp_tab <- fread(ocr_snp_file, verbose = F, showProgress = F) %>%
        dplyr::mutate(condition = cond)
    } else {
      ocr_snp_tab <- fread(ocr_snp_file, verbose = F, showProgress = F) %>%
        dplyr::mutate(condition = cond)
      rbind(ocr_snp_tab)
    }
  }

  # Integrating SNP information
  ocr_annot_tab <- heter_snp_tab %>%
    dplyr::select(contig, position, variantID, celltype, condition) %>%
    dplyr::left_join(ocr_snp_tab %>% dplyr::select(RSid, Peak_pos, peak2gene, SNP_in_motif),
      by = c("variantID" = "RSid"), na_matches = "never"
    ) %>%
    dplyr::left_join(asoc_snp_tab %>% dplyr::select(variantID, celltype, condition, p_value_adj),
      by = c("variantID", "celltype", "condition"), na_matches = "never"
    ) %>%
    dplyr::mutate(
      SNP_in_motif = if_else(is.na(SNP_in_motif), "no_motif", SNP_in_motif),
      peak2gene = if_else(is.na(peak2gene), "no_gene", peak2gene),
      Peak_pos = if_else(is.na(Peak_pos), "no_peak", Peak_pos),
      p_value_adj = if_else(is.na(p_value_adj), 1, p_value_adj)
    )

  ocr_enr_tab <- ocr_annot_tab %>%
    dplyr::group_by(condition, celltype) %>%
    dplyr::summarise(enr = {
      ct_tbl <- cur_data() %>%
        dplyr::mutate(
          is_asoc = factor(p_value_adj < p_val_th, levels = c(T, F)),
          is_peak = factor(Peak_pos != "no_peak", levels = c(T, F))
        ) %>%
        dplyr::select(is_peak, is_asoc) %>%
        base::table()

      est_summ <- data.frame(
        p_value = 1, odds_ratio = 1, ci_left = 1, ci_right = 1,
        is_peak.is_asoc = ct_tbl[1, 1], # is_asoc and is_xxx
        is_peak.nt_asoc = ct_tbl[1, 2], # not_asoc and is_xxx
        nt_peak.is_asoc = ct_tbl[2, 1], # is_asoc and not_xxx
        nt_peak.nt_asoc = ct_tbl[2, 2] # not_asoc and not_xxx
      )

      if (all(ct_tbl > 0)) {
        est <- fisher.test(ct_tbl, simulate.p.value = T)
        est_summ <- data.frame(
          p_value = est$p.value,
          odds_ratio = est$estimate,
          ci_left = min(est$conf.int),
          ci_right = max(est$conf.int),
          is_peak.is_asoc = ct_tbl[1, 1], # is_asoc and is_xxx
          is_peak.nt_asoc = ct_tbl[1, 2], # not_asoc and is_xxx
          nt_peak.is_asoc = ct_tbl[2, 1], # is_asoc and not_xxx
          nt_peak.nt_asoc = ct_tbl[2, 2] # not_asoc and not_xxx
        )
      }
      est_summ
    }) %>%
    as.data.table() %>%
    dplyr::rename_with(starts_with("enr."), .fn = function(x) str_remove(x, "enr.")) %>%
    dplyr::filter(celltype %in% celltypes) %>%
    dplyr::mutate(
      p_value_adj = p.adjust(p_value, method = "fdr"),
      condition = factor(condition, levels = disease_group),
      celltype = factor(celltype, levels = celltypes),
      Label = case_when(
        p_value_adj < 0.001 ~ format(p_value_adj, digits = 3, scientific = T),
        p_value_adj < 0.01 ~ format(round(p_value_adj, 3), scientific = F),
        T ~ format(round(p_value_adj, 1), scientific = F)
      )
    )

  g <- ocr_enr_tab %>%
    ggplot(aes(x = celltype, fill = condition)) +
    geom_col(aes(y = odds_ratio),
      position = position_dodge(width = 0.9), color = "white"
    ) +
    geom_text(aes(y = 0, label = Label),
      angle = 90, size = 9, hjust = 0, position = position_dodge(width = 0.9)
    ) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    scale_fill_discrete(name = "Condition", labels = condition_map) +
    labs(y = NULL, x = NULL) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 15),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
  save_fig_to <- str_glue("{wkdir}/ocr_enrichment_pcpd.pdf")
  ggsave(save_fig_to, plot = g, width = 9, height = 6)


  # atac_motif_enr_lst <- c(
  #   "SPIC", "SPIB", "POU2F1", "POU2F3", "POU5F1B", "POU2F2", "POU3F4", "ARID3A",
  #   "NR2F1", "BCL11A", "BCL11B", "STAT2", "IRF4", "SPI1", "RUNX1", "RUNX2", "SP2",
  #   "KLF5", "KLF4", "ZNF148", "SP1", "SP3", "WT1", "FOS", "FOSL2", "JUN", "JUNB",
  #   "FOSL1", "JUND", "CEBPD", "CEBPA"
  # )
  # atac_motif_lst <- ocr_snp_tab %$% SNP_in_motif %>%
  #   paste0(collapse = ";") %>%
  #   str_split(pattern = ";", simplify = T) %>%
  #   str_remove(pattern = "_.*") %>%
  #   unlist() %>%
  #   unique() %>%
  #   purrr::discard(.p = function(x) x %in% c("", "no", "T"))

  # motif_enr_file <- str_glue("{wkdir}/motif_enr_est.csv")
  # if (file.exists(motif_enr_file)) {
  #   motif_enr_tab <- fread(motif_enr_file, verbose = F, showProgress = F)
  # } else {
  #   motif_enr_tab <- NULL
  #   for (mt_gene in atac_motif_lst) {
  #     .mt_gene <- paste0("\\<", mt_gene, "\\>")
  #     motif_enr_tab <- ocr_annot_tab %>%
  #       dplyr::group_by(condition, celltype) %>%
  #       dplyr::summarise(enr = {
  #         ct_tbl <- cur_data() %>%
  #           dplyr::mutate(
  #             is_asoc = factor(p_value_adj < p_val_th, levels = c(T, F)),
  #             is_motif = factor(base::grepl(mt_gene, SNP_in_motif, fixed = T),
  #               levels = c(T, F)
  #             )
  #           ) %>%
  #           dplyr::select(is_asoc, is_motif) %>%
  #           base::table()

  #         est_summ <- data.frame(
  #           p_value = NA,
  #           odds_ratio = NA,
  #           motif = mt_gene,
  #           is_motif.is_asoc = ct_tbl[1, 1], # is_asoc and is_xxx
  #           is_motif.nt_asoc = ct_tbl[1, 2], # not_asoc and is_xxx
  #           nt_motif.is_asoc = ct_tbl[2, 1], # is_asoc and not_xxx
  #           nt_motif.nt_asoc = ct_tbl[2, 2]
  #         ) # not_asoc and not_xxx
  #         if (all(ct_tbl >= 3)) {
  #           est <- fisher.test(ct_tbl)
  #           est_summ <- data.frame(
  #             p_value = est$p.value,
  #             odds_ratio = est$estimate,
  #             motif = mt_gene,
  #             is_motif.is_asoc = ct_tbl[1, 1], # is_asoc and is_xxx
  #             is_motif.nt_asoc = ct_tbl[1, 2], # not_asoc and is_xxx
  #             nt_motif.is_asoc = ct_tbl[2, 1], # is_asoc and not_xxx
  #             nt_motif.nt_asoc = ct_tbl[2, 2] # not_asoc and not_xxx
  #           )
  #         }
  #         est_summ
  #       }) %>%
  #       as.data.table() %>%
  #       dplyr::rename_with(starts_with("enr."),
  #         .fn = function(x) str_remove(x, pattern = "enr.")
  #       ) %>%
  #       rbind(motif_enr_tab)
  #   }
  #   motif_enr_tab %>% fwrite(motif_enr_file, sep = ",", quote = F, row.names = F)
  # }

  # cur_ctype <- "cMono"
  # cond_colors <- c(
  #   "Other" = "black", "active" = "#F8766D", "post" = "#7CAE00"
  # )
  # for (cur_ctype in celltypes) {
  #   motif_enr_for_plot_tab <- motif_enr_tab %>%
  #     dplyr::filter(celltype == cur_ctype & !is.na(p_value)) %>%
  #     dplyr::mutate(
  #       fdr = p.adjust(p_value, method = "fdr"),
  #       condition = if_else(fdr < 0.05, condition, "Other"),
  #       condition = factor(condition, levels = c(disease_group, "Other")),
  #       label = if_else(fdr < 0.05 & motif %in% atac_motif_enr_lst, motif, "")
  #     )

  #   g <- motif_enr_for_plot_tab %>%
  #     ggplot(aes(x = -log10(p_value), y = log2(odds_ratio))) +
  #     geom_point(aes(fill = condition, color = condition),
  #       size = 1.5,
  #       alpha = 0.7, shape = 21
  #     ) +
  #     geom_text_repel(aes(label = label), min.segment.length = 0, max.overlaps = 25) +
  #     theme_classic() +
  #     labs(
  #       x = quote(~ -Log[10] ~ "(p-value) by Fisher's exact test"),
  #       y = quote(~ -Log[2] ~ "(odds ratio)")
  #     ) +
  #     scale_fill_manual(name = "Condition", values = cond_colors, labels = condition_map) +
  #     scale_color_manual(name = "Condition", values = cond_colors, labels = condition_map)

  #   save_to <- str_glue("{wkdir}/motif_enrichment_{cur_ctype}.pdf")
  #   ggsave(save_to, plot = g, width = 7.5, height = 6.5)
  # }
}



###############################################################################
# Enrichment of VEP annotation
###############################################################################
enrich_for_vep_annot <- FALSE
if (enrich_for_vep_annot) {
  annot_enrich_file <- file.path(pjdir, "outputs/scATAC-seq/vep/regulatory_feature_enrichment.csv")
  if (file.exists(annot_enrich_file)) {
    annot_enrich_tab <- fread(annot_enrich_file, verbose = F, showProgress = F)
  } else {
    vep_annot_file <- file.path(pjdir, "outputs/scATAC-seq/vep/heter_snps_vep.txt.gz")
    vep_annot_tab <- data.table::fread(vep_annot_file, showProgress = F, verbose = F) %>%
      dplyr::rename("variantID" = "#Uploaded_variation") %>%
      dplyr::group_by(variantID)

    tf_sets <- vep_annot_tab %>%
      dplyr::filter(TRANSCRIPTION_FACTORS != "-") %$%
      TRANSCRIPTION_FACTORS %>%
      stringr::str_split(pattern = ",") %>%
      unlist() %>%
      unique()

    vep_annot_map <- vep_annot_tab %>%
      dplyr::summarise(
        is_enh = any(base::grepl("enhancer", BIOTYPE)),
        is_ocr = any(base::grepl("open_chromatin_region", BIOTYPE)),
        is_prom = any(base::grepl("promoter", BIOTYPE)),
        CTCFBS = any(base::grepl("CTCF_binding_site", BIOTYPE)),
        TFBS = any(base::grepl("TF_binding_site_variant", Consequence))
      )

    # NOTE: 1. Dimers will not be counted when only try to fetch a specific TF.
    for (tf_sym in tf_sets) {
      vep_annot_map <- vep_annot_tab %>%
        dplyr::summarise("{tf_sym}" := any(base::grepl(TRANSCRIPTION_FACTORS,
          pattern = paste0("(^|,)", tf_sym, "(,|$)")
        ))) %>%
        dplyr::full_join(vep_annot_map, by = "variantID")
    }

    annot_enrich_tab <- NULL
    all_celltype <- c("cMono", "ncMono", "B", "CD4T", "CD8T", "NK")
    for (cur_cond in disease_group) {
      asoc_snp_pcd_tab <- asoc_snp_tab %>%
        select(variantID, condition, celltype, p_value_adj) %>%
        filter(condition == cur_cond)

      cur_cond_snps <- heter_snp_tab %>%
        dplyr::filter(condition == cur_cond) %$%
        variantID %>%
        unique()

      for (cur_celltype in all_celltype) {
        aim_celltype <- cur_celltype
        if (cur_celltype == "ALL") aim_celltype <- all_celltype

        promenh_tab <- heter_snp_tab %>%
          dplyr::filter(condition == cur_cond & celltype %in% aim_celltype) %>%
          dplyr::mutate(
            Promoter = base::grepl("\\<promoter\\>", chrom_annot_type),
            Enhancer = base::grepl("\\<enhancer\\>", chrom_annot_type)
          )

        prom_snps <- promenh_tab %>% dplyr::filter(Promoter) %$% variantID
        enh_snps <- promenh_tab %>% dplyr::filter(Enhancer) %$% variantID

        cur_ct_asoc_snps <- asoc_snp_pcd_tab %>%
          dplyr::filter(p_value_adj < p_val_th & celltype %in% aim_celltype) %>%
          dplyr::select(variantID) %>%
          unlist()

        all_annots <- unique(c("Promoter", "Enhancer", colnames(vep_annot_map)))
        for (annot in all_annots) {
          if (annot == "variantID") next

          annot_enrich_tab <- vep_annot_map %>%
            (function(x) {
              if (annot %in% c("Promoter", "Enhancer")) {
                x %<>%
                  dplyr::mutate(Promoter = F, Enhancer = F) %>%
                  dplyr::mutate(
                    Promoter = (variantID %in% prom_snps),
                    Enhancer = (variantID %in% enh_snps)
                  )
              }
              x
            }) %>%
            dplyr::select(variantID, !!rlang::sym(annot)) %>%
            dplyr::mutate(
              is_asoc = factor((variantID %in% cur_ct_asoc_snps), levels = c(T, F)),
              "{annot}" := factor(!!rlang::sym(annot), levels = c(T, F))
            ) %>%
            dplyr::filter(variantID %in% cur_cond_snps) %>%
            dplyr::select(-variantID) %>%
            base::table() %>%
            (function(x) {
              est <- stats::fisher.test(x)

              if (any(x == 0)) {
                return(NULL)
              }

              se <- sqrt(sum(1 / x))
              log2or <- log2(est$estimate)
              ci_left <- log2or - 1.96 * se
              ci_right <- log2or + 1.96 * se
              p_value <- est$p.value

              data.table(
                annotation = annot,
                condition = cur_cond,
                celltype = cur_celltype,
                p_value = p_value,
                log2or = log2or,
                ci_left = ci_left,
                ci_right = ci_right,
                is_xxx.is_asoc = x[1, 1], # is_asoc and is_xxx
                is_xxx.nt_asoc = x[1, 2], # not_asoc and is_xxx
                nt_xxx.is_asoc = x[2, 1], # is_asoc and not_xxx
                nt_xxx.nt_asoc = x[2, 2] # not_asoc and not_xxx
              )
            }) %>%
            base::rbind(annot_enrich_tab)
        }
      }
    }

    annot_enrich_tab %>%
      fwrite(annot_enrich_file, sep = ",", row.names = F, quote = F)
  }

  # What enrichment of annotation to plot.
  cell2plot <- c("cMono", "ncMono", "CD4T", "CD8T", "NK", "B")
  anno2plot <- c("Promoter", "Enhancer", "TFBS")
  place_holder <- expand.grid(disease_group, cell2plot, anno2plot) %>%
    as.data.frame() %>%
    set_colnames(c("condition", "celltype", "annotation"))

  annot_enrich_for_plot <- annot_enrich_tab %>%
    dplyr::filter(annotation %in% anno2plot & celltype %in% cell2plot) %>%
    dplyr::select(condition, celltype, annotation, log2or, p_value) %>%
    dplyr::full_join(place_holder, by = c("condition", "celltype", "annotation")) %>%
    dplyr::mutate(
      condition = factor(condition, levels = disease_group),
      celltype = factor(celltype, levels = cell2plot),
      annotation = factor(annotation, levels = anno2plot),
      fdr = p.adjust(p_value, method = "fdr"),
      label = case_when(
        is.na(fdr) ~ ".",
        fdr < 0.001 ~ format(fdr, scientific = T, digits = 2),
        fdr < 0.1 ~ format(round(fdr, 3), scientific = F),
        T ~ format(round(fdr, 1), scientific = F)
      ),
      log2or = if_else(is.na(log2or), 0, log2or)
    )

  width <- length(cell2plot) * length(anno2plot) * 0.6
  g_bar <- annot_enrich_for_plot %>%
    ggplot(aes(x = annotation, fill = condition)) +
    geom_bar(aes(y = log2or), stat = "identity", position = position_dodge2(0.5)) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.1) +
    geom_text(aes(label = label, y = 0.1),
      angle = 90, hjust = 0, size = 5,
      position = position_dodge2(0.9)
    ) +
    facet_grid(. ~ celltype) +
    scale_fill_discrete(name = "Condition", labels = condition_map) +
    # labs(y = quote(~ -Log[2] ~ "(Odds ratio)"), x = "Annotation") +
    labs(y = NULL, x = NULL) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.title = element_text(size = 15),
      strip.text = element_text(size = 15)
    )

  save_fig_to <- str_glue("{wkdir}/annotation_enrichment_prom_enh_tfbs.pdf")
  ggsave(save_fig_to, plot = g_bar, width = width, height = 6)
}



###############################################################################
#  Enrichment of ASoC SNPs for epigenomic annotation
###############################################################################
enrich_for_epimark <- FALSE
if (enrich_for_epimark) {
  enrichment_annot <- NULL
  for (cond in disease_group) {
    for (annot in c("promoter", "enhancer")) {
      enrichment_annot <- asoc_snp_tab %>%
        dplyr::filter(condition == cond) %>%
        dplyr::mutate(
          is_asoc = if_else(p_value_adj < p_val_th, "ASoC", "NonASoC"),
          is_annot = if_else(grepl(annot, chrom_annot_type), "Annotation", "NonAnno")
        ) %>%
        dplyr::select(is_asoc, is_annot, variantID) %>%
        dplyr::distinct() %>%
        dplyr::select(is_asoc, is_annot) %>%
        table() %>%
        fisher.test() %>%
        as.data.frame.htest() %>%
        dplyr::mutate(condition = cond, annotation = annot) %>%
        rbind(enrichment_annot)
    }
  }

  g <- enrichment_annot %>%
    dplyr::mutate(label = if_else(
      p.value < 0.001,
      format(p.value, digits = 3),
      format(round(p.value, 2), digits = 2)
    )) %>%
    ggplot(aes(y = odds.ratio, x = condition, fill = annotation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_text(aes(label = label), vjust = -0.1, position = position_dodge2(width = 0.9)) +
    labs(y = "Odds ratio", x = "Condition") +
    theme_classic()
  save_fig_to <- str_glue("{wkdir}/enrich_epimarkers.pdf")
  ggsave(save_fig_to, plot = g, width = 3.5, height = 5)

  enrichment_annot <- NULL
  for (cond in disease_group) {
    for (annot in c("promoter", "enhancer")) {
      for (ct in c("Mono", "cMono", "CD4T", "CD8T", "NK", "B")) {
        enrichment_annot <- asoc_snp_tab %>%
          dplyr::filter(condition == cond, grepl(paste0("\\<", ct, "\\>"), celltype)) %>%
          dplyr::mutate(
            is_asoc = if_else(p_value_adj < p_val_th, "ASoC", "NonASoC"),
            is_annot = if_else(grepl(annot, chrom_annot_type), "Annotation", "NonAnno")
          ) %>%
          dplyr::select(is_asoc, is_annot, variantID) %>%
          dplyr::distinct() %>%
          dplyr::select(is_asoc, is_annot) %>%
          table() %>%
          (function(tbl) {
            if (base::nrow(tbl) == 2) {
              return(tbl %>%
                fisher.test() %>%
                as.data.frame.htest() %>%
                dplyr::mutate(condition = cond, annotation = annot, celltype = ct))
            }

            return(NULL)
          }) %>%
          rbind(enrichment_annot)
      }
    }
  }

  g <- enrichment_annot %>%
    dplyr::filter(is.finite(odds.ratio)) %>%
    dplyr::mutate(label = if_else(
      p.value < 0.001,
      format(p.value, digits = 3, scientific = TRUE),
      format(round(p.value, 3), digits = 2)
    )) %>%
    ggplot(aes(y = odds.ratio, x = condition, fill = annotation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_text(aes(label = label), vjust = -0.1, position = position_dodge2(width = 0.9)) +
    facet_grid(~celltype) +
    labs(y = "Odds ratio", x = "Condition") +
    theme_classic()

  save_fig_to <- str_glue("{wkdir}/enrich_epimarkers_pct.pdf")
  ggsave(save_fig_to, plot = g, width = 9, height = 5)
}



###############################################################################
#  Enrichment for Cytokine QTLs
###############################################################################
enrich_for_cqtl <- FALSE
if (enrich_for_cqtl) {
  cqtl_snp_dir <- str_glue("{pjdir}/inputs/cQTL_500FG")
  cqtl_snp_files <- list.files(path = cqtl_snp_dir, pattern = "*.qtls.gz", full.names = FALSE)

  for (per_file in cqtl_snp_files) {
    cytokine <- str_remove(per_file, ".qtls.gz")
    cqtl_snp_tab <- fread(file.path(cqtl_snp_dir, per_file),
      tmpdir = dirname(cqtl_snp_dir), verbose = F,
      showProgress = F
    )

    asoc_vs_cqtl_tab <- cqtl_snp_tab %>%
      dplyr::select(-c(rsId2, chr, ps)) %>%
      dplyr::inner_join(asoc_snp_tab, by = c("rsId" = "variantID")) %>%
      dplyr::select(-c(chrom_annot_evi)) %>%
      dplyr::rename("cytokine" = "gene", "variantID" = "rsId") %>%
      dplyr::relocate(cytokine, .after = "condition") %>%
      dplyr::group_by(condition, celltype) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::mutate(
        observed = -log10(pvalue),
        expected = -log10(ppoints(n())),
        ci_upper = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = seq(n()), shape2 = rev(seq(n())))),
        ci_lower = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = seq(n()), shape2 = rev(seq(n()))))
      )

    g <- asoc_vs_cqtl_tab %>%
      ggplot(aes(x = expected, y = observed)) +
      geom_point(aes(color = condition), size = 1) +
      geom_line(aes(color = condition), size = 0.05) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
      geom_text_repel(aes(label = variantID), . %>% filter(observed >= -log10(5e-5)),
        min.segment.length = 0, nudge_x = 0.2
      ) +
      facet_wrap(~celltype, nrow = 2) +
      labs(
        x = quote("Expected" ~ -Log[10] ~ "(p-value)"),
        y = quote("Observed" ~ -Log[10] ~ "(p-value)"),
        title = cytokine
      ) +
      scale_color_discrete(name = "Condition", labels = condition_map) +
      theme_classic()

    ggsave(str_glue("{wkdir}/qqplot.asoc_vs_cqtl.{cytokine}.pdf"),
      plot = g, width = 10, height = 6
    )
  }
}



###############################################################################
#  Enrichment for public available GWAS SNPs
###############################################################################
enrich_for_pub_gwas <- FALSE
if (enrich_for_pub_gwas) {
  gwas_ver <- "B2"
  gwas_title <- "hg_b2_2021_covid19"
  gwas_snp_file <- str_glue("{pjdir}/inputs/covid19_gwas/covid19hg.org/COVID19_HGI_{gwas_ver}_ALL_leave_23andme_20210607.txt.gz")
  gwas_snp_tab <- fread(
    gwas_snp_file,
    tmpdir = dirname(as.character(gwas_snp_file))
  ) %>%
    dplyr::rename("CHR" = "#CHR")

  snps_38 <- SNPlocs.Hsapiens.dbSNP144.GRCh38
  covid19_asoc_vs_gwas <- assign_gwas_snp(asoc_snp_tab, gwas_snp_tab, alpha = 5e-5, rsid_col = "rsid", pval_col = "all_inv_var_meta_p") %>%
    dplyr::filter(p_value_adj <= 0.05 & all_inv_var_meta_p <= 5e-6) %>%
    dplyr::select(rsid:condition, p_value_adj, all_inv_var_meta_p, all_inv_var_meta_beta, chrom_annot, by_promoter, by_degene) %>%
    dplyr::rename("p_value_gwas" = all_inv_var_meta_p, "p_value_asoc" = p_value_adj) %T>%
    fwrite(str_glue("{wkdir}/asoc_snps_vs_gwas_snps_{gwas_ver}.csv"), sep = ",", row.names = F)

  # Nr. of SNPs that are both COVID19 SNPs and detected heterozygous SNPs
  covid_gwas_snps <- gwas_snp_tab %>%
    filter(all_inv_var_meta_p < 5e-8) %$% rsid %>%
    unique()
  # Nr. of detected heterozygous SNPs
  heter_snps <- heter_snp_tab %>%
    filter(condition == "Severe" & celltype == "cMono") %$% variantID %>%
    unique()
  # Nr. of ASoC SNPs
  covid_asoc_snps <- asoc_snp_tab %>%
    filter(p_value_adj < 0.1 & condition == "Severe" & celltype == "cMono") %$% variantID %>%
    unique()


  # Nature critical vs population. GRCh37
  # gwas_title <- "nature_pario_2020_covid19"
  # gwas_snp_file <- "../../inputs/covid19_gwas/2020_nature_Pairo-Castineira-etal_Genetic/genomicc.EUR.PLINK2.txt.gz"
  # gwas_snp_tab <- fread(gwas_snp_file, tmpdir=dirname(gwas_snp_file))

  # NEJM severe (835) vs population (1,255). No overlpped SNPs. GRCh38
  # gwas_title <- "nejm_ellinghaus_2020_covid19"
  # gwas_snp_file <- "../../inputs/covid19_gwas/2020_NEJM_Ellinghaus-etal_Genomewide/GCST90000255_GRCh38.tsv.gz"
  # gwas_snp_tab <- fread(gwas_snp_file, tmpdir=dirname(gwas_snp_file))

  # Other non-COVID19 GWAS summary statistics
  snp_enrichment <- list()
  snps_37 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  vcffile <- "~/Documents/projects/wd_plink_data/1kg_phase3/1kg_p3_eur.vcf.gz"
  snpbin_dir <- str_glue("{pjdir}/outputs/scATAC-seq/enrichment/pub_gwas_clump")
  n_snps_per_study <- fread(file.path(snpbin_dir, "n_snps_per_study.csv"))

  for (snpbin_file in list.files(snpbin_dir, pattern = "*.clumped")) {
    snp_clump <- fread(file.path(snpbin_dir, snpbin_file)) %>%
      dplyr::mutate(SP2 = if_else(SP2 == "NONE", SNP, SP2))

    gwas_title <- str_remove(snpbin_file, ".clumped")
    n_bk_snps <- n_snps_per_study %>% dplyr::filter(study == snpbin_file) %$% n_snps

    for (cond in disease_group) {
      cat(snpbin_file, cond, "\n")
      snp_enrichment[[gwas_title]][[cond]] <- asoc_snp_tab %>%
        dplyr::filter(p_value_adj < p_val_th & condition == cond) %>%
        dplyr::select(variantID, contig) %>%
        dplyr::rename("rsid" = variantID, "chrom" = contig) %>%
        dplyr::distinct() %>%
        enrich_snps(snp_clump, vcffile = vcffile, snp_locs = snps_37)
    }
  }

  enrdf <- snp_enrichment %>%
    (function(snp_enrich_list) {
      studies <- names(snp_enrich_list)
      dtfm <- NULL
      for (sd in studies) {
        conditions <- names(snp_enrich_list[[sd]])
        for (cd in conditions) {
          p_value <- snp_enrich_list[[sd]][[cd]][["p_value"]]
          ratio <- snp_enrich_list[[sd]][[cd]][["ratio"]]

          p_value <- ifelse(is.null(p_value), 1, p_value)
          ratio <- ifelse(is.null(ratio), 0, ratio)

          dtfm <- data.frame(
            study = sd, condition = cd, p_value = p_value,
            ratio = ratio
          ) %>%
            rbind(dtfm)
        }
      }
      dtfm
    }) %>%
    dplyr::mutate(
      FDR = p.adjust(p_value, method = "fdr"),
      condition = factor(condition, levels = disease_group)
    ) %T>%
    fwrite(str_glue("{wkdir}/asoc_snps_enrich_in_gwas.csv"), row.names = F)


  g <- enrdf %>%
    ggplot(aes(x = study, y = condition)) +
    geom_point(aes(size = -log10(FDR), color = log2(ratio))) +
    scale_color_gradient2(
      low = "darkblue", mid = "gray", high = "darkred",
      name = bquote(~ Log[2] ~ "(ratio)")
    ) +
    scale_size(name = bquote(~ -Log[10] ~ "(FDR)")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1)) +
    labs(y = "ASoC SNPs by condition", x = "Study")

  save_fig_to <- str_glue("{wkdir}/gwas_snp_enrichment.pdf")
  ggsave(save_fig_to, plot = g, width = 7, height = 6)
}



###############################################################################
#  Enrichment of genes having promoter ASoC SNPs
###############################################################################
enrich_for_gsea <- TRUE
if (enrich_for_gsea) {
  # Enrichment of genes by promoter
  gmt_fmap <- c(
    "c1.all.v7.4.symbols.gmt" = "positional",
    "c3.all.v7.4.symbols.gmt" = "regulatory-target",
    "c4.all.v7.4.symbols.gmt" = "computational",
    "c6.all.v7.4.symbols.gmt" = "oncogenic-sig",
    "c7.all.v7.4.symbols.gmt" = "immunologic-sig",
    "c8.all.v7.4.symbols.gmt" = "cell-type-sig"
  )

  for (cond in disease_group) {
    gene_by_prom <- asoc_snp_tab %>%
      dplyr::filter(condition == cond & by_promoter != "" & p_value_adj < p_val_th) %$%
      by_promoter %>%
      base::paste(collapse = "|") %>%
      str_split(pattern = "\\|") %>%
      base::unlist() %>%
      base::unique() %>%
      as.vector()

    cat(cond, ": ", length(gene_by_prom), "\n", sep = "")

    gene_by_prom <- data.frame(
      SYMBOL = gene_by_prom,
      ENTREZID = mapIds(org.Hs.eg.db, gene_by_prom, keytype = "SYMBOL", column = "ENTREZID")
    )
    gene_by_prom <- gene_by_prom %>%
      dplyr::arrange(-as.integer(ENTREZID)) %>%
      as.data.frame()

    ekegg <- enrichKEGG(gene = gene_by_prom$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
    if (nrow(ekegg) > 0) {
      g <- dotplot(ekegg, showCategory = 20)
      ggsave(str_glue("{wkdir}/{cond}_KEGG.pdf"), plot = g, width = 10, height = 7)
    } else {
      cat("No enriched terms found!\n")
    }

    for (ont in c("CC", "BP", "MF")) {
      ego <- enrichGO(gene = gene_by_prom$ENTREZID, OrgDb = org.Hs.eg.db, ont = ont, readable = T)

      if (nrow(ego)) {
        g <- dotplot(ego, showCategory = 20)
        ggsave(str_glue("{wkdir}/{cond}_GO_{ont}.pdf"), plot = g, width = 10, height = 7)
      } else {
        cat("No enriched terms found!\n")
      }
    }

    for (gmtfile in list.files(str_glue("{pjdir}/inputs/MSigDB"))) {
      hmname <- gmt_fmap[gmtfile]
      gmtobj <- read.gmt(str_glue("../../inputs/MSigDB/{gmtfile}"))

      ehm <- enricher(gene_by_prom$SYMBOL, TERM2GENE = gmtobj)
      if (nrow(ehm)) {
        g <- dotplot(ehm, showCategory = 20)
        ggsave(str_glue("{wkdir}/{cond}_GSEA_{hmname}.pdf"), plot = g, width = 10, height = 7)
      } else {
        cat("No enriched terms found!\n")
      }
    }
  }
}

sessionInfo()
