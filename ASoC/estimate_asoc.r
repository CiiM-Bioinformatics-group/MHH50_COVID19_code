#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 02, 2021
# Updated: Jan 28, 2022

options(stringsAsFactors = F, error = traceback)

library(ggrepel)
library(ggupset)
library(magrittr)
library(tidyverse)
library(data.table)
library(ggVennDiagram)

pjdir <- "~/Documents/projects/wp_covid19_mhh50"
wkdir <- str_glue("{pjdir}/outputs/scATAC-seq/summary")

#
## Allele-specific read counts
#
# Multiple test adjustment
adj_mthd <- "fdr"

# Minimal per allele read counts
min_pa_rc <- 3
min_gr_rc <- 15

# Wheter merge Severe and Mild into a group called 'active'.
merge_active <- TRUE

# Mapping from donor ID to condition
idmap_file <- str_glue("{pjdir}/inputs/idmapping/metainfo_20220128.csv")
idmap_tab <- fread(idmap_file) %>%
  dplyr::filter(remove_for_uniq_sample == "n")

group_by_col <- "New_Severity"
cond_map <- idmap_tab %>%
  dplyr::select(one_of("scATACID", group_by_col)) %>%
  tibble::deframe()

geno_map <- idmap_tab %>%
  dplyr::select("scATACID", genoID) %>%
  tibble::deframe()

cand_donor_id <- names(cond_map)


# Allelic read counts at heterozygous sites.
# asoc_fpath <- c("../../outputs/scATAC-seq/run_v1/optdir/readcounts/asoc_estimation_run_v1_dup.txt")
# asoc_fpath <- c("{wkdir}/outputs/scATAC-seq/readcounts/part_one/optdir/readcounts/readcounts_run_v2_dup.csv")
asoc_fpath <- c(
  file.path(
    pjdir, "outputs", "scATAC-seq", "readcounts", "part_one", "optdir",
    "readcounts", "readcounts_run_v2_nondup.csv"
  )
)

.readcounts_tab <- NULL
for (asocfp in asoc_fpath) {
  .readcounts_tab <- fread(asocfp) %>% rbind(.readcounts_tab)
}

# Keep only individuals we need.
readcounts_tab <- .readcounts_tab %>%
  dplyr::filter(donor_id %in% cand_donor_id) %>%
  dplyr::mutate(
    condition = factor(cond_map[donor_id], levels = c("severe", "mild", "post")),
  ) %>%
  dplyr::select(-c(genotype, lowBaseQDepth, lowMAPQDepth, rawDepth, otherBases, improperPairs))

if (merge_active) {
  disease_group <- c("active", "post")
  readcounts_tab %<>%
    dplyr::mutate(condition = factor(if_else(condition == "post", "post", "active"),
      levels = disease_group
    ))
}


# SNP information per cell type per individual
# Each heterozygous SNPs will be annotated according to Roadmap project.
# Check the scripts/bash/intersect_asocsnps_and_epimarks.sh for details
heter_snps_bed_file <- file.path(wkdir, "heter-snps.bed")
if (!file.exists(heter_snps_bed_file)) {
  readcounts_tab %>%
    dplyr::select(
      c(
        contig, position, variantID, refAllele, altAllele, donor_id, celltype,
        condition, refCount, altCount, totalCount
      )
    ) %>%
    dplyr::arrange(contig, position) %T>%
    data.table::fwrite(str_glue("{wkdir}/heter-snps_metainfo.csv"),
      sep = ",", row.names = F, quote = F
    ) %>%
    dplyr::select(c(contig, position, variantID)) %>%
    dplyr::distinct() %>%
    dplyr::summarise(
      chrom = paste0("chr", contig),
      chromStart = position - 1, chromEnd = position, name = variantID
    ) %>%
    dplyr::rename("#chrom" = "chrom") %>%
    data.table::fwrite(heter_snps_bed_file, sep = "\t", row.names = F, quote = F)
}


# Add condtion information, keep only heterozygous sites.
#  condition n_donors
#  active          21
#  post            18
readcounts_tab %>%
  dplyr::group_by(condition) %>%
  dplyr::select(donor_id) %>%
  dplyr::distinct() %>%
  dplyr::summarize(n_donors = n())

# Total heterozygous SNPs.
# n_het_snps
#     175082
readcounts_tab %>%
  dplyr::select(variantID) %>%
  dplyr::distinct() %>%
  dplyr::summarise(n_het_snps = n())

# Heterozygous SNPs per condition.
# condition n_snps
# active    128893
# post       79143
readcounts_tab %>%
  dplyr::group_by(condition) %>%
  dplyr::select(variantID) %>%
  dplyr::distinct() %>%
  dplyr::summarize(n_snps = n())

# Only keep heterozygous SNPs exists at least in three individuals
heter_snps <- readcounts_tab %>%
  dplyr::select(donor_id, variantID) %>%
  dplyr::distinct() %>%
  dplyr::group_by(variantID) %>%
  dplyr::filter(n() >= 3) %$%
  variantID %>%
  unique()

# Heterozygous SNPs shared by conditions
# shared     n
# <chr>  <int>
# no      3189
# yes    26270
readcounts_tab %>%
  dplyr::select(variantID, condition) %>%
  dplyr::distinct() %>%
  dplyr::filter(variantID %in% heter_snps) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(variantID) %>%
  dplyr::summarise(shared = if_else(n() != 1, "yes", "no")) %>%
  dplyr::group_by(shared) %>%
  dplyr::summarise(n = n())

# n_snps
# 423824
readcounts_tab %>%
  dplyr::select(variantID) %>%
  dplyr::filter(variantID %in% heter_snps) %>%
  dplyr::summarise(n_snps = n())

# condition n_snps
# active     29106
# post       26623
readcounts_tab %>%
  dplyr::group_by(condition) %>%
  dplyr::select(variantID) %>%
  dplyr::filter(variantID %in% heter_snps) %>%
  dplyr::distinct() %>%
  dplyr::summarize(n_snps = n())

# condition celltype    n_snps
# active    B            13946
# active    CD4T         42831
# active    CD8T         26935
# active    cMono        69633
# active    ncMono        6270
# active    NK           30814
# active    pDC            567
# active    Plasmablast   3138
# post      B            12435
# post      CD4T         42484
# post      CD8T         22710
# post      cMono        27967
# post      ncMono        6061
# post      NK           11726
# post      pDC            284
# post      Plasmablast    396
readcounts_tab %>%
  group_by(condition, celltype) %>%
  select(variantID) %>%
  distinct() %>%
  summarize(n_snps = n()) %>%
  print(n = 25)



###############################################################################
## Summary
###############################################################################
# ASoC per sample, read counts were merged across cell types per sample
ps_asoc_summ <- readcounts_tab %>%
  select(donor_id, condition, altCount, totalCount) %>%
  group_by(donor_id, condition) %>%
  summarize(
    altCount = sum(altCount), totalCount = sum(totalCount),
    mle_theta = altCount / totalCount
  ) %>%
  mutate(
    condition = factor(condition, levels = disease_group),
    labels = paste0(
      format(totalCount, big.mark = ","),
      "(", format(mle_theta, digits = 3), ")"
    )
  )

# Sort donor id by totalCount
ps_donor_id_order <- ps_asoc_summ %>%
  arrange(totalCount) %>%
  select(donor_id) %>%
  unlist() %>%
  as.vector()

# Bar plot of read counts per sample
# Will be in the supplementary figure
save_to <- file.path(wkdir, "read-counts_per-sample.pdf")
if (!file.exists(save_to)) {
  g <- ps_asoc_summ %>%
    ggplot() +
    geom_bar(aes(y = donor_id, x = totalCount, fill = condition), stat = "identity") +
    geom_text_repel(aes(y = donor_id, x = -250000, label = labels), direction = "x", nudge_x = 1) +
    theme_classic() +
    theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5)) +
    labs(y = "Donor ID", x = "Total read counts at heterozygous sites") +
    scale_y_discrete(limits = ps_donor_id_order)

  ggsave(save_to, plot = g, width = 7, height = 12)
} else {
  cat("[W]: File exists, skip it\n")
}


###############################################################################
## Pseudo-bulk. Regardless of samples.
###############################################################################
# Estimate ASoC SNPs
pbpc_asoc_snp <- readcounts_tab %>%
  group_by(contig, position, variantID, refAllele, altAllele, condition, celltype) %>%
  summarise(refCount = sum(refCount), altCount = sum(altCount), totalCount = sum(totalCount)) %>%
  ungroup() %>%
  relocate(refCount, altCount, totalCount, .before = "condition") %>%
  (function(dtfm) {
    optpath <- str_glue("{wkdir}/pseudo-bulk_per-celltype_heter-snps.csv")
    if (!file.exists(optpath)) {
      fwrite(dtfm, optpath, sep = ",", row.names = F, quote = F)
    } else {
      cat("[W]: File exists, skip it\n")
    }

    dtfm
  }) %>%
  filter(refCount >= min_pa_rc & altCount >= min_pa_rc & totalCount >= min_gr_rc) %>%
  mutate(
    p_value = mapply(function(x, n) binom.test(x, n)$p.value, altCount, totalCount),
    p_value_adj = p.adjust(p_value, method = adj_mthd)
  ) %>%
  arrange(condition, celltype, variantID)

readcounts_source_file <- file.path(
  wkdir,
  "pseudo-bulk_per-celltype_heter-snps_reads-source.csv"
)
if (!file.exists(readcounts_source_file)) {
  readcounts_tab %>%
    dplyr::filter(variantID %in% unique(pbpc_asoc_snp %$% variantID)) %>%
    dplyr::select(-c(refCount, altCount, totalCount, geno_id)) %>%
    dplyr::arrange(contig, position, condition, celltype) %>%
    dplyr::distinct() %>%
    data.table::fwrite(readcounts_source_file, sep = ",", row.names = F, quote = F)
} else {
  cat("[W]: File exists, skip it.\n")
}

# Number of identified ASoC SNPs per cell type per condition.
# condition celltype n_snps
# active    B             2
# active    CD4T         92
# active    CD8T         31
# active    cMono       157
# active    NK           47
# post      CD4T         62
# post      CD8T         19
# post      cMono        11
# post      NK            1
pbpc_asoc_snp %>%
  filter(p_value_adj <= 0.05) %>%
  group_by(condition, celltype) %>%
  summarise(n_snps = n())

# ASoC SNPs (p_value_adj<0.05 by FDR)
pbpc_cand_asoc_snp <- pbpc_asoc_snp %>%
  filter(p_value_adj <= 0.05) %>%
  arrange(p_value_adj, contig)

#
pbpc_asoc_file <- str_glue("{wkdir}/pseudo-bulk_per-celltype_asoc-snps.csv")
if (file.exists(pbpc_asoc_file)) {
  pbpc_asoc_snp <- fread(pbpc_asoc_file, verbose = F, showProgress = F)
} else {
  pbpc_asoc_snp %>%
    dplyr::arrange(contig, position, p_value) %>%
    data.table::fwrite(pbpc_asoc_file, sep = ",", row.names = F, quote = F)
}


###############################################################################
# Venn plot by conditions
###############################################################################
if (merge_active) {
  data_gross <- list(
    Hospitalized = pbpc_cand_asoc_snp %>% filter(condition == "active") %>% select(variantID) %>% distinct() %>% unlist() %>% as.vector(),
    Convalescent = pbpc_cand_asoc_snp %>% filter(condition == "post") %>% select(variantID) %>% distinct() %>% unlist() %>% as.vector()
  )
} else {
  data_gross <- list(
    Seve. = pbpc_cand_asoc_snp %>% filter(condition == "Severe") %>% select(variantID) %>% distinct() %>% unlist() %>% as.vector(),
    Mild = pbpc_cand_asoc_snp %>% filter(condition == "Mild") %>% select(variantID) %>% distinct() %>% unlist() %>% as.vector(),
    Convalescent = pbpc_cand_asoc_snp %>% filter(condition == "Post") %>% select(variantID) %>% distinct() %>% unlist() %>% as.vector()
  )
}

g_venn <- ggVennDiagram(data_gross,
  edge_size = 3, label_size = 12, set_size = 8, label_percent_digit = 2
) +
  scale_fill_gradient(low = "white", high = "white") +
  guides(fill = "none")

save_fig_to <- file.path(wkdir, "venn_diagram_asoc_snp.pdf")
ggsave(save_fig_to, plot = g_venn, width = 8, height = 8)

# Venn plot by cell type
for (ct in c("Mono", "B", "CD4T", "CD8T", "NK")) {
  # Venn diagram
  if (ct == "Mono") {
    ct_comb <- c("cMono", "ncMono")
  } else {
    ct_comb <- ct
  }

  if (merge_active) {
    data_pct <- list(
      Hospitalized = pbpc_cand_asoc_snp %>% filter(celltype %in% ct_comb & condition == "active") %>% select(variantID) %>% distinct() %>% unlist() %>% as.vector(),
      Convalescent = pbpc_cand_asoc_snp %>% filter(celltype %in% ct_comb & condition == "post") %>% select(variantID) %>% distinct() %>% unlist() %>% as.vector()
    )
  } else {
    data_pct <- list(
      Seve. = pbpc_cand_asoc_snp %>% filter(celltype %in% ct_comb & condition == "Severe") %>% select(variantID) %>% unlist(),
      Mild = pbpc_cand_asoc_snp %>% filter(celltype %in% ct_comb & condition == "Mild") %>% select(variantID) %>% unlist(),
      convalescent = pbpc_cand_asoc_snp %>% filter(celltype %in% ct_comb & condition == "Post") %>% select(variantID) %>% unlist()
    )
  }

  g_venn <- ggVennDiagram(data_pct, edge_size = 3, label_size = 12, set_size = 8, label_percent_digit = 2) +
    scale_fill_gradient(low = "white", high = "grey", name = "Counts")

  save_to <- str_glue("{wkdir}/venn_diagram_asoc_snp_{ct}.pdf")
  ggsave(save_to, plot = g_venn, width = 8, height = 8)

  # Upset plots by ggupset
  # g_upset <- pbpc_cand_asoc_snp %>%
  #   filter(celltype %in% ct_comb) %>%
  #   select(variantID, condition) %>%
  #   distinct() %>%
  #   group_by(variantID) %>%
  #   summarize(SNPs = list(condition)) %>%
  #   ggplot(aes(x = SNPs)) +
  #   geom_bar() +
  #   geom_text(aes(label = after_stat(count)), stat = "count", vjust = -1) +
  #   scale_x_upset() +
  #   theme_classic()
  # save_to <- str_glue("{wkdir}/upset_plot_asoc_snp_{ct}.pdf")
  # ggsave(save_to, plot = g_upset, height = 8, width = 8)
}



for (cond in disease_group) {
  # Venn diagram
  data_pcond <- list(
    Mono = pbpc_cand_asoc_snp %>% filter(condition == cond & celltype %in% c("cMono", "ncMono")) %>% select(variantID) %>% unlist(),
    B = pbpc_cand_asoc_snp %>% filter(condition == cond & celltype == "B") %>% select(variantID) %>% unlist(),
    CD4T = pbpc_cand_asoc_snp %>% filter(condition == cond & celltype == "CD4T") %>% select(variantID) %>% unlist(),
    CD8T = pbpc_cand_asoc_snp %>% filter(condition == cond & celltype == "CD8T") %>% select(variantID) %>% unlist(),
    NK = pbpc_cand_asoc_snp %>% filter(condition == cond & celltype == "NK") %>% select(variantID) %>% unlist()
  )

  g_venn <- ggVennDiagram(data_pcond,
    edge_size = 3, label_size = 8, set_size = 10, label_percent_digit = 2
  ) +
    scale_fill_gradient(low = "white", high = "grey", name = "Counts")

  save_to <- str_glue("{wkdir}/venn_diagram_asoc_snp_{cond}.pdf")
  ggsave(save_to, plot = g_venn, width = 10, height = 10)

  # # Upset plots by ggupset
  # g_upset <- pbpc_cand_asoc_snp %>%
  #   mutate(celltype = case_when(
  #     celltype %in% c("cMono", "ncMono") ~ "Mono",
  #     celltype %in% c("Plasmablast") ~ "B",
  #     T ~ celltype
  #   )) %>%
  #   filter(condition == cond) %>%
  #   select(variantID, celltype) %>%
  #   distinct() %>%
  #   group_by(variantID) %>%
  #   summarize(SNPs = list(celltype)) %>%
  #   ggplot(aes(x = SNPs)) +
  #   geom_bar() +
  #   geom_text(aes(label = after_stat(count)), stat = "count", vjust = -1) +
  #   scale_x_upset() +
  #   theme_classic()

  # save_to <- str_glue("{wkdir}/upset_plot_asoc_snp_{cond}.pdf")
  # ggsave(save_to, plot = g_upset, height = 8, width = 8)
}


# Volcano plot to show mapping bias
pval_vs_ratio_tab <- pbpc_asoc_snp %>%
  select(refCount, totalCount, p_value, p_value_adj) %>%
  mutate(
    ratio = refCount / totalCount,
    log10_p_value = -log10(p_value),
    fdr_le_0.05 = p_value_adj < 0.05
  )

g_dot <- pval_vs_ratio_tab %>%
  ggplot() +
  theme_classic() +
  geom_point(aes(x = ratio, y = log10_p_value, color = fdr_le_0.05), size = 0.2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "FDR<0.05") +
  labs(x = "Reference allele ratio", y = bquote(~ -Log[10] ~ "(p-value)"))
ggsave(str_glue("{wkdir}/check_mapping_bias.pdf"),
  plot = g_dot, width = 6, height = 3
)



# TODO: The following code should be removed when submitting the manuscript.
###############################################################################
# The rest analyzes were exploratory and not reported in main text.           #
###############################################################################

if (F) {
  # Check whether the condition specific ASoC SNPs only exists in the conditions (singleton)
  # (i.e. only exist in cohort of a specific condition)

  # Nr. of conditional specific SNPs
  cond_spec_asoc_snp <- pbpc_cand_asoc_snp %>%
    select(variantID, condition) %>%
    distinct() %>%
    group_by(variantID) %>%
    summarize(n_condition = n(), condition = paste0(condition, collapse = "|")) %>%
    filter(n_condition == 1) %>%
    group_by(condition) %>%
    summarize(n_specific = n())

  specific_snps_per_celltype <- pbpc_cand_asoc_snp %>%
    select(variantID, condition) %>%
    distinct() %>%
    group_by(variantID) %>%
    summarize(n_condition = n()) %>%
    filter(n_condition == 1) %>%
    select(variantID) %>%
    unlist()

  # Number of singleton condition specific ASoC SNPs.
  singleton_cond_spec_asoc_snp <- pbpc_asoc_snp %>%
    select(variantID, condition) %>%
    distinct() %>%
    group_by(variantID) %>%
    summarize(n_condition = n(), condition = paste0(condition, collapse = "|")) %>%
    filter(n_condition == 1 & variantID %in% specific_snps_per_celltype) %>%
    group_by(condition) %>%
    summarize(n_singleton = n())

  # Plot a bar plot
  g <- cond_spec_asoc_snp %>%
    inner_join(singleton_cond_spec_asoc_snp) %>%
    mutate(n_non_singleton = n_specific - n_singleton) %>%
    pivot_longer(cols = c(n_non_singleton, n_singleton), names_to = "groups", values_to = "n_snps") %>%
    mutate(
      groups = factor(groups, levels = c("n_singleton", "n_non_singleton")),
      condition = factor(condition, levels = disease_group)
    ) %>%
    ggplot(aes(x = condition, y = n_snps, fill = groups)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    geom_text(aes(label = n_snps), position = position_fill(vjust = 0.5)) +
    theme_classic() +
    labs(x = "Condition", y = "Proportion") +
    scale_fill_discrete(name = "Condition specific SNPs", labels = c("n_non_singleton" = "Non-singleton", "n_singleton" = "Singleton"))
  ggsave(str_glue("{wkdir}/check_cell_type_specific_snps.pdf"), width = 5, height = 7)
}


# Cell type specific and shared in >= 3 cell type ASoC
# As the number of ASoC SNPs per cell grows with number of the corresponding
# cell type, it's not fair to identify cell type specific ASoC SNPs for
# scATAC-seq data, but shared ASoC SNPs are still meaningful.
# pbpc_asoc_snp <- pbpc_asoc_snp %>%
#   group_by(variantID, condition) %>%
#   mutate(specificity=case_when(sum(p_value_adj<0.05)==1 & sum(p_value>0.05)>=0 ~ "Specific",
#                                sum(p_value_adj<0.05)>=3                        ~ "Shared",
#                                T                                               ~ "Other"),
#          n_celltype=n()) %>%
#   ungroup()

# Write the estimation resutls to disk
if (F) {
  # NOTE: The cell counts should be updated, or remove this code block.
  # Check the correlation between cell counts and nr. of identified ASoC SNPs
  cell_counts <- data.frame(
    condition = c("Mild", "Post", "Severe"),
    CD4T = c(3729, 4015, 2485),
    CD8T = c(2053, 2508, 2223),
    cMono = c(3185, 2514, 3935),
    ncMono = c(322, 472, 323),
    NK = c(1490, 1473, 2139),
    B = c(1094, 1330, 992)
  )
  cell_counts <- cell_counts %>%
    pivot_longer(
      cols = c(CD4T, CD8T, cMono, ncMono, NK, B), names_to = "celltype",
      values_to = "cell_counts"
    )

  ccas_tab <- pbpc_asoc_snp %>%
    filter(p_value_adj < 0.05) %>%
    group_by(condition, celltype) %>%
    select(variantID) %>%
    distinct() %>%
    summarize(n_snps = n()) %>%
    ungroup() %>%
    full_join(cell_counts, by = c("condition", "celltype")) %>%
    mutate(condition = factor(condition, levels = c("Severe", "Mild", "Post"))) %>%
    filter(!is.na(cell_counts))

  g <- ccas_tab %>%
    ggplot(aes(x = cell_counts, y = n_snps, color = condition)) +
    geom_point() +
    geom_smooth(method = "lm", se = F, size = 0.1) +
    geom_text_repel(aes(label = celltype)) +
    labs(x = "Cell counts", y = "Nr. of identified ASoC SNPs (FDR<0.05)") +
    scale_color_discrete(name = "Condition") +
    theme_classic()
  ggsave(str_glue("{wkdir}/check_cell_counts_vs_asoc_snps.pdf"), plot = g)

  hsas_tab <- readcounts_tab %>%
    group_by(condition, celltype) %>%
    select(variantID) %>%
    distinct() %>%
    summarise(n_snps = n()) %>%
    ungroup() %>%
    full_join(cell_counts, by = c("condition", "celltype")) %>%
    mutate(condition = factor(condition, levels = c("Severe", "Mild", "Post"))) %>%
    filter(!is.na(cell_counts))

  g <- hsas_tab %>%
    ggplot(aes(x = cell_counts, y = n_snps, color = condition)) +
    geom_point() +
    geom_smooth(method = "lm", se = F, size = 0.1) +
    geom_text_repel(aes(label = celltype)) +
    labs(x = "Cell counts", y = "Nr. of heterozygous SNPs") +
    scale_color_discrete(name = "Condition") +
    theme_classic()
  ggsave(str_glue("{wkdir}/check_cell_counts_vs_heter_snps.pdf"), plot = g)

  # Number of ASoC SNPs. Cell type specific and shared in >= 3 cell type
  # pbpc_asoc_snp %>%
  #   group_by(specificity, condition) %>%
  #   summarize(n_asoc_snps=n()) %>%
  #   filter(specificity %in% c("Shared", "Specific"))

  # Assign celltype spcificity to candidate SNPs
  # celltype_spcificity <- pbpc_asoc_snp %>%
  #   filter(variantID %in% pbpc_cand_asoc_snp$variantID) %>%
  #   select(variantID, specificity) %>%
  #   distinct() %>%
  #   deframe()

  # pbpc_cand_asoc_snp <- pbpc_cand_asoc_snp %>%
  #   mutate(specificity=celltype_spcificity[variantID])


  #
  ## Pseudo-bulk. Regardless of either sample or cell type.
  #
  pb_asoc_snp <- readcounts_tab %>%
    select(-c(donor_id, geno_id, condition, celltype)) %>%
    filter(refCount >= min_pa_rc & altCount >= min_pa_rc & totalCount >= min_gr_rc) %>%
    group_by(variantID) %>%
    mutate(refCount = sum(refCount), altCount = sum(altCount), totalCount = sum(totalCount)) %>%
    ungroup() %>%
    distinct() %>%
    mutate(
      p_value = mapply(function(x, n) {
        binom.test(x, n)$p.value
      }, altCount, totalCount),
      p_value_adj = p.adjust(p_value, method = adj_mthd)
    )

  pb_cand_asoc_snp <- pb_asoc_snp %>%
    filter(p_value_adj < 0.05) %>%
    arrange(p_value_adj, contig)
  pb_cand_asoc_snp %>%
    write.table(str_glue("{wkdir}/pseudo-bulk.csv"), sep = ",", row.names = F, quote = F)
}

if (F) {
  #
  ## Pseudo-bulk. Regardless of cell types
  #
  # Estimate ASoC per sample
  pbps_asoc_snp <- readcounts_tab %>%
    select(-c(celltype, geno_id)) %>%
    filter(refCount >= min_pa_rc & altCount >= min_pa_rc & totalCount >= min_gr_rc) %>%
    group_by(donor_id, variantID) %>%
    mutate(refCount = sum(refCount), altCount = sum(altCount), totalCount = sum(totalCount)) %>%
    ungroup() %>%
    distinct() %>%
    mutate(
      p_value = mapply(function(x, n, i) {
        binom.test(x, n)$p.value
      }, altCount, totalCount),
      p_value_adj = p.adjust(p_value, method = adj_mthd)
    )

  # Keep variants with at least min_pa_rc REF and min_pa_rc ALT read counts
  pbps_cand_asoc_snp <- pbps_asoc_snp %>%
    filter(p_value_adj < 0.05) %>%
    arrange(p_value_adj)

  pbps_cand_asoc_snp %>%
    write.table(str_glue("{wkdir}/pseudo-bulk_per-sample.csv"), sep = ",", row.names = F, quote = F)

  # Candidate variants
  pbps_cand_asoc_snp_id <- pbps_asoc_snp %>%
    group_by(variantID) %>%
    mutate(n_snps = n()) %>%
    filter(n_snps >= 3 & p_value_adj < 0.05 & sum(p_value_adj < 0.05) > 1) %>%
    select(variantID) %>%
    distinct() %>%
    unlist() %>%
    as.vector()

  # Pseudo-bulk per sample ASoC plot
  g <- pbps_asoc_snp %>%
    filter(variantID %in% pbps_cand_asoc_snp_id) %>%
    ggplot(aes(x = refCount, y = altCount)) +
    facet_wrap(~variantID, nrow = 5) +
    geom_rect(xmin = min_pa_rc, xmax = Inf, ymin = min_pa_rc, ymax = Inf, fill = "green", alpha = 0.01) +
    geom_point(aes(color = p_value_adj < 0.05, shape = condition), size = 3, alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_text_repel(aes(label = donor_id), min.segment.length = 0) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "P-value-adj<0.05") +
    scale_shape_discrete(name = "Condtion") +
    theme_classic() +
    theme(text = element_text(size = 15)) +
    labs(x = "Reference read counts", y = "Alternative read counts")

  ggsave(str_glue("{wkdir}/pseudo-bulk_per-sample_asoc.pdf"), plot = g, width = 20, height = 20)


  #
  ## ASoC specific
  #
  for (snp_id in pbps_cand_asoc_snp_id) {
    next
    pbps_cand_asoc <- readcounts_tab %>%
      filter(donor_id %in% c("2d_3") & variantID == snp_id) %>%
      arrange(variantID) %>%
      mutate(
        p_value_adj = p.adjust(p_value, method = adj_mthd),
        celltype = factor(celltype, levels = c("NK", "CD4T", "cMono", "B", "CD8T", "ncMono")),
        success_prob = altCount / totalCount,
        failed_prob = refCount / totalCount,
        errbar_up = mapply(function(x, n) binom.test(x, n)$conf.int[1], x = altCount, n = totalCount),
        errbar_dw = mapply(function(x, n) binom.test(x, n)$conf.int[2], x = altCount, n = totalCount)
      )

    # Will be used in the figure.
    g <- pbps_cand_asoc %>%
      ggplot(aes(x = celltype, y = success_prob)) +
      geom_bar(aes(y = 1), stat = "identity", width = 0.6, fill = "lightgrey") +
      geom_bar(stat = "identity", width = 0.6, fill = "#00AFBB") +
      geom_errorbar(aes(ymin = errbar_dw, ymax = errbar_up), width = 0.1) +
      geom_hline(yintercept = 0.5, linetype = "dotted") +
      geom_text(aes(y = 0.02, label = paste0(altCount, "/", totalCount))) +
      ylim(c(0, 1)) +
      labs(x = "Cell type", y = "Proportion of alternative read counts") +
      theme_classic() +
      theme(text = element_text(size = 15))

    ggsave(str_glue("{wkdir}/pseudo-bulk_per-sample_asoc_{snp_id}.pdf"), plot = g, width = 5, height = 7)
  }



  #
  ## Per cell type per sample
  #
  # Only heterozygous sites with at least 5 refCount and 5 altCount.
  pcps_asoc_snp <- readcounts_tab %>%
    filter(refCount >= min_pa_rc & altCount >= min_pa_rc & totalCount >= min_gr_rc) %>%
    group_by(donor_id, celltype, condition) %>%
    mutate(
      p_value = mapply(function(x, n) {
        binom.test(x, n)$p.value
      }, altCount, totalCount),
      p_value_adj = p.adjust(p_value, method = adj_mthd)
    ) %>%
    ungroup()


  # SNP of ASoC per cell type per sample.
  pcps_cand_asoc_snp <- pcps_asoc_snp %>%
    filter(p_value_adj < 0.05) %>%
    arrange(p_value_adj, contig)

  pcps_cand_asoc_snp %>%
    write.table(str_glue("{wkdir}/per-celltype_per-sample_asoc-snp.csv"), sep = ",", row.names = F, quote = F)

  # The order of donor ID sorted by disease condition.
  pcps_donor_id_order <- pcps_asoc_snp %>%
    select(condition, donor_id) %>%
    distinct() %>%
    arrange(condition) %>%
    select(donor_id) %>%
    unlist() %>%
    as.vector()

  # The summary of ASoC per cell type per sample.
  pcps_asoc_snp_summary <- pcps_asoc_snp %>%
    group_by(donor_id, condition, celltype) %>%
    summarize(
      "5e-2" = sum(p_value < 0.05),
      "1e-2" = sum(p_value < 0.01),
      "5e-3" = sum(p_value < 0.005),
      "1e-3" = sum(p_value < 0.001),
      "5e-4" = sum(p_value < 0.0005),
      "5e-5" = sum(p_value < 0.00005),
      "0.05(FDR)" = sum(p_value_adj < 0.05)
    ) %>%
    ungroup() %>%
    reshape2::melt(id.vars = c("donor_id", "condition", "celltype"), variable.name = "alpha", value.name = "count") %>%
    as.data.frame()

  # A bar plot to show the summary
  # filter(! celltype %in% c("ncMono", "Plasmablast")) %>%
  g <- pcps_asoc_snp_summary %>%
    ggplot() +
    theme_classic() +
    facet_grid(alpha ~ celltype, scale = "free_x", switch = "y") +
    geom_bar(aes(x = donor_id, y = count, fill = condition), stat = "identity") +
    scale_fill_manual(name = "Disease condition", values = c("Mild" = "blue", "Severe" = "red", "Post" = "green")) +
    labs(x = "Donor ID", y = "Number of ASoC SNPs") +
    scale_x_discrete(limits = pcps_donor_id_order) +
    theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))

  ggsave(str_glue("{wkdir}/per-celltype_per-sample_asoc-snp_summary.pdf"), plot = g, width = 20, height = 8)

  # SNP IDs of ASoC per cell type per sample.
  pcps_cand_asoc_snp_id <- pcps_cand_asoc_snp %>%
    ungroup() %>%
    group_by(variantID) %>%
    select(variantID) %>%
    mutate(n_snps = n()) %>%
    distinct() %>%
    filter(n_snps > 1) %>%
    select(variantID) %>%
    unlist()

  # Dot plot to show each chosen SNP.
  for (snpid in pcps_cand_asoc_snp_id) {
    next
    pcps_tar_asoc_snp <- readcounts_tab %>%
      mutate() %>%
      ungroup() %>%
      filter(variantID == snpid)

    lim_max <- pcps_tar_asoc_snp %>%
      select(refCount, altCount) %>%
      max()
    g <- pcps_tar_asoc_snp %>%
      ggplot(aes(x = refCount, y = altCount)) +
      theme_classic() +
      facet_wrap(~celltype, nrow = 2) +
      geom_point(aes(color = p_value_adj < 0.05, shape = condition), size = 2) +
      geom_text_repel(aes(label = donor_id), min.segment.length = 0) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Significance") +
      scale_shape_discrete(name = "Condtion") +
      theme(text = element_text(size = 20)) +
      lims(x = c(1, lim_max), y = c(1, lim_max)) +
      labs(x = "Reference read counts", y = "Alternative read counts")

    ggsave(str_glue("{wkdir}/per-celltype_per-sample_asoc_{snpid}.pdf"), plot = g, width = 8, height = 7)
  }
}


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
# [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#  [1] ggVennDiagram_1.2.0 data.table_1.14.2   forcats_0.5.1
#  [4] stringr_1.4.0       dplyr_1.0.7         purrr_0.3.4
#  [7] readr_2.1.1         tidyr_1.1.4         tibble_3.1.6
# [10] tidyverse_1.3.1     magrittr_2.0.1      ggupset_0.3.0
# [13] ggrepel_0.9.1       ggplot2_3.3.5       devtools_2.4.3
# [16] usethis_2.1.5
#
# loaded via a namespace (and not attached):
#  [1] httr_1.4.2         pkgload_1.2.4      jsonlite_1.7.3     modelr_0.1.8
#  [5] brio_1.1.3         assertthat_0.2.1   cellranger_1.1.0   remotes_2.4.2
#  [9] sessioninfo_1.2.2  pillar_1.6.4       backports_1.4.1    glue_1.6.1
# [13] digest_0.6.29      rvest_1.0.2        RVenn_1.1.0        colorspace_2.0-2
# [17] pkgconfig_2.0.3    broom_0.7.11       haven_2.4.3        scales_1.1.1
# [21] processx_3.5.2     tzdb_0.2.0         proxy_0.4-26       generics_0.1.1
# [25] farver_2.1.0       ellipsis_0.3.2     cachem_1.0.6       withr_2.4.3
# [29] cli_3.1.1          crayon_1.4.2       readxl_1.3.1       memoise_2.0.1
# [33] ps_1.6.0           fs_1.5.2           fansi_1.0.2        class_7.3-20
# [37] xml2_1.3.3         pkgbuild_1.3.1     tools_4.1.2        prettyunits_1.1.1
# [41] hms_1.1.1          lifecycle_1.0.1    munsell_0.5.0      reprex_2.0.1
# [45] callr_3.7.0        e1071_1.7-9        compiler_4.1.2     rlang_0.4.12
# [49] units_0.7-2        classInt_0.4-3     grid_4.1.2         rstudioapi_0.13
# [53] labeling_0.4.2     testthat_3.1.2     gtable_0.3.0       DBI_1.1.2
# [57] R6_2.5.1           lubridate_1.8.0    fastmap_1.1.0      utf8_1.2.2
# [61] rprojroot_2.0.2    KernSmooth_2.23-20 desc_1.4.0         stringi_1.7.6
# [65] Rcpp_1.0.8         vctrs_0.3.8        sf_1.0-5           dbplyr_2.1.1
# [69] tidyselect_1.1.1
