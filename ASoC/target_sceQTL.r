#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
library(Seurat)
library(ggpubr)
library(ggrepel)
library(magrittr)
library(tidyverse)
library(data.table)
library(preprocessCore)

pjdir <- "~/Documents/projects/wp_covid19_mhh50"

tar_donor <- c(
  "10", "11", "12", "14", "15", "17", "18", "19", "2", "20", "23", "24", "25",
  "26", "27", "29", "3", "31", "34", "35", "36", "38", "39", "4", "40", "41",
  "42", "44", "45", "46", "5", "8"
)



# Genotype to samples
id_mapping_file <- str_glue("{pjdir}/inputs/scRNA-seq/id_mapping_rnaseq.txt")
if (!exists("idmap_tab")) {
  idmap_tab <- fread(id_mapping_file)
}

# Other meta information
metainfo_file <- str_glue("{pjdir}/inputs/scRNA-seq/metainfo_scRNA-seq.csv")
if (!exists("metainfo_scrna")) {
  metainfo_scrna <- fread(metainfo_file) %>%
    dplyr::filter(remove_for_uniq_sample == "n" | sampleID == "27a")
}

# Genotypes
tar_snp <- "rs6800484"
genotype_file <- str_glue("{pjdir}/outputs/genotypes/{tar_snp}.vcf.gz")
gntp_tab <- fread(genotype_file, showProgress = FALSE, verbose = FALSE, skip = "CHROM") %>%
  dplyr::select(-one_of(c("#CHROM", "POS", "QUAL", "FILTER", "INFO", "FORMAT"))) %>%
  dplyr::rename_with(contains("-"), .fn = function(e) paste0("XXXX_", e)) %>%
  dplyr::filter(ID %in% tar_snp) %>%
  dplyr::mutate(across(starts_with("XXXX_"), ~ str_extract(.x, "[01]\\|[01]"))) %>%
  tidyr::pivot_longer(cols = c(starts_with("XXXX_"))) %>%
  tidyr::pivot_wider(names_from = "ID") %>%
  dplyr::mutate(
    name = str_remove(name, "^XXXX_"),
    SNP = case_when(
      rs6800484 %in% c("0|0") ~ paste0(REF, REF),
      rs6800484 %in% c("0|1", "1|0") ~ paste0(REF, ALT),
      rs6800484 %in% c("1|1") ~ paste0(ALT, ALT)
    ),
    dosage = case_when(
      rs6800484 %in% c("0|1", "1|0") ~ 1,
      rs6800484 %in% c("0|0") ~ 0,
      rs6800484 %in% c("1|1") ~ 2
    ),
    SNP_merge = case_when(
      rs6800484 %in% c("0|0", "0|1", "1|0") ~ paste0(REF, REF, "/", REF, ALT),
      rs6800484 %in% c("1|1") ~ paste0(ALT, ALT)
    )
  ) %>%
  dplyr::select(-c(REF, ALT)) %>%
  dplyr::arrange(name)


# scRNA-seq results
covid_obj_file <- str_glue("{pjdir}/inputs/scRNA-seq/pbmc.filter.annotev5.rds")
if (!exists("covid_obj")) {
  covid_obj <- readRDS(covid_obj_file)
}

meta_info <- covid_obj@meta.data %>%
  base::as.data.frame() %>%
  dplyr::select(sampleID, Age, gender, Severity) %>%
  base::unique() %>%
  dplyr::inner_join(metainfo_scrna, by = c("sampleID")) %>%
  dplyr::inner_join(idmap_tab, by = c("sampleID")) %>%
  dplyr::inner_join(gntp_tab, by = c("genoID" = "name"))

new_meta_info <- covid_obj@meta.data %>%
  dplyr::left_join(dplyr::select(meta_info, sampleID, SNP, SNP_merge),
    by = "sampleID"
  ) %>%
  (function(d) {
    rownames(d) <- d$cellbarcodes
    d
  })

covid_obj <- Seurat::AddMetaData(covid_obj, metadata = new_meta_info)


#
## Check expression per sample
#
check_avg_diff <- FALSE
if (check_avg_diff) {
  mono_obj <- covid_obj[, covid_obj@meta.data$celltypeL0 == "cMono"]
  tar_gene <- c("FYCO1", "XCR1", "CCRL2", "CXCR6", "CCR1", "CCR2", "CCR3", "CCR5", "CCR9")

  # Averaged (mean) expression
  mean_exp_tab <- AverageExpression(mono_obj, "RNA", tar_gene, group.by = "SampleID") %$%
    RNA %>%
    as.data.frame() %>%
    dplyr::mutate(gene_symbol = rownames(.)) %>%
    tidyr::pivot_longer(-c(gene_symbol), values_to = "avg_exp", names_to = "sampleID") %>%
    tidyr::pivot_wider(sampleID, names_from = "gene_symbol", values_from = "avg_exp") %>%
    dplyr::rename_with(-c(sampleID), .fn = function(x) paste0("avg.", x))

  rm("mono_obj")
  gc()

  # agg_tab informations
  agg_tab <- meta_info %>%
    dplyr::inner_join(mean_exp_tab, by = "sampleID") %>%
    dplyr::select(-c(imputeID, loadID))

  # Regression
  fml <- formula("value ~ Age + gender + Genotype")
  reg_tab <- agg_tab %>%
    dplyr::filter(sampleID %in% tar_donor, Severity %in% c("mild", "severe")) %>%
    dplyr::select(c(
      sampleID, Age, gender, Severity, WHO_score, SNP, dosage, SNP_merge, starts_with("avg.")
    )) %>%
    tidyr::pivot_longer(c(starts_with("avg."))) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(qtl = {
      snp_reg <- cur_data() %>%
        dplyr::rename("Genotype" = "SNP") %>%
        lm(fml, data = .) %>%
        anova() %>%
        as.data.frame() %>%
        dplyr::mutate(Parameters = rownames(.)) %>%
        dplyr::filter(Parameters != "Residuals") %>%
        dplyr::select(-Df) %>%
        dplyr::rename(
          "P_value_F" = "Pr(>F)", "F_value" = "F value",
          "Sum_sq" = "Sum Sq", "Mean_sq" = "Mean Sq"
        )

      snp_merge_reg <- cur_data() %>%
        dplyr::rename("Genotype" = "SNP_merge") %>%
        lm(fml, data = .) %>%
        anova() %>%
        as.data.frame() %>%
        dplyr::mutate(Parameters = rownames(.)) %>%
        dplyr::filter(Parameters != "Residuals") %>%
        dplyr::select(-Df) %>%
        dplyr::rename(
          "P_value_F" = "Pr(>F)", "F_value" = "F value",
          "Sum_sq" = "Sum Sq", "Mean_sq" = "Mean Sq"
        )

      dosage_reg <- cur_data() %>%
        dplyr::rename("Genotype" = "dosage") %>%
        lm(fml, data = .) %>%
        anova() %>%
        as.data.frame() %>%
        dplyr::mutate(Parameters = rownames(.)) %>%
        dplyr::filter(Parameters != "Residuals") %>%
        dplyr::select(-Df) %>%
        dplyr::rename(
          "P_value_F" = "Pr(>F)", "F_value" = "F value",
          "Sum_sq" = "Sum Sq", "Mean_sq" = "Mean Sq"
        )

      full_join(snp_reg, snp_merge_reg, by = "Parameters", suffix = c(".snp", ".snp_merge")) %>%
        full_join(dosage_reg, by = "Parameters", suffix = c("", ".dosage")) %>%
        dplyr::relocate("Parameters")
    }) %>%
    as.data.table() %>%
    dplyr::rename_with(starts_with("qtl."), .fn = function(e) str_remove_all(e, "qtl."))

  for (per_gene in tar_gene) {
    my_comps <- agg_tab %$% SNP_merge %>%
      unique() %>%
      combn(2) %>%
      as.data.frame() %>%
      as.list()

    g_box <- agg_tab %>%
      # dplyr::filter(Severity %in% c("mild", "severe")) %>%
      ggboxplot(x = "SNP_merge", y = str_glue("avg.{per_gene}"), color = "SNP_merge", add = "jitter") +
      stat_compare_means(comparisons = my_comps, method = "t.test") +
      labs(x = "Genotype", y = "Expression (average)")

    ggsave(str_glue("{pjdir}/outputs/scRNA-seq/tar_sceqtl/boxplot_{per_gene}.pdf"),
      plot = g_box
    )
  }
}



#
## DE genes by rs6800484
#
fc_thr <- 0.15 # avg_log2FC threshold
pv_thr <- 0.01 # p_val_adj threshold

color_names <- c(
  paste0(
    "abs(log2FC)", c(">", ">", "<"), fc_thr, ",P-adj", c("<", ">", "<"), pv_thr
  ),
  "Others"
)
color_map <- c("darkred", "darkblue", "darkgreen", "gray")
names(color_map) <- color_names


# Target gene
kept_gene <- covid_obj %>%
  rownames() %>%
  purrr::discard(~ str_detect(.x, "^MT") | str_detect(.x, "^RP"))

tar_gene <- "CCR2"


check_severe <- TRUE
if (check_severe) { # CC vs TT in severe
  ident_1 <- "CC"
  group_by <- "SNP"
  tar_cell <- "cMono"
  tar_cohort <- "severe"
  comp <- paste0(ident_1, ".vs.TT")

  kept_cells <- covid_obj@meta.data %>%
    dplyr::mutate(
      kept_cells = !is.na(SNP) &
        celltypeL0 == tar_cell &
        Severity == "severe" &
        SNP %in% c("CC", "TT")
    ) %$%
    kept_cells
} else { # Severe CC vs Post TT
  ident_1 <- "CC"
  group_by <- "SNP"
  tar_cell <- "cMono"
  tar_cohort <- "severeCC.postTT"
  comp <- paste0(ident_1, ".vs.TT")

  kept_cells <- covid_obj@meta.data %>%
    dplyr::mutate(
      kept_cells = !is.na(SNP) &
        celltypeL0 == tar_cell &
        ((Severity == "severe" & SNP == "CC") |
          (Severity == "post" & SNP == "TT"))
    ) %$%
    kept_cells
}

# Cell subset
sub_obj <- covid_obj[, kept_cells]


# Box and violin plot
ccr2_expr <- sub_obj@assays$RNA@data[tar_gene, ] %>%
  as.data.frame() %>%
  purrr::set_names(tar_gene) %>%
  dplyr::mutate(cellbarcodes = row.names(.)) %>%
  dplyr::inner_join(sub_obj@meta.data, by = "cellbarcodes")

g_box_save_to <- file.path(
  pjdir, "outputs/scRNA-seq/de_genes",
  str_glue("deg_boxplot-{tar_cohort}-{tar_cell}-{comp}-{tar_snp}.pdf")
)
g_box <- ggplot(ccr2_expr, aes_string(x = group_by, y = tar_gene)) +
  geom_violin(width = 0.4) +
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  labs(
    x = str_glue("Genotype ({tar_snp})"),
    y = str_glue("{tar_gene} avg expression")
  ) +
  theme_classic() +
  theme(text = element_text(size = 15))
ggsave(g_box_save_to, plot = g_box, width = 5, height = 7)


# DEG
de_tab_save_to <- file.path(
  pjdir, "outputs/scRNA-seq/de_genes",
  str_glue("deg_table-{tar_cohort}-{tar_cell}-{comp}-{tar_snp}.csv")
)

de_tab <- Seurat::FindMarkers(sub_obj,
  ident.1 = ident_1, group.by = group_by, method = "wilcox",
  min.pct = 0.1, features = kept_gene, logfc.threshold = 0.01
) %>%
  dplyr::mutate(
    cell = tar_cell, comp = comp, gene = rownames(.), method = "wilcox",
    direction = dplyr::if_else(avg_log2FC < 0, "Down", "Up"),
    dot_color = case_when(
      p_val_adj < pv_thr & abs(avg_log2FC) > fc_thr ~ names(color_map)[1],
      p_val_adj > pv_thr & abs(avg_log2FC) > fc_thr ~ names(color_map)[2],
      p_val_adj < pv_thr & abs(avg_log2FC) < fc_thr ~ names(color_map)[3],
      T ~ names(color_map)[4]
    )
  )
data.table::fwrite(de_tab, de_tab_save_to)


# Volvano plot
ext_genes <- c("CCR1", "CCR2")
g_pnt <- ggplot() +
  geom_point(
    mapping = aes(x = avg_log2FC, y = -log10(p_val_adj), color = dot_color),
    data = de_tab, alpha = 0.5
  ) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dotted") +
  geom_hline(yintercept = -log10(pv_thr), linetype = "dotted") +
  geom_label_repel(
    aes(
      label = gene, x = avg_log2FC, y = -log10(p_val_adj), fill = gene_color
    ),
    de_tab %>%
      dplyr::mutate(gene_color = if_else(gene %in% ext_genes, "COVID-19", "Others")) %>%
      dplyr::filter(dplyr::row_number() <= 25 | gene %in% ext_genes),
    max.overlaps = 25, min.segment.length = 0, show.legend = FALSE, alpha = 0.7
  ) +
  scale_color_manual(name = NULL, values = color_map) +
  scale_fill_manual(values = c("COVID-19" = "pink", "Others" = "white")) +
  labs(
    x = quote(~ Log[2] ~ "(fold-change)"), y = quote(~ -Log[10] ~ "(P-adj)")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )

g_pnt_save_to <- file.path(
  pjdir, "outputs/scRNA-seq/de_genes",
  str_glue("deg_volcano-{tar_cohort}-{tar_cell}-{comp}-{tar_snp}.pdf")
)
ggsave(g_pnt_save_to, plot = g_pnt, width = 10, height = 7)
