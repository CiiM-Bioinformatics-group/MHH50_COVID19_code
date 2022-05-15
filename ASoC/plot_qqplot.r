#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 05, 2021
# Updated: Apr 05, 2022

options(stringsAsFactors = FALSE, error = traceback)

library(ggrepel)
library(magrittr)
library(tidyverse)
library(data.table)

# Parameters
label_min_p <- 5e-5
wkdir <- file.path("~/Documents/projects/wp_covid19_mhh50", "outputs/scATAC-seq/summary")


# TFs per subcluster of cMono (severe specific)
# "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", # C4.Severe
# "NFE2L2", "BATF", "BACH2", "FOSB", "BACH1", "JDP2", "JUND", "SMARCC1", "FOSL1", "JUNB", "JUN", "FOSL2", "FOS", "NFE2", # C3.Severe
# "DNMT1", "KLF15", "EGR1", "SP4", "ZNF263", "ZFX", "PURA", "KLF7", "ZBTB7A", "KLF16", "WT1", "SP3", "SP1", "ZNF148", "KLF4", "KLF5", "SP2" # C2.Severe
severe_specific <- c("CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG")


# COVID19 GWAS SNPs in around the target TFs (tar_tf)
if (F) {
  snp_around_tfgene_file <- "../../outputs/scATAC-seq/tfmotif/gwas_snp_around_tf_gene.tbl"
  snp_around_tfgene_tab <- fread(snp_around_tfgene_file, sep = "\t") %>%
    group_by(TFname) %>%
    arrange(SNPPval) %>%
    mutate(
      obs = -log10(SNPPval),
      exp = -log10(ppoints(n())),
      ciup = -log10(qbeta(p = (1 - 0.95) / 2, seq(n()), rev(seq(n())))),
      cilw = -log10(qbeta(p = (1 + 0.95) / 2, seq(n()), rev(seq(n())))),
      condition = if_else(TFname %in% severe_specific, "Severe", "Post")
    ) %>%
    ungroup()

  snp_around_tfgene_tab %>%
    group_by(TFname) %>%
    summarise(plot = {
      cdata <- cur_data() %>%
        select(obs, exp, SNPId) %>%
        group_by(SNPId) %>%
        summarise(obs = max(obs), exp = max(exp))

      g <- cdata %>%
        ggplot(aes(x = exp, y = obs)) +
        geom_point(size = 1) +
        geom_line(size = 0.05) +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
        geom_text_repel(aes(label = SNPId), . %>% filter(obs >= -log10(label_min_p)),
          min.segment.length = 0, nudge_x = 0.2
        ) +
        labs(
          x = quote("Expected" ~ -Log[10] ~ "(p-value)"),
          y = quote("Observed" ~ -Log[10] ~ "(p-value)")
        ) +
        theme_classic()

      tf_name <- cur_group()
      save_name <- str_glue("../../outputs/scATAC-seq/tfmotif/gwas_snp_around_tf_gene_{tf_name}.pdf")
      ggsave(save_name, plot = g, width = 5.5, height = 5)

      cdata %>% dplyr::filter(obs >= -log10(label_min_p))
    }) %>%
    as.data.table()
}

# COVID19 GWAS SNPs located in TF motifs
asoc_snp_file <- "../../outputs/scATAC-seq/summary/ver_beta/pseudo-bulk_per-celltype_asoc-snps_annotated.csv"
asoc_snp_tab <- fread(asoc_snp_file)

snp_in_tfmotifs_file <- "../../outputs/scATAC-seq/tfmotif/gwas_snp_in_tf_motifs.tbl"
snp_in_tfmotifs_tab <- fread(snp_in_tfmotifs_file, sep = "\t") %>%
  group_by(MotifTFName) %>%
  separate(col = MotifTFName, into = c("MotifTFName", "CellType", "Condition"), sep = "@|\\.") %>%
  filter(MotifTFName != "no_motif") %>%
  arrange(SNPPval) %>%
  mutate(
    obs = -log10(SNPPval),
    exp = -log10(ppoints(n())),
    ciup = -log10(qbeta(p = (1 - 0.95) / 2, seq(n()), rev(seq(n())))),
    cilw = -log10(qbeta(p = (1 + 0.95) / 2, seq(n()), rev(seq(n())))),
    Specific = if_else(MotifTFName %in% severe_specific,
      "SevereSpecific", "PostSpecific"
    )
  ) %>%
  ungroup()

# QQ plot per motif
snp_in_tfmotifs_tab %>%
  group_by(MotifTFName, Condition) %>%
  summarise(plot = {
    cdata <- cur_data() %>%
      select(obs, exp, SNPId) %>%
      group_by(SNPId) %>%
      summarise(obs = max(obs), exp = max(exp))

    g <- cdata %>%
      ggplot(aes(x = exp, y = obs)) +
      geom_point(size = 1) +
      geom_line(size = 0.05) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
      geom_text_repel(aes(label = SNPId),
        cdata %>% filter(obs >= -log10(label_min_p)),
        min.segment.length = 0, nudge_x = 0.2
      ) +
      labs(
        x = quote("Expected" ~ -Log[10] ~ "(p-value)"),
        y = quote("Observed" ~ -Log[10] ~ "(p-value)")
      ) +
      theme_classic()

    cgroup <- cur_group()
    tfname <- cgroup["MotifTFName"]
    cond <- cgroup["Condition"]
    save_name <- str_glue("../../outputs/scATAC-seq/tfmotif/gwas_snp_in_tf_motifs_{tfname}.{cond}.pdf")
    ggsave(save_name, plot = g, width = 5.5, height = 5)

    cdata %>% dplyr::filter(obs >= -log10(label_min_p))
  }) %>%
  as.data.table()

#     MotifTFName Condition plot.SNPId plot.obs plot.exp
#  1:       CEBPD      Mild  rs2004260   7.0339   3.4115
#  2:       CEBPD    Severe  rs2004260   7.0339   3.3389
#  3:       RUNX1    Severe rs12626952   4.4249   2.9343
#  4:       RUNX2      Mild  rs1491959   5.5104   3.2224
#  5:       RUNX2      Mild rs55886855   4.8226   2.9905
#  6:       RUNX2      Mild   rs981792   4.3855   2.9088
#  7:       RUNX2      Post  rs1491959   5.5104   3.1741
#  8:       RUNX2      Post   rs981792   4.3855   2.8847
#  9:       RUNX2    Severe  rs1491959   5.5104   3.1306
# 10:       RUNX2    Severe rs55886855   4.8226   2.9615
# 11:       RUNX2    Severe   rs981792   4.3855   2.8618
# 12:        SPI1      Mild  rs2154567   4.9220   3.0911
# 13:        SPI1      Mild rs41289622  43.2528   4.4529
# 14:        SPI1      Post  rs2154567   4.9220   3.0549
# 15:        SPI1    Severe  rs2154567   4.9220   3.0215
# 16:        SPI1    Severe rs41289622  43.2528   3.9757
# 17:       STAT2      Mild rs35951367  29.8415   3.7539
# 18:       STAT2      Post rs35951367  29.8415   3.6078
# 19:       STAT2    Severe rs34953890   5.8879   3.2768
# 20:       STAT2    Severe rs35951367  29.8415   3.4986

# QQ plot combined
motif_order <- c(
  "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG",
  "BCL11A", "BCL11B", "IRF4", "RUNX1", "RUNX2", "SPI1", "STAT2"
)
cond_order <- c("Severe", "Mild", "Post")
spec_order <- c("SevereSpecific", "PostSpecific")

snp_in_tfmotifs_plot_tab <- snp_in_tfmotifs_tab %>%
  group_by(SNPId, MotifTFName, Condition, Specific) %>%
  summarise(obs = max(obs), exp = max(exp)) %>%
  ungroup() %>%
  mutate(
    MotifTFName = factor(MotifTFName, levels = motif_order),
    Condition = factor(Condition, levels = cond_order),
    Specific = factor(Specific, levels = spec_order)
  ) %>%
  filter(Specific == "SevereSpecific")

g <- snp_in_tfmotifs_plot_tab %>%
  ggplot(aes(x = exp, y = obs)) +
  geom_point(aes(color = MotifTFName), size = 1) +
  geom_line(aes(color = MotifTFName), size = 0.05) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  geom_text_repel(aes(label = SNPId),
    snp_in_tfmotifs_plot_tab %>% filter(obs >= -log10(label_min_p)),
    min.segment.length = 0, nudge_x = 0.2, max.overlaps = 20
  ) +
  facet_grid(~Condition) +
  labs(
    x = quote("Expected" ~ -Log[10] ~ "(p-value)"),
    y = quote("Observed" ~ -Log[10] ~ "(p-value)")
  ) +
  theme_classic()

# Save to PNG for quick check (low resolution)
save_name <- str_glue("../../outputs/scATAC-seq/tfmotif/pngs/gwas_snp_in_tf_motifs.png")
ggsave(save_name, plot = g, width = 11, height = 5, dpi = 300)

# Save to PDF for high resolution
save_name <- str_glue("../../outputs/scATAC-seq/tfmotif/gwas_snp_in_tf_motifs.pdf")
ggsave(save_name, plot = g, width = 11, height = 5)


inner_join(asoc_snp_tab, snp_in_tfmotifs_tab, by = c("variantID" = "SNPId")) %>%
  filter(p_value_adj < 0.05, SNPPval < 5e-5) %>%
  select(MotifTFName, SNPPval, variantID, p_value_adj)
