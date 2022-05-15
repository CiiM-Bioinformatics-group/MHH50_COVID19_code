#!/bin/bash
#
## Utility functions
#
# Error message
err() {
    echo -e "$(tput bold; tput setaf 1)[E]: $(tput sgr0)$@ Exit..." >&2
    exit
}

# Warning message
wrn() {
    echo -e "$(tput bold; tput setaf 3)[W]: $(tput sgr0)$@" >&2
}

# Information message
msg() {
    echo -e "$(tput bold; tput setaf 7)[I]: $(tput sgr0)$@"
}

# Save my ass
set -Ee -o pipefail

projdir=~/Documents/projects/wp_covid19_mhh50

# Fetch motifs of target TFs
# tar_tf="CEBPA CEBPB CEBPD CEBPE CEBPG" # S4.Severe
# tar_tf=$tar_tf" NFE2L2 BATF BACH2 FOSB BACH1 JDP2 JUND SMARCC1 FOSL1 JUNB JUN FOSL2 FOS NFE2" # C3.Severe
# tar_tf=$tar_tf" DNMT1 KLF15 EGR1 SP4 ZNF263 ZFX PURA KLF7 ZBTB7A KLF16 WT1 SP3 SP1 ZNF148 KLF4 KLF5 SP2" # C2.Severe

#
tar_tf="CEBPA CEBPB CEBPD CEBPE CEBPG" # C2
tar_tf=$tar_tf" RUNX1 RUNX2 SPI1 IRF4 STAT2 BCL11A BCL11B"

~/tools/R/4.0.0/bin/Rscript -e '
options(stringsAsFactors=FALSE)

library(biomaRt)
library(magrittr)
library(tidyverse)
library(data.table)

winsize <- 2.5e4

# C2 severe, C3 severe, and C4 severe.
# 
tar_tf <- c("CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", # S4.Severe
            "NFE2L2", "BATF", "BACH2", "FOSB", "BACH1", "JDP2", "JUND", "SMARCC1", "FOSL1", "JUNB", "JUN", "FOSL2", "FOS", "NFE2", # C3.Severe
            "DNMT1", "KLF15", "EGR1", "SP4", "ZNF263", "ZFX", "PURA", "KLF7", "ZBTB7A", "KLF16", "WT1", "SP3", "SP1", "ZNF148", "KLF4", "KLF5", "SP2" # C2.Severe
)
tar_tf <- c("CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", # C4.severe
            "RUNX1", "RUNX2", "SPI1", "IRF4", "STAT2", "BCL11A", "BCL11B") # C2.post-covid


# Fetch TF position
if (FALSE) listEnsembl()
ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

tar_attrs <- c("chromosome_name", "transcript_start", "transcript_end",
               "strand", "ensembl_gene_id", "external_gene_name")
tar_flkey <- c("chromosome_name", "biotype", "external_gene_name",
              "transcript_gencode_basic", "transcript_is_canonical")
tar_flval <- list(as.character(1:22), c("protein_coding"), tar_tf, TRUE, TRUE)
tfpos_tab <- getBM(attributes=tar_attrs, filters=tar_flkey, value=tar_flval,
                   mart=ensembl)

if (!all(c("upper", "lower") %in% colnames(tfpos_tab)))
  tfpos_tab %>%
    dplyr::mutate(upper=transcript_start - winsize/2,
                  lower=transcript_end + winsize/2) %>%
    dplyr::rowwise() %>%
    dplyr::rename("chrom"="chromosome_name", "start"="upper", "stop"="lower",
                  "name"="external_gene_name") %>%
    dplyr::select(chrom, start, stop, name) %>%
    dplyr::arrange(chrom, start, stop) %>%
    fwrite("../../outputs/scATAC-seq/tfmotif/tf_gene_pos.bed", sep="\t",
           col.names=FALSE)'

echo $tar_tf \
  | tr " " "\n" \
  | xargs -n 1 -I '{}' zgrep {}_ $projdir/inputs/scATAC-seq/all_motif.txt.gz \
  | awk '/^chr[0-9]{1,2}/{split($NF, A, "_"); gsub("chr", "", $1); print $1"\t"$2"\t"$3"\t"A[1]}' \
  | sort -k1,1g -k2,2g -k3,3g \
  > $projdir/outputs/scATAC-seq/tfmotif/tf_motif_pos.bed

# Enrich_motif.txt header line
awk -v mincount=50 '/PeakID/ { print "#CHROM\tSTART\tEND\tNAME"; next}
/chrX|chrMT/ {next} 
{
  split($1, POS, ":|-|\\|");
  gsub("chr", "", POS[1]);
  gsub("_[0-9]+", "", POS[4]);
  record=POS[1]"\t"POS[2]"\t"POS[3]"\t"POS[4];

  if($2>=mincount) {print record"@cMono.Mild"}
  if($13>=mincount) {print record"@cMono.Post"}
  if($14>=mincount) {print record"@cMono.Severe"}
}' $projdir/inputs/scATAC-seq/Enrich_motif.txt \
  | sort -k1,1g -k2,2g -k3,3g \
  > $projdir/outputs/scATAC-seq/tfmotif/atacseq_peak_pos.bed

zcat $projdir/inputs/covid19_gwas/covid19hg.org/COVID19_HGI_A2_ALL_leave_23andme_20210607.txt.gz \
  | awk '$1<23 {name=$3","$4","$15","$7","$9","$14; if(NR>1){start=$2-1; end=$2+1} else{start="START"; end="END"}; print $1"\t"start"\t"end"\t"name}' \
  | sort -k1,1g -k2,2g -k3,3g \
  > $projdir/outputs/scATAC-seq/tfmotif/gwas_snp.bed

# Interset ATACseq peak with motif
bedtools intersect \
  -a $projdir/outputs/scATAC-seq/tfmotif/tf_motif_pos.bed \
  -b $projdir/outputs/scATAC-seq/tfmotif/atacseq_peak_pos.bed \
  -wa -wb \
  | awk '{split($8, INF, "@"); if($4==INF[1]) {print $1"\t"$2"\t"$3"\t"$8}}' \
  > $projdir/outputs/scATAC-seq/tfmotif/tf_motif_pos_in_atacseq_peak.bed

{
  echo "CHROM\tMotifStart\tMotifEnd\tMotifTFName\tSNPPos\tSNPRef\tSNPAlt\tSNPId\tSNPBeta\tSNPPval\tSNPAF";
  bedtools intersect \
    -a $projdir/outputs/scATAC-seq/tfmotif/tf_motif_pos_in_atacseq_peak.bed \
    -b $projdir/outputs/scATAC-seq/tfmotif/gwas_snp.bed \
    -wa -wb \
    | tr "," "\t" \
    | cut -f1-4,7-13
} > $projdir/outputs/scATAC-seq/tfmotif/gwas_snp_in_tf_motifs.tbl

{
  echo "CHROM\tTFUpStart\tTFDwEnd\tTFname\tSNPPos\tSNPRef\tSNPAlt\tSNPId\tSNPBeta\tSNPPval\tSNPAF";
  bedtools intersect \
    -a $projdir/outputs/scATAC-seq/tfmotif/tf_motif_pos_in_atacseq_peak.bed \
    -b $projdir/outputs/scATAC-seq/tfmotif/gwas_snp.bed \
    -wa -wb \
    | tr "," "\t" \
    | cut -f1-4,7-13
} > $projdir/outputs/scATAC-seq/tfmotif/gwas_snp_around_tf_gene.tbl

msg Job finished!
