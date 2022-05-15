#!/bin/bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 27, 2021
# Updated: Jan 29, 2022

# XXX: Using bash, zsh doesn't function sometime.

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

pjdir=~/Documents/projects/wp_covid19_mhh50

# Epigenomic anotation with corresponding description.
# declare -A celltype=( [E029]="Mono" [E032]="B" [E034]="T" [E046]="NK" [E062]="PBMC" )
declare -A celltype=(
  [E062]="Primary_mononuclear_cells_from_peripheral_blood"
  [E034]="Primary_T_cells_from_peripheral_blood"
  [E045]="Primary_T_cells_effector_memory_enriched_from_peripheral_blood"
  [E044]="Primary_T_regulatory_cells_from_peripheral_blood"
  [E043]="Primary_T_helper_cells_from_peripheral_blood"
  [E039]="Primary_T_helper_naive_cells_from_peripheral_blood"
  [E040]="Primary_T_helper_memory_cells_from_peripheral_blood_1"
  [E037]="Primary_T_helper_memory_cells_from_peripheral_blood_2"
  [E048]="Primary_T_CD8+_memory_cells_from_peripheral_blood"
  [E038]="Primary_T_helper_naive_cells_from_peripheral_blood"
  [E047]="Primary_T_CD8+_naive_cells_from_peripheral_blood"
  [E029]="Primary_monocytes_from_peripheral_blood"
  # [E035]="Primary_hematopoietic_stem_cells"
  # [E051]="Primary_hematopoietic_stem_cells_G-CSF-mobilized_Male"
  # [E050]="Primary_hematopoietic_stem_cells_G-CSF-mobilized_Female"
  # [E036]="Primary_hematopoietic_stem_cells_short_term_culture"
  [E032]="Primary_B_cells_from_peripheral_blood"
  [E046]="Primary_Natural_Killer_cells_from_peripheral_blood"
  # [E030]="Primary_neutrophils_from_peripheral_blood"
  # [E115]="Dnd41_TCell_Leukemia_Cell_Line"
  # [E116]="GM12878_Lymphoblastoid_Cells"
  # [E123]="K562_Leukemia_Cells"
  [E124]="Monocytes-CD14+_RO01746_Primary_Cells")

base_url=https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/bed_hg38_lifted_over/
save_path=~/Documents/projects/wp_reference/Roadmap
for sample_id in ${!celltype[@]}; do
  save_file=${sample_id}_25_imputed12marks_hg38lift_dense.bed.gz
  if [[ ! -e $save_path/$save_file ]]; then
    wget -c -P $save_path $base_url/$save_file
    ln -s $save_path/$save_file $pjdir/inputs/Roadmap
  else
    msg "The epigenomic database for sample $sample_id was ready."
  fi
done

# Annotate ASoC SNPs by intersecting them with Roadmap chromatin states
echo "contig,position,variantID,chrom_start,chrom_stop,chrom_state,celltype" \
  > $pjdir/outputs/scATAC-seq/summary/heter-snps_epiannot.csv
for x in $pjdir/inputs/Roadmap/*.bed.gz; do
  sample_id=$(basename $x)
  sample_id=${sample_id%%_*}
  if [[ -z ${celltype[$sample_id]} ]]; then
    continue
  fi

  bedtools intersect \
    -a $pjdir/outputs/scATAC-seq/summary/heter-snps.bed \
    -b $x \
    -wb \
    | sed 's/chr//g' \
    | awk -F $'\t' -v CT="${celltype[$sample_id]}" \
      '{print $1","$3","$4","$6","$7","$8","CT}' \
    >> $pjdir/outputs/scATAC-seq/summary/heter-snps_epiannot.csv
done

# This script generates a file as big as 2.1G, be carefule when loading it into R.
# bedtools is much faster than R scripts code using dplyr pipes.

