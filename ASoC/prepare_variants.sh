#!/bin/bash
#$ -N anno
#$ -cwd
#$ -o /vol/projects/wli/projects/mhh50/log/
#$ -e /vol/projects/wli/projects/mhh50/log/
#$ -l arch=linux-x64
#$ -b n
#$ -l vf=100G
#$ -pe multislot 2
#$ -q all.q
#$ -i /dev/null

# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 18, 2021
# Updated: May 18, 2021

#
## A short script to annotate and query inputed variants (DS field is required)
#

# Usage: bash annotate_variants.sh vcf_file id_mapping_file

iptvcf=/vol/projects/wli/projects/mhh50/inputs/genotypes/covid.merge.flt.recode.vcf.gz
id_map_g2p=/vol/projects/wli/projects/mhh50/inputs/idmapping/id_mapping_g2p.txt
chrom_map=/vol/projects/wli/projects/mhh50/inputs/idmapping/chrom_mapping.txt

# export bcftools=/vol/biotools/bin/bcftools
# export PATH=$PATH:/vol/projects/wli/projects/mhh50/code/tools/bin
# export BCFTOOLS_PLUGINS=/vol/projects/wli/projects/mhh50/code/tools/libexec/bcftools
bcftools=~/mhh50/zzh/tools/bin/bcftools

iptdir=/vol/projects/wli/projects/mhh50/inputs/genotypes
optdir=/vol/projects/wli/projects/mhh50/outputs/genotypes
tmpdir=/vol/projects/wli/projects/mhh50/outputs/genotypes/tmpdir
cpus=10

# Prepare workspace
mkdir -p $optdir $tmpdir

# Download dbSNP SNPs (b151)
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz
wget -cP $iptdir \
    -o $iptdir/wget.log \
    https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz

wget -cP $iptdir \
    -o $iptdir/wget.log \
    https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz.tbi

# Index the input VCF file
# $bcftools index -f --threads $cpus $iptvcf

# Reheader the VCF: change genotype ID to scRNA-seq ID
$bcftools reheader \
    -s $id_map_g2p \
    -o $tmpdir/variants_nh.vcf.gz \
    $iptvcf

#Index the reheadered vcf file
$bcftools index --threads $cpus $tmpdir/variants_nh.vcf.gz

# Rename chrs
$bcftools annotate \
    -O z \
    -o $tmpdir/variants_nc.vcf.gz \
    --rename-chrs $chrom_map \
    --threads $cpus \
    $tmpdir/variants_nh.vcf.gz

# Index
$bcftools index --threads $cpus $tmpdir/variants_nc.vcf.gz

# Annotate the genotypes
# Add MAF and HWE
# Filter variants. NO MAF, R2, or ER2 in FROMAT, this is unusual!
# -i '%ID!="." & MAF[0]>=0.1 & (R2[0]>=0.3 || ER2[0]>=0.3)' \
$bcftools view \
    -S $iptdir/sample_ids.txt \
    -o $tmpdir/variants_ss.vcf.gz \
    --threads $cpus \
    $tmpdir/variants_nc.vcf.gz \
    {{1..22},X,Y}

$bcftools index --threads $cpus $tmpdir/variants_ss.vcf.gz

$bcftools annotate \
	-c ID \
        -a $iptdir/00-common_all.vcf.gz \
        --threads 4 \
	$tmpdir/variants_ss.vcf.gz \
    | $bcftools norm -d snps \
    | $bcftools view \
	    -O z \
	    -v snps \
	    -i '%ID!="." && COUNT(GT=="het")>=3 && COUNT(GT=="RR") >=3 && COUNT(GT=="AA") >=3' \
	    --threads 4 \
	    -o $optdir/variants_fl.vcf.gz

# Index the annotated and filtered VCF file
$bcftools index --threads $cpus $optdir/variants_fl.vcf.gz

# Query variant information
$bcftools query \
    -Hf '%ID %CHROM %POS %REF %ALT\n' \
    $optdir/variants_fl.vcf.gz \
    | awk '/^#/ {print "ID CHROM POS REF ALT"; next} {print}' \
    | gzip -f > $optdir/variants_info.txt.gz

# Query dosage of the genotypes
$bcftools query \
    -Hf '%ID[ %DS]\n' \
    $optdir/variants_fl.vcf.gz \
    | awk -F' ' '/^#/ {for(i=2;i<=NF;i++){nm=gensub(/\[[0-9]+\](.+)($|:DS)/,"\\1","g",$i); if(i<NF){printf nm" "} else{print nm}}; next} {print}' \
    | gzip -f > $optdir/variants_doage.txt.gz

# rm -f $tmpdir/variants*
echo "Check $optdir for annotated and splitted genotypes (dosage + variant information)"

# vim: set nocompatible number nowrap:
