#!/bin/bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 12, 2021
# Updated: May 26, 2021

set -E -o pipefail

# A script to estimate allele-specific open chromatin.


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



#
## Resource, workspace, tools, and input files
#
cpus=10

# Workspace
projdir=/vol/projects/wli/projects/mhh50
seqtech=scATAC-seq
runver=${1:?run_v1 or run_v2 should be given!} # First round of ATAC-seq data, remove duplicated reads.

iptdir=$projdir/inputs
optdir=$projdir/outputs/$seqtech/$runver/optdir
tmpdir=$projdir/outputs/$seqtech/$runver/tmpdir
rfrdir=$projdir/outputs/rfrdir

# Tools
tooldir=~/mhh50/zzh/tools
py3env=$tooldir/py3env/bin/activate

# Barcode to cell type table: a table presenting which cell type the barcode blongs to.
# Among different pool, the barcode could be duplicated, therefore, it should be splited by pool.
b2c_file=$iptdir/idmapping/metadata.Dm2.with_donor_id.txt
imp_file=$iptdir/idmapping/id_mapping_atac.csv

# Reference genome
# The reference genome which used in WASP pipeline and SNPsplit.
# Reference genome from Ensembl, the toplevel contigs, unmasked, GRCh38
# ref_genome_url=http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
ref_genome=$iptdir/ensembl/reference_genome/Homo_sapiens.GRCh38.dna.toplevel.nonpatch.fa.gz

# The length of each chromosome in the reference genome
chrom_info_file=$iptdir/ensembl/reference_genome/chrom_length_grch38.txt

# Genotypes with a fake FI flag which required by SNPsplit
gntp_file=$projdir/outputs/genotypes/variants_fl_addFI.vcf.gz

# Which operation?
# cmd=snpsplit # Split SNPs by allele.
cmd=aserc    # Count allelic reads.


#
## Set up workspace, update ENV, etc.
mkdir -p $tmpdir $optdir $rfrdir

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$tooldir/lib
[[ ! -z $PYTHONPATH ]] && export PYTHONPATH=$py3env/lib/site-packages
[[ -e $py3env ]] && source $py3env || wrn "$py3env NOT found!"


#
## HDF5 database for WASP pipeline. Only autosomes were included
#
if [[ ! -e $rfrdir/snph5db_dir/snps_tab.h5 ]]; then
    mkdir -p $tmpdir/snph5db_dir $rfrdir/snph5db_dir
    
    $tooldir/bin/bcftools index --force --threads $cpus $gntp_file

    seq 1 22 | xargs -n 1 -I '{}' \
        $tooldir/bin/bcftools view -Oz -o $tmpdir/snph5db_dir/chr{}.vcf.gz $gntp_file {}
    
    $tooldir/bin/snp2h5 \
        --chrom $chrom_info_file \
        --format vcf \
        --snp_tab $rfrdir/snph5db_dir/snps_tab.h5 \
        --snp_index $rfrdir/snph5db_dir/snps_index.h5 \
        --haplotype $rfrdir/snph5db_dir/haplotype.h5 \
        $tmpdir/snph5db_dir/chr*.vcf.gz

    rm -fr $tmpdir/snph5db_dir
fi


#
## Build general Bowtie2 genome index
#
bt2idx=$rfrdir/bt2idx/grch38_unmasked
if [[ ! -e $bt2idx.1.bt2 ]]; then
    mkdir -p $(dirname $bt2idx)
    $tooldir/bin/bowtie2-build --threads $cpus $ref_genome $bt2idx
fi


count_allelic_reads=false
keep_calibrated_bam=true
bam_perct_available=true
case $cmd in
    aserc)
        bam_perct=$optdir/bam_perct
        [[ ! -d $bam_perct ]] && mkdir -p $bam_perct
        [[ ! -d $optdir/gatk ]] && mkdir -p $optdir/gatk

        cur_poolid=NULL
        while read -r record; do
            pool_id=$(cut -f1 -d$'\t' <<<$record | sed 's/pool//g')
            gntp_id=$(cut -f2 -d$'\t' <<<$record)
            donor_id=$(cut -f3 -d$'\t' <<<$record)

            if [[ $cur_poolid != $pool_id ]]; then
                if [[ ! -h $iptdir/bamfiles/$runver/atac_$pool_id ]]; then
                    wrn "Not found $iptdir/bamfiles/$runver/atac_$pool_id. Skipping it..."
                    continue
                fi

                if [[ $bam_perct_available != true ]]; then
                    python3 $tooldir/scripts/split_reads_by_cell_type.py \
                        -b $iptdir/bamfiles/$runver/atac_$pool_id/outs/possorted_bam.bam \
                        -c <(grep -e pool$pool_id -e Pool $b2c_file) \
                        --celltype-col Celltype \
                        --barcode-col Cellbarcode \
                        --donor-col DonorID \
                        --threads $cpus \
                        -o $bam_perct
                else
                    wrn "BAM files splitted by cell type is available, skipping splitting step..."
                fi

                cur_poolid=$pool_id
            fi

            # Donor genotype ID, donor temporary directory
            donor_tmpdir=$tmpdir/$donor_id
            [[ ! -d $donor_tmpdir ]] && mkdir -p $donor_tmpdir

            for celltype in B cMono ncMono CD4T CD8T NK pDC Plasmablast; do
                #
                ## WASP to fix reference bias
                #
                # The input BAM file per donor per cell type
                donor_bamid=${donor_id}-${celltype}
                donor_bamfile=$optdir/bam_perct/$donor_bamid.bam
                [[ ! -e $donor_bamfile ]] && wrn "$donor_bamfile NOT found!" && continue

                # Index
                $tooldir/bin/samtools index -@ $cpus $donor_bamfile

                donor_tmp_pref=$donor_tmpdir/$donor_bamid
                # Remove chr notation
                $tooldir/bin/samtools view -h -@ $cpus $donor_bamfile chr{1..22} \
                    | sed 's/chr//g' \
                    | $tooldir/bin/samtools view -h -@ $cpus -bo $donor_tmp_pref.bam

                # Index
                $tooldir/bin/samtools index -@ $cpus $donor_tmp_pref.bam

                # Find intersecting SNPs
                python3 $tooldir/WASP/mapping/find_intersecting_snps.py \
                    --is_sorted \
                    --is_paired_end \
                    --output_dir $donor_tmpdir \
                    --snp_tab $rfrdir/snph5db_dir/snps_tab.h5 \
                    --snp_index $rfrdir/snph5db_dir/snps_index.h5 \
                    --haplotype $rfrdir/snph5db_dir/haplotype.h5 \
                    --samples $gntp_id \
                    $donor_tmp_pref.bam

                # Rename fake BAM to SAM, and compress SAM
                mv -f $donor_tmp_pref.keep.bam $donor_tmp_pref.keep.sam
                $tooldir/bin/samtools view \
                    -@ $cpus \
                    -hbo $donor_tmp_pref.keep.bam \
                    $donor_tmp_pref.keep.sam

                # Rename fake BAM to SAM, and compress SAM
                mv -f $donor_tmp_pref.to.remap.bam $donor_tmp_pref.to.remap.sam
                $tooldir/bin/samtools view \
                    -@ $cpus \
                    -hbo $donor_tmp_pref.to.remap.bam \
                    $donor_tmp_pref.to.remap.sam

                # Remove SAM
                rm -f $donor_tmp_pref.keep.sam $donor_tmp_pref.to.remap.sam


                # If no reads were piped to *.fq[12].gz, then we skip next steps.
                n_biased_reads=$(zgrep -c . $donor_tmp_pref.remap.fq1.gz)
                if [[ $n_biased_reads -ge 4 ]]; then
                    # Remapping by Bowtie2, sort and compress by SAMtools.
                    # Here use $donor_bamid as read group, including donor ID and cell type.
                    $tooldir/bin/bowtie2 \
                        -X 2000 \
                        -p $cpus \
                        -x $bt2idx \
                        -1 $donor_tmp_pref.remap.fq1.gz \
                        -2 $donor_tmp_pref.remap.fq2.gz \
                        --rg-id $donor_bamid \
                        | $tooldir/bin/samtools sort \
                        | $tooldir/bin/samtools view \
                            -hbo $donor_tmp_pref.remap.bam

                    # Index the remapped reads
                    $tooldir/bin/samtools index -@ $cpus $donor_tmp_pref.remap.bam

                    # Sort *.to.remap.bam
                    $tooldir/bin/samtools sort \
                        -@ $cpus \
                        -o $donor_tmp_pref.to.remap.sort.bam \
                        $donor_tmp_pref.to.remap.bam

                    # Index *.to.remap.bam
                    $tooldir/bin/samtools index -@ $cpus $donor_tmp_pref.to.remap.sort.bam

                    # Filter rempped reads
                    python $tooldir/WASP/mapping/filter_remapped_reads.py \
                        $donor_tmp_pref.to.remap.sort.bam \
                        $donor_tmp_pref.remap.bam \
                        $donor_tmp_pref.remap.keep.bam

                    # Merge
                    $tooldir/bin/samtools merge -f \
                        $donor_tmp_pref.keep.merged.bam \
                        $donor_tmp_pref.remap.keep.bam \
                        $donor_tmp_pref.keep.bam
                else
                    mv -f $donor_tmp_pref.keep.bam $donor_tmp_pref.keep.merged.bam
                fi

                # Sort
                $tooldir/bin/samtools sort \
                    -@ $cpus \
                    -o $donor_tmp_pref.keep.merged.sorted.bam \
                    $donor_tmp_pref.keep.merged.bam

                # Index
                $tooldir/bin/samtools index -@ $cpus $donor_tmp_pref.keep.merged.sorted.bam


                #
                ## GATK/ASEReadCounter to count allelic reads
                #
                if [[ $count_allelic_reads == true ]]; then
                    $tooldir/bin/gatk ASEReadCounter \
                        --disable-read-filter NotDuplicateReadFilter \
                        --variant $gntp_file \
                        --output-format CSV \
                        --reference $ref_genome \
                        --input $donor_tmp_pref.keep.merged.sorted.bam \
                        --output $optdir/gatk/$donor_bamid.gatkAlleleReadCounts.csv \
                        2>&1 | grep -v WARN
                fi

                if [[ $keep_calibrated_bam == true ]]; then
                    if [[ ! -d $optdir/bam_perct_calibrated ]]; then
                        mkdir -p $optdir/bam_perct_calibrated
                    fi

                    mv -f $donor_tmp_pref.keep.merged.sorted.bam $optdir/bam_perct_calibrated
                fi

                rm -f $donor_tmp_pref.*
            done
        done < <(cut -f2,5,11 -d$'\t' $b2c_file | sort -k3,3g | uniq | grep -v DonorID)
        ;;

    snpsplit)
        [[ ! -d $optdir/snpsplit ]] && mkdir -p $optdir/snpsplit
        while read -r record; do
            pool_id=$(cut -f1 -d$'\t' <<<$record)
            gntp_id=$(cut -f2 -d$'\t' <<<$record)
            donor_id=$(cut -f3 -d$'\t' <<<$record)

            # Build Bowtie2 index, the genome should be processed by SNPsplit_genome_preparation from SNPsplit.
            donor_bt2idx=$rfrdir/bt2idx/grch38_nm_$gntp_id
            donor_tmpdir=$tmpdir/$donor_id

            if [[ ! -e $donor_bt2idx.1.bt2 ]]; then
                cd $donor_tmpdir
                $tooldir/bin/SNPsplit_genome_preparation \
                    --vcf_file $gntp_file \
                    --strain $gntp_id \
                    --reference_genome $(dirname $ref_genome) \
                    --genome_build GRCh38
                cd -
                
                cat $donor_tmpdir/$gntp_id"_N-masked"/chr{1..22}.N-masked.fa \
                    | gzip -c > $donor_tmpdir/${gntp_id}_grch38_N_masked.fa.gz
                
                $tooldir/bin/bowtie2-build \
                    --threads $cpus \
                    $donor_tmpdir/${gntp_id}_grch38_N_masked.fa.gz \
                    $donor_bt2idx
                
                rm -fr $donor_tmpdir/{${gntp_id}_N-masked,SNPs_$gntp_id}
                mv -f $donor_tmpdir/all_SNPs_${gntp_id}_GRCh38.txt.gz $optdir/snpsplit
            fi

            
            # Split reads by SNPsplit per cell type
            for celltype in B cMono ncMono CD4T CD8T NK pDC Plasmablast; do
                donor_bamid=${donor_id}-${celltype}
                donor_bamfile=$optdir/bam_perct/$donor_bamid.bam
                [[ ! -d $optdir/snpsplit ]] && mkdir -p $optdir/snpsplit/$donor_bamid

                donor_tmp_pref=$donor_tmpdir/$donor_bamid

                # Create fastq
                $tooldir/bin/samtools fastq \
                    -N \
                    -@ $cpus \
                    -1 $donor_tmp_pref"_R1.fq.gz" \
                    -2 $donor_tmp_pref"_R2.fq.gz" \
                    $donor_bamfile

                # Align
                $tooldir/bin/bowtie2 \
                    -X 2000 \
                    -p $cpus \
                    -x $donor_bt2idx \
                    -1 $donor_tmp_pref"_R1.fq.gz" \
                    -2 $donor_tmp_pref"_R2.fq.gz" \
                    --rg-id $donor_bamid \
                    | $tooldir/bin/samtools sort \
                    | $tooldir/bin/samtools view \
                        -hbo $donor_tmp_pref.nm.bam
                
                # Index
                $tooldir/bin/samtools index -@ $cpus $donor_tmp_pref.nm.bam
                
                # SNPsplit
                $tooldir/bin/SNPsplit \
                    -o $optdir/snpsplit/$donor_bamid \
                    --paired \
                    --no_sort \
                    --singletons \
                    --snp_file $optdir/snpsplit/all_SNPs_${gntp_id}_GRCh38.txt.gz \
                    --samtools_path $tooldir/bin \
                    $donor_tmp_pref.nm.bam

                # Clean up
                rm -f $donor_tmp_pref*
            done
        done < <(cut -f2,5,11 -d$'\t' $b2c_file | sort -k3,3g | uniq | grep -v DonorID)
        ;;

    *)
        err "Bad subcmd"
esac

# vim: set nowrap number relativenumber expandtab cindent tw=500 ts=4:
