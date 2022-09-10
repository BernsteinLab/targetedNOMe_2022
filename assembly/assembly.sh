#!/bin/bash

#==============================================================================
# Code for reassembling H19/IGF2 locus
# Pools reads from 5 separate experiments (aligned to hg38)
# Use NECAT to create de novo assembly
#==============================================================================

### Dotkits ###
source /broad/software/scripts/useuse
use .samtools-1.7
use Java-1.8
use -q .perl-5.28.0
###############

### Paths ###
k562_fq="/seq/epiprod02/Battaglia/NanoNOMe/K562merged/workspace/K562.fastq.gz"
h9_fq="/seq/epiprod02/Battaglia/NanoNOMe/H9_ESC/H9_merged/H19_SegDup/workspace/H9.fastq.gz"
gm_fq="/seq/epiprod02/Battaglia/NanoNOMe/190330_GpC/190330.fastq.gz"
h9_fq_2="/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9.fastq.gz"
hsmm_fq="/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/HSMM/workspace/HSMM.fastq.gz"

k562_bam="/seq/epiprod02/Battaglia/NanoNOMe/K562merged/workspace/K562.bam"
h9_bam="/seq/epiprod02/Battaglia/NanoNOMe/H9_ESC/H9_merged/H19_SegDup/workspace/H9.bam"
gm_bam="/seq/epiprod02/Battaglia/NanoNOMe/190330_GpC/190330.bam"
h9_bam_2="/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9.bam"
hsmm_bam="/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/HSMM/workspace/HSMM.bam"

bbmap="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/bbmap"
necat="/seq/epiprod02/kdong/software/NECAT/Linux-amd64/bin/"
#############

export PATH=$PATH:${necat}

# Get the read names for all reads that aligned to locus
locus="chr11:1696314-2276795"

samtools view -h ${k562_bam} ${locus} > combined.sam
samtools view ${h9_bam} ${locus} >> combined.sam
samtools view ${gm_bam} ${locus} >> combined.sam
samtools view ${h9_bam_2} ${locus} >> combined.sam
samtools view ${hsmm_bam} ${locus} >> combined.sam

# Generate combined BAM to visualize reads before fixing
samtools sort combined.sam | samtools view -bh > combined.bam && rm combined.sam
samtools index combined.bam
samtools view combined.bam | cut -f1 | sort | uniq > all_names.txt

${bbmap}/filterbyname.sh in=${k562_fq} ignorebadquality ow=t include=t names=all_names.txt out=K562.fastq
${bbmap}/filterbyname.sh in=${h9_fq} ignorebadquality ow=t include=t names=all_names.txt out=H9.fastq
${bbmap}/filterbyname.sh in=${gm_fq} ignorebadquality ow=t include=t names=all_names.txt out=GM12878.fastq
${bbmap}/filterbyname.sh in=${h9_fq_2} ignorebadquality ow=t include=t names=all_names.txt out=H9_2.fastq
${bbmap}/filterbyname.sh in=${hsmm_fq} ignorebadquality ow=t include=t names=all_names.txt out=HSMM.fastq

# # Old search function (very slow)
# searchByName () {
#     original=$1
#     bam=$2
#     locus=$3
#     output=$4

#     zcat $original | awk 'NR % 4 == 1' | cut -d ' ' -f 1 | sort | uniq > names_tmp.txt
#     echo -n '' > tmp
#     samtools view $bam $locus | cut -f1 | xargs -I {} grep -n {} names_tmp.txt |\
# 		cut -d: -f1 | awk '{OFS=","};{print ($1-1)*4+1, ($1-1)*4+4}' >> tmp

#     echo -n '' > $output
#     cat tmp | xargs -I {} sed -n '{}p' <(zcat $original) >> $output && rm tmp names_tmp.txt
# }

# searchByName $k562_fq $k562_bam $locus K562.fastq
# searchByName $h9_fq $h9_bam $locus H9.fastq
# searchByName $gn_fq $gm_bam $locus GM12878.fastq
# searchByName $h9_fq_2 $h9_bam_2 $locus H9_2.fastq
# searchByName $hsmm_fq $hsmm_bam $locus HSMM.fastq

# Write fastqs to read_list for NECAT
ls -1 $PWD/*.fastq > read_list.txt

# Create config file for NECAT prior to running
necat.pl correct h19_config.txt
necat.pl assemble h19_config.txt
necat.pl bridge h19_config.txt

# Generate reverse complement in case of opposite orientation
cat h19/6-bridge_contigs/polished_contigs.fasta | \
	while read L; do  echo $L; read L; echo "$L" | \
	rev | tr "ATGC" "TACG" ; done > h19/reverse_complement.fasta