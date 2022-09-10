#!/bin/bash

#==============================================================================
# Code to perform simple realignment and visualization
#==============================================================================

### Dotkits ###
source /broad/software/scripts/useuse
use .samtools-1.7
use .minimap2-2.11
use BWA
use Bowtie2
use Picard-Tools
use .macs2-2.1.2
###############

### Paths ###
original='/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/hg38.chrXYM.fa'
assembly='h19/6-bridge_contigs/polished_contigs.fasta'
repaired='hg38_repair/hg38.chrXYM.repaired.fa'
#############

# Use BWA
bwa mem -t 8 $repaired $assembly | samtools sort | \
	samtools view -b > assembly.bam && samtools index assembly.bam

# Use minimap2
minimap2 -a -t 8 -x map-ont $repaired $assembly | \
	samtools sort -T tmp -o assembly2.bam  && samtools index assembly2.bam

# Realign other nanopore reads
minimap2 -a -t 8 -x map-ont $repaired <(cat *.fastq) |\
	samtools sort -T tmp -o realigned.bam  && samtools index realigned.bam

# Use Bowtie2 for illumina reads (e.g. ENCODE)
preproc () {
	fastq=$1
	prefix=${fastq%_*}
	ref=$2

	bowtie2 -k 1 --mm --threads 8 -x ${ref} -U $fastq | samtools view -b -F 1804 -q 30 | samtools sort > ${prefix}_filt.bam

	java -Xmx30g -jar $PICARD MarkDuplicates I=${prefix}_filt.bam O=${prefix}_dup.bam M=${prefix}_dupMet.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

	samtools view -F 1804 -b ${prefix}_dup.bam > ${prefix}_dedup.bam && samtools index ${prefix}_dedup.bam
}

# Perform realignment on all fastqs in folder (assumed single-end)
export -f preproc
ls -1 *fastq* | xargs -I {} bash -c "preproc {} $repaired"
# macs2 callpeak -t *_dedup.bam -f BAM -n $prefix -g hs -p 0.01 --nomodel --shift 0 --extsize 200 --keep-dup all --SPMR

preproc_pe () {
	ref=$1
	fastq1=$(find . -maxdepth 1 -name "*R1*fastq.gz")
	fastq2=$(find . -maxdepth 1 -name "*R2*fastq.gz")
	prefix=${fastq1%_*}

	bowtie2 -k 1 -X2000 --mm --threads 8 -x ${ref} -1 $fastq1 -2 $fastq2 | samtools view -b -F 1804 -q 30 | samtools sort > tmp.filtered1.bam
	samtools fixmate -r tmp.filtered1.bam tmp.fixed.bam
	samtools view -F 1804 -f 2 -u tmp.fixed.bam | samtools sort > ${prefix%_*}_filt.bam

	java -Xmx30g -jar $PICARD MarkDuplicates I=${prefix}_filt.bam O=${prefix}_dup.bam M=${prefix}_dupMet.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

	samtools view -F 1804 -b ${prefix}_dup.bam > ${prefix}_dedup.bam && samtools index ${prefix}_dedup.bam
}

# Perform realignment on all fastqs in folder (assumed paired-end and only one set)
preproc_pe $repaired
# macs2 callpeak -t *_dedup.bam -f BAM -n $prefix -g hs -p 0.01 --nomodel --shift 0 --extsize 200 --keep-dup all --SPMR

# Code for merging BAMs and generating TDFs
# # Slower version
# mergeBam () {
#     file1=$(find . -name "*dedup.bam" | head -n 1)
#     prefix=${file1%_*}

#     samtools view -h $file1 > ${prefix%_*}_combined.sam

#     find . -name "*dedup.bam" | tail -n +2 | xargs -I {} samtools view {} >> ${prefix%_*}_combined.sam

#     samtools sort ${prefix%_*}_combined.sam | samtools view -b > ${prefix%_*}_combined.bam && samtools index ${prefix%_*}_combined.bam
# }
mergeBam () {
	file1=$(find . -name "*dedup.bam" | head -n 1)
	files=$(find . -name "*dedup.bam")
	arg=$(echo $files | sed 's/ / I=/g')
	prefix=$(basename ${file1%_*})

	java -Xmx30g -jar $PICARD MergeSamFiles AS=true I=$arg O=${prefix%_*}_combined.bam
	samtools index ${prefix%_*}_combined.bam
}

toTDF () {
	bam=$(find . -name "*combined.bam")
	prefix=${bam%.*}
	genome=${1/fa/chrom.sizes}
	locus=$2
	
	igvtools count --query $locus $bam $prefix.tdf $genome
}

locus="chr11:1400000-2400000"
mergeBam
toTDF $repaired $locus