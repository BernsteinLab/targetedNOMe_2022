#!/bin/bash
# 1_primary_analysis.sh executes the first pipeline for analyzing NanoNOMe data
# Performs basecalling with Guppy and flip-flop (--bcall)
# Performs FastQ indexing with nanopolish (--index)
# Aligns FastQ to Hg38 using minimap2 (--align)
# Filters out poorly mapped alignments (--filter)
# Performs CpG and GpC methylation detection (--meth)
# User must provide the following preceded by the appropriate flag:
# Working directory (-wd) with the QSUB command
# Optionally provide a prefix for file names (-n)
# If using --bcall or --index, path to Fast5 directory (-f)
# If using --filter, path to loci BED file (-L) for filtering 

## qsub arguments
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
#$ -l os=RedHat7
#$ -pe smp 8
#$ -binding linear:8
#$ -r y
#$ -j y
#$ -b y
#$ -m e
#$ -o pipeline.log
#$ -e pipeline.error

# Dotkits
source /broad/software/scripts/useuse
reuse -q .minimap2-2.11
reuse -q Samtools
reuse -q Anaconda3
reuse -q .htslib-1.8
reuse -q GCC-5.2
reuse -q .jdk-10-2018.3.20
reuse -q .gatk-4.2.0.0
reuse -q Picard-Tools

# Infer pipeline directory
SCRIPT_DIR=$(readlink -f $(dirname "$0"))

## Additional software/paths
ref="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/hg38.chrXYM.fa"
guppy="/seq/epiprod02/Sofia/software/ont-guppy-cpu/bin/guppy_basecaller"
nanopolish="/seq/epiprod02/kdong/software/nanopolish-13.2/nanopolish"
nanopolish_meth="/seq/epiprod02/kdong/software/nanopolish-cpggpc/nanopolish"
np_frequency_cpg="${SCRIPT_DIR}/scripts/calculate_methylation_frequency_cpg_strands.py"
np_frequency_gpc="${SCRIPT_DIR}/scripts/calculate_methylation_frequency_gpc_strands.py"

# Defaults
PREFIX="output"
loci=""
bcall=false
index=false
align=false
filter=false
realign=false
meth=false

# Parse inputs and store parmeters in appropriate variables
while [ $# -gt 0 ]
do
	case $1 in
		-f)
			shift
			FILE=$1
			shift
			;;
		-n)
			shift
			PREFIX=$1
			shift
			;;
		-L)
			shift
			loci=$1
			shift
			;;
		--bcall)
			bcall=true
			shift
			;;
		--index)
			index=true
			shift
			;;
		--align)
			align=true
			shift
			;;
		--filter)
			filter=true
			shift
			;;
		--meth)
			meth=true
			shift
			;;
		*)
			echo "ERROR: Unknown flag"
			exit 1
			;;
	esac
done

# Check for valid inputs
if [ "$bcall" = true -o "$index" = true ]; then
	# Check if Fast5 is given
	if [ ! -f "$FILE" -a ! -d "$FILE" ]; then
		echo "ERROR: Missing Fast5 directory"
		exit 1
	fi

	# Unzip
	if [ "${FILE:(-3)}" = "zip" ]; then
		DIR=$(echo "${FILE}" | sed 's/....$//')
		if [ ! -d "$DIR" ]; then
			unzip $FILE
		fi
		FILE=$DIR
	fi
fi

if [ ! -z "$loci" ]; then
	# If specified, check for existence of loci BED file
	if [ ! -f "$loci" ]; then
		echo "ERROR: Missing valid loci file"
		exit 2
	fi
fi

D=`date "+%Y%m%d"`

echo "*** Analysis Started ***"
echo "* Date: ${D}"

# Perform basecalling with Guppy plus Flip-Flop 
if [ "$bcall" = true ]; then
	echo "Basecalling with Guppy plus Flip-Flop"

	echo -e "Running command:\n\
	[$guppy -i $FILE -r -s ./fastq -c dna_r9.4.1_450bps_fast.cfg --cpu_threads_per_caller 8]"
	$guppy -i $FILE -r -s ./fastq -c dna_r9.4.1_450bps_fast.cfg --cpu_threads_per_caller 8

	# Concatenate output and gzip
	mkdir -p ./workspace
	cat ./fastq/*.fastq > ./workspace/${PREFIX}.fastq
fi

# Nanopolish index
if [ "$index" = true ]; then
	echo "Indexing FastQ with Nanopolish..."

	# Check working directory
	DIR=""
	if [ -d ./workspace/ ]; then
		DIR="./workspace/"
	fi

	# Get absolute directory of fast5
	fast5_dir=$(readlink -f $FILE)
	echo -e "Running command:\n\
	[$nanopolish index -d $fast5_dir ${DIR}${PREFIX}.fastq]"
	$nanopolish index -d $fast5_dir ${DIR}${PREFIX}.fastq
fi

# Minimap alignment
if [ "$align" = true ]; then
	echo "Aligning reads with minimap2..."

	# Check working directory
	DIR=""
	if [ -d ./workspace/ ]; then
		DIR="./workspace/"
	fi
	echo -e "Running command:\n\
	[minimap2 -a -t 8 -x map-ont $ref ${DIR}${PREFIX}.fastq | samtools sort -T tmp -o ${DIR}${PREFIX}.bam && samtools index ${DIR}${PREFIX}.bam]"
	minimap2 -a -t 8 -x map-ont $ref ${DIR}${PREFIX}.fastq | \
		samtools sort -T tmp -o ${DIR}${PREFIX}.bam && samtools index ${DIR}${PREFIX}.bam
fi

# Filter
if [ "$filter" = true ]; then
	echo "Filtering alignments for MapQ > 50..."	

	# Check working directory
	DIR=""
	if [ -d ./workspace/ ]; then
		DIR="./workspace/"
	fi

	if [ -z "$loci" ]; then
		echo -e "Running command:\n\
	[samtools view -bh -q 50 ${DIR}${PREFIX}.bam > ${DIR}${PREFIX}.filtered.bam]"
		samtools view -bh -q 50 ${DIR}${PREFIX}.bam > ${DIR}${PREFIX}.filtered.bam
	else
		echo -e "Running command:\n\
	[samtools view -bh -q 50 -L $loci ${DIR}${PREFIX}.bam > ${DIR}${PREFIX}.filtered.bam]"
		samtools view -bh -q 50 -L $loci ${DIR}${PREFIX}.bam > ${DIR}${PREFIX}.filtered.bam
	fi
	
	samtools index ${DIR}${PREFIX}.filtered.bam
fi

# Call CpG and GpC methylation
if [ "$meth" = true ]; then
	echo "Calling methylation with Nanopolish..."

	# Check working directory
	DIR=""
	if [ -d ./workspace/ ]; then
		DIR="./workspace/"
	fi

	meth_bam=${DIR}${PREFIX}.filtered.bam
	# Check if filter step has run previously
	if [ ! -f "$meth_bam" ]; then
		meth_bam = ${DIR}${PREFIX}.bam
	fi

	# Call both methylation types simultaneously
	echo -e "Running command:\n\
	[$nanopolish_meth call-methylation --methylation=cpggpc -t 8 -r ${DIR}${PREFIX}.fastq -b $meth_bam -g $ref > ${DIR}${PREFIX}_methylation_calls.tsv]"
	$nanopolish_meth call-methylation --methylation=cpggpc -t 8 -r ${DIR}${PREFIX}.fastq \
		-b $meth_bam -g $ref > ${DIR}${PREFIX}_methylation_calls.tsv

	echo "Tabulating methylation frequencies..."
	# Create separate file for CpG
	awk 'BEGIN{FS="\t";OFS="\t"} NR == 1 || $11 == "CG" {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $12}' \
		${DIR}${PREFIX}_methylation_calls.tsv > ${DIR}${PREFIX}_cpg_methylation_calls.tsv
	echo -e "Running command:\n\
	[$np_frequency_cpg -i ${DIR}${PREFIX}_cpg_methylation_calls.tsv > ${DIR}${PREFIX}_cpg_methylation.bed]"
	$np_frequency_cpg -i ${DIR}${PREFIX}_cpg_methylation_calls.tsv > \
		${DIR}${PREFIX}_cpg_methylation.bed
	awk '{total += $5}; END {print "Average CpG Methylation: " total / NR}' \
		${DIR}${PREFIX}_cpg_methylation.bed

	# Create separate file for GpC
	awk 'BEGIN{FS="\t";OFS="\t"} NR == 1 || $11 == "GC" {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $12}' \
		${DIR}${PREFIX}_methylation_calls.tsv > ${DIR}${PREFIX}_gpc_methylation_calls.tsv
	echo -e "Running command:\n\
	[$np_frequency_gpc -i ${DIR}${PREFIX}_gpc_methylation_calls.tsv > ${DIR}${PREFIX}_gpc_methylation.bed]"
	$np_frequency_gpc -i ${DIR}${PREFIX}_gpc_methylation_calls.tsv > \
		${DIR}${PREFIX}_gpc_methylation.bed
	awk '{total += $5}; END {print "Average GpC Methylation: " total / NR}' \
		${DIR}${PREFIX}_gpc_methylation.bed
fi

echo "*** Analysis Finished ***"