#!/bin/bash

# 3_r_analysis.sh generates rdata files with methylation data
# User must provide the following preceded by the appropriate flag:
# Working directory (-wd) with the QSUB command
# Optionally provide a prefix for file names (-n)
# Provide a locus file (-L)

## qsub arguments
#$ -l h_vmem=16G
#$ -l h_rt=6:00:00
#$ -l os=RedHat7
#$ -r y
#$ -j y
#$ -b y
#$ -m e
#$ -o r_analysis.log
#$ -e r_analysis.error

# Dotkits (RedHat7)
source /broad/software/scripts/useuse
reuse -q .r-3.6.0-bioconductor

# Infer pipeline directory
SCRIPT_DIR=$(readlink -f $(dirname "$0"))

# Additional software/paths
analysis="${SCRIPT_DIR}/scripts/r_analysis.R"
np_frequency_cpg="${SCRIPT_DIR}/scripts/calculate_methylation_frequency_cpg_reads.py"
np_frequency_gpc="${SCRIPT_DIR}/scripts/calculate_methylation_frequency_gpc_reads.py"

# Defaults
PREFIX="output"
loci=""

# Parse inputs and store parmeters in appropriate variables
while [ $# -gt 0 ]
do
	case $1 in
		-n)
			shift
			PREFIX=$1
			shift
			;;
		-L)
			shift
			loci=$(readlink -f $1)
			shift
			;;
		*)
			echo "ERROR: Unknown flag"
			exit 1
			;;
	esac
done

# Check working directory
if [ -d ./workspace/ ]; then
	cd ./workspace
fi

# Check for valid inputs
if [ ! -f "$loci" ]; then
	echo "ERROR: Missing loci BED file"
	exit 2
fi

D=`date "+%Y%m%d"`
echo "*** Analysis Started ***"
echo "* Date: ${D}"

mkdir -p rdata

echo "Generating BEDs..."

echo -e "Running command:\n\
	[$np_frequency_cpg -s -i ${PREFIX}_cpg_methylation_calls.tsv > ${PREFIX}_cpg_methylation_reads.bed]"
$np_frequency_cpg -s -i ${PREFIX}_cpg_methylation_calls.tsv > \
	${PREFIX}_cpg_methylation_reads.bed
echo -e "Running command:\n\
	[$np_frequency_gpc -s -i ${PREFIX}_gpc_methylation_calls.tsv > ${PREFIX}_gpc_methylation_reads.bed]"
$np_frequency_gpc -s -i ${PREFIX}_gpc_methylation_calls.tsv > \
	${PREFIX}_gpc_methylation_reads.bed

echo "Running R script..."
echo -e "Running command:\n\
	[Rscript $analysis $loci ${PREFIX}_cpg_methylation_reads.bed ${PREFIX}_gpc_methylation_reads.bed]"
Rscript $analysis $loci ${PREFIX}_cpg_methylation_reads.bed ${PREFIX}_gpc_methylation_reads.bed

echo "*** Analysis Finished ***"