#!/bin/bash
# 2_visualize.sh generates files for viewing the analyzed data in IGV
# Performs GCG exclusion and generates a TDF file of average methylation (--avg)
# Produces smoothed TDF file (--smoo) with optional window size (-w) [DEFAULT: 100]
# Modifies single reads for bisulfite mode viewing (--bsm)
# User must provide the following preceded by the appropriate flag:
# Working directory (-wd) with the QSUB command
# Optionally provide a prefix for file names (-n)

## qsub arguments
#$ -l h_vmem=4G
#$ -l h_rt=6:00:00
#$ -l os=RedHat7
#$ -pe smp 8
#$ -binding linear:8
#$ -r y
#$ -j y
#$ -b y
#$ -m e
#$ -o visualize.log
#$ -e visualize.error

# Dotkits (RedHat7)
source /broad/software/scripts/useuse
reuse -q .anaconda3-5.3.1
eval "$(conda shell.bash hook)"
reuse -q .samtools-1.9
reuse -q .htslib-1.8
reuse -q .r-3.6.0-bioconductor
reuse -q .igvtools-2.4.16

# Infer pipeline directory
SCRIPT_DIR=$(readlink -f $(dirname "$0"))

# Additional software/paths
ref="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/hg38.chrXYM.fa"
exclude_gcg="${SCRIPT_DIR}/scripts/exclude_GCG.R"
mtsv2bed="${SCRIPT_DIR}/scripts/mtsv2bed.py"
convert="${SCRIPT_DIR}/scripts/convertBam.py"
convert_env="/seq/epiprod02/kdong/environments/convert_bam"

# Defaults
PREFIX="output"
avg=false
bsm=false
con=false
smoo=false
window=100

# Parse inputs and store parmeters in appropriate variables
while [ $# -gt 0 ]
do
	case $1 in
		-n)
			shift
			PREFIX=$1
			shift
			;;
		-w)
			shift
			window=$1
			shift
			;;
		--avg)
			avg=true
			shift
			;;
		--smoo)
			smoo=true
			shift
			;;
		--bsm)
			bsm=true
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
mtsv_cpg="${PREFIX}_cpg_methylation_calls.tsv"
mtsv_gpc="${PREFIX}_gpc_methylation_calls.tsv"

if [ "$bsm" = true ]; then
	if [ ! -f "${PREFIX}_cpg_methylation.bed" ]; then
		echo "ERROR: Missing CpG methylation calls"
		exit 2
	fi
	if [ ! -f "${PREFIX}_gpc_methylation.bed" ]; then
		echo "ERROR: Missing GpC methylation calls"
		exit 3
	fi
fi

if [ "$avg" = true -o "$smoo" = true ]; then
	if [ ! -f "$mtsv_cpg" ]; then
		echo "ERROR: Missing CpG methylation BED file"
		exit 4
	fi
	if [ ! -f "$mtsv_gpc" ]; then
		echo "ERROR: Missing GpC methylation BED file"
		exit 5
	fi
fi

D=`date "+%Y%m%d"`

echo "*** Analysis Started ***"
echo "* Date: ${D}"

if [ "$avg" = true -o "$smoo" = true ]; then
	echo "Generating TDFs..."

	flag=""
	if [ "$avg" = true ]; then
		flag="${flag}avg"
	fi
	if [ "$smoo" = true ]; then
		flag="${flag}Smooth"
	fi

	echo -e "Running command:\n\
	[Rscript $exclude_gcg ${PREFIX}_cpg_methylation.bed CG $flag $window]"
	Rscript $exclude_gcg ${PREFIX}_cpg_methylation.bed CG $flag $window
	echo -e "Running command:\n\
	[Rscript $exclude_gcg ${PREFIX}_gpc_methylation.bed GC $flag $window]"
	Rscript $exclude_gcg ${PREFIX}_gpc_methylation.bed GC $flag $window

	makeTDF () {
		local filename=$1
		local igv=${filename/.bedgraph/.igv}
		local tdf=${igv/.igv/.tdf}

		echo -e 'Chromosome\tStart\tEnd\tFeature\tMethylation' > $igv
		awk '{FS="\t"; OFS="\t"}; {print $1, $2, $3, ".", $4}' $filename >> $igv

		igvtools toTDF $igv $tdf hg38 && rm $igv
	}

	if [ "$avg" = true ]; then
		echo "Making averaged TDFs..."
		echo -e "Running command:\n\
	[makeTDF ${PREFIX}_cpg_methylation.bedgraph && rm ${PREFIX}_cpg_methylation.bedgraph]"
		makeTDF ${PREFIX}_cpg_methylation.bedgraph && rm ${PREFIX}_cpg_methylation.bedgraph
		
		echo -e "Running command:\n\
	[makeTDF ${PREFIX}_gpc_methylation.bedgraph && rm ${PREFIX}_gpc_methylation.bedgraph]"
		makeTDF ${PREFIX}_gpc_methylation.bedgraph && rm ${PREFIX}_gpc_methylation.bedgraph
	fi
	if [ "$smoo" = true ]; then
		echo "Making smoothed TDFs..."
		echo -e "Running command:\n\
	[makeTDF ${PREFIX}_cpg_methylation_smooth.bedgraph && rm ${PREFIX}_cpg_methylation_smooth.bedgraph]"
		makeTDF ${PREFIX}_cpg_methylation_smooth.bedgraph \
			&& rm ${PREFIX}_cpg_methylation_smooth.bedgraph
		
		echo -e "Running command:\n\
	[makeTDF ${PREFIX}_gpc_methylation_smooth.bedgraph && rm ${PREFIX}_gpc_methylation_smooth.bedgraph]"
		makeTDF ${PREFIX}_gpc_methylation_smooth.bedgraph \
			&& rm ${PREFIX}_gpc_methylation_smooth.bedgraph
	fi
fi

if [ "$bsm" = true ]; then
	echo "Converting BAM..."
	conda activate $convert_env

	# Set file names
	mbed_cpg=${mtsv_cpg/.tsv/.bed.gz}
	mbed_gpc=${mtsv_gpc/.tsv/.bed.gz}
	bam="${PREFIX}.filtered.bam"

	# Check if filter step has run previously
	if [ ! -f "$bam" ]; then
		bam = ${PREFIX}.bam
	fi

	echo "Reading CpG data..."
	echo -e "Running command:\n\
	[python $mtsv2bed -m cpg -c 1.5 -i $mtsv_cpg | sort -k1,1 -k2,2n | bgzip > $mbed_cpg && tabix -p bed $mbed_cpg]"
	python $mtsv2bed -m cpg -c 1.5 -i $mtsv_cpg | sort -k1,1 -k2,2n | \
		bgzip > $mbed_cpg && tabix -p bed $mbed_cpg
	echo "Reading GpC data..."
	echo -e "Running command:\n\
	[python $mtsv2bed -m gpc -c 1 -i $mtsv_gpc | sort -k1,1 -k2,2n | bgzip > $mbed_gpc && tabix -p bed $mbed_gpc]"
	python $mtsv2bed -m gpc -c 1 -i $mtsv_gpc | sort -k1,1 -k2,2n | \
		bgzip > $mbed_gpc && tabix -p bed $mbed_gpc
	echo "Modifying Bam..."
	echo -e "Running command:\n\
	[python -u $convert -t 8 -f $ref -b $bam -c $mbed_cpg -g $mbed_gpc | samtools sort -o ${PREFIX}.converted.bam && samtools index ${PREFIX}.converted.bam]"
	python -u $convert -t 8 -f $ref -b $bam -c $mbed_cpg -g $mbed_gpc | \
		samtools sort -o ${PREFIX}.converted.bam && \
		samtools index ${PREFIX}.converted.bam
    conda deactivate
fi

echo "*** Analysis Finished ***"