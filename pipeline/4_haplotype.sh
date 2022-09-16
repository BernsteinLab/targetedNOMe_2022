#!/bin/bash
# haplotype.sh calls variants and tags the reads of a bam file
# Calls variants with Nanopolish (--var)
# Modifies single reads for bisulfite mode viewing (--bsm)
# Generates a TDF file of average methylation (--avg)
# Produces smoothed TDF file (--smoo) with optional window size (-w) [DEFAULT: 100]
# Generates run-based bigwigs for each haplotype (--bw)
# User must provide the following preceded by the appropriate flag:
# Working directory (-wd) with the QSUB command
# Provide a file with loci of interest (-L)
# Optionally provide a prefix for file names (-n)

## qsub arguments
#$ -l h_vmem=4G
#$ -l h_rt=4:00:00
#$ -l os=RedHat7
#$ -pe smp 8
#$ -binding linear:8
#$ -r y
#$ -j y
#$ -b y
#$ -m e
#$ -o haplotype.log
#$ -e haplotype.error

# Dotkits (RedHat7)
source /broad/software/scripts/useuse
source /broad/software/free/Linux/redhat_7_x86_64/pkgs/anaconda3_5.3.1/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
reuse -q Samtools
reuse -q .anaconda3-5.3.1
reuse -q .htslib-1.8
reuse -q GCC-5.2
reuse -q Bcftools
reuse -q .r-3.6.0-bioconductor
reuse -q .igvtools-2.4.16

# Infer pipeline directory
SCRIPT_DIR=$(readlink -f $(dirname "$0"))

# Additional software/paths
ref="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/hg38.chrXYM.fa"
nanopolish="/seq/epiprod02/kdong/software/nanopolish-13.2/nanopolish"
whatshap_env="/seq/epiprod02/kdong/environments/whatshap"
convert="${SCRIPT_DIR}/scripts/convertBam.py"
convert_env="/seq/epiprod02/kdong/environments/convert_bam"
tdf="${SCRIPT_DIR}/scripts/hap_tdf.R"
bigwig="${SCRIPT_DIR}/scripts/hap_bigwig.R"

# Defaults
PREFIX="output"
window=100
var=false
bsm=false
bw=false
loci=""
switched=false

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
		--avg)
			avg=true
			shift
			;;
		--smoo)
			smoo=true
			shift
			;;
		-w)
			shift
			window=$1
			shift
			;;
		--var)
			var=true
			shift
			;;
		--bsm)
			bsm=true
			shift
			;;
		--bw)
			bw=true
			shift
			;;
		*)
			echo "ERROR: Unknown flag"
			exit 1
			;;
	esac
done

# Check working directory
DIR=""
if [ -d ./workspace/ ]; then
	DIR="./workspace/"
fi

# Check for valid inputs
if [ "$var" = true -o "$avg" = true -o "$smoo" = true ]; then
	if [ ! -f "${loci}" ]; then
		echo "ERROR: Missing loci file"
		exit 2
	fi
fi

if [ "$bsm" = true ]; then
	mbed_cpg="${DIR}${PREFIX}_cpg_methylation_calls.bed.gz"
	mbed_gpc="${DIR}${PREFIX}_gpc_methylation_calls.bed.gz"
	if [ ! -f "${mbed_cpg}" ]; then
		echo "ERROR: Missing CpG methylation BED file (gzipped)."
		exit 3
	fi
	if [ ! -f "${mbed_gpc}" ]; then
		echo "ERROR: Missing GpC methylation BED file (gzipped)"
		exit 4
	fi
fi

if [ "$avg" = true -o "$smoo" = true -o "$bw" = true ]; then
	rds_hcg="${DIR}${PREFIX}_hcg_reduc.rds"
	rds_gch="${DIR}${PREFIX}_gch_reduc.rds"
	phased_bam="${DIR}${PREFIX}.phased.bam"

	if [ ! -f "${rds_hcg}" ]; then
		echo "ERROR: Missing HCG methylation RDS file. Please run r_analysis.sh first."
		exit 5
	fi
	if [ ! -f "${rds_gch}" ]; then
		echo "ERROR: Missing GCH methylation RDS file. Please run r_analysis.sh first."
		exit 6
	fi
	if [ ! -f "${phased_bam}" ]; then
		if [ "$var" != true ]; then
			echo "ERROR: Missing phased BAM file. Please rerun this script with --var."
			exit 7
		fi
	fi
fi

D=`date "+%Y%m%d"`
echo "*** Analysis Started ***"
echo "* Date: ${D}"

if [ "$var" = true ]; then

	bam=${DIR}${PREFIX}.filtered.bam
	# Check if filter step has run previously
	if [ ! -f "$bam" ]; then
		bam = ${DIR}${PREFIX}.reduced.bam
		echo "Filtering BAM..."
		echo -e "Running command:\n\
		[samtools view -bh -L $loci ${DIR}${PREFIX}.bam > $bam && samtools index $bam]"
		samtools view -bh -L $loci ${DIR}${PREFIX}.bam > $bam && samtools index $bam
	fi

	# Create combined VCF file
	# VCF Header
	echo '##fileformat=VCFv4.2' > ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=TotalReads,Number=1,Type=Integer,Description="The number of event-space reads used to call the variant">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=SupportFraction,Number=1,Type=Float,Description="The fraction of event-space reads that support the variant">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=SupportFractionByStrand,Number=2,Type=Float,Description="Fraction of event-space reads that support the variant for each strand">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=BaseCalledReadsWithVariant,Number=1,Type=Integer,Description="The number of base-space reads that support the variant">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=BaseCalledFraction,Number=1,Type=Float,Description="The fraction of base-space reads that support the variant">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=AlleleCount,Number=1,Type=Integer,Description="The inferred number of copies of the allele">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=StrandSupport,Number=4,Type=Integer,Description="Number of reads supporting the REF and ALT allele, by strand">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=StrandFisherTest,Number=1,Type=Integer,Description="Strand bias fisher test">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=SOR,Number=1,Type=Float,Description="StrandOddsRatio test from GATK">' >> ${DIR}${PREFIX}.vcf
	echo '##INFO=<ID=RefContext,Number=1,Type=String,Description="The reference sequence context surrounding the variant call">' >> ${DIR}${PREFIX}.vcf
	echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${DIR}${PREFIX}.vcf
	echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample' >> ${DIR}${PREFIX}.vcf

	echo "Calling Variants..."
	echo -e "Running command:\n\
		[cat $loci | awk 'BEGIN{FS="\t"}; {print $1":"$2"-"$3}' | xargs -I {} sh -c \"$nanopolish variants -x 5000 --snps -w {} -p 2 -q cpg,gpc -t 8 \
-r ${DIR}${PREFIX}.fastq -b $bam -g $ref | grep -v '#'\" | sort -k1,1 -k2,2n >> ${DIR}${PREFIX}.vcf]"
	cat $loci | awk 'BEGIN{FS="\t"}; {print $1":"$2"-"$3}' | xargs -I {} sh -c \
		"$nanopolish variants -x 5000 --snps -w {} -p 2 -q cpg,gpc -t 8 \
		-r ${DIR}${PREFIX}.fastq -b $bam -g $ref | grep -v '#'" | \
		sort -k1,1 -k2,2n >> ${DIR}${PREFIX}.vcf

	echo "Filtering Variants..."
	echo -e "Running command:\n\
		[bcftools view -i \"SOR<3\" ${DIR}${PREFIX}.vcf > ${DIR}${PREFIX}.filtered.vcf]"
	# Filter by SOR (Strand Odds Ratio) < 3. Can also use StrandFisherTest < 30 
	bcftools view -i "SOR<3" ${DIR}${PREFIX}.vcf > ${DIR}${PREFIX}.filtered.vcf

	# Tag BAM file with haplotypes
	conda activate $whatshap_env
	echo "Phasing Variants..."
	echo -e "Running command:\n\
		[whatshap phase --distrust-genotypes  --ignore-read-groups ${DIR}${PREFIX}.filtered.vcf $bam > ${DIR}${PREFIX}.phased.vcf]"
	whatshap phase --distrust-genotypes --ignore-read-groups ${DIR}${PREFIX}.filtered.vcf \
		$bam > ${DIR}${PREFIX}.phased.vcf

	echo -e "Running command:\n\
		[bgzip -c ${DIR}${PREFIX}.phased.vcf > ${DIR}${PREFIX}.phased.vcf.gz && tabix -p vcf ${DIR}${PREFIX}.phased.vcf.gz]"
	bgzip -c ${DIR}${PREFIX}.phased.vcf > ${DIR}${PREFIX}.phased.vcf.gz && \
		tabix -p vcf ${DIR}${PREFIX}.phased.vcf.gz

	echo "Tagging BAM reads..."
	echo -e "Running command:\n\
		[whatshap haplotag --ignore-read-groups ${DIR}${PREFIX}.phased.vcf.gz $bam > ${DIR}${PREFIX}.phased.bam && samtools index ${DIR}${PREFIX}.phased.bam]"
	whatshap haplotag --ignore-read-groups ${DIR}${PREFIX}.phased.vcf.gz $bam > \
		${DIR}${PREFIX}.phased.bam && samtools index ${DIR}${PREFIX}.phased.bam
	conda deactivate
fi

if [ "$bsm" = true ]; then
	# Create converted BAM after haplotype calling
	conda activate $convert_env

	echo "Converting phased BAM..."
	echo -e "Running command:\n\
	[python -u $convert -t 8 -f $ref -b ${DIR}${PREFIX}.phased.bam -c $mbed_cpg -g $mbed_gpc | samtools sort -o ${DIR}${PREFIX}.phased_converted.bam && samtools index ${DIR}${PREFIX}.phased_converted.bam]"
	python -u $convert -t 8 -f $ref -b ${DIR}${PREFIX}.phased.bam -c $mbed_cpg -g $mbed_gpc |\
		samtools sort -o ${DIR}${PREFIX}.phased_converted.bam \
		&& samtools index ${DIR}${PREFIX}.phased_converted.bam

	conda deactivate
fi

if [ "$avg" = true -o "$smoo" = true ]; then
	echo "Generating TDFs..."

	flag=""
	if [ "$avg" = true ]; then
		flag="${flag}avg"
	fi
	if [ "$smoo" = true ]; then
		flag="${flag}Smooth"
	fi

	# Read names
	samtools view -L $loci $phased_bam | grep HP:i:1 | cut -f 1 > ${DIR}${PREFIX}_hap1.txt
	samtools view -L $loci $phased_bam | grep HP:i:2 | cut -f 1 > ${DIR}${PREFIX}_hap2.txt

	if [ -d ./workspace/ ]; then
		cd ./workspace
		switched=true
	fi

	echo -e "Running command:\n\
	[Rscript $tdf ${PREFIX}_hap1.txt ${PREFIX}_hap2.txt $loci $flag $window]"
	Rscript $tdf ${PREFIX}_hap1.txt ${PREFIX}_hap2.txt $loci $flag $window

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
	
		makeTDF ${PREFIX}_hap1_cpg.bedgraph && rm ${PREFIX}_hap1_cpg.bedgraph
		makeTDF ${PREFIX}_hap2_cpg.bedgraph && rm ${PREFIX}_hap2_cpg.bedgraph
		makeTDF ${PREFIX}_hap1_gpc.bedgraph && rm ${PREFIX}_hap1_gpc.bedgraph
		makeTDF ${PREFIX}_hap2_gpc.bedgraph && rm ${PREFIX}_hap2_gpc.bedgraph
	fi

	if [ "$smoo" = true ]; then
		echo "Making averaged TDFs..."
	
		makeTDF ${PREFIX}_hap1_cpg_smooth_${window}.bedgraph \
			&& rm ${PREFIX}_hap1_cpg_smooth_${window}.bedgraph
		makeTDF ${PREFIX}_hap2_cpg_smooth_${window}.bedgraph \
			&& rm ${PREFIX}_hap2_cpg_smooth_${window}.bedgraph
		makeTDF ${PREFIX}_hap1_gpc_smooth_${window}.bedgraph \
			&& rm ${PREFIX}_hap1_gpc_smooth_${window}.bedgraph
		makeTDF ${PREFIX}_hap2_gpc_smooth_${window}.bedgraph \
			&& rm ${PREFIX}_hap2_gpc_smooth_${window}.bedgraph
	fi

	if [ "$switched" = true ]; then
		cd ..
		switched=false
	fi
fi

if [ "$bw" = true ]; then
	echo "Generating bigwigs..."

	# Read names
	samtools view -L $loci $phased_bam | grep HP:i:1 | cut -f 1 > ${DIR}${PREFIX}_hap1.txt
	samtools view -L $loci $phased_bam | grep HP:i:2 | cut -f 1 > ${DIR}${PREFIX}_hap2.txt

	if [ -d ./workspace/ ]; then
		cd ./workspace
		switched=true
	fi

	echo -e "Running command:\n\
	[Rscript $bigwig ${PREFIX}_hap1.txt ${PREFIX}_hap2.txt $loci $window]"
	Rscript $bigwig ${PREFIX}_hap1.txt ${PREFIX}_hap2.txt $loci $window

	if [ "$switched" = true ]; then
		cd ..
		switched=false
	fi
fi

echo "*** Analysis Finished ***"