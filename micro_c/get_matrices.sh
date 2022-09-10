#!/bin/bash

#==============================================================================
# Code to extract matrices from Micro-C coolers
#==============================================================================

### Dotkits ###
source /broad/software/scripts/useuse
use Anaconda3
###############

### Paths ###
repair_cool="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/microC_merged/results_repair/coolers_library_group/mergedLib.hg38.mapq_30.1000.mcool"
chm13_cool="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/microC_merged/results_chm13/coolers_library_group/mergedLib.hg38.mapq_30.1000.mcool"
hg38_cool="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/microC/4DNFI9GMP2J8.mcool"

repair_coord="chr11:1768505-2167285"
chm13_coord="chr11:1854166-2253688"
hg38_coord="chr11:1768505-2164321"

repair_locus="chr11:1774757-1987668"
chm13_locus="chr11:2034394-2247771"

repair_out="coolMatrix_repaired_merged_all.txt"
repair_out_locus="coolMatrix_repaired_merged.txt"
chm13_out="coolMatrix_chm13_merged_all.txt"
chm13_out_locus="coolMatrix_chm13_merged.txt"
hg38_out="coolMatrix_paper_all.txt"

cooltools="/seq/epiprod02/kdong/environments/cooltools"
#############

conda activate cooltools

# Extract entire reassembled region
python cooltools_extract.py --input ${repair_cool} --output ${repair_out} \
	--region ${repair_coord} --flip
python cooltools_extract.py --input ${chm13_cool} --output ${chm13_out} \
	--region ${chm13_coord}
python cooltools_extract.py --input ${hg38_cool} --output ${hg38_out} \
	--region ${hg38_coord}

# Extract just the H19/IGF2 region
python cooltools_extract.py --input ${repair_cool} --output ${repair_out_locus} \
	--region ${repair_locus} --flip
python cooltools_extract.py --input ${chm13_cool} --output ${chm13_out_locus} \
	--region ${chm13_locus}
