#!/bin/bash

#==============================================================================
# Code for repairing hg38 assembly
# Replaces H19/IGF2 locus with new h19 assembly from NECAT 
# Use NECAT to create de novo assembly
#==============================================================================

### Dotkits ###
source /broad/software/scripts/useuse
use .samtools-1.7
use BWA
###############

### Paths ###
hg38="/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/hg38.chrXYM.fa"
#############

# Create new directory
mkdir -p hg38_repair
cd hg38_repair

# Split original genome by chromosome
csplit -s -z ${hg38} '/>/' '{*}'

# Split chromosome 11 (exact coordinates determined by manual inspection)
mv xx10 xx10_old

samtools faidx xx10_old
samtools faidx xx10_old "chr11:1-1768505" > chr11_1.fa
samtools faidx xx10_old "chr11:2164321-135086622" > chr11_3.fa

# # Insert h19 assembly into chromosome 11
# cat chr11_1.fa ../h19/6-bridge_contigs/polished_contigs.fasta chr11_3.fa | sed -e '1!{/^>.*/d;}' | \
# 	sed ':a;N;$!ba;s/\n//2g' | sed '1!s/.\{60\}/&\n/g' > xx10

# Insert h19 assembly (reversed) into chromosome 11
cat chr11_1.fa ../h19/reverse_complement.fasta chr11_3.fa | sed -e '1!{/^>.*/d;}' | \
 	sed ':a;N;$!ba;s/\n//2g' | sed '1!s/.\{60\}/&\n/g' > xx10

# Manually change name inside xx10 to "chr11" before concatenation
cat xx*[0-9] > hg38.chrXYM.repaired.fa

# Generate indices and chrom sizes
samtools faidx hg38.chrXYM.repaired.fa
cut -f1,2 hg38.chrXYM.repaired.fa.fai > hg38.chrXYM.repaired.chrom.sizes

bwa index hg38.chrXYM.repaired.fa