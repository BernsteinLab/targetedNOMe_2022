# Call this script with the first and second haplotype names file (args[1] and args[2])
# Pass in the loci BED file which specify the enriched regions (args[3])
# Specify smoothing window (args[4]), half the width will be added on each side
# Optionally pass the genome (args[5]). Defaults to hg38 if left empty

args <- commandArgs(trailingOnly=F)

# Infer scripts director
scripts.dir <- dirname(gsub('--file=', '', args[4]))
utils.dir <- paste0(scripts.dir, '/../../util/')

args <- commandArgs(trailingOnly=TRUE)

# Load libraries and functions
source(paste0(utils.dir, 'helper_functions.r'))
source(paste0(utils.dir, 'pipeline_functions.r'))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
library(rtracklayer)
genome <- BSgenome.Hsapiens.UCSC.hg38 

hap1.names <- read.table(args[1], stringsAsFactors=F)$V1
hap2.names <- read.table(args[2], stringsAsFactors=F)$V1

loci <- import.bed(args[3])

# Identify haplotypes
# Extract correct rows from runMat
# Write the scores to same buffer, different metadata fields
# Create final buffer with correct metadata for writing out

# Generate BigWig Tracks
message("Generating BigWigs...")

prefix1 <- gsub('.txt', '', args[1])
prefix2 <- gsub('.txt', '', args[2])

write_bw(loci, 'open', prefix1, smooth=F, hap1.names)
write_bw(loci, 'open', prefix2, smooth=F, hap2.names)
write_bw(loci, 'open', prefix1,  smooth=T, hap1.names)
write_bw(loci, 'open', prefix2,  smooth=T, hap2.names)
write_bw(loci, 'short', prefix1,  smooth=F, hap1.names)
write_bw(loci, 'short', prefix2,  smooth=F, hap2.names)