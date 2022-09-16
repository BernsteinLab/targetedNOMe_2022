# Call this script with the first and second haplotype names file (args[1] and args[2])
# Pass in the loci BED file which specify the enriched regions (args[3])
# Specify the mode (args[4]) [One of either "avg" or "smoo"]
# Specify smoothing window (args[5]), half the width will be added on each side

args <- commandArgs(trailingOnly=F)

# Infer scripts director
scripts.dir <- dirname(gsub('--file=', '', args[4]))
utils.dir <- paste0(scripts.dir, '/../../util/')

args <- commandArgs(trailingOnly=TRUE)

# Load libraries and functions
source(paste0(utils.dir, 'helper_functions.r'))
source(paste0(utils.dir, 'pipeline_functions.r'))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
genome <- BSgenome.Hsapiens.UCSC.hg38

# References
gpc_tiles.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_gpc_tiles.rds'
cpg_tiles.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_cpg_tiles.rds'

hap1.names <- read.table(args[1], stringsAsFactors=F)$V1
hap2.names <- read.table(args[2], stringsAsFactors=F)$V1

loci <- import.bed(args[3])

message('Loading CpG data...')
cpg.list <- loadCpG(loci, 'rdata', hap1.names, hap2.names)
cpgs <- do.call(c, cpg.list)

message('Loading GpC data...')
gpc.list <- loadGpC(loci, 'rdata', hap1.names, hap2.names)
gpcs <- do.call(c, gpc.list)

if(grepl('avg', args[4])){
	message('Averaging...')

    cpgs$score <- cpgs$hap1
    export.bedGraph(sort(cpgs[!is.na(cpgs$score)]),gsub('.txt', '_cpg.bedgraph', args[1]))
    message(sprintf("Output successfully written to %s", gsub('.txt', '_cpg.bedgraph', args[1])))

    cpgs$score <- cpgs$hap2
    export.bedGraph(sort(cpgs[!is.na(cpgs$score)]),gsub('.txt', '_cpg.bedgraph', args[2]))
    message(sprintf("Output successfully written to %s", gsub('.txt', '_cpg.bedgraph', args[2])))

    gpcs$score <- gpcs$hap1
    export.bedGraph(sort(gpcs[!is.na(gpcs$score)]),gsub('.txt', '_gpc.bedgraph', args[1]))
    message(sprintf("Output successfully written to %s", gsub('.txt', '_gpc.bedgraph', args[1])))

    gpcs$score <- gpcs$hap2
    export.bedGraph(sort(gpcs[!is.na(gpcs$score)]),gsub('.txt', '_gpc.bedgraph', args[2]))
    message(sprintf("Output successfully written to %s", gsub('.txt', '_gpc.bedgraph', args[2])))
}

if(grepl('Smooth', args[4])){
    message('Smoothing...')

	cpgs$score <- cpgs$hap1
	write_smooth_tdf(cpgs, args[5], cpg_tiles.path,
					 gsub('.txt', sprintf('_cpg_smooth_%s.bedgraph', args[5]), args[1]))

	cpgs$score <- cpgs$hap2
	write_smooth_tdf(cpgs, args[5], cpg_tiles.path,
					 gsub('.txt', sprintf('_cpg_smooth_%s.bedgraph', args[5]), args[2]))

	gpcs$score <- gpcs$hap1
	write_smooth_tdf(gpcs, args[5], gpc_tiles.path,
					 gsub('.txt', sprintf('_gpc_smooth_%s.bedgraph', args[5]), args[1]))

	gpcs$score <- gpcs$hap2
	write_smooth_tdf(gpcs, args[5], gpc_tiles.path,
					 gsub('.txt', sprintf('_gpc_smooth_%s.bedgraph', args[5]), args[2]))
}