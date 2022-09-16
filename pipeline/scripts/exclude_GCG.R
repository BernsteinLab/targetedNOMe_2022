# Call this script with BED filename (args[1]) followed with methylation type (args[2])
# Available outputs determined by flag (args[3]): avg, Smooth, or avgSmooth
# Specify smoothing window (args[4]), half the width will be added on each side

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
gpc.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_gpc.rds'
cpg.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_cpg.rds'
gcg.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/gcg_strands.rds'
gpc_tiles.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_gpc_tiles.rds'
cpg_tiles.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_cpg_tiles.rds'

gcg <- readRDS(gcg.path)

message('Loading methylation BED file...')
gr <- import.bed(args[1], extraCols=c(total='integer', meth='integer'),
				 colnames=c('chrom', 'start', 'end', 'strand', 'total', 'meth'))
gr$score <- gr$meth / gr$total

# Tally GCG hits and CG/GC hits
if (args[2] == 'CG'){
	cpg <- readRDS(cpg.path)
	cpg_rev <- shift(cpg, 1)
	strand(cpg) <- '+'
	strand(cpg_rev) <- '-'
	
	ref <- c(cpg, cpg_rev)
	rm(cpg, cpg_rev)
} else if (args[2] == 'GC'){
	gpc <- readRDS(gpc.path)
	gpc <- keepSeqlevels(gpc, seqlevels(gr), pruning='coarse')
	gpc_rev <- shift(gpc, -1)
	strand(gpc) <- '+'
	strand(gpc_rev) <- '-'
	
	ref <- c(gpc, gpc_rev)
	rm(gpc, gpc_rev)
}

# Resize GCG to width of 3
message('Removing groups with > 50% GCG...')
gcg_3 <- resize(gcg, width=3, fix='center')

gr$tot <- tallyHits(gr, ref)
rm(ref)

gr$gcg <- tallyHits(gr, gcg_3)
gr$gcg[is.na(gr$gcg)] <- 0

# Exclude bad groups with > 50% GCG, then weave
message('Combining strands...')
gr.c <- weave(gr[gr$gcg/gr$tot < 0.5], args[2], c('total', 'meth'))
gr.c$score <- gr.c$meth / gr.c$total

if(grepl('avg', args[3])){
	message('Averaging...')

	export.bedGraph(sort(gr.c), paste0(args[1], 'graph'))
	message(sprintf("Output successfully written to %s", paste0(args[1], 'graph')))
}

# Smoothing
if(grepl('Smooth', args[3])){
	message('Smoothing...')

	write_smooth_tdf(gr.c, args[4], ifelse(args[2]=='CG', cpg_tiles.path, gpc_tiles.path),
					 gsub('.bed', '_smooth.bedgraph', args[1]))
}
