# This code examines the run signal around CTCF motifs.
# It uses GM12878 data and examines CTCF motifs from HOMER.
# CTCF ChIP-seq data from ENCODE is used to verify binding sites.

# libraries and functions
library(RColorBrewer)
library(rtracklayer)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
rdata.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/200913_GM12878/rdata'

motifs.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/gm_k562_motifs.RDS'
ctcf.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/ENCODE/gm12878_ctcf_hg38_idr_ENCFF796WRU.bed'

# Load in data
loci <- import.bed(loci.path)
motifs <- readRDS(motifs.path)

# CTCF
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
						  qValue = "numeric", peak = "integer")
ctcf <- import.bed(ctcf.path, extraCols=extraCols_narrowPeak)

over  <- motifs[motifs %over% loci & grepl('\\<ctcf', motifs$name, ignore.case=T) 
				& motifs %over% ctcf]
message(sprintf('%s CTCF motifs overlapping ChIP-Seq peak', length(over)))
out <- motifs[motifs %over% loci & grepl('\\<ctcf', motifs$name, ignore.case=T) 
			  & motifs %outside% ctcf]
message(sprintf('%s CTCF motifs not overlapping ChIP-Seq peak', length(out)))

# Get plot matrices and scores for motifs overlapping ChIP peaks
over.list <- lapply(seq_along(over), function(i, roi.width=1000){
	motif <- over[i]
	locus <- loci[loci %over% motif]
	message(locus$name)
	roi <- resize(motif, width=roi.width, fix='center')
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))
	
	# Filter reads by open run quality and coverage
	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% roi]))
	select <- readCov >= 0.7 & open_run.qc < 2
	
	# Calculate short closed run score over motif
	peak <- resize(roi, width=50, fix='center')
	score <- mergeMat(runMat.short[select,], tiles, peak)
	rev_idx <- which(grepl('M-', rownames(meMat.gc.f[select,])))
	score[rev_idx] <- -score[rev_idx]

	# Calculate CpG score over motif
	cgIdx <- match(rownames(meMat.gc.f[select,]), rownames(meMat.cg.f))
	cpg_score <- mergeMat(meMat.cg.f[cgIdx,], CpG_fil, peak)

	# Calculate open run score
	open_score <- rowMeans(runMat.short[select,tiles %over% peak] | 
										runMat.open[select,tiles %over% peak], na.rm=T)

	# Get plot matrix
	results <- plotReadMat(meMat.gc.f, GpC_fil, meMat.cg.f, CpG_fil, open_runs,
						   closed_runs, subset=select, gap=5, region=roi, open=T, 
						   linker=F, short=T, cpg=F, peaks=NULL, plot=F)
	return(list(results, rownames(meMat.gc.f)[select], score, cpg_score, open_score))
}, roi.width=1000)

# Get plot matrices and scores for motifs not overlapping ChIP peaks
out.list <- lapply(seq_along(out), function(i, roi.width=1000){
	motif <- out[i]
	locus <- loci[loci %over% motif]
	message(locus$name)
	roi <- resize(motif, width=roi.width, fix='center')
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))
	
	# Filter reads by open run quality and coverage
	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% roi]))
	select <- readCov >= 0.7 & open_run.qc < 2
	
	# Calculate short closed run score over motif
	peak <- resize(roi, width=50, fix='center')
	score <- mergeMat(runMat.short[select,], tiles, peak)
	rev_idx <- which(grepl('M-', rownames(meMat.gc.f[select,])))
	score[rev_idx] <- -score[rev_idx]

	# Calculate CpG score over motif
	cgIdx <- match(rownames(meMat.gc.f[select,]), rownames(meMat.cg.f))
	cpg_score <- mergeMat(meMat.cg.f[cgIdx,], CpG_fil, peak)

	# Calculate open run score
	open_score <- rowMeans(runMat.short[select,tiles %over% peak] | 
										runMat.open[select,tiles %over% peak], na.rm=T)

	# Get plot matrix
	results <- plotReadMat(meMat.gc.f, GpC_fil, meMat.cg.f, CpG_fil, open_runs,
						   closed_runs, subset=select, gap=5, region=roi, open=T, 
						   linker=F, short=T, cpg=F, peaks=NULL, plot=F)
	return(list(results, rownames(meMat.gc.f)[select], score, cpg_score, open_score))
}, roi.width=1000)

plot_ctcf <- function(list){
	n <- sum(sapply(list, function(el){length(el[[3]])}))
	message(sprintf('%s total reads', n))

	results.list <- lapply(list, function(el){
		el[[1]]
	})
	results2 <- do.call(rbind, results.list)

	scores.list <- lapply(list, function(el){
		return(el[[3]])
	})
	scores <- do.call(c, scores.list)

	cpg_scores.list <- lapply(list, function(el){
		return(el[[4]])
	})
	cpg_scores <- do.call(c, cpg_scores.list)

	open_scores.list <- lapply(list, function(el){
		return(el[[5]])
	})
	open_scores <- do.call(c, open_scores.list)
	
	# Assign colors to scores
	idx <- which(cpg_scores >= 0.5)
	cpg_col <- vector('character', length(cpg_scores))
	cpg_col[idx] <- '#cb181d'
	cpg_col[-idx] <- '#DCDCDC'

	col_fun <- colorRamp(brewer.pal(9, 'Greens'))
	open_col<- rgb(col_fun(open_scores), maxColorValue=255)

	col_fun <- colorRamp(c('#ffffff', '#f0eaf3', '#e2d5e7', '#d3c0db', '#c5abcf', '#b697c3', '#a884b7', '#9970ab'))
	score_col<- rgb(col_fun(abs(scores)), maxColorValue=255)

	# Reorder the reads in the results plot
	groups <- split(seq(nrow(results2)), ceiling(seq(nrow(results2))/4))
	order <- order(abs(scores), decreasing=T)
	mat_order <- do.call(c, groups[order])

	# # Test pattern
	# grid.newpage()

	# pushViewport(viewport(layout = grid.layout(3, 9, widths=c(1,15,0.25,0.75,0.25,0.75,0.25,0.75,1), 
	# 										   heights=c(1.5,17,1.5))))
	# for(i in seq(3)){
	# 	for(j in seq(9)){
	# 		pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	# 		grid.rect()
	# 		grid.circle()
	# 		popViewport()
	# 	}
	# }
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(3, 9, widths=c(1,15,0.25,0.75,0.25,0.75,0.25,0.75,1), 
											heights=c(1.5,17,1.5))))
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=2))
	grid.raster(results2[mat_order,], interpolate=F, width=unit(1,'npc'), height=unit(1,'npc'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=4))
	grid.rect()
	sideColorBar(cpg_col[order], rep=3, width=1)
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=6))
	grid.rect()
	sideColorBar(open_col[order], rep=3, width=1)
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=8))
	grid.rect()
	sideColorBar(score_col[order], rep=3, width=1)
	popViewport()
}

pdf('ctcf_ordered_reads_meth.pdf', width=20, height=12)
plot_ctcf(over.list)
plot_ctcf(out.list)
dev.off()

ctcf_on.list <- lapply(seq_along(over), function(i, roi.width){
	motif <- over[i]
	locus <- loci[loci %over% motif]
	message(locus$name)
	roi <- resize(motif, width=roi.width, fix='center')
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	# Filter reads by open run quality and coverage
	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% roi]))
	select <- readCov >= 0.7 & open_run.qc < 2

	# Calculate CpG score over motif
	cgIdx <- match(rownames(meMat.gc.f[select,]), rownames(meMat.cg.f))

	tile_idx <- tiles %over% roi
	tiles <- resize(tiles, width=25, fix='center')

	return(mergeMat(meMat.cg.f, CpG_fil, tiles[tile_idx]))
}, roi.width=1000)
cpg_ctcf_on <- do.call(rbind, ctcf_on.list)

ctcf_off.list <- lapply(seq_along(out), function(i, roi.width){
	motif <- out[i]
	locus <- loci[loci %over% motif]
	message(locus$name)
	roi <- resize(motif, width=roi.width, fix='center')
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	# Filter reads by open run quality and coverage
	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% roi]))
	select <- readCov >= 0.7 & open_run.qc < 2

	# Calculate CpG score over motif
	cgIdx <- match(rownames(meMat.gc.f[select,]), rownames(meMat.cg.f))

	tile_idx <- tiles %over% roi
	tiles <- resize(tiles, width=25, fix='center')

	return(mergeMat(meMat.cg.f, CpG_fil, tiles[tile_idx]))
}, roi.width=1000)
cpg_ctcf_off <- do.call(rbind, ctcf_off.list)

pdf('ctcf_meta_plot.pdf', width=18, height=6)
plot(colMeans(cpg_ctcf_on, na.rm=T), ylim=c(0,1), type='b', lwd=2, col='blue', xaxt='n', ylab='CpG Methylation', xlab='Distance from CTCF Motif')
axis(1, at=seq(0,200,50), labels=c(-500, -250, 0, 250, 500))
lines(colMeans(cpg_ctcf_off, na.rm=T), col='red', lwd=2, type='b')

legend(-5, 1, legend=c("No ChIP Peak", "ChIP Peak"),
	   col=c("red", "blue"), lwd=1, cex=0.8)
dev.off()