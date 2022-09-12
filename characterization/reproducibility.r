# libraries and functions
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
# H9
loci_h9_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/210603_loci.bed'
rdata_h9_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/rdata/'

loci_h9_2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/H9_ESC/191231_H9/191231_loci.bed'
rdata_h9_2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/H9_ESC/191231_H9/workspace/rdata/'

# HSMM 
loci_hsmm_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/210603_loci.bed'
rdata_hsmm_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/HSMM/workspace/rdata/'
loci_hsmm_2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/210603_loci.bed'
rdata_hsmm_2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/HSMM_rep2/workspace_1/rdata/'

# GM12878
loci_gm_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/190522_GM12878/190522_loci.bed.txt'
rdata_gm_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/190522_GM12878/workspace/rdata/'
loci_gm_2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
rdata_gm_2.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/200913_GM12878/rdata/'

# T-cell
loci_tcell_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_march2020/Tstim/donor73/200326_loci.bed'
rdata_tcell_1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_march2020/Tres/donor72/workspace/rdata/'
loci_tcell_2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/TcellNov2020_loci.bed'
rdata_tcell_2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/rdata/'

get_tile_averages <- function(loci.path, rdata.path, cov, min, type){
	loci <- import.bed(loci.path)
	chrs <- unique(as.character(runValue(seqnames(loci))))
	# 250 bp tiles of each relevant chromosome
	tiles <- tileGenome(seqlengths=seqinfo(BSgenome.Hsapiens.UCSC.hg38)[chrs], 
						 tilewidth=250, cut.last.tile.in.chrom=T)
	tiles <- tiles[tiles %within% loci]

	results.list <- lapply(seq_along(loci), function(i, rdata.path, tiles, min){
		locus <- loci[i]
		message(locus$name)
		load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
		
		if (type == 'CG'){
			idx <- which(tallyHits(tiles, CpG_fil) > min)
			binMat <- mergeMat(meMat.cg.f, CpG_fil, tiles[idx])
		} else if (type == 'GC'){
			idx <- which(tallyHits(tiles, GpC_fil) > min)
			binMat <- mergeMat(meMat.gc.f, GpC_fil, tiles[idx])
		}

		final <- tiles[idx]
		final$avg <- colMeans(binMat, na.rm=T)
		final$cov <- colSums(!is.na(binMat))
		return(final)
	}, rdata.path=rdata.path, tiles=tiles, min=min)

	results <- do.call(c, results.list)
	results <- results[results$cov > cov]
	return(results)
}

correlate <- function(loci1.path, rdata1.path, loci2.path, rdata2.path, title, cov, min=0){
	# Load CpG and GpC data
	cpgs.1 <- get_tile_averages(loci1.path, rdata1.path, cov, min, type='CG')
	gpcs.1 <- get_tile_averages(loci1.path, rdata1.path, cov, min, type='GC')

	cpgs.2 <- get_tile_averages(loci2.path, rdata2.path, cov, min, type='CG')
	gpcs.2 <- get_tile_averages(loci2.path, rdata2.path, cov, min, type='GC')

	hits1 <- findOverlaps(cpgs.1, cpgs.2)
	cor1 <- cor(cpgs.1$avg[queryHits(hits1)], cpgs.2$avg[subjectHits(hits1)], 
				use='complete.obs')
	color <- rgb(0, 0, 0, alpha=0.33)
	plot(cpgs.1$avg[queryHits(hits1)], cpgs.2$avg[subjectHits(hits1)], col=color, pch=19,
		 xlab=sprintf('%s Rep 1', title), ylab=sprintf('%s Rep 2', title), xlim=c(0,1), 
		 ylim=c(0,1), main=sprintf('CpG Reproducibility (r=%s)', round(cor1,3)))
	abline(a=0, b=1, lty=2, col='red', lwd=2)
	
	cpg <- data.frame(cpg1=cpgs.1$avg[queryHits(hits1)], cpg2=cpgs.2$avg[subjectHits(hits1)])

	hits2 <- findOverlaps(gpcs.1, gpcs.2)
	cor2 <- cor(gpcs.1$avg[queryHits(hits2)], gpcs.2$avg[subjectHits(hits2)], 
				use='complete.obs')
	plot(gpcs.1$avg[queryHits(hits2)], gpcs.2$avg[subjectHits(hits2)], col=color, pch=19,
		 xlab=sprintf('%s Rep 1', title), ylab=sprintf('%s Rep 2', title), xlim=c(0,1), 
		 ylim=c(0,1), main=sprintf('GpC Reproducibility (r=%s)', round(cor2,3)))
	abline(a=0, b=1, lty=2, col='red', lwd=2)
	gpc <- data.frame(gpc1=gpcs.1$avg[queryHits(hits2)], gpc2=gpcs.2$avg[subjectHits(hits2)])

	return(list(cpg, gpc))
}

# Thresholds
min <- 1
cov <- 20

pdf('reproducibility_plots.pdf', width=12, height=6)
par(mfrow=c(1,2))
h9.list <- correlate(loci_h9_1.path, rdata_h9_1.path, loci_h9_2.path, rdata_h9_2.path, 'H9', cov, min)
hsmm.list <- correlate(loci_hsmm_1.path, rdata_hsmm_1.path, loci_hsmm_2.path, rdata_hsmm_2.path, 'HSMM', cov, min)
gm.list <- correlate(loci_gm_1.path, rdata_gm_1.path, loci_gm_2.path, rdata_gm_2.path, 'GM12878', cov, min)
tcell.list <- correlate(loci_tcell_1.path, rdata_tcell_1.path, loci_tcell_2.path, rdata_tcell_2.path, 'T-cell', cov, min)
dev.off()
