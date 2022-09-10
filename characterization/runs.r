# This code examines open and closed runs.
# It visualizes linker runs as well as the distribution
# of open and closed run lengths.

# libraries and functions
library(Gviz)
library(org.Hs.eg.db)
library(rtracklayer)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
rdata.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/200913_GM12878/rdata'

dnase.1.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/ENCODE/gm12878_dnase_hg38_rep1_ENCFF598KWZ.bed'
dnase.2.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/ENCODE/gm12878_dnase_hg38_rep2_ENCFF073ORT.bed'

txdb.path <- '/seq/epiprod02/kdong/references/txdb/hg38.knownGene'
ideo.path <- '/seq/epiprod02/kdong/references/gviz/cytoBandIdeo.txt'

# Load in data
loci <- import.bed(loci.path)

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
						  qValue = "numeric", peak = "integer")
dnase.1 <- import(dnase.1.path, extraCols=extraCols_narrowPeak)
dnase.2 <- import(dnase.1.path, extraCols=extraCols_narrowPeak)
dnase.gm <- mergePeaks(dnase.1, dnase.2)
dnase.gm$idx <- seq(length(dnase.gm))

idx <- which(loci$name == 'IL2RA_RBM17')
locus <- loci[idx]

# Define Gviz plotting region and tracks
chr <- as.character(seqnames(locus))
from <- start(locus)
to <- end(locus)

# Ideogram
bands <- read.table(ideo.path, sep='\t', header=F, stringsAsFactors=F, 
					col.names=c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
itrack <- IdeogramTrack(genome='hg38', bands=bands, chromosome=chr, cex=1)

# Genome axis
gtrack <- GenomeAxisTrack(add53=T, add35=T, littleTicks=F, labelPos='below', cex=1)

# Genes
txdb <- loadDb(txdb.path)
txTr <- GeneRegionTrack(txdb, chromosome=chr, start=from, end=to, cex.title=1, 
						rot.title=0, name='Genes', cex.group=1.2, 
						collapseTranscripts ="meta", transcriptAnnotation='symbol')

# Get correct gene symbols
if(length(gene(txTr)) > 0){
	symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
	symbol(txTr) <- symbols[gene(txTr)]
}


# Plot Linker track
pdf('GM12878_linker_track.pdf', width=12, height=8)
lapply(idx, function(i){
	locus <- loci[i]
	message(locus$name)	
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	meanTrack1 <- DataTrack(range=tiles, data=colMeans(runMat.open, na.rm=T), type='polygon', 
							name='Non-linker', cex.title=1, cex.axis=1, ylim=c(0,1))
	meanTrack2 <- DataTrack(range=tiles, data=colMeans(runMat.linker, na.rm=T), type='polygon', 
							name='Linker', cex.title=1, cex.axis=1, ylim=c(0,1))
	peaks <- AnnotationTrack(dnase.gm[dnase.gm %over% locus], chromosome=chr, start=from, 
							 end=to, name='DNase', rot.title=0, cex.title=1, stacking='dense')

	plotTracks(c(meanTrack1, meanTrack2, peaks, txTr, gtrack, itrack), type='polygon', 
			   chromosome=chr, from=from, to=to, collapseTranscripts ="meta", 
			   transcriptAnnotation='symbol', title.width=1)
	
	return(NULL)
})
dev.off()

# Make histograms of length distributions
all_loci.list <- lapply(seq_along(loci), function(i){
	locus <- loci[i]
	message(locus$name)	
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))

	return(return(do.call(rbind,open_runs)))
})
open.df <- do.call(rbind, all_loci.list)

all_loci.list <- lapply(seq_along(loci), function(i){
	locus <- loci[i]
	message(locus$name)	
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))

	return(return(do.call(rbind, closed_runs)))
})
closed.df <- do.call(rbind, all_loci.list)

pdf('GM12878_all_run_lengths.pdf', width=9, height=6)
hist(open.df$dist, main='Distribution of Open Run Lengths', xlab='Length', ylab='Count', 
	 breaks=1000, xlim=c(0,1000), col=rgb(27, 158, 119, maxColorValue=255))
hist(closed.df$dist, main='Distribution of Closed Run Lengths', xlab='Length', ylab='Count', 
	 breaks=1000, xlim=c(0,1000), col='darkred')
# Actual Peaks
# abline(v=c(25,130,310,500), col='red', lty='dashed', lwd=2)
# abline(v=147*2 + 45, col='red', lty='dashed', lwd=2)
# abline(v=147*3 + 45*2, col='red', lty='dashed', lwd=2)
# # Expected nucleosomes
# abline(v=147, col='red', lty='dashed', lwd=2)
# abline(v=147*2 + 45, col='red', lty='dashed', lwd=2)
# abline(v=147*3 + 45*2, col='red', lty='dashed', lwd=2)
dev.off()

# Histograms of 3 methods
pdf('GM12878_closed_run_min_mid_max.pdf', width=9, height=9)
par(mfrow=c(3,1))

dens <- density(closed.df$min, n=1024)
hist(closed.df$min, main='Distribution of Closed Run Lengths (Minimum)', xlab='Length', 
	 ylab='Density', freq=F, breaks=seq(0,12000, by=10), xlim=c(0,1000), col='darkred')
lines(dens)
abline(v=dens$x[dens$x > 0 & dens$x < 80][which.max(dens$y[dens$x > 0 & dens$x < 80])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 80 & dens$x < 200][which.max(dens$y[dens$x > 80 & dens$x < 200])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 200 & dens$x < 400][which.max(dens$y[dens$x > 200 & dens$x < 400])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 400 & dens$x < 600][which.max(dens$y[dens$x > 400 & dens$x < 600])], 
	   col='red', lty='dashed', lwd=2)

dens <- density(closed.df$dist, n=1024)
hist(closed.df$dist, main='Distribution of Closed Run Lengths (Midpoint)', xlab='Length', 
	 ylab='Density', freq=F, breaks=seq(0,12000, by=10), xlim=c(0,1000), col='darkred')
lines(dens)
abline(v=dens$x[dens$x > 0 & dens$x < 80][which.max(dens$y[dens$x > 0 & dens$x < 80])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 80 & dens$x < 200][which.max(dens$y[dens$x > 80 & dens$x < 200])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 200 & dens$x < 400][which.max(dens$y[dens$x > 200 & dens$x < 400])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 400 & dens$x < 600][which.max(dens$y[dens$x > 400 & dens$x < 600])], 
	   col='red', lty='dashed', lwd=2)

dens <- density(closed.df$max, n=1024)
hist(closed.df$max, main='Distribution of Closed Run Lengths (Maximum)', xlab='Length', 
	 ylab='Density', freq=F, breaks=seq(0,12000, by=10), xlim=c(0,1000), col='darkred')
lines(dens)
abline(v=dens$x[dens$x > 0 & dens$x < 80][which.max(dens$y[dens$x > 0 & dens$x < 80])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 80 & dens$x < 250][which.max(dens$y[dens$x > 80 & dens$x < 250])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 250 & dens$x < 450][which.max(dens$y[dens$x > 250 & dens$x < 450])], 
	   col='red', lty='dashed', lwd=2)
abline(v=dens$x[dens$x > 430 & dens$x < 600][which.max(dens$y[dens$x > 450 & dens$x < 600])], 
	   col='red', lty='dashed', lwd=2)
dev.off()