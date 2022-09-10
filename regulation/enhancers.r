# This code examines open run signal at enhancer sites in T-cells.
# Enhancers are inferred based on open run peaks with
# the exclusion of transcription start sites. A heatmap
# is produced to visualize the changes across time points.

# libraries and functions
library(grid)
library(RColorBrewer)
library(rtracklayer)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/TcellNov2020_loci.bed'
rdata.0.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/rdata/'
rdata.24.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t24h/workspace/rdata/'
rdata.48.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/rdata/'

peaks.0.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/t0h_open_peaks.bed'
peaks.24.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t24h/workspace/t24h_open_peaks.bed'
peaks.48.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/t48h_open_peaks.bed'

tss_tcell.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/tss_tcell.rds'

# Load in data
loci <- import.bed(loci.path)

peaks.0 <- import.bed(peaks.0.path)
peaks.24 <- import.bed(peaks.24.path)
peaks.48 <- import.bed(peaks.48.path)

tss.l <- readRDS(tss_tcell.path)

gpc.list <- loadGpC(loci, rdata.0.path)
gpcs <- do.call(c, gpc.list)

# Filter and merge called peaks
peaks.0$tot <- tallyHits(peaks.0, gpcs)
peaks.24$tot <- tallyHits(peaks.24, gpcs)
peaks.48$tot <- tallyHits(peaks.48, gpcs)

peaks.0.f <- peaks.0[which(peaks.0$tot > 1)]
peaks.24.f <- peaks.24[which(peaks.24$tot > 1)]
peaks.48.f <- peaks.48[which(peaks.48$tot > 1)]

peaks.all <- unique(mergePeaks(mergePeaks(peaks.0.f, peaks.24.f), peaks.48.f))
peaks.all$zero <- peaks.all %over% peaks.0
peaks.all$day <- peaks.all %over% peaks.24
peaks.all$day2 <- peaks.all %over% peaks.48
peaks.all$diff <- (!peaks.all$zero & peaks.all$day2)

# Enhancer Heatmap
avoid <- resize(tss.l, width=4000, fix='center')

get_enhancers <- function(loci, rdata.path, peaks, tss){

	result.list <- lapply(seq_along(loci), function(i, rdata.path, peaks, tss){
		locus <- loci[i]
		message(locus$name)
		load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
		load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))
		
		buffer <- peaks[peaks %over% locus & peaks %outside% tss]
		buffer$score <- colMeans(mergeMat(runMat.open[open_run.qc < 2,], tiles, buffer), 
								 na.rm=T)
		return(buffer)

	}, rdata.path=rdata.path, peaks=peaks, tss=tss)

	result <- do.call(c, result.list)
}

enh.0 <- get_enhancers(loci, rdata.0.path, peaks.all, avoid)
enh.24 <- get_enhancers(loci, rdata.24.path, peaks.all, avoid)
enh.48 <- get_enhancers(loci, rdata.48.path, peaks.all, avoid)

hm <- rbind(enh.48$score, enh.24$score, enh.0$score)
rownames(hm) <- c('t48h', 't24h', 't0h')
colnames(hm) <- NULL

pdf('enhancher_heatmap.pdf', width=6, height=6)

# Skew scale towards high open run signal
scale <- c(rev(brewer.pal(9, 'Blues')[-c(1,3,5,7)]), brewer.pal(9, 'Reds'))
heatmap(hm, Rowv=NA, labCol=NA, col=scale)

# Plot color scale on new page
tot <- 20
tot_h <- 30
margin <- 1
width <- tot - (4*margin)
height_hm <- width / 2
height_tracks <- 5
height_reads <- 12

grid.newpage()

heights <- c(margin*1.5,height_hm,1.5*margin,height_tracks,1.5*margin, height_reads,margin)
widths <- c(3*margin, width, margin)
pushViewport(viewport(layout = grid.layout(7,3, heights=heights, widths=widths, respect=T)))
pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
pushViewport(viewport(layout = grid.layout(3,1, heights=c(1.5,2,1))))
pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
grid.raster(rev(scale))
grid.text(label=paste0('-  ', seq(0,0.8,l=5)), x=unit(0.5+1/3*3*margin/tot, 'npc'), 
		  y=seq(0,1,l=5), just='left', gp=gpar(fontsize=10))
popViewport()
dev.off()

