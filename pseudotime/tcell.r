# This code examines the CD28 locus in activated T-cells.
# It generates a correlation heatmap and focuses on 
# four regions of interest. In addition, it performs PCA
# and pseudotime analysis and orders the reads accordingly.

# libraries and functions
library(cowplot)
library(ggfortify)
library(Gviz)
library(org.Hs.eg.db)
library(reshape2)
library(rtracklayer)
library(TSCAN)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/TcellNov2020_loci.bed'
peaks.0.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/t0h_open_peaks.bed'
peaks.0.c.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/t0h_closed_peaks.bed'
rdata.0.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/rdata'

peaks.24.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t24h/workspace/t24h_open_peaks.bed'
peaks.24.c.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t24h/workspace/t24h_closed_peaks.bed'
rdata.24.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t24h/workspace/rdata'

peaks.48.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/t48h_open_peaks.bed'
peaks.48.c.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/t48h_closed_peaks.bed'
rdata.48.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/rdata'

txdb.path <- '/seq/epiprod02/kdong/references/txdb/hg38.knownGene'
ideo.path <- '/seq/epiprod02/kdong/references/gviz/cytoBandIdeo.txt'

# Load in data
loci <- import.bed(loci.path)

peaks.0 <- import.bed(peaks.0.path)
peaks.0.c <- import.bed(peaks.0.c.path)

peaks.24 <- import.bed(peaks.24.path)
peaks.24.c <- import.bed(peaks.24.c.path)

peaks.48 <- import.bed(peaks.48.path)
peaks.48.c <- import.bed(peaks.48.c.path)

gpc.list <- loadGpC(loci, rdata.0.path)
gpcs <- do.call(c, gpc.list)

# Count GpCs in peaks
hits <- findOverlaps(peaks.0, gpcs)
hits.df <- aggregate(rep(1, length(hits)), by=list(queryHits(hits)), FUN=sum)
peaks.0$tot <- NA
peaks.0$tot[hits.df$Group.1] <- hits.df$x

hits <- findOverlaps(peaks.24, gpcs)
hits.df <- aggregate(rep(1, length(hits)), by=list(queryHits(hits)), FUN=sum)
peaks.24$tot <- NA
peaks.24$tot[hits.df$Group.1] <- hits.df$x  

hits <- findOverlaps(peaks.48, gpcs)
hits.df <- aggregate(rep(1, length(hits)), by=list(queryHits(hits)), FUN=sum)
peaks.48$tot <- NA
peaks.48$tot[hits.df$Group.1] <- hits.df$x

peaks.all <- mergePeaks(mergePeaks(peaks.0[which(peaks.0$tot > 1)], 
								   peaks.24[which(peaks.24$tot > 1)]), peaks.48[which(peaks.48$tot > 1)])
peaks.all <- unique(peaks.all)
peaks.all$zero <- peaks.all %over% peaks.0
peaks.all$day <- peaks.all %over% peaks.24
peaks.all$day2 <- peaks.all %over% peaks.48
peaks.all$diff <- (!peaks.all$zero & peaks.all$day2)

peaks.all.c <- unique(mergePeaks(mergePeaks(peaks.0.c, peaks.24.c), peaks.48.c))

idx <- which(loci$name == 'CD28')
locus <- loci[idx]

# Define Gviz plotting region and tracks
chr <- as.character(seqnames(locus))
from <- start(locus)
to <- end(locus)

# Ideogram
bands <- read.table(ideo.path, sep='\t', header=F, stringsAsFactors=F, 
					col.names=c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
itrack <- IdeogramTrack(genome='hg38', bands=bands, chromosome=chr, cex=1.5)

# Genome axis
gtrack <- GenomeAxisTrack(add53=T, add35=T, littleTicks=F, labelPos='below', cex=1.5)

# Genes
txdb <- loadDb(txdb.path)
txTr <- GeneRegionTrack(txdb, chromosome=chr, start=from, end=to, cex.title=2, 
						rot.title=0, name='Genes', cex.group=1.2, 
						collapseTranscripts ="meta", transcriptAnnotation='symbol')

# Get correct gene symbols
if(length(gene(txTr)) > 0){
	symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
	symbol(txTr) <- symbols[gene(txTr)]
}

# Example accessibility correlation
pdf('cd28_combined.pdf', width=12, height=18)
all_loci.list <- lapply(idx, function(i, rdata.0.path, rdata.24.path, rdata.48.path,
									  roi.width, peaks.all, itrack, gtrack, txTr){
	locus <- loci[i]
	message(locus$name)	
	locus <- resize(locus, width=36000, fix='start')
	
	win.width=1000
	windows <- tile(locus, width=win.width)[[1]]

	# Handpicked peaks
	roi <- peaks.all[peaks.all %over% locus][c(4,5,11,14)]

	load(sprintf('%s/%s-run.RData', rdata.0.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.0.path, locus$name))
	meMat.gc.0 <- meMat.gc.f
	meMat.cg.0 <- meMat.cg.f
	open_runs.0 <- open_runs
	closed_runs.0 <- closed_runs
	open_run.qc.0 <- open_run.qc
	runMat.open.0 <- runMat.open
	load(sprintf('%s/%s-run.RData', rdata.24.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.24.path, locus$name))
	meMat.gc.24 <- meMat.gc.f
	meMat.cg.24 <- meMat.cg.f
	open_runs.24 <- open_runs
	closed_runs.24 <- closed_runs
	open_run.qc.24 <- open_run.qc
	runMat.open.24 <- runMat.open
	load(sprintf('%s/%s-run.RData', rdata.48.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.48.path, locus$name))
	meMat.gc.48 <- meMat.gc.f
	meMat.cg.48 <- meMat.cg.f
	open_runs.48 <- open_runs
	closed_runs.48 <- closed_runs
	open_run.qc.48 <- open_run.qc
	runMat.open.48 <- runMat.open

	meMat.gc.f <- rbind(meMat.gc.0, meMat.gc.24, meMat.gc.48)
	meMat.cg.f <- rbind(meMat.cg.0, meMat.cg.24, meMat.cg.48)
	open_runs <- c(open_runs.0, open_runs.24, open_runs.48)
	closed_runs <- c(closed_runs.0, closed_runs.24, closed_runs.48)
	open_run.qc <- c(open_run.qc.0, open_run.qc.24, open_run.qc.48)
	runMat.open <- rbind(runMat.open.0, runMat.open.24, runMat.open.48)

	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% roi]))
	select <- readCov >= 0.7 & open_run.qc < 2

	peaks <- peaks.all[peaks.all %over% locus]
	buffer <- mergeMat(runMat.open[open_run.qc < 2,], tiles, windows)

	corMat <- cor(buffer, method = 'pearson', use='pairwise.complete.obs')

	# Consider margins of unit 1
	# If right and bottom margins are equal, left margin equals 3X regular margin, 
	# top equals 1.5X, middle panel equals 1.5X as well
	tot <- 20
	tot_h <- 30
	margin <- 1
	width <- tot - (4*margin)
	height_hm <- width / 2
	height_tracks <- 5
	height_reads <- 12
	margin_top <- tot_h - height_hm - 1.5*margin - height_tracks - 1.5*margin - height_reads - margin

	# # Test Pattern
	# grid.newpage()
	# pushViewport(viewport(width=unit(2/3, 'snpc'), height=unit(1,'snpc')))
	# grid.rect()
	# pushViewport(viewport(x=(3*margin+width/2)/tot, y=(4*margin+height_reads+height_tracks)/tot_h, width=unit(width/tot / sqrt(2),'snpc'), height=unit(width/tot / sqrt(2),'snpc'), angle=45))
	# grid.rect()
	# grid.circle()
	# popViewport()

	# pushViewport(viewport(layout = grid.layout(7,3, heights=c(margin_top,height_hm,1.5*margin,height_tracks,1.5*margin, height_reads,margin), 
	# 						widths=c(3*margin, width, margin), respect=F)))
	# for (i in seq(7)){
	# 	for (j in seq(3)){
	# 		if(i==2 & j == 2){
	# 			next
	# 		}
	# 		pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	# 		grid.rect()
	# 		grid.circle()
	# 		popViewport()
	# 	}
	# }

	grid.newpage()
	# # Fix 2:3 ratio
	# pushViewport(viewport(width=unit(2/3, 'snpc'), height=unit(1,'snpc')))
	pushViewport(viewport(x=(3*margin+width/2)/tot, y=(4*margin+height_reads+height_tracks)/tot_h, width=unit(width/tot / sqrt(2),'npc'), height=unit(width/tot_h / sqrt(2),'npc'), angle=45))
	
	col_pal <- c(colorRampPalette(colors = c("blue"))(5), colorRampPalette(colors = c("blue","yellow"))(5), 
		colorRampPalette(colors = c("yellow"))(10), colorRampPalette(colors = c("yellow","red"))(5), colorRampPalette(colors = c("red"))(5))
	scale <- plotTriMatrix((corMat+1)/2, col_pal)
	n <- nrow(corMat)
	# if((n %% 2) == 0){
	#     grid.polyline(unit(c(rep(c(0, 1,1, 0), n/2),0,1), 'npc'), unit(rep(0:n, rep(2,n+1))/n, 'npc'), gp=gpar(col='white', lwd=2))
	#     grid.polyline(unit(rep(0:n, rep(2,n+1))/n, 'npc'), unit(c(rep(c(0, 1,1, 0), n/2),0,1), 'npc'), gp=gpar(col='white', lwd=2))
	# } else {
	#     grid.polyline(unit(rep(c(0, 1,1, 0), floor(n/2)), 'snpc'), unit(rep(0:n, rep(2,n+1))/n, 'snpc'), gp=gpar(col='white', lwd=2))
	#     grid.polyline(unit(rep(0:n, rep(2,n+1))/n, 'snpc'), unit(rep(c(0, 1,1, 0), floor(n/2)), 'snpc'), gp=gpar(col='white', lwd=2))
	# }
	# highlight(windows, peaks)

	# popViewport()
	popViewport()

	pushViewport(viewport(layout = grid.layout(7,3, heights=c(margin*1.5,height_hm,1.5*margin,height_tracks,1.5*margin, height_reads,margin), 
						  widths=c(3*margin, width, margin), respect=T)))
	# Color scale
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(sprintf('%s - Open Run Correlation', locus$name), gp=gpar(fontsize=20))
	popViewport()

	pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
	pushViewport(viewport(layout = grid.layout(3,1, heights=c(1.5,2,1))))
	
	pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
	grid.raster(scale)
	grid.text(label=paste0('-  ', seq(-1,1,l=5)), x=unit(0.5+1/3*3*margin/tot, 'npc'), y=seq(0,1,l=5), just='left', gp=gpar(fontsize=10))
	popViewport()

	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	grid.text('Correlation', y=1/6, just='bottom',  gp=gpar(fontsize=12))
	popViewport()
	
	popViewport()
	popViewport()

	# Arcs
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	arcs(corMat, n=7, peak.ind=which(windows %over% peaks), min.height=0, hl.ratio=1.5*margin/width)
	popViewport()

	# Gviz
	pushViewport(viewport(layout.pos.row = c(4), layout.pos.col=c(2)))
	chr <- as.character(seqnames(locus))
	from <- start(locus)
	to <- end(locus)
	
	meanTrack1 <- DataTrack(range=tiles, data=colMeans(runMat.open.0, na.rm=T), name='0h', cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)
	meanTrack2 <- DataTrack(range=tiles, data=colMeans(runMat.open.24, na.rm=T), name='24h', cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)
	meanTrack3 <- DataTrack(range=tiles, data=colMeans(runMat.open.48, na.rm=T), name='48h', cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)

	peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Peaks', rot.title=0, cex.title=1, stacking='dense')

	sizes <- c(1,2,2,2)
	plotTracks(c(peakCalls, meanTrack1, meanTrack2, meanTrack3), sizes=sizes, chromosome=chr, from=from, to=to, max.height=1, panel.only=T, 
			   innerMargin=0, collapseTranscripts ="meta", transcriptAnnotation='symbol', add=T, margin=0, title.width=1)

	popViewport()

	# Custom axes and titles for Gviz plots
	pushViewport(viewport(layout.pos.row = c(4), layout.pos.col=c(1,2,3)))
	pushViewport(viewport(layout = grid.layout(length(sizes),3, heights=sizes, widths=c(3*margin, width, margin))))

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
	grid.text('Called Peaks',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
	grid.text('0h\nOpen Run',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
	grid.text('24h\nOpen Run',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=1))
	grid.text('48h\nOpen Run',  gp=gpar(fontsize=12))
	popViewport()

	popViewport()
	popViewport()

	# Legend
	pushViewport(viewport(layout=grid.layout(1,2, widths=c(14,2)), layout.pos.row=6, layout.pos.col=1))
	
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	grid.rect(x=unit(1/6,'npc'), y=unit(8/9, 'npc'), width=unit(0.5/10, 'snpc'), height=unit(0.5/10, 'snpc'), gp=gpar(fill='#1B9E77'))

	grid.rect(x=unit(1/6,'npc'), y=unit(6/9, 'npc'), width=unit(1/10, 'snpc'), height=unit(0.25/10, 'snpc'), gp=gpar(col='#992C0E', fill='#992C0E'))
	grid.rect(x=unit(1/6,'npc'), y=unit(5/9, 'npc'), width=unit(1/10, 'snpc'), height=unit(0.25/10, 'snpc'), gp=gpar(col='#254BD9', fill='#254BD9'))
	grid.rect(x=unit(1/6,'npc'), y=unit(4/9, 'npc'), width=unit(1/10, 'snpc'), height=unit(0.25/10, 'snpc'), gp=gpar(col='#529906', fill='#529906'))
	grid.text(label=c('Open Run', '', 'Resting T-cells', '24h Activation', '48h Activation'), x=unit(3/12, 'npc'), y=rev((4:8)/9), just='left', gp=gpar(fontsize=12))
	popViewport()

	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	colors <- rep(c('#992C0E', '#254BD9', '#529906'), c(sum(select[1:nrow(meMat.gc.0)]), sum(select[(nrow(meMat.gc.0)+1):(nrow(meMat.gc.0)+nrow(meMat.gc.24))]), 
														sum(select[(nrow(meMat.gc.0)+nrow(meMat.gc.24)+1):length(select)])))
	sideColorBar(colors, 3)
	grid.lines(x=c(1,1), y=c(0,1), gp=gpar(lwd=2))
	popViewport()
	
	popViewport()

	pushViewport(viewport(layout=grid.layout(1,length(roi)), layout.pos.row=c(6), layout.pos.col=2))
	for(i in seq_along(roi)){
		peak <- resize(roi[i], width=roi.width, fix='center')
		cat(sprintf('Region of Interest %s\n',i))
		pushViewport(viewport(layout.pos.row=1, layout.pos.col=i))
		plotReadMat(meMat.gc.f, GpC_fil, meMat.cg.f, CpG_fil, open_runs, closed_runs, subset=select, gap=5, region=peak, open=T, linker=F, short=F, cpg=F, peaks=NULL)
		if(i > 1){
			grid.lines(x=c(0,0), y=c(0,1), gp=gpar(lwd=2.5))
		}
		popViewport()
		cat('\n')
	}
	popViewport()

	# Show context
	pushViewport(viewport(layout.pos.row = 5, layout.pos.col=2))
	context(locus, roi, roi.width, reverse=T)
	popViewport()

	popViewport()

	return()
}, rdata.0.path=rdata.0.path, rdata.24.path=rdata.24.path, rdata.48.path=rdata.48.path, 
	roi.width=500, peaks.all=peaks.all, itrack=itrack, gtrack=gtrack, txTr=txTr)
dev.off()

# PCA
meMat.list <- lapply(idx, function(i, rdata.0.path, rdata.24.path, rdata.48.path, peaks.all){
	locus <- loci[i]

	# Load and save run matrices
	message(locus$name)
	load(sprintf('%s/%s-runMat.RData', rdata.0.path, locus$name))
	runMat.open.0 <- runMat.open[open_run.qc < 2,]
	load(sprintf('%s/%s-runMat.RData', rdata.24.path, locus$name))
	runMat.open.24 <- runMat.open[open_run.qc < 2,]
	load(sprintf('%s/%s-runMat.RData', rdata.48.path, locus$name))
	runMat.open.48 <- runMat.open[open_run.qc < 2,]

	# Focus on first part of locus
	locus <- resize(locus, width=36500, fix='start')
	peaks <- peaks.all[peaks.all %over% locus]
	message(sprintf('Using %s peaks', length(peaks)))

	runMat.open <- rbind(runMat.open.0, runMat.open.24, runMat.open.48)

	# Create data frame
	buffer <- mergeMat(runMat.open, tiles, peaks)
	buffer <- as.data.frame(buffer)
	cond <- rep(c('0h', '24h', '48h'), c(nrow(runMat.open.0), nrow(runMat.open.24), nrow(runMat.open.48)))
	buffer$Stimulation <- cond

	return(buffer)
}, rdata.0.path=rdata.0.path, rdata.24.path=rdata.24.path, 
   rdata.48.path=rdata.48.path, peaks.all=peaks.all)

data <- meMat.list[[1]]
data.pca <- prcomp(data[rowSums(!is.na(data)) == ncol(data), 1:(ncol(data)-1)], scale=T)

# Add random seed setting to TSCAN's exprmclust function
exprmclust2 <- function (data, clusternum = 2:9, modelNames = "VVV", reduce = T) {
	set.seed(1234)
	if (reduce) {
		sdev <- prcomp(t(data), scale = T)$sdev[1:20]
		x <- 1:20
		optpoint <- which.min(sapply(2:10, function(i) {
		x2 <- pmax(0, x - i)
		sum(lm(sdev ~ x + x2)$residuals^2)
		}))
		pcadim = optpoint + 1
		tmpdata <- t(apply(data, 1, scale))
		colnames(tmpdata) <- colnames(data)
		tmppc <- prcomp(t(tmpdata), scale = T)
		pcareduceres <- t(tmpdata) %*% tmppc$rotation[, 1:pcadim]
  	}
  	else {
		pcareduceres <- t(data)
  	}
  	clusternum <- clusternum[clusternum > 1]
  	res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, 
	modelNames = modelNames))
	clusterid <- apply(res$z, 1, which.max)
	clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = res$G)
	for (cid in 1:res$G) {
		clucenter[cid, ] <- colMeans(pcareduceres[names(clusterid[clusterid == 
		cid]), , drop = F])
	}
	dp <- as.matrix(dist(clucenter))
	gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
	dp_mst <- minimum.spanning.tree(gp)
	list(pcareduceres = pcareduceres, MSTtree = dp_mst, clusterid = clusterid, 
		 clucenter = clucenter)
}
environment(exprmclust2) <- asNamespace('TSCAN')

data.clust <- exprmclust2(t(data[rowSums(!is.na(data)) == ncol(data), 1:(ncol(data)-1)]), modelNames='EVV')
data.order <- TSCANorder(data.clust)
order <- match(data.order, seq(nrow(data)))

# Visualize PCA and clustering
pdf('cd28_pca.pdf', width=10, height=14)
p1 <- autoplot(data.pca, data=data[rowSums(!is.na(data)) == ncol(data),], colour='Stimulation', size=4) + 
		theme_bw(base_size=18) + theme(legend.position='top') 
p2 <- plotmclust(data.clust, show_cell_names=F) + xlab('Principal Component 1') + ylab('Principal Component 2') + ylim(-3,3) + xlim(-4,4)
plot_grid(p1, p2, align = "vh", ncol=1, nrow=2)
dev.off()

# Reorder and smooth data
data.ordered <- data[order, -ncol(data)]
colnames(data.ordered) <- seq(ncol(data.ordered))
smooth <- function(x, n=5){
	stats::filter(x, rep(1/n, n), sides=2)
}
data.ordered.s <- as.data.frame(apply(data.ordered, 2, function(col){
	smooth(col, 51)
}))
data.ordered.s$x <- seq(nrow(data.ordered))

data.melt <- melt(data.ordered.s, id='x')

# Plot peak accessibility over pseudotime
pdf('cd28_peaks.pdf', width=10, height=6)
ggplot(data.melt, aes(x=x, y=value, col=variable)) + geom_line() +
		xlab('Pseudotime') + ylab('Open Run Signal') + ylim(0,1) + scale_color_discrete(name='Peak ID') + 
		theme_bw(base_size=20) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
dev.off()

# Order reads by pseudotime
pdf('cd28_ordered_reads.pdf', width=16, height=12)
all_loci.list <- lapply(idx, function(i, rdata.0.path, rdata.24.path, 
									  rdata.48.path, peaks.all, itrack, gtrack, txTr){
	locus <- loci[i]
	message(locus$name)	
	locus <- resize(locus, width=38000, fix='start')

	load(sprintf('%s/%s-run.RData', rdata.0.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.0.path, locus$name))
	meMat.gc.0 <- meMat.gc.f
	meMat.cg.0 <- meMat.cg.f
	open_runs.0 <- open_runs
	closed_runs.0 <- closed_runs
	open_run.qc.0 <- open_run.qc
	runMat.open.0 <- runMat.open
	load(sprintf('%s/%s-run.RData', rdata.24.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.24.path, locus$name))
	meMat.gc.24 <- meMat.gc.f
	meMat.cg.24 <- meMat.cg.f
	open_runs.24 <- open_runs
	closed_runs.24 <- closed_runs
	open_run.qc.24 <- open_run.qc
	runMat.open.24 <- runMat.open
	load(sprintf('%s/%s-run.RData', rdata.48.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.48.path, locus$name))
	meMat.gc.48 <- meMat.gc.f
	meMat.cg.48 <- meMat.cg.f
	open_runs.48 <- open_runs
	closed_runs.48 <- closed_runs
	open_run.qc.48 <- open_run.qc
	runMat.open.48 <- runMat.open

	meMat.gc.f <- rbind(meMat.gc.0, meMat.gc.24, meMat.gc.48)
	meMat.cg.f <- rbind(meMat.cg.0, meMat.cg.24, meMat.cg.48)
	open_runs <- c(open_runs.0, open_runs.24, open_runs.48)
	closed_runs <- c(closed_runs.0, closed_runs.24, closed_runs.48)
	open_run.qc <- c(open_run.qc.0, open_run.qc.24, open_run.qc.48)
	runMat.open <- rbind(runMat.open.0, runMat.open.24, runMat.open.48)

	peaks <- peaks.all[peaks.all %over% locus]
	peaks <- peaks[-length(peaks)]

	buffer <- mergeMat(runMat.open, tiles, peaks)
	buffer <- as.data.frame(buffer)
	cond <- rep(c('0h', '24h', '48h'), c(nrow(runMat.open.0), nrow(runMat.open.24), nrow(runMat.open.48)))
	buffer$cond <- cond

	select <- rowSums(!is.na(buffer)) == ncol(buffer) & open_run.qc < 2
	message(sprintf('%s reads selected', sum(select)))
	data.clust <- exprmclust2(t(buffer[select, 1:(ncol(buffer)-1)]), modelNames='EVV')
	plotmclust(data.clust)
	data.order <- TSCANorder(data.clust)
	order <- match(data.order, seq(nrow(buffer))[select])

	# Consider margins out of 20
	# If right and bottom margins are equal, left margin equals 3X regular margin, top equals 1.5X, middle panel equals 1.5X as well
	tot <- 20
	margin <- 1
	width <- tot - (4*margin)
	height <- width / 2

	# # Test pattern
	# grid.newpage()

	# pushViewport(viewport(layout = grid.layout(4,3, heights=c(margin,(tot-margin)*2/3,(tot-margin)/3,margin), widths=c(3*margin, width, margin))))
	# for(i in seq(4)){
	#     for(j in seq(3)){
	#         pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	#         grid.rect()
	#         grid.circle()
	#         popViewport()
	#     }
	# }

	grid.newpage()
	# Create layout
	pushViewport(viewport(layout = grid.layout(4,3, heights=c(margin,(tot-margin)*2/3,(tot-margin)/3,margin), widths=c(3*margin, width, margin))))
	
	# Title
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(sprintf('%s - Read Visualization', locus$name), gp=gpar(fontsize=20))
	popViewport()
	
	# Reads
	pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))
	plotReadMat(meMat.gc.f, GpC_fil, meMat.cg.f, CpG_fil, open_runs, closed_runs, subset=select, order=order,  gap=5, region=locus, open=T, peaks=NULL)
	popViewport()

	# Legend
	pushViewport(viewport(layout=grid.layout(1,2, widths=c(14,2)), layout.pos.row=2, layout.pos.col=1))

	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	grid.rect(x=unit(1/6,'npc'), y=unit(8/9, 'npc'), width=unit(0.5/10, 'snpc'), height=unit(0.5/10, 'snpc'), gp=gpar(fill='#1B9E77'))
	grid.rect(x=unit(1/6,'npc'), y=unit(3/9, 'npc'), width=unit(1/10, 'snpc'), height=unit(0.25/10, 'snpc'), gp=gpar(col='#992C0E', fill='#992C0E'))
	grid.rect(x=unit(1/6,'npc'), y=unit(2/9, 'npc'), width=unit(1/10, 'snpc'), height=unit(0.25/10, 'snpc'), gp=gpar(col='#254BD9', fill='#254BD9'))
	grid.rect(x=unit(1/6,'npc'), y=unit(1/9, 'npc'), width=unit(1/10, 'snpc'), height=unit(0.25/10, 'snpc'), gp=gpar(col='#529906', fill='#529906'))
	grid.text(label=c('Open Run', '', '', '', '', 'Resting T-cells', '24h Activation', '48h Activation'), x=unit(3/12, 'npc'), y=rev(seq(1:8)/9), just='left', gp=gpar(fontsize=12))
	popViewport()

	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	colors <- rep(c('#992C0E', '#254BD9', '#529906'), c(nrow(runMat.open.0), nrow(runMat.open.24), nrow(runMat.open.48)))
	sideColorBar(colors[select][order], 3)
	grid.lines(x=c(1,1), y=c(0,1), gp=gpar(lwd=2))
	popViewport()
	
	popViewport()

	# Gviz Tracks
	pushViewport(viewport(layout.pos.row = c(3), layout.pos.col=c(2)))
	chr <- as.character(seqnames(locus))
	from <- start(locus)
	to <- end(locus)

	meanTrack1 <- DataTrack(range=tiles, data=colMeans(runMat.open.0, na.rm=T), name='0h', cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)
	meanTrack2 <- DataTrack(range=tiles, data=colMeans(runMat.open.24, na.rm=T), name='24h', cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)
	meanTrack3 <- DataTrack(range=tiles, data=colMeans(runMat.open.48, na.rm=T), name='48h', cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)

	peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Peaks', rot.title=0, cex.title=1, stacking='dense')

	sizes <- c(1, 2, 2, 2, 1.5, 1.5, 1)
	plotTracks(c(peakCalls, meanTrack1, meanTrack2, meanTrack3, txTr, gtrack, itrack), 
			   sizes=sizes, chromosome=chr, from=from, to=to, max.height=1, 
			   innerMargin=0, collapseTranscripts ="meta", transcriptAnnotation='symbol', 
			   add=T, panel.only=T, margin=0, title.width=1)

	popViewport()

	# Custom Titles for Gviz plots
	pushViewport(viewport(layout = grid.layout(length(sizes),3, heights=sizes, widths=c(3*margin, width, margin)), layout.pos.row = c(3), layout.pos.col=c(1,2,3)))

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
	grid.text('Called Peaks',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
	grid.text('0h\nOpen Run',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
	grid.text('24h\nOpen Run',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=1))
	grid.text('48h\nOpen Run',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 5, layout.pos.col=1))
	grid.text('Genes',  gp=gpar(fontsize=12))
	popViewport()

	popViewport()

	return()
}, rdata.0.path=rdata.0.path, rdata.24.path=rdata.24.path, rdata.48.path=rdata.48.path, 
   peaks.all=peaks.all, itrack=itrack, gtrack=gtrack, txTr=txTr)
dev.off()