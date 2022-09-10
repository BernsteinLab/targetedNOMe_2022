# This code examines the GNAS locus in H9 cells.
# It calls DMRs and uses PCA to order the reads accordingly.
# It also generates a correlation heatmap of the locus.

# libraries and functions
library(rtracklayer)
library(GenomicFeatures)
library(ggfortify)
library(Gviz)
library(org.Hs.eg.db)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/210603_loci.bed'
rdata.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/rdata'
hap1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9_hap1.txt'
hap2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9_hap2.txt'

txdb.path <- '/seq/epiprod02/kdong/references/txdb/hg38.knownGene'
ideo.path <- '/seq/epiprod02/kdong/references/gviz/cytoBandIdeo.txt'

# Load in data
loci <- import.bed(loci.path)
hap1.names <- read.table(hap1.path, stringsAsFactors=F)$V1
hap2.names <- read.table(hap2.path, stringsAsFactors=F)$V1

cpg.list <- loadCpG(loci, rdata.path, hap1.names, hap2.names)

# Define Gviz plotting region and tracks
idx <- which(loci$name == 'GNAS')
locus <- loci[idx]
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

cpg <- cpg.list[[idx]]

# Force regions with ~50% methylation to have highest signal
cpg.avg <- -abs(cpg$avg - 0.5) + 0.5

pdf('dmr_calling_gnas.pdf', width=24, height=16)
# Reduce smoothing thresholds because of higher standard deviation
dmrpeaks <- findPeaks(cpg.avg, cpg, thresh=25, restrict=locus, smoothing=T, 
					  mult=1, smooth_mult=0.75, mean=median(cpg.avg), sd=sd(cpg.avg), 
					  erosion=4, dilation=4, plot=T, chr=chr, from=from, to=to, 
					  itrack=itrack, gtrack=gtrack, txTr=txTr)
dev.off()

# Retrieve methylation data as a data frame
meMat.list <- lapply(idx, function(i, rdata.path, hap1.names, hap2.names, peaks){
	locus <- loci[i]
	message(locus$name)
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))

	# Replace zeros with -1
	meMat.cg.f[meMat.cg.f == 0] <- -1
	buffer <- mergeMat(meMat.cg.f, CpG_fil, peaks)
	buffer <- as.data.frame(buffer)

	# Add haplotype information
	hap1.idx <- chopStr(rownames(meMat.cg.f), 2) %in% hap1.names
	hap2.idx <- chopStr(rownames(meMat.cg.f), 2) %in% hap2.names

	buffer$Haplotype <- NA
	buffer$Haplotype[hap1.idx] <- '1'
	buffer$Haplotype[hap2.idx] <- '2'

	return(buffer)
}, rdata.path=rdata.path, hap1.names=hap1.names, hap2.names=hap2.names, peaks=dmrpeaks)

data <- meMat.list[[1]]

# Perform PCA on reads without NAs
data.pca <- prcomp(data[rowSums(!is.na(data)) == ncol(data), 1:(ncol(data)-1)])

pdf('gnas_pca.pdf', width=10, height=6)
autoplot(data.pca, data=data[rowSums(!is.na(data)) == ncol(data),], colour='Haplotype', size=4)

# Look at PC weights to determine 
pc1 <- data.pca$rotation[,1]
names(pc1) <- paste0('Peak ', seq(length(dmrpeaks)))
barplot(pc1, main='Peak Importance for Principal Component 1')
dmrpeaks$pc1 <- pc1

# Get indices of true DMRs 
peak_idx <- which(abs(pc1) > sd(pc1))

# Project all data onto PC1 and PC2
pc_scores <- as.matrix(data[,-ncol(data)])
pc_scores[is.na(pc_scores)] <- 0
pc_scores <- pc_scores %*% data.pca$rotation[,c(1,2)]
pc_scores <- as.data.frame(pc_scores)
pc_scores$Haplotype <- data$Haplotype

# Plot all reads that cover at least one DMR
ggplot(pc_scores[rowSums(!is.na(data[,peak_idx])) > 0,], aes(x=PC1, y=PC2, color=Haplotype)) + 
	geom_point(size=4) + theme_bw(base_size=24)
dev.off()

# Order Reads by PCA results
pdf('gnas_ordered_reads.pdf', width=20, height=12)
all_loci.list <- lapply(idx, function(i, rdata.path,  
									  peaks, pc_scores, itrack, gtrack, txTr){
	locus <- loci[i]
	CpG <- cpg.list[[i]]

	# Load pipeline outputs
	message(sprintf('Loading %s data', locus$name))
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	# Select reads with both CpG and GpC information that cover 50% of the locus
	buffer <- mergeMat(meMat.cg.f, CpG, peaks)
	good_reads <- rownames(meMat.cg.f)[rowSums(!is.na(buffer)) > 0]
	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% locus]))
	select <- rownames(meMat.gc.f) %in% good_reads & (readCov > 0.5)
	message(sprintf('%s reads selected', sum(select)))

	# Order reads along PC1
	match_idx <- match(rownames(meMat.gc.f)[select], good_reads)
	order <- order(pc_scores$PC1[match_idx])
	haps <- pc_scores$Haplotype[match_idx][order]

	# Assign colors to haplotypes
	colors <- haps
	colors <- replace(colors, colors=='1', '#f8766d')
	colors <- replace(colors, colors=='2', '#00bfc4')

	# Generate plot
	# Consider margins out of 20
	# If right and bottom margins are equal, left margin equals 3X regular margin, 
	# top equals 1.5X, middle panel equals 1.5X as well
	tot <- 20
	margin <- 1
	width <- tot - (6 * margin)
	height <- width / 2
	# # Test pattern
	# grid.newpage()

	# pushViewport(viewport(layout=grid.layout(4, 3, heights=c(margin, (tot-2*margin)*2/3, (tot-2*margin)/3, margin), 
	# 										 widths=c(2.5*margin, width, 2.5*margin))))
	# for(i in seq(4)){
	#     for(j in seq(3)){
	#         pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	#         grid.rect()
	#         grid.circle()
	#         popViewport()
	#     }
	# }
	plot.new()
	grid.newpage()
	# Create layout
	pushViewport(viewport(layout=grid.layout(4, 3, heights=c(margin, (tot-2*margin)*2/3, (tot-2*margin)/3, margin), 
											 widths=c(2.5*margin, width, 2.5*margin))))
	
	# Title
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(sprintf('%s - Read Visualization', locus$name), gp=gpar(fontsize=20))
	popViewport()
	
	# Reads
	pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))
	plotReadMat(meMat.gc.f, GpC_fil, meMat.cg.f, CpG, open_runs, closed_runs, subset=select, 
				order=order, gap=5, region=locus, open=F, linker=F, short=F, cpg=T, peaks=NULL)
	popViewport()

	# Dendrogram & Haplotype colors
	pushViewport(viewport(layout=grid.layout(1, 2, widths=c(14, 2)), layout.pos.row=2, layout.pos.col=1))
	
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	dend <- as.dendrogram(hclust(dist(pc_scores$PC1[match_idx])))
	par(plt=c(0.75*margin/tot, 2.5*margin/tot, (margin+(tot-2*margin)/3)/tot, (tot-margin)/tot), new=T)
	plot(dend, leaflab='none', horiz=T, yaxt='n', yaxs='i')
	popViewport()

	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	sideColorBar(colors, 3)
	grid.lines(x=c(1,1), y=c(0,1), gp=gpar(lwd=2))
	popViewport()

	popViewport()

	# Legend
	pushViewport(viewport(layout.pos.row=2, layout.pos.col=3))
	grid.points(x=unit(1/6, 'npc'), y=unit(8/9, 'npc'), pch=21, size=unit(0.9, 'char'), 
				gp=gpar(fill='#0b529d', col='#0b529d'))
	grid.points(x=unit(1/6, 'npc'), y=unit(7/9, 'npc'), pch=21, size=unit(0.9, 'char'), 
				gp=gpar(fill='#ed2127', col='#ed2127'))
	grid.rect(x=unit(1/6, 'npc'), y=unit(3/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#f8766d', fill='#f8766d'))
	grid.rect(x=unit(1/6, 'npc'), y=unit(2/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#00bfc4', fill='#00bfc4'))
	grid.text(label=c('Unmethylated CpG', 'Methylated CpG', '', '', '', 'Haplotype 1', 'Haplotype 2', ''), 
			  x=unit(3/12, 'npc'), y=rev(seq(1:8)/9), just='left', gp=gpar(fontsize=12))
	popViewport()
	
	# Gviz Tracks
	pushViewport(viewport(layout.pos.row = c(3), layout.pos.col=c(2)))
	chr <- as.character(seqnames(locus))
	from <- start(locus)
	to <- end(locus)

	meanTrack1 <- DataTrack(range=CpG, data=CpG$hap1, type='mountain', name='0h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0, 1), size=3)
	meanTrack2 <- DataTrack(range=CpG, data=CpG$hap2, type='mountain', name='24h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0, 1), size=3)

	peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Peaks', 
								 rot.title=0, cex.title=1, stacking='dense')

	# Relative sizes
	sizes <- c(1, 2, 2, 1.5, 1.5, 1)
	plotTracks(c(peakCalls, meanTrack1, meanTrack2, txTr, gtrack, itrack), type='polygon', 
			   sizes=sizes, chromosome=chr, from=from, to=to, max.height=1, innerMargin=0,
			   collapseTranscripts ="meta", transcriptAnnotation='symbol', add=T, panel.only=T, 
			   margin=0, title.width=1)

	popViewport()

	# Custom Title for Gviz plots
	pushViewport(viewport(layout=grid.layout(length(sizes), 3, heights=sizes, widths=c(2.5*margin, width, 2.5*margin)), 
						  layout.pos.row = c(3), layout.pos.col=c(1, 2, 3)))

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA, col='#A9A9A9'))
	grid.yaxis(at=seq(0, 1, l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10, col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0, 1, l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10, col='#A9A9A9'))
	popViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
	grid.text('Called DMRs', gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
	grid.text('Haplotype 1\nCpG Methylation', gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
	grid.text('Haplotype 2\nCpG Methylation', gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=1))
	grid.text('Genes', gp=gpar(fontsize=12))
	popViewport()

	popViewport()

	return()

}, rdata.path=rdata.path, 
   peaks=dmrpeaks[peak_idx], pc_scores, itrack=itrack, gtrack=gtrack, txTr=txTr)
dev.off()

# Bin-based correlation
pdf('gnas_200.pdf', width=12, height=12)
all_loci.list <- lapply(idx, function(i, rdata.path, 
									  peaks, pc_scores, itrack, gtrack, txTr){
	locus <- loci[i]
	CpG <- cpg.list[[i]]

	# Load pipeline outputs
	message(sprintf('Loading %s data', locus$name))
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	windows <- tile(locus, width=200)[[1]]

	buffer <- mergeMat(meMat.cg.f, CpG, windows)
	corMat <- cor(buffer, method = 'pearson', use='pairwise.complete.obs')

	# Consider margins out of 20
	# If right and bottom margins are equal, left margin equals 3X regular margin, top equals 1.5X, middle panel equals 1.5X as well
	tot <- 20
	margin <- 1
	width <- tot - (4*margin)
	height <- width / 2

	# # Test pattern
	# grid.newpage()
	# pushViewport(viewport(width=unit(1, 'snpc'), height=unit(1,'snpc')))
	# pushViewport(viewport(x=(3*margin+width/2)/tot, y=(2.5*margin+height)/tot, width=height/tot* 2/sqrt(2), height=height/tot * 2/sqrt(2), angle=45))
	# grid.rect()
	# grid.circle()
	# popViewport()

	# pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin*1.5,height,1.5*margin,height,margin), widths=c(3*margin, width, margin))))
	# for(i in seq(5)){
	#     for(j in seq(3)){

	#         if(i== 2 & j ==2 ){
	#             next
	#         }
	#         pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	#         grid.rect()
	#         grid.circle()
	#         popViewport()
	#     }
	# }

	grid.newpage()
	pushViewport(viewport(width=unit(1, 'snpc'), height=unit(1,'snpc')))

	# Plot Trianglular Matrix
	pushViewport(viewport(x=(3*margin+width/2)/tot, y=(2.5*margin+height)/tot, width=unit(height/tot* 2/sqrt(2),'snpc'), height=unit(height/tot * 2/sqrt(2),'snpc'), angle=45))
	col_pal <- colorRampPalette(colors=rev(c('#b85c42', '#b95c42', '#c0755d', '#c47f6b', '#c98c7a', '#ce9989', '#d1a497', '#d8b2a7', '#ddbfb6', '#e1cbc5', '#e8d8d3', '#ede5e2', '#f2f1f1', '#f4f4f4', '#e1e9eb', '#d5e0e3', '#c8d6db', '#b9cad2', '#abc1ca', '#9fb8c3', '#91aeb9', '#81a4b0', '#759aa9', '#6791a0', '#5a8799', '#4e7e91')))(26)
	scale <- plotTriMatrix((corMat+1)/2, col_pal)  
	# highlight(windows, peaks)
	popViewport()

	pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin*1.5,height,1.5*margin,height,margin), widths=c(3*margin, width, margin))))
	
	# Title
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(sprintf('%s - CpG Methylation Correlation', locus$name), gp=gpar(fontsize=20))
	popViewport()
	
	# Color scale
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

	# Gviz Tracks
	pushViewport(viewport(layout.pos.row = c(3,4), layout.pos.col=c(2)))
	chr <- as.character(seqnames(locus))
	from <- start(locus)
	to <- end(locus)

	meanTrack1 <- DataTrack(range=CpG, data=CpG$hap1, type='mountain', name='0h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0, 1), size=3)
	meanTrack2 <- DataTrack(range=CpG, data=CpG$hap2, type='mountain', name='24h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0, 1), size=3)

	peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Peaks', 
								 rot.title=0, cex.title=1, stacking='dense')

	# Relative sizes
	sizes <- c(1, 2, 2, 1.5, 1.5, 1)
	plotTracks(c(peakCalls, meanTrack1, meanTrack2, txTr, gtrack, itrack), type='polygon', 
			   sizes=sizes, chromosome=chr, from=from, to=to, max.height=1, innerMargin=0,
			   collapseTranscripts ="meta", transcriptAnnotation='symbol', add=T, panel.only=T, 
			   margin=0, title.width=1)

	popViewport()

	# Custom Titles
	pushViewport(viewport(layout = grid.layout(length(sizes),3, heights=sizes, widths=c(3*margin, width, margin)), layout.pos.row = c(3,4), layout.pos.col=c(1,2,3)))

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA, col='#A9A9A9'))
	grid.yaxis(at=seq(0, 1, l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10, col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0, 1, l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), gp=gpar(fontsize=10, col='#A9A9A9'))
	popViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
	grid.text('Called DMRs', gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
	grid.text('Haplotype 1\nCpG Methylation', gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
	grid.text('Haplotype 2\nCpG Methylation', gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=1))
	grid.text('Genes', gp=gpar(fontsize=12))
	popViewport()

	popViewport()

	return()
}, rdata.path=rdata.path, 
   peaks=dmrpeaks[peak_idx], pc_scores, itrack=itrack, gtrack=gtrack, txTr=txTr)
dev.off()
