# This code examines the imprinted SNRPN locus in H9 cells.
# It focuses on two selected regions and visualizes the 
# open run signal, CpG methylation, and SNPs between the alleles.

# libraries and functions
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(GenomicAlignments)
library(GenomicFeatures)
library(Gviz)
library(org.Hs.eg.db)
library(VariantAnnotation)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/210603_loci.bed'
bam.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9.phased.bam'
vcf.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9.phased.vcf'
rdata.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/rdata'
hap1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9_hap1.txt'
hap2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9_hap2.txt'

txdb.path <- '/seq/epiprod02/kdong/references/txdb/hg38.knownGene'
ideo.path <- '/seq/epiprod02/kdong/references/gviz/cytoBandIdeo.txt'

refLoc2queryLoc <- function(pos, read){
	# Converts reference positions to read positions
	#
	# Args:
	#	pos: vector of reference positions
	#	read: single GenomicAlignments read
	#
	# Returns:
	#	Vector of query positions
	cigar <- read@cigar
	refRanges <- cigarRangesAlongReferenceSpace(cigar)[[1]]
	qRanges <- cigarRangesAlongQuerySpace(cigar)[[1]]
	snp.loc <- pos - start(read) + 1
	return.list <- lapply(snp.loc, function(loc){
		id <- refRanges %over% IRanges(loc,loc) & width(refRanges) != 0
		q <- qRanges[id]
		if(width(q) == 0){
			return(NA)
		}
		diff <- loc - start(refRanges[id])
		return(start(q)+diff)
	})
	return(do.call(c, return.list))
}

# Load in data
loci <- import.bed(loci.path)
bam <- readGAlignments(bam.path, use.names=T, param=ScanBamParam(what='seq'))
genome <- BSgenome.Hsapiens.UCSC.hg38
snps <-readVcf(vcf.path, seqinfo(genome))

hap1.names <- read.table(hap1.path, stringsAsFactors=F)$V1
hap2.names <- read.table(hap2.path, stringsAsFactors=F)$V1

gpc.list <- loadGpC(loci, rdata.path, hap1.names, hap2.names)

# Determine SNP type
snps.het <- snps@rowRanges[grepl('\\|', geno(snps)@listData$GT[,'sample'])]
snps.het$ref <- substr(names(snps.het), nchar(names(snps.het))-2, nchar(names(snps.het))-2)
snps.het$alt <- substr(names(snps.het), nchar(names(snps.het)), nchar(names(snps.het)))
snps.het$geno <- 'hetero'
snps.hom <- snps@rowRanges[grepl('1/1', geno(snps)@listData$GT[,'sample'])]
snps.hom$ref <- substr(names(snps.hom), nchar(names(snps.hom))-2, nchar(names(snps.hom))-2)
snps.hom$alt <- substr(names(snps.hom), nchar(names(snps.hom)), nchar(names(snps.hom)))
snps.hom$geno <- 'homo'
snps.all <- sort(c(snps.hom, snps.het))

idx <- which(loci$name == 'SNRPN')
locus <- loci[idx]

# Create SNP matrix
reads <- bam[bam %over% locus]
snps <- snps.all[snps.all %over% locus]

snp.mat <- matrix(NA, nrow=length(reads), ncol=length(snps))
rownames(snp.mat) <- mapply(function(strand, name){
							return(sprintf('%s-%s', switch(strand, '+'='P', '-'='M'), name))
						}, as.character(strand(reads)), names(reads))

for(i in seq_along(reads)){
	read <- reads[i]

	# Find overlapping SNPs
	snp.id <- snps %over% read
	if(sum(snp.id) == 0){
		next
	}    
	
	# Get query positions from reference positions
	seq.id <- refLoc2queryLoc(start(snps[snp.id]), read)
	good.seq.id <- !is.na(seq.id)
	
	# Determine if sequenced base is the same as alternate allele
	nts <- as.character(as.vector(mcols(read)@listData$seq[[1]][seq.id[good.seq.id]]))
	bools <- snps$alt[snp.id][good.seq.id] == nts
	snp.mat[i, which(snp.id)[good.seq.id]] <- bools
}

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

# Regions of interest
peaks <- c(GRanges('chr15', IRanges(24770250, 24775250)), 
		   GRanges('chr15', IRanges(24821350, 24826350)))

pdf('snrpn_plots.pdf', width=20, height=12)
# Open Run
all_loci.list <- lapply(idx, function(i, rdata.path, hap1.names, hap2.names,
									  peaks, pc_scores, itrack, gtrack, txTr){
	locus <- loci[i]
	GpC <- gpc.list[[i]]

	# Load pipeline outputs
	message(sprintf('Loading %s data', locus$name))
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% peaks]))
	select <- readCov >= 0.7

	hap1.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap1.names
	hap2.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap2.names

	order <- c(which(hap1.idx), which(hap2.idx))
	colors <- c(rep('#f8766d', sum(hap1.idx)), rep('#00bfc4', sum(hap2.idx)))
	# Consider margins out of 20
	# If right and bottom margins are equal, left margin equals 3X regular margin, 
	# top equals 1.5X, middle panel equals 1.5X as well
	tot <- 20
	margin <- 1
	width <- tot - (6*margin)
	height <- width / 2

	# # Test pattern
	# grid.newpage()

	# pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin,(tot-margin)*2/3,margin, (tot-margin)/3,margin), 
	# 										   widths=c(3*margin, width, margin))))
	# for(i in seq(5)){
	#     for(j in seq(3)){
	#         pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	#         grid.rect()
	#         grid.circle()
	#         popViewport()
	#     }
	# }

	grid.newpage()
	
	# Create layout
	pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin,(tot-margin)*2/3,margin, (tot-margin)/3,margin), 
											   widths=c(3*margin, width, margin))))
	
	# Title
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(sprintf('%s - Read Visualization', locus$name), gp=gpar(fontsize=20))
	popViewport()
	
	pushViewport(viewport(layout=grid.layout(1,2, widths=c(14,2)), 
						  layout.pos.row=2, layout.pos.col=1))

	# Legend
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	grid.rect(x=unit(1/6,'npc'), y=unit(8/9, 'npc'), width=unit(.2/10, 'snpc'), 
			  height=unit(0.8/10, 'snpc'), gp=gpar(fill='#1b9e77', col='#1b9e77'))
	grid.rect(x=unit(1/6,'npc'), y=unit(3/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#f8766d', fill='#f8766d'))
	grid.rect(x=unit(1/6,'npc'), y=unit(2/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#00bfc4', fill='#00bfc4'))
	grid.text(label=c('Open Run', '', '', '', '', 'Haplotype 1', 'Haplotype 2', ''), 
			  x=unit(3/12, 'npc'), y=rev(seq(1:8)/9), just='left', gp=gpar(fontsize=12))
	popViewport()

	# Haplotype Colors
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	sideColorBar(colors, 3)
	grid.lines(x=c(1,1), y=c(0,1), gp=gpar(lwd=2))
	popViewport()

	popViewport()

	# Reads
	pushViewport(viewport(layout=grid.layout(1,length(peaks)), layout.pos.row=c(2), 
						  layout.pos.col=2))
	for (i in seq_along(peaks)){
		peak <- peaks[i]
		cat(sprintf('Region of Interest %s\n',i))
		pushViewport(viewport(layout.pos.row=1, layout.pos.col=i))
		plotReadMat(meMat.gc.f, GpC, meMat.cg.f, CpG_fil, open_runs, closed_runs, 
					snp.mat, snps, subset=select, order=order, gap=5, region=peak, 
					open=T, peaks=NULL)
		if (i > 1){
			grid.lines(x=c(0, 0), y=c(0, 1), gp=gpar(lwd=2.5))
		}
		popViewport()
		cat('\n')
	}
	popViewport()

	# Show context
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=2))
	context(locus, peaks, 5000)
	popViewport()

	# Gviz Tracks
	pushViewport(viewport(layout.pos.row = c(4), layout.pos.col=c(2)))
	chr <- as.character(seqnames(locus))
	from <- start(locus)
	to <- end(locus)

	meanTrack1 <- DataTrack(range=GpC, data=GpC$hap1, type='mountain', name='0h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)
	meanTrack2 <- DataTrack(range=GpC, data=GpC$hap2, type='mountain', name='24h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)

	peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Peaks', 
								 rot.title=0, cex.title=1, stacking='dense')

	sizes <- c(1, 2, 2, 1.5, 1.5, 1)
	plotTracks(c(peakCalls, meanTrack1, meanTrack2, txTr, gtrack, itrack), 
				 type='polygon', sizes=sizes, chromosome=chr, from=from, 
				 to=to, max.height=1, innerMargin=0, collapseTranscripts ="meta", 
				 transcriptAnnotation='symbol', add=T, panel.only=T, margin=0, title.width=1)

	popViewport()

	# Custom Titles for Gviz plots
	pushViewport(viewport(layout = grid.layout(length(sizes), 3, heights=c(sizes), widths=c(3*margin, width, margin)), 
											   layout.pos.row = c(4), layout.pos.col=c(1,2,3)))

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), 
			   gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), 
			   gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
	grid.text('ROIs',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
	grid.text('Haplotype 1\nGpC Methylation',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
	grid.text('Haplotype 2\nGpC Methylation',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=1))
	grid.text('Genes',  gp=gpar(fontsize=12))
	popViewport()

	popViewport()

	return()
}, rdata.path=rdata.path, hap1.names=hap1.names, hap2.names=hap2.names,
   peaks=peaks, pc_scores, itrack=itrack, gtrack=gtrack, txTr=txTr)

# CpG
all_loci.list <- lapply(idx, function(i, rdata.path, hap1.names, hap2.names,
									  peaks, pc_scores, itrack, gtrack, txTr){
	locus <- loci[i]
	GpC <- gpc.list[[i]]

	# Load pipeline outputs
	message(sprintf('Loading %s data', locus$name))
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% peaks]))
	select <- readCov >= 0.7

	hap1.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap1.names
	hap2.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap2.names

	order <- c(which(hap1.idx), which(hap2.idx))
	colors <- c(rep('#f8766d', sum(hap1.idx)), rep('#00bfc4', sum(hap2.idx)))
	# Consider margins out of 20
	# If right and bottom margins are equal, left margin equals 3X regular margin, 
	# top equals 1.5X, middle panel equals 1.5X as well
	tot <- 20
	margin <- 1
	width <- tot - (6*margin)
	height <- width / 2

	# # Test pattern
	# grid.newpage()

	# pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin,(tot-margin)*2/3,margin, (tot-margin)/3,margin), 
	# 										   widths=c(3*margin, width, margin))))
	# for(i in seq(5)){
	#     for(j in seq(3)){
	#         pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	#         grid.rect()
	#         grid.circle()
	#         popViewport()
	#     }
	# }

	grid.newpage()
	
	# Create layout
	pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin,(tot-margin)*2/3,margin, (tot-margin)/3,margin), 
											   widths=c(3*margin, width, margin))))
	
	# Title
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(sprintf('%s - Read Visualization', locus$name), gp=gpar(fontsize=20))
	popViewport()
	
	pushViewport(viewport(layout=grid.layout(1,2, widths=c(14,2)), 
						  layout.pos.row=2, layout.pos.col=1))

	# Legend
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	grid.points(x=unit(1/6,'npc'), y=unit(8/9, 'npc'), pch=21, 
				size=unit(0.9, 'char'), gp=gpar(fill='#0b529d', col='#0b529d'))
	grid.points(x=unit(1/6,'npc'), y=unit(7/9, 'npc'), pch=21, 
				size=unit(0.9, 'char'), gp=gpar(fill='#ed2127', col='#ed2127'))
	grid.rect(x=unit(1/6,'npc'), y=unit(3/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#f8766d', fill='#f8766d'))
	grid.rect(x=unit(1/6,'npc'), y=unit(2/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#00bfc4', fill='#00bfc4'))
	grid.text(label=c('Unmethylated CpG', 'Methylated CpG', '', '', '', 'Haplotype 1', 'Haplotype 2', ''), 
			  x=unit(3/12, 'npc'), y=rev(seq(1:8)/9), just='left', gp=gpar(fontsize=12))
	popViewport()

	# Haplotype Colors
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	sideColorBar(colors, 3)
	grid.lines(x=c(1,1), y=c(0,1), gp=gpar(lwd=2))
	popViewport()

	popViewport()

	# Reads
	pushViewport(viewport(layout=grid.layout(1,length(peaks)), layout.pos.row=c(2), 
						  layout.pos.col=2))
	for (i in seq_along(peaks)){
		peak <- peaks[i]
		cat(sprintf('Region of Interest %s\n',i))
		pushViewport(viewport(layout.pos.row=1, layout.pos.col=i))
		plotReadMat(meMat.gc.f, GpC, meMat.cg.f, CpG_fil, open_runs, closed_runs, 
					snp.mat, snps, subset=select, order=order, gap=5, region=peak, 
					cpg=T, peaks=NULL)
		if (i > 1){
			grid.lines(x=c(0, 0), y=c(0, 1), gp=gpar(lwd=2.5))
		}
		popViewport()
		cat('\n')
	}
	popViewport()

	# Show context
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=2))
	context(locus, peaks, 5000)
	popViewport()

	# Gviz Tracks
	pushViewport(viewport(layout.pos.row = c(4), layout.pos.col=c(2)))
	chr <- as.character(seqnames(locus))
	from <- start(locus)
	to <- end(locus)

	meanTrack1 <- DataTrack(range=GpC, data=GpC$hap1, type='mountain', name='0h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)
	meanTrack2 <- DataTrack(range=GpC, data=GpC$hap2, type='mountain', name='24h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)

	peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Peaks', 
								 rot.title=0, cex.title=1, stacking='dense')

	sizes <- c(1, 2, 2, 1.5, 1.5, 1)
	plotTracks(c(peakCalls, meanTrack1, meanTrack2, txTr, gtrack, itrack), 
				 type='polygon', sizes=sizes, chromosome=chr, from=from, 
				 to=to, max.height=1, innerMargin=0, collapseTranscripts ="meta", 
				 transcriptAnnotation='symbol', add=T, panel.only=T, margin=0, title.width=1)

	popViewport()

	# Custom Titles for Gviz plots
	pushViewport(viewport(layout = grid.layout(length(sizes), 3, heights=c(sizes), widths=c(3*margin, width, margin)), 
											   layout.pos.row = c(4), layout.pos.col=c(1,2,3)))

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), 
			   gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), 
			   gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
	grid.text('ROIs',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
	grid.text('Haplotype 1\nGpC Methylation',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
	grid.text('Haplotype 2\nGpC Methylation',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=1))
	grid.text('Genes',  gp=gpar(fontsize=12))
	popViewport()

	popViewport()

	return()
}, rdata.path=rdata.path, hap1.names=hap1.names, hap2.names=hap2.names,
   peaks=peaks, pc_scores, itrack=itrack, gtrack=gtrack, txTr=txTr)

# SNPs
all_loci.list <- lapply(idx, function(i, rdata.path, hap1.names, hap2.names,
									  peaks, pc_scores, itrack, gtrack, txTr){
	locus <- loci[i]
	GpC <- gpc.list[[i]]

	# Load pipeline outputs
	message(sprintf('Loading %s data', locus$name))
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	readCov <- rowMeans(!is.na(runMat.open[,tiles %over% peaks]))
	select <- readCov >= 0.7

	hap1.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap1.names
	hap2.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap2.names

	order <- c(which(hap1.idx), which(hap2.idx))
	colors <- c(rep('#f8766d', sum(hap1.idx)), rep('#00bfc4', sum(hap2.idx)))
	# Consider margins out of 20
	# If right and bottom margins are equal, left margin equals 3X regular margin, 
	# top equals 1.5X, middle panel equals 1.5X as well
	tot <- 20
	margin <- 1
	width <- tot - (6*margin)
	height <- width / 2

	# # Test pattern
	# grid.newpage()

	# pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin,(tot-margin)*2/3,margin, (tot-margin)/3,margin), 
	# 										   widths=c(3*margin, width, margin))))
	# for(i in seq(5)){
	#     for(j in seq(3)){
	#         pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
	#         grid.rect()
	#         grid.circle()
	#         popViewport()
	#     }
	# }

	grid.newpage()
	
	# Create layout
	pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin,(tot-margin)*2/3,margin, (tot-margin)/3,margin), 
											   widths=c(3*margin, width, margin))))
	
	# Title
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(sprintf('%s - Read Visualization', locus$name), gp=gpar(fontsize=20))
	popViewport()
	
	pushViewport(viewport(layout=grid.layout(1,2, widths=c(14,2)), 
						  layout.pos.row=2, layout.pos.col=1))

	# Legend
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	grid.rect(x=unit(1/6,'npc'), y=unit(8/9, 'npc'), width=unit(.2/10, 'snpc'), 
			  height=unit(0.8/10, 'snpc'), gp=gpar(fill='#b25a27', col='#b25a27'))
	grid.rect(x=unit(1/6,'npc'), y=unit(7/9, 'npc'), width=unit(.2/10, 'snpc'), 
			  height=unit(0.8/10, 'snpc'), gp=gpar(fill='#542d87', col='#542d87'))
	grid.rect(x=unit(1/6,'npc'), y=unit(3/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#f8766d', fill='#f8766d'))
	grid.rect(x=unit(1/6,'npc'), y=unit(2/9, 'npc'), width=unit(1/10, 'snpc'), 
			  height=unit(0.25/10, 'snpc'), gp=gpar(col='#00bfc4', fill='#00bfc4'))
	grid.text(label=c('Homozygous SNP', 'Heterozygous SNP', '', '', '', 'Haplotype 1', 'Haplotype 2', ''), 
			  x=unit(3/12, 'npc'), y=rev(seq(1:8)/9), just='left', gp=gpar(fontsize=12))
	popViewport()

	# Haplotype Colors
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	sideColorBar(colors, 3)
	grid.lines(x=c(1,1), y=c(0,1), gp=gpar(lwd=2))
	popViewport()

	popViewport()

	# Reads
	pushViewport(viewport(layout=grid.layout(1,length(peaks)), layout.pos.row=c(2), 
						  layout.pos.col=2))
	for (i in seq_along(peaks)){
		peak <- peaks[i]
		cat(sprintf('Region of Interest %s\n',i))
		pushViewport(viewport(layout.pos.row=1, layout.pos.col=i))
		plotReadMat(meMat.gc.f, GpC, meMat.cg.f, CpG_fil, open_runs, closed_runs, 
					snp.mat, snps, subset=select, order=order, gap=5, region=peak, 
					snp=T, peaks=NULL)
		if (i > 1){
			grid.lines(x=c(0, 0), y=c(0, 1), gp=gpar(lwd=2.5))
		}
		popViewport()
		cat('\n')
	}
	popViewport()

	# Show context
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=2))
	context(locus, peaks, 5000)
	popViewport()

	# Gviz Tracks
	pushViewport(viewport(layout.pos.row = c(4), layout.pos.col=c(2)))
	chr <- as.character(seqnames(locus))
	from <- start(locus)
	to <- end(locus)

	meanTrack1 <- DataTrack(range=GpC, data=GpC$hap1, type='mountain', name='0h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)
	meanTrack2 <- DataTrack(range=GpC, data=GpC$hap2, type='mountain', name='24h', 
							cex.title=1, cex.axis=1, type='l', ylim=c(0,1), size=3)

	peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Peaks', 
								 rot.title=0, cex.title=1, stacking='dense')

	sizes <- c(1, 2, 2, 1.5, 1.5, 1)
	plotTracks(c(peakCalls, meanTrack1, meanTrack2, txTr, gtrack, itrack), 
				 type='polygon', sizes=sizes, chromosome=chr, from=from, 
				 to=to, max.height=1, innerMargin=0, collapseTranscripts ="meta", 
				 transcriptAnnotation='symbol', add=T, panel.only=T, margin=0, title.width=1)

	popViewport()

	# Custom Titles for Gviz plots
	pushViewport(viewport(layout = grid.layout(length(sizes), 3, heights=c(sizes), widths=c(3*margin, width, margin)), 
											   layout.pos.row = c(4), layout.pos.col=c(1,2,3)))

	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), 
			   gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=c(2)))
	grid.rect(gp=gpar(fill=NA,col='#A9A9A9'))
	grid.yaxis(at=seq(0,1,l=5), label=c(' ', '0.25', '0.5', '0.75', ' '), 
			   gp=gpar(fontsize=10,col='#A9A9A9'))
	popViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
	grid.text('ROIs',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
	grid.text('Haplotype 1\nGpC Methylation',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
	grid.text('Haplotype 2\nGpC Methylation',  gp=gpar(fontsize=12))
	popViewport()
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col=1))
	grid.text('Genes',  gp=gpar(fontsize=12))
	popViewport()

	popViewport()

	return()
}, rdata.path=rdata.path, hap1.names=hap1.names, hap2.names=hap2.names,
   peaks=peaks, pc_scores, itrack=itrack, gtrack=gtrack, txTr=txTr)
dev.off()