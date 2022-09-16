# Call this script with loci BED filename (args[1]),
# CpG BED file (args[2]) and GpC BED file (args[3])

args <- commandArgs(trailingOnly=F)

# Infer scripts director
scripts.dir <- dirname(gsub('--file=', '', args[4]))
utils.dir <- paste0(scripts.dir, '/../../util/')

args <- commandArgs(trailingOnly=TRUE)

# Load libraries and functions
source(paste0(utils.dir, 'helper_functions.r'))
source(paste0(utils.dir, 'pipeline_functions.r'))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("Gviz"))
suppressMessages(library("org.Hs.eg.db"))

genome <- BSgenome.Hsapiens.UCSC.hg38

# References
gcg.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/gcg_strands.rds'
txdb.path <- '/seq/epiprod02/kdong/references/txdb/hg38.knownGene'
ideo.path <- '/seq/epiprod02/kdong/references/gviz/cytoBandIdeo.txt'

# Load files and write RDS files
gcg <- readRDS(gcg.path)
loci <- import.bed(args[1])
cpg.bed <- args[2]
cpg.out <- gsub("_cpg_methylation_reads.bed", "_hcg_reduc.rds", cpg.bed)
gpc.bed <- args[3]
gpc.out <- gsub("_gpc_methylation_reads.bed", "_gch_reduc.rds", gpc.bed)
prefix <- gsub("_cpg_methylation_reads.bed", "", cpg.bed)

message("Importing CpG files")
cpg <- readMethyl(cpg.bed, 'CG')
hcg <- exclude_gcg(cpg[cpg %over% loci], gcg, 'CG')

# Remove for memory
rm(cpg)
saveRDS(hcg, file=cpg.out, version=2)

message("Importing GpC files...")
gpc <- readMethyl(gpc.bed, 'GC')
gch <- exclude_gcg(gpc[gpc %over% loci], gcg, 'GC')

# Remove for memory
rm(gpc)
saveRDS(gch, file=gpc.out, version=2)

# Save 
short_c <- 80
long_c <- 80
linker <- 80

message("Saving methylation matrices...")
for(i in seq(length(loci))){
	message(loci$name[i])

	CpG <- diNT(loci[i], 'CG')
	GpC <- diNT(loci[i], 'GC')

	GCG <- CpG[CpG %over% GpC]
	CGC <- findCGC(loci[i])

	meMat.cg <- meMat(loci[i], hcg, CpG)
	meMat.gc <- meMat(loci[i], gch, GpC)

	groups <- findGroups(loci[i])

	CpG_mer <- groups[groups %over% CpG]
	GpC_mer <- shift(groups[shift(groups,1) %over% GpC],1)

	# Tally GCG
	CpG_mer$tot <- tallyHits(CpG_mer, CpG)

	CpG_mer$gcg_p <- tallyHits(CpG_mer, GCG)
	CpG_mer$gcg_p[is.na(CpG_mer$gcg_p)] <- 0

	CpG_mer$gcg_m <- tallyHits(CpG_mer, resize(CGC, width=1, fix='start'))
	CpG_mer$gcg_m[is.na(CpG_mer$gcg_m)] <- 0

	GpC_mer$tot <- tallyHits(GpC_mer, GpC)

	GpC_mer$gcg_p <- tallyHits(GpC_mer, GCG)
	GpC_mer$gcg_p[is.na(GpC_mer$gcg_p)] <- 0

	GpC_mer$gcg_m <- tallyHits(GpC_mer, resize(CGC, width=1, fix='start'))
	GpC_mer$gcg_m[is.na(GpC_mer$gcg_m)] <- 0
	
	CpG_mer$p <- (CpG_mer$gcg_p / CpG_mer$tot) < 0.5
	CpG_mer$m <- (CpG_mer$gcg_m / CpG_mer$tot) < 0.5
	GpC_mer$p <- (GpC_mer$gcg_p / GpC_mer$tot) < 0.5
	GpC_mer$m <- (GpC_mer$gcg_m / GpC_mer$tot) < 0.5
	CpG_fil <- CpG_mer[(CpG_mer$gcg_p / CpG_mer$tot) < 0.5 | 
					   (CpG_mer$gcg_m / CpG_mer$tot) < 0.5]
	GpC_fil <- GpC_mer[(GpC_mer$gcg_p / GpC_mer$tot) < 0.5 | 
					   (GpC_mer$gcg_m / GpC_mer$tot) < 0.5]
	
	meMat.cg.f <- mergeMat(meMat.cg, CpG, CpG_fil, strands=T)
	meMat.gc.f <- mergeMat(meMat.gc[dev(meMat.gc) < 0.3 & rowMeans(meMat.gc, na.rm=T) > 0.02,], 
						   GpC, GpC_fil, strands=T)
	
	if(nrow(meMat.gc.f) == 0){
		message("SKIP")
		next
	}
	
	# Determine closed runs in NA aware manner
	closed_runs <- apply(meMat.gc.f, 1, function(row){
		nas <- findRuns(as.logical(row), val=NA)
		closed <- findRuns(as.logical(row), val=0)
		if(nrow(closed) == 0){
			return(data.frame())
		}
		stitched <- stitchRuns(closed, nas)
		stitchedDist(stitched, GpC_fil)  
	})

	# Infer open runs from closed runs
	open_runs <- lapply(closed_runs, function(closed_run){
		if(nrow(closed_run) > 1){
			starts <- closed_run$ends[-nrow(closed_run)]
			starts[starts %% 1 == 0] <- starts[starts %% 1 == 0] + 1
			ends <- closed_run$starts[-1]
			ends[ends %% 1 == 0] <- ends[ends %% 1 == 0] - 1
			temp.df <- data.frame(starts=starts, ends=ends)
			minDist <- runDistMin(temp.df, GpC_fil)
			maxDist <- runDistMax(temp.df, GpC_fil)
			cov <- ends - starts
			start <- closed_run$end[-nrow(closed_run)]
			end <- closed_run$start[-1]
			dist <- end-start + 1
			left <- closed_run$dist[-nrow(closed_run)]
			right <- closed_run$dist[-1]
			
			df <- data.frame(starts=starts, ends=ends, cov=cov, start=start, end=end, 
							 dist=dist, minDist=minDist, maxDist=maxDist, left=left, 
							 right=right)
			return(df)
		}
		return(data.frame())
	})

	save(CpG_fil, GpC_fil, meMat.cg.f, meMat.gc.f, closed_runs, open_runs, 
		 file=sprintf('rdata/%s-run.RData', loci$name[i]), version=2)
}

message("Saving run matrices")
for(i in seq(length(loci))){
	locus <- loci[i]
	message(locus$name)
	load(sprintf('rdata/%s-run.RData', locus$name))
	
	open_run.qc <- sapply(open_runs, function(el){sum(el$dist > 500)})

	# Create GRanges from run data frames
	closed.df <- do.call(rbind, closed_runs[open_run.qc < 2])

	closed.gr <- GRanges(seqnames(locus), IRanges(closed.df$start, closed.df$end))
	closed.gr$dist <- closed.df$dist
	closed.gr$minDist <- closed.df$minDist
	closed.gr$maxDist <- closed.df$maxDist

	open.df <- do.call(rbind, open_runs[open_run.qc < 2])

	open.gr <- GRanges(seqnames(locus), IRanges(open.df$start, open.df$end))
	open.gr$dist <- open.df$dist
	open.gr$minDist <- open.df$minDist
	open.gr$maxDist <- open.df$maxDist

	open.gr$left <- open.df$left
	open.gr$right <- open.df$right

	# Tiles
	tiles <- tile(locus, width=5)[[1]]
	tiles <- resize(tiles, width=1, fix='center')

	# Calculate total run coverage
	tiles$tot <- tallyHits(tiles, c(open.gr, closed.gr, ignore.mcols=T))
	
	# Calculate filtered open run coverage
	tiles$open <- tallyHits(tiles, open.gr[!(open.gr$dist < linker & 
											 open.gr$right > long_c & 
											 open.gr$left > long_c)])
	
	tiles$closed <- tallyHits(tiles, closed.gr[closed.gr$dist <= short_c])

	runMat.open <- runMat(meMat.gc.f, GpC_fil, open_runs, locus, 5, 'open')
	runMat.linker <- runMat(meMat.gc.f, GpC_fil, open_runs, locus, 5, 'linker')

	runMat.short <- runMat(meMat.gc.f, GpC_fil, closed_runs, locus, 5, 'short')
	runMat.mono <- runMat(meMat.gc.f, GpC_fil, closed_runs, locus, 5, 'mono')
	runMat.di <- runMat(meMat.gc.f, GpC_fil, closed_runs, locus, 5, 'di')
	runMat.tri <- runMat(meMat.gc.f, GpC_fil, closed_runs, locus, 5, 'tri')

	save(runMat.open, runMat.linker, runMat.short, runMat.mono, runMat.di, runMat.tri, 
		 tiles, open_run.qc, file=sprintf('rdata/%s-runMat.RData', locus$name), version=2)
}

# Generate BigWig Tracks
message("Generating BigWigs...")

write_bw(loci, 'open', prefix, smooth=F)
write_bw(loci, 'open', prefix,  smooth=T)
write_bw(loci, 'short', prefix,  smooth=F)

# Peak Calling
message("Calling peaks...")
open_tiles.list <- sapply(seq_along(loci), function(i){
	locus <- loci[i]

	message(locus$name)
	load(sprintf('rdata/%s-runMat.RData', locus$name))

	return(tiles$open/tiles$tot)
})
open_tiles <- do.call(c, open_tiles.list)
open_mean <- mean(open_tiles, na.rm=T)
open_sd <- sd(open_tiles, na.rm=T)

pdf('peak_calling.pdf', width=24, height=16)
txdb <- loadDb(txdb.path)
open_peak.list <- lapply(seq_along(loci), function(i){
	locus <- loci[i]
	message(locus$name)
	load(sprintf('rdata/%s-run.RData', locus$name))
	load(sprintf('rdata/%s-runMat.RData', locus$name))

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
	txTr <- GeneRegionTrack(txdb, chromosome=chr, start=from, end=to, cex.title=2, 
							rot.title=0, name='Genes', cex.group=1.2,
							collapseTranscripts ="meta", transcriptAnnotation='symbol')
	
	if(length(gene(txTr)) > 0){
		symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
		symbol(txTr) <- symbols[gene(txTr)]
	}
	
	tiles$cov <- tiles$tot 
	openPeaks1 <- findPeaks(tiles$open/tiles$tot, tiles, thresh=30, restrict=locus, 
							smoothing=T, mean=open_mean, sd=open_sd, erosion=4, 
							dilation=4, plot=F, chr=chr, from=from, to=to, 
							itrack=itrack, gtrack=gtrack, txTr=txTr)
	if (length(openPeaks1) == 0){
		message('No Peaks')
		return(openPeaks1)
	}
	
	openPeaks2 <- findPeaks(tiles$open/tiles$tot, tiles, thresh=0, restrict=openPeaks1, 
							smoothing=F, mean=open_mean, sd=open_sd, erosion=2, 
							dilation=2, plot=T, chr=chr, from=from, to=to, 
							itrack=itrack, gtrack=gtrack, txTr=txTr)
	
	openPeaks2$tot <- tallyHits(openPeaks2, GpC_fil)

	return(openPeaks2[which(openPeaks2$tot > 1)])
})
dev.off()

open_peaks <- suppressWarnings(do.call(c, open_peak.list))
export.bed(open_peaks, paste0(prefix, '_open_peaks.bed'))

# Closed Peak Calling
message("Calling closed peaks...")
closed_tiles.list <- sapply(seq_along(loci), function(i){
	locus <- loci[i]

	message(locus$name)
	load(sprintf('rdata/%s-runMat.RData', locus$name))

	return(tiles$closed/tiles$tot)
})
closed_tiles <- do.call(c, closed_tiles.list)
closed_mean <- mean(closed_tiles, na.rm=T)
closed_sd <- sd(closed_tiles, na.rm=T)

pdf('peak_calling_closed.pdf', width=24, height=16)
closed_peak.list <- lapply(seq_along(loci), function(i){
	locus <- loci[i]
	message(locus$name)
	load(sprintf('rdata/%s-run.RData', locus$name))
	load(sprintf('rdata/%s-runMat.RData', locus$name))

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
	txTr <- GeneRegionTrack(txdb, chromosome=chr, start=from, end=to, cex.title=2, rot.title=0, name='Genes', cex.group=1.2,
		collapseTranscripts ="meta", transcriptAnnotation='symbol')
	
	if(length(gene(txTr)) > 0){
		symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
		symbol(txTr) <- symbols[gene(txTr)]
	}
	tiles$cov <- tiles$tot
	closedPeaks1 <- findPeaks(tiles$closed/tiles$tot, tiles, thresh=0, restrict=locus, 
							  smoothing=F, mean=closed_mean*2, sd=closed_sd*2, erosion=1, 
							  dilation=1, plot=T, chr=chr, from=from, to=to, 
							itrack=itrack, gtrack=gtrack, txTr=txTr)

	return(closedPeaks1)
})
dev.off()

closed_peaks <- suppressWarnings(do.call(c, closed_peak.list))
export.bed(closed_peaks, paste0(prefix, '_closed_peaks.bed'))