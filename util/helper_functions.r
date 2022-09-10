#******************************************************************************
# Collection of helper functions commonly used for analysis
#******************************************************************************

chopStr <- function(x, n){
	# Extracts the end of a string
	#
	# Args:
	#	x: string to subset
	#	n: how many characters to skip
	#
	# Returns:
	#	Subsetted string
	substr(x, n + 1, nchar(x))
}

loadCpG <- function(loci, rdata.path, hap1=NULL, hap2=NULL){
	# Loads RData files generated by pipeline for the specified loci
	#
	# Args:
	#	loci: a GRanges object of named loci
	#	rdata.path: the path to the directory of RData files
	#	hap1: optional array of read names assigned to haplotype 1
	#	hap2: optional array of read names assigned to haplotype 2
	#
	# Returns:
	#	List of GRanges objects with CpG positions
	#	Average methylation and coverage are provided as metadata columns
	list <- lapply(seq_along(loci), function(i, hap1.names, hap2.names){
		locus <- loci[i]
		message(sprintf('Loading %s data', locus$name))
		load(sprintf('%s/%s-run.RData', rdata.path, locus$name))

		if (!is.null(hap1) & !is.null(hap2)){
			hap1.idx <- chopStr(rownames(meMat.cg.f), 2) %in% hap1.names
			hap2.idx <- chopStr(rownames(meMat.cg.f), 2) %in% hap2.names

			CpG_fil$hap1 <- colMeans(meMat.cg.f[hap1.idx,], na.rm=T)
			CpG_fil$hap2 <- colMeans(meMat.cg.f[hap2.idx,], na.rm=T)
		}

		CpG_fil$avg <- colMeans(meMat.cg.f, na.rm=T)
		CpG_fil$cov <- colSums(is.na(meMat.cg.f))
		return(CpG_fil)
	}, hap1.names=hap1, hap2.names=hap2)
	return(list)
}

loadGpC <- function(loci, rdata.path, hap1=NULL, hap2=NULL){
	# Loads RData files generated by pipeline for the specified loci
	#
	# Args:
	#	loci: a GRanges object of named loci
	#	rdata.path: the path to the directory of RData files
	#	hap1: optional array of read names assigned to haplotype 1
	#	hap2: optional array of read names assigned to haplotype 2
	#
	# Returns:
	#	List of GRanges objects with GpC positions
	#	Average methylation and coverage are provided as metadata columns
	list <- lapply(seq_along(loci), function(i, hap1.names, hap2.names){
		locus <- loci[i]
		message(sprintf('Loading %s data', locus$name))
		load(sprintf('%s/%s-run.RData', rdata.path, locus$name))

		if (!is.null(hap1) & !is.null(hap2)){
			hap1.idx <- chopStr(rownames(meMat.gc.f), 2) %in% hap1.names
			hap2.idx <- chopStr(rownames(meMat.gc.f), 2) %in% hap2.names

			GpC_fil$hap1 <- colMeans(meMat.gc.f[hap1.idx,], na.rm=T)
			GpC_fil$hap2 <- colMeans(meMat.gc.f[hap2.idx,], na.rm=T)
		}

		GpC_fil$avg <- colMeans(meMat.gc.f, na.rm=T)
		GpC_fil$cov <- colSums(is.na(meMat.gc.f))
		return(GpC_fil)
	}, hap1.names=hap1, hap2.names=hap2)
	return(list)
}

mergeMat <- function(reads, diNT, mergeNT, strands=F){
	# Merges columns of methylation matrix according to GRanges objects
	#
	# Args:
	#	reads: methylation matrix (read x position)
	#	diNT: GRanges object corresponding to columns of reads
	#	mergeNT: GRanges object corresponding to desired output positions
	#	strands: a boolean indicating strand-awareness
	#
	# Returns:
	#	A matrix with aggregated columns (equal in length to mergeNT)

	# Get overlaps between diNT and mergeNT
	hits <- findOverlaps(diNT, mergeNT)

	ids <- unique(subjectHits(hits))

	# Initialize output matrix
	mat <- matrix(NA, nrow=nrow(reads), ncol=length(mergeNT))

	# Aggregate the columns
	for (i in seq_along(ids)){
		cols <- queryHits(hits)[subjectHits(hits) == ids[i]]
		if (length(cols) == 1){
			mat[,ids[i]] <- reads[,cols]
		} else {
			mat[,ids[i]] <- rowMeans(reads[,cols], na.rm=T)
		}
	}

	if (strands){
		# Get strand information from read names
		read_strands <- sapply(rownames(reads), function(name){
			ifelse(substr(name, 1, 1)=='P', '+' ,'-')
		})

		# Set NAs if mergeNT requests to ignore the strand
		mat[read_strands=='+',!mergeNT$p] <- NA
		mat[read_strands=='-',!mergeNT$m] <- NA
	}
	
	# Name output matrix rows and columns
	colnames(mat) <- (start(mergeNT) + end(mergeNT)) / 2
	rownames(mat) <- rownames(reads)

	return(mat)
}

findPeaks <- function(summary, diNT, thresh, restrict, smoothing=T, 
					  mult=2, smooth_mult=1, mean=0.15, sd=0.1, 
					  erosion=4, dilation=4, plot=T, chr=NULL, from=NULL,
					  to=NULL, itrack=NULL, gtrack=NULL, txTr=NULL){
	# Identifies peaks based on scores tied to GRanges object
	#
	# Args:
	#	summary: an array of scores equal in length to diNT
	#	diNT: a GRanges object of genomic positions
	#	thresh: number of basepairs to resize diNT for smoothing (will be doubled)
	#	restrict: a GRanges object limiting where a peak can be called
	# 	smoothing: a boolean determining whether or not to use smoothing
	#	mult: upper cutoff (will be multiplied by standard deviation)
	#	smooth_mult: lower cutoff for smoothing (will be multiplied by standard deviation)
	#	mean: mean value of scores
	#	sd: standard deviation of scores
	#	erosion: how many indices to erode the putative peaks by
	#	dilation: how many indices to dilate the putative peaks by
	#	plot: boolean indicating whether or not to plot the scores and called peaks
	#	chr: if plotting, the chromosome name of the locus to plot
	#	from: if plotting, the start position of the locus to plot
	#	to: if plotting, the end position of the locus to plot
	#	itack: if plotting, the ideogram track to be plotted
	#	gtrack: if plotting, the genomix axis track to be plotted
	#	txTr: if plotting, the gene region track to be plotted
	#
	# Returns:
	# 	GRanges object of identified peaks
	#	If specified, Gviz plot showing the identified peaks

	# Smooth an array of scores tied to GRanges object
	.smooth <- function(means, diNT, thresh){
		# Create smoothing windows and find overlaps
		windows <- resize(diNT, width=2*thresh, fix='center')
		diNT <- resize(diNT, width=1, fix='center')
		hits <- findOverlaps(diNT, windows)

		# Average score per window
		avg.df <- aggregate(means[queryHits(hits)], by=list(subjectHits(hits)), FUN=mean, na.rm=T)
		colnames(avg.df) <- c('new_idx', 'x')

		windows$score <- NA
		windows$score[avg.df$new_idx] <- avg.df$x 

		# Average score across windows
		rev_avg.df <- aggregate(windows$score[subjectHits(hits)], by=list(queryHits(hits)), FUN=mean, na.rm=T)
		colnames(rev_avg.df) <- c('new_idx', 'x')
		rev_avg.df$cov <- aggregate(rep(1,length(hits)), by=list(queryHits(hits)), FUN=sum)$x
		diNT$score <- NA
		diNT$score[rev_avg.df$new_idx] <- rev_avg.df$x
		diNT$cov <- 1
		diNT$cov[rev_avg.df$new_idx] <- rev_avg.df$cov

		return(diNT)
	}
	
	# Identify positive runs in array of numbers
	.findContig <- function(bools, max){
		starts <- vector('integer')
		ends <- vector('integer')
		peak <- FALSE
		i <- 1

		while(i <= length(bools)){
			if (bools[i] == max){
				if(!peak){
					j <- i - 1
					# Search for indices below max but above zero
					while(j > 0){
						if(bools[j] > 0){
							j <- j - 1
							next
						} else {
							break
						}
					}
					# Declare start of peak
					starts <- c(starts, j + 1)
					peak <- TRUE
				} 
			} else if (bools[i] == 0){
				# Declare end of peak
				if(peak){
					peak <- FALSE
					ends <- c(ends, i - 1)
				} 
			} 
			i <- i + 1
		} 
		# End last peak if applicable
		if (peak){
			ends <- c(ends, i)
		}

		return(data.frame(starts=starts, ends=ends))
	}

	# Expand peak boundaries (by index)
	.dilate <- function(bools, max, num){
		result <- bools

		idx <- .findContig(bools, max)
		# Skip dilation if it exceeds the total number of indices
		idx <- idx[(idx$starts - num > 0) & (idx$ends + num < length(bools)),]

		for (i in seq(num)){
			result[idx$starts - i] <- TRUE
			result[idx$ends + i] <- TRUE
		}
		return(result)
	}
	
	# Shrink peak boundaries (by index)
	.erode <- function(bools, max, num){
		result <- bools

		idx <- .findContig(bools, max)
		# Skip erosion if it exceeds the total number of indices
		idx <- idx[(idx$starts - num > 0) & (idx$ends + num < length(bools)),]

		for (i in seq(num)){
			result[idx$starts - 1 + i] <- FALSE
			result[idx$ends - i + 1] <- FALSE
		}
		return(result)
	}

	if (smoothing){
		smoothed <- .smooth(summary, diNT, thresh)

		# Override baseline if more active
		if (median(smoothed$score, na.rm=T) > mean){
			mean <- median(smoothed$score, na.rm=T)
			sd <- sd(smoothed$score, na.rm=T)
		}

		max <- 2
		hi <- (smoothed$score > (mean + mult * sd)) + (smoothed$score > (mean + smooth_mult * sd))
		hi[is.na(hi)] <- 0

		# Dilate & erode then erode & dilate
		dil <- .erode(.dilate(hi, max, dilation), max, erosion)
		ero <- .dilate(.erode(dil, max, erosion), max, dilation)
	} else {
		smoothed <- diNT
		smoothed$score <- summary

		# Override baseline if more active
		if (median(smoothed$score, na.rm=T) > mean){
			mean <- median(smoothed$score, na.rm=T)
			sd <- sd(smoothed$score, na.rm=T)
		}
		
		max <- 1
		hi <- (smoothed$score > (mean + mult *sd))
		hi[is.na(hi)] <- 0

		# Erode then dilate
		ero <- .dilate(.erode(hi, max, erosion), max, dilation)
	}

	# Define peaks
	idx <- .findContig(ero, max)

	if(nrow(idx) == 0){
		message('No Peaks found')
		return(GRanges())
	}

	peaks <- GRanges(as.character(unique(seqnames(smoothed))), IRanges(start(smoothed[idx$starts]), end(smoothed[idx$ends])))
	peaks <- peaks[peaks %over% restrict]

	# # Find First and Second Derivatives
	# dist1 <- end(smoothed[-1]) - start(smoothed[-length(smoothed)])
	# dx <- shift(smoothed[-length(smoothed)], dist1 / 2)
	# d1 <- diff(smoothed$score) / dist
	
	# dist2 <- end(dx[-1]) - start(dx[-length(dx)])
	# dx2 <- shift(dx[-length(dx)], dist2 / 2)
	# d2 <- diff(d1) / dist2
	if(plot){
		meanTrack <- DataTrack(range=diNT, data=summary, name='Mean GCH\n Methylation', cex.title=2, cex.axis=1.2, type='l', ylim=c(0,1), size=3)
		smoothTrack <- DataTrack(range=smoothed[smoothed$cov > 1], data=smoothed$score[smoothed$cov > 1], name='Smoothed GCH\n Methylation', cex.title=2, cex.axis=1.2, type='l', ylim=c(0,1), size=3)

		peakCalls <- AnnotationTrack(peaks, chromosome=chr, start=from, end=to, name='Called Peaks', rot.title=0, cex.title=2)

		ht <- HighlightTrack(trackList=list(meanTrack, smoothTrack, peakCalls), range=peaks)
		plotTracks(list(itrack, gtrack, ht, txTr), chromosome=chr, from=from, to=to, max.height=1,
		collapseTranscripts ="meta", transcriptAnnotation='symbol')
	}

	return(peaks)
}

mergePeaks <- function(peak1, peak2){
	# Merge two sets of peaks together
	#
	# Args:
	#	peak1: Granges object of first set of peaks
	#	peak2: Granges object of second set of peaks
	#
	# Returns:
	#	GRanges object with the union of both peak sets
	hits <- suppressWarnings(findOverlaps(peak1, peak2))

	newStart <- pmin(start(peak2)[subjectHits(hits)], start(peak1)[queryHits(hits)])
	newEnd <- pmax(end(peak2)[subjectHits(hits)], end(peak1)[queryHits(hits)])
	
	result <- peak1

	start(result[queryHits(hits)]) <- newStart
	end(result[queryHits(hits)]) <- newEnd

	result <- c(result, peak2[-subjectHits(hits)], ignore.mcols=T)

	return(sort(result))
}

tallyHits <- function(gr1, gr2){
	# Tally overlaps of two sets of genomic ranges
	#
	# Args:
	#	gr1: Granges object of interest
	#	gr2: Granges object to count
	#
	# Returns:
	#	A vector of integers equal in length to gr1 counting hits of gr2
	hits <- findOverlaps(gr1, gr2)
	hits.df <- aggregate(rep(1, length(hits)), by=list(queryHits(hits)), FUN=sum)
	result <- rep(NA, length(gr1))
	result[hits.df$Group.1] <- hits.df$x
	return(result)
}

filterClosedPeaks <- function(closedPeaks, openPeaks, loci, peakDist=10, lociDist=100){
	# Filters closed peaks
	#
	# Args:
	#	closedPeaks: Granges object of closed run peaks
	#	openPeaks: Granges object of open run peaks
	#	loci: Granges object of targeted loci
	#	peakDist: distance to next closest open run peak
	#	lociDist: minimum distance to loci boundary
	#
	# Returns:
	#	GRanges object of closed run peaks that are flanked by open runs and
	#	far from the loci boundaries

	# Function to check flanks of peaks
	.checkFlank <- function(peak1, peak2, dist=10){
		left <- shift(resize(peak1, width=1, fix='start'), -dist)
		right <- shift(resize(peak1, width=1, fix='end'), dist)

		peak1$left <- left %over% peak2
		peak1$right <- right %over% peak2

		return(peak1)
	}

	buffer <- .checkFlank(closedPeaks, openPeaks, dist=peakDist)
	buffer2 <- .checkFlank(closedPeaks, loci, dist=lociDist)

	return(buffer[(buffer$left | buffer$right) & (buffer2$left & buffer2$right)])
}

smooth <- function(x, n=5){
	# Smooth a vector with sliding windos
	#
	# Args:
	#	x: a numeric vector to smooth
	#	n: an odd integer specifying the total positions to smooth over
	#
	# Returns:
	#	A numeric vector that has been averaged with (n-1)/2 positions 
	#	to the left and to the right

	stats::filter(x, rep(1/n, n), sides=2)
}