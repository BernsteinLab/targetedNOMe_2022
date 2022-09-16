#******************************************************************************
# Collection of functions commonly used in the pipeline
#******************************************************************************

## Import ##
readMethyl <- function(filename, diNT){
	# Imports the 'methylation_reads.bed' file which is split by read
	#
	# Args:
	#	filename: string path to 'methylation_read.bed' file
	#	diNT: string type of methylation ('CG' or 'GC')
	#
	# Returns:
	#	A GRanges object of methylation calls
	#	Names capture read and strand information

	fix_start <- ifelse(diNT == 'CG', TRUE, FALSE)
	
	buffer <- import.bed(filename)
	buffer[strand(buffer)=='-'] <- flank(buffer[strand(buffer)=='-'], 1, start=fix_start)
	buffer$name <- mapply(function(strand, name){
		return(sprintf('%s-%s', switch(strand, '+'='P', '-'='M'), name))
	}, as.character(strand(buffer)), buffer$name)
	
	return(buffer)
}

## Nucleotides ##
exclude_gcg <- function(gr, gcg, diNT){
	# Removes any dinucleotides that overlap with GCGs
	#
	# Args:
	#	gr: GRanges object of interest
	#	gcg: GRanges object containing positions of all GCGs
	#
	# Returns:
	#	An unstranded GRanges object with GCGs removed
	fix_start <- ifelse(diNT == 'CG', FALSE, TRUE)
	
	gcg_gr <- findOverlaps(gcg, gr)
	gr.f <- gr[-subjectHits(gcg_gr)]
	gr.f[strand(gr.f)=='-'] <- flank(gr.f[strand(gr.f)=='-'], 1, start=fix_start)
	return(unstrand(gr.f))
}

diNT <- function(locus, diNT){
	# Identify dinucleotides in locus of interest
	#
	# Args:
	#	locus: GRanges object of locus
	#	diNT: string type of methylation ('CG' or 'GC')
	#
	# Returns:
	#	A GRanges object with all positions of diNT within the locus
	ntSeq <- getSeq(genome, locus)
	matches <- vmatchPattern(diNT, ntSeq, fixed='subject')[[1]]
	fix <- ifelse(diNT == 'CG', 'start', 'end')
	mC <- GRanges(seqnames(locus), shift(resize(matches, width=1, fix=fix), start(locus)-1))
	return(mC)
}

findCGC <- function(locus){
	# Identify CGCs (reverse complement of GCG) in locus of interest
	#
	# Args:
	#	locus: GRanges object of locus
	#
	# Returns:
	#	A GRanges object with all positions of CGC within the locus
	#	on the positive strand
	ntSeq <- getSeq(genome, locus)
	matches <- vmatchPattern('CGC', ntSeq, fixed='subject')[[1]]
	mC <- GRanges(seqnames(locus), shift(matches, start(locus)-1))
	return(mC)
}

mergeNT <- function(diNT, thresh=5){
	# Merge together nearby dinucleotides
	#
	# Args:
	#	diNT: GRanges object of dinucleotides of interest
	#	thresh: integer specifying maximum base pair distance
	#
	# Returns:
	#	A GRanges object with ranges merged that are within
	#	the threshold distance.
	dist <- start(diNT)[-1] - end(diNT)[-length(diNT)]
	merge_idx <- which(dist <= thresh & dist > 0)

	result <- diNT

	for(i in merge_idx){
	  start(result[i+1]) <- start(result[i])
	}

	result <- result[-merge_idx]

	return(result)
}

findGroups <- function(locus, thresh=5){
	# Group together nearby dinucleotides
	#
	# Args:
	#	locus: GRanges object of locus
	#	thresh: integer specifying maximum base pair distance
	#
	# Returns:
	#	A GRanges object of all CpG and GpC positionswith ranges merged 
	#	that are within the threshold distance.
	ntSeq <- getSeq(genome, locus)
	CG_hits <- vmatchPattern('CG', ntSeq, fixed='subject')[[1]]
	GC_hits <- vmatchPattern('GC', ntSeq, fixed='subject')[[1]]
	CG <- GRanges(seqnames(locus), shift(resize(CG_hits, width=1, fix='start'), start(locus)-1))
	GC <- GRanges(seqnames(locus), shift(resize(GC_hits, width=1, fix='start'), start(locus)-1))

	mC <- mergeNT(sort(c(CG, GC)), thresh)
	return(mC)
}

## Methylation Matrices ##
meMat <- function(locus, calls, diNT){
	# Create matrix of read methylation values
	#
	# Args:
	#	locus: GRanges object of locus
	#	calls: GRanges object of methylation calls
	#	diNT: GRanges object of dinucleotides of interest
	#
	# Returns:
	#	A binary matrix (with NAs) where each row is a read
	#	and each column is a cytosine position within the locus.
	#	1 mean methylated, 0 means unmethylated, and 
	#	NA means undetermined or missing

	# Find overlaps with locus
	overlaps <- findOverlaps(calls, locus)
	buffer <- calls[queryHits(overlaps)]

	# Identify unique reads
	reads <- unique(buffer$name)
	
	# # Print some interesting statistics
	# cat(sprintf('Locus: %s\n', loci$name[i]))
	# cat(sprintf('Loci size: %s Kb, Unique GCHs: %s\n', width(loci[i])/ 1000, length(mC)))
	# cat(sprintf('Unique reads: %s, Average GCH coverage: %.2f\n', length(reads), length(overlaps) / length(mC)))

	# Initialize empty matrix for results  
	meMat <- matrix(data=NA, nrow=length(reads), ncol=length(diNT)) 
	rownames(meMat) <- reads
	colnames(meMat) <- start(diNT)

	# Input methylation values in matrix
	for(i in seq(length(reads))){
		read <- buffer[buffer$name == reads[i]]
		matches <- findOverlaps(read, diNT)
		if (length(matches) > 0){
		meMat[i, subjectHits(matches)] <- read$score[queryHits(matches)] 
		}
	}

	row_id <- order(rowSums(!is.na(meMat)), decreasing=T)
	return(meMat[row_id,])
}

dev <- function(meMat){
	# Calculate row-wise deviation from mean
	#
	# Args:
	#	meMat: matrix of methylation values
	#
	# Returns:
	#	A numeric vector equal in length to the number of rows
	#	with average squared deviation values per row.
	
	exp <- colMeans(meMat, na.rm=T)
	
	# Subtract mean from each row and get average squared deviation
	dev <- rowMeans((sweep(meMat, 2, exp, '-'))^2, na.rm=T)

	return(dev)
}

## Runs ##
findRuns <- function(bools, val){
	# Identifies runs of value of interest in boolean vector
	#
	# Args:
	#	bools: vector of booleans (can include NAs)
	#	val: value to generate runs for
	#
	# Returns:
	#	A data frame of start and stop indices for each run

	# Identify contiguous segments
	starts <- vector('integer')
	ends <- vector('integer')
	peak <- FALSE
	
	# Identify length of read
	call_idx <- which(!is.na(bools))
	i <- min(call_idx)

	# Replace NA
	bools[is.na(bools)] <- 99
	if(is.na(val)){ val <- 99 }

	while(i < max(call_idx)){
		if (bools[i] == val){
			if(!peak){
				starts <- c(starts, i)
				peak <- TRUE
			} 
		} else {
			if(peak){
				peak <- FALSE
				ends <- c(ends, i-1)
			} 
		} 
		i <- i + 1
	} 
	if (peak){
		ends <- c(ends, i)
	}

	return(data.frame(starts=starts, ends=ends))
}

runDistMax <- function(contigs, diNT){
	# Calculates maximum possible run length in base pairs
	#
	# Args:
	#	contigs: data frame of merged closed run start and stop indices
	#	diNT: GRanges object of GpC dinucleotides
	#
	# Returns:
	#	A numeric vector equal in length to the rows of contigs.
	#	Each value corresponds to the maximum length of the respective run.

	# Hack to account for first and last index
	contigs$starts[contigs$starts == 1] <- 2
	contigs$ends[contigs$ends == length(diNT)] <- length(diNT) - 1
	
	# Hack to account for intermediate idx
	int_starts <- contigs$starts %% 1 != 0
	int_ends <- contigs$ends %% 1 != 0

	contigs$start <- NA
	contigs$end <- NA

	contigs$start[!int_starts] <- end(diNT[contigs$starts[!int_starts]-1]) + 1
	contigs$start[int_starts] <- (start(diNT[floor(contigs$starts[int_starts])]) + 
								  end(diNT[floor(contigs$starts[int_starts])])) / 2
	contigs$end[!int_ends] <- start(diNT[contigs$ends[!int_ends]+1]) - 1
	contigs$end[int_ends] <- (start(diNT[ceiling(contigs$ends[int_ends])]) + 
							  end(diNT[ceiling(contigs$ends[int_ends])])) / 2

	return(contigs$end - contigs$start + 1)
}

runDistMin <- function(contigs, diNT){
	# Calculates minimum possible run length in base pairs
	#
	# Args:
	#	contigs: data frame of merged closed run start and stop indices
	#	diNT: GRanges object of GpC dinucleotides
	#
	# Returns:
	#	A numeric vector equal in length to the rows of contigs.
	#	Each value corresponds to the minimum length of the respective run.

	# Hack to account for intermediate idx
	int_starts <- contigs$starts %% 1 != 0
	int_ends <- contigs$ends %% 1 != 0

	contigs$start <- NA
	contigs$end <- NA

	contigs$start[!int_starts] <- start(diNT[contigs$starts[!int_starts]])
	contigs$start[int_starts] <- round((start(diNT[floor(contigs$starts[int_starts])]) + 
										end(diNT[floor(contigs$starts[int_starts])])) / 2)
	contigs$end[!int_ends] <- end(diNT[contigs$ends[!int_ends]])
	contigs$end[int_ends] <- round((start(diNT[ceiling(contigs$ends[int_ends])]) + 
									end(diNT[ceiling(contigs$ends[int_ends])])) / 2)

	return(contigs$end - contigs$start + 1)
}

stitchRuns <- function(closed, nas){
	# Combines two data frames of runs
	#
	# Args:
	#	closed: data frame of closed run starts and stops
	#	nas: data frame of NA run starts and stops
	#
	# Returns:
	#	A merged data frame of start and stop indices for closed runs.
	#	NA runs are fully incorporated if flanked by closed runs.
	#	Otherwise only half of the NA run is incorporated if flanked by 
	#	only one closed run.

	if(nrow(nas) > 0){
		for (j in seq(nrow(nas))){
			start_idx <- nas$starts[j]
			end_idx <- nas$ends[j]

			pre <- (start_idx - 1) %in% closed$ends
			post <- (end_idx + 1) %in% closed$starts

			if(pre & post){
				closed_pre <- which(closed$ends %in% (start_idx-1))
				closed_post <- which(closed$starts %in% (end_idx + 1))

				closed$ends[closed_pre] <- closed$ends[closed_post]
				closed <- closed[-closed_post, ]
			} else if(pre){
				closed_pre <- which(closed$ends %in% (start_idx-1))
				mid <- median(c(start_idx, end_idx)) - 0.5
				closed$ends[closed_pre] <- mid
			} else if(post){
				closed_post <- which(closed$starts %in% (end_idx + 1))
				mid <- median(c(start_idx, end_idx)) + 0.5
				closed$starts[closed_post] <- mid
			}
		}
	}
	return(closed)
}

stitchedDist <- function(contigs, diNT){
	# Calculates actual run lengths in base pairs
	#
	# Args:
	#	contigs: data frame of merged closed run start and stop indices
	#	diNT: GRanges object of GpC dinucleotides
	#
	# Returns:
	#	The same input data frame with additional meta data columns:
	#	start, end, dist, minDist, and maxDist

	contigs$maxDist <- runDistMax(contigs, diNT)
	contigs$minDist <- runDistMin(contigs, diNT)
	
	contigs$cov <- contigs$ends - contigs$starts + 1

	# Hack to account for first and last index
	contigs$starts[contigs$starts == 1] <- 2
	contigs$ends[contigs$ends == length(diNT)] <- length(diNT) - 1

	# Hack to account for intermediate idx
	int_starts <- contigs$starts %% 1 != 0
	int_ends <- contigs$ends %% 1 != 0

	starts <- start(diNT[contigs$starts[!int_starts]])
	ends <- end(diNT[contigs$ends[!int_ends]])
		
	prev_ends <- end(diNT[contigs$starts[!int_starts] - 1])
	next_start <- start(diNT[contigs$ends[!int_ends] + 1])
	
	contigs$start[!int_starts] <- (starts + prev_ends)/2
	contigs$start[int_starts] <- (start(diNT[floor(contigs$starts[int_starts])]) + 
								  end(diNT[floor(contigs$starts[int_starts])])) / 2
	contigs$end[!int_ends] <- (ends + next_start)/2
	contigs$end[int_ends] <- (start(diNT[ceiling(contigs$ends[int_ends])]) + 
							  end(diNT[ceiling(contigs$ends[int_ends])])) / 2

	contigs$dist <- round(contigs$end - contigs$start) + 1
	contigs$start <- round(contigs$start)
	contigs$end <- round(contigs$end)

	return(contigs)
}

runMat <- function(gcReads, GpC, runs, locus, gap, mode){
	# Generates run matrices with pre-determined cutoffs
	#
	# Args:
	#	gcReads: matrix of GpC methylation values
	#	GpC: GRanges object of GpC dinucleotides
	#	runs: data frame of run start and stop indices
	#	locus: GRanges object of locus
	#	gap: integer base pair spacing of tiles
	#	mode: string type of run matrix to generate
	#
	# Returns:
	#	A binary matrix (with NAs) where each row is a read
	#	and each column is a tile position within the locus.

	chr <- as.character(seqnames(locus))
	tiles <- tile(locus, width=gap)[[1]]
	tiles <- resize(tiles, width=1, fix='center')

	results <- matrix(data=NA, nrow=nrow(gcReads), ncol=length(tiles))

	for(i in seq(nrow(gcReads))){
		gcRow <- gcReads[i, ]
		gcIdx <- !is.na(gcRow)
		gcFirst <- min(which(gcIdx))
		gcLast <- max(which(gcIdx))

		read <- GRanges(chr, IRanges(start(GpC[gcFirst]), end(GpC[gcLast])))

		results[i, tiles %over% read] <- 0

		run.df <- runs[[i]]
		if(nrow(run.df) > 0){
			run.gr <- GRanges(chr, IRanges(run.df$start, run.df$end))
			run.gr$dist <- run.df$dist
			run.gr$left <- run.df$left
			run.gr$right <- run.df$right
			run.gr$minDist <- run.df$minDist
			run.gr$maxDist <- run.df$maxDist
		} else {
			run.gr <- GRanges()
		}

		if (mode == 'open'){
			results[i, tiles %over% run.gr[!(run.gr$dist < 80 &
											 run.gr$left > 80 & run.gr$right > 80)]] <-  1
		} else if (mode == 'short'){
			results[i, tiles %over% run.gr[run.gr$dist <= 80]] <-  1
		} else if (mode == 'mono'){
			results[i, tiles %over% run.gr[(run.gr$dist > 80 & run.gr$dist <= 200)]] <-  1
		} else if (mode == 'di'){
			results[i, tiles %over% run.gr[(run.gr$dist > 200 & run.gr$dist <= 390)]] <-  1
		} else if (mode == 'tri'){
			results[i, tiles %over% run.gr[(run.gr$dist > 390 & run.gr$dist < 590)]] <-  1
		} else if (mode == 'linker'){
			results[i, tiles %over% run.gr[(run.gr$dist < 80 &
											run.gr$left > 80 & run.gr$right > 80)]] <-  1
		}
	}

	return(results)
}

## Visualization ##
write_smooth_tdf <- function(gr, win, tiles.path, filename){
	# Writes out a smoothed TDF file
	#
	# Args:
	#	gr: GRanges object with data in 'scores' column
	#	win: integer number of basepairs to extend gr
	#	tiles.path: path to dinucleotide tiles RDS object
	#	filename: string output name of file
	#
	# Returns:
	#	Expands input gr and averages all overlapping ranges.
	#	Outputs the gr to a bedgraph to be converted to a TDF.

	# Expand windows by specified window size
	gr.w <- resize(gr, width(gr) + as.integer(win), fix='center')

	gr$mean <- avgHits(gr.w, gr)
	
	xx <- readRDS(tiles.path)

	tiles <- keepSeqlevels(xx, seqlevels(gr.w), pruning='coarse')
	rm(xx)

	# Find the mean of all overlapping methylation events
	hits <- findOverlaps(tiles, gr)
	tiles$score <- NA
	tiles$score[queryHits(hits)] <- gr$mean[subjectHits(hits)]

	export.bedGraph(sort(tiles[!is.na(tiles$score)]), filename)
	message(sprintf("Output successfully written to %s", filename))
}

write_bw <- function(loci, type, prefix, smooth=F, hap.names=NULL){
	# Writes out a bigwig file from all run matrices for given loci
	#
	# Args:
	#	loci: GRanges object of all loci of interest
	#	type: string type of run matrix to use ('open' or 'short')
	#	prefix: string prefix of file name
	#	smooth: boolean specifying whether to apply smoothing
	#	hap.names: string vector of read names assigned to a haplotype
	#
	# Returns:
	#	Gets the desired run matrix and applies optional smoothing.
	#	Writes the output to a bigwig file

	bw.list <- lapply(seq_along(loci), function(i, type){
		locus <- loci[i]
		message(locus$name)
		load(sprintf('rdata/%s-run.RData', locus$name))
		load(sprintf('rdata/%s-runMat.RData', locus$name))
		
		if (is.null(hap.names)) {
			idx <- rep(TRUE, nrow(meMat.gc.f))
		} else {
			idx <- chopStr(rownames(meMat.gc.f), 2) %in% hap.names
		}
		select <- idx & open_run.qc < 2
		buffer <- tiles

		if (type == 'open'){
			runMat <- runMat.open
		} else if(type == 'short'){
			runMat <- runMat.short
		}
		buffer$tot <- colSums(!is.na(runMat[select,]))
		buffer$run <- colSums(runMat[select,], na.rm=T)
		buffer$score <- buffer$run / buffer$tot

		if (smooth){
			buffer <- smooth_gr(buffer$score, buffer, thresh=12.5)
		}
		buffer <- resize(buffer, width=5, fix='center')

		return(buffer)
	}, type=type)

	bw <- do.call(c, bw.list)
	seqinfo(bw) <- keepSeqlevels(seqinfo(genome), seqlevels(bw))
	# Fill in gaps
	idx <- queryHits(findOverlaps(bw, drop.self=T))[c(TRUE, FALSE)]
	bw[idx] <- resize(bw[idx], width=4, fix='start')
	if(type == 'open'){
		if(smooth){
			filename <- paste0(prefix, '_open_runs_smooth.bw')
			export.bw(bw[!is.na(bw$score)], filename)
			message(sprintf("Output successfully written to %s", filename))
		} else {
			filename <- paste0(prefix, '_open_runs.bw')
			export.bw(bw[!is.na(bw$score)], filename)
			message(sprintf("Output successfully written to %s", filename))
		}
	} else if (type == 'short'){
		filename <- paste0(prefix, '_short_closed_runs.bw')
		export.bw(bw[!is.na(bw$score)], filename)
		message(sprintf("Output successfully written to %s", filename))
	}
}