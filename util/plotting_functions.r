#******************************************************************************
# Collection of plotting functions used for figures
#******************************************************************************

plotReadMat <- function(gcReads, GpC, cgReads=NULL, CpG=NULL, open_runs, 
						closed_runs, snpReads=NULL, snps=NULL, subset=NULL, 
						order=NULL, gap, region, rep=3, open=F, linker=F, 
						short=F, cpg=F, gpc=F, snp=F, peaks=NULL, plot=T){  
	# Plots single reads with color annotations
	#
	# Args:
	#	gcReads: GpC methylation matrix (read x position)
	#	GpC: GRanges object corresponding to columns of gcReads
	#	cgReads: optional CpG methylation matrix (read x position)
	#	CpG: GRanges object corresponding to columns of cgReads
	#	open_runs: list of open runs corresponding to gcReads
	#	closed_runs: list of closed runs corresponding to gcReads
	#	subset: vector for subsetting gcReads, open_runs, and closed_runs
	#	order: vector for ordering gcReads, open_runs, and closed_runs (after subsetting)
	#	gap: number of basepairs to tile the locus by (reduces the resolution)
	#	region: GRanges object corresponding to locus that is being plotted
	#	rep: the number lines (height) of each read
	#	open: boolean indicating whether to color open runs
	#	linker: boolean indicating whether to color linker runs
	#	short: boolean indicating whether to color short runs
	#	cpg: boolean indicating whether to color CpG methylation (cannot combine with gpc)
	#	gpc: boolean indicating whether to color GpC methylation (cannot combine with cpg)
	#	peaks: GRanges object of peak regions to highlight within plot
	#	plot: boolean indicating whether to plot or to just return plot matrix
	#
	# Returns:
	# 	Plots single reads with color annotations in current viewport

	# Split region of interest into tiles of specified gap
	tiles <- tile(region, width=gap)[[1]]
	tiles <- resize(tiles, width=1, fix='center')

	chr <- unique(seqnames(GpC))

	# Subset
	if(!is.null(subset)){
		gcReads <- gcReads[subset,]
		open_runs <- open_runs[subset]
		closed_runs <- closed_runs[subset]
	}

	# Set order
	if(!is.null(order)){
		gcReads <- gcReads[order,]
		open_runs <- open_runs[order]
		closed_runs <- closed_runs[order]
	}

	# Determine extent of reads
	gcBounds <- apply(gcReads, 1, function(row){
		buffer <- which(!is.na(row))
		return(cbind(min(buffer), max(buffer)))
	})
	gcStarts <- start(GpC[gcBounds[1,]])
	gcEnds <- end(GpC[gcBounds[2,]])

	# Look at CpG matrix if applicable
	if(!is.null(cgReads)){
		cgIdx <- match(rownames(gcReads), rownames(cgReads))
		cgBounds <- apply(cgReads[cgIdx,], 1, function(row){
			buffer <- which(!is.na(row))
			if(length(buffer) == 0){
				return(cbind(ncol(cgReads), 0))
			}
			return(cbind(min(buffer), max(buffer)))
		})
		cgStarts <- start(CpG[cgBounds[1,]])
		cgEnds <- end(CpG[cgBounds[2,]])
		reads <- GRanges(chr, IRanges(pmin(gcStarts, cgStarts), pmax(gcEnds, cgEnds))) 

		# Expand CpG
		cgTemplate <- CpG
		cgTemplate <- resize(CpG, width=gap, fix='center') 
	} else {
		reads <- GRanges(chr, IRanges(gcStarts, gcEnds))
	}

	# SNPs
	if(!is.null(snps)){
		snpIdx <- match(rownames(gcReads), rownames(snpReads))
		expand <- gap * 5
		snpTemplate <- resize(snps, width=expand, fix='center')
		snpDist <- start(snpTemplate)[-1] - end(snpTemplate)[-length(snpTemplate)]
		end(snpTemplate[which(snpDist <= 0 )]) <- end(snpTemplate[which(snpDist<= 0)]) - floor(abs(snpDist)[snpDist<= 0] / 2)
		start(snpTemplate[which(snpDist <= 0) + 1]) <- start(snpTemplate[which(snpDist<= 0)+1]) + ceiling(abs(snpDist)[snpDist<= 0] / 2)
	}

	gcTemplate <- GpC
	gcTemplate <- resize(GpC, width=gap, fix='center') 

	# Create GRanges for open and closed runs
	# TODO: precalculate this
	open_runs.gr <- lapply(open_runs, function(open.df){
		if(nrow(open.df) > 0){
			open.gr <- GRanges(chr, IRanges(open.df$start, open.df$end))
			open.gr$dist <- open.df$dist
			open.gr$minDist <- open.df$minDist
			open.gr$maxDist <- open.df$maxDist

			open.gr$left <- open.df$left
			open.gr$right <- open.df$right
		} else {
			open.gr <- GRanges()
		}
		return(open.gr)
	})

	closed_runs.gr <- lapply(closed_runs, function(closed.df){
		if(nrow(closed.df) > 0){
			closed.gr <- GRanges(chr, IRanges(closed.df$start, closed.df$end))
			closed.gr$dist <- closed.df$dist
			closed.gr$minDist <- closed.df$minDist
			closed.gr$maxDist <- closed.df$maxDist
		} else {
			closed.gr <- GRanges()
		}
		return(closed.gr)
	})

	# Initialize data structures
	results <- matrix(data=NA, nrow=nrow(gcReads), ncol=length(tiles))

	if(short){
		x.points <- vector('list', nrow(gcReads))
		y.points <- vector('list', nrow(gcReads))
	}
	if(cpg | gpc){
		x2.points <- vector('list', nrow(gcReads))
		y2.points <- vector('list', nrow(gcReads))

		x2.5.points <- vector('list', nrow(gcReads))
		y2.5.points <- vector('list', nrow(gcReads))
	}

	n <- nrow(gcReads)
	for(i in seq(n)){
		cat(sprintf('\rRendering %s out of %s reads', i, n))
		
		# TODO: Remove hardcoded colors and cutoffs
		# Paint Gray with reads
		result <- results[i,]
		result[tiles %over% reads[i]] <- rgb(224,224,224, 150, maxColorValue=255)

		if (open){
			open.gr <- open_runs.gr[[i]]
			result[tiles %over% open.gr[!(open.gr$dist < 80 & open.gr$right > 80 & open.gr$left > 80)]] <- rgb(27, 158, 119, maxColorValue=255)
		}
		if (linker){
			open.gr <- open_runs.gr[[i]]
			result[tiles %over% open.gr[(open.gr$dist < 80 & open.gr$right > 80 & open.gr$left > 80)]] <- rgb(255, 178, 102, maxColorValue=255)
		}

		if (short){
			closed.gr <- closed_runs.gr[[i]]
			result[tiles %over% closed.gr[closed.gr$dist <= 80]] <- rgb(153, 112, 171, maxColorValue=255)
		}

		if (cpg){
			cgBuffer <- cgTemplate
			cgBuffer$score <- cgReads[cgIdx[i],]
			
			x2.5.points[[i]] <- which(tiles %over% cgBuffer[which(cgBuffer$score == 0)])
			y2.5.points[[i]] <- rep((rep / 2) + (rep + 1) * (n - i), length(x2.5.points[[i]]))
			
			x2.points[[i]] <- which(tiles %over% cgBuffer[which(cgBuffer$score == 1)])
			y2.points[[i]] <- rep((rep / 2) + (rep + 1) * (n - i), length(x2.points[[i]]))
		}

		if (gpc){
			gcBuffer <- gcTemplate
			gcBuffer$score <- gcReads[i,]
			x2.points[[i]] <- which(tiles %over% gcBuffer[which(gcBuffer$score == 1)])
			y2.points[[i]] <- rep((1 + rep / 2) + (rep + 1) * (n - i), length(x2.points[[i]]))

			x2.5.points[[i]] <- which(tiles %over% gcBuffer[which(gcBuffer$score == 0)])
			y2.5.points[[i]] <- rep((1 + rep / 2) + (rep + 1) * (n - i), length(x2.5.points[[i]]))
		}

		if (snp){
			snpBuffer <- snpTemplate
			snpBuffer$score <- snpReads[snpIdx[i],]
			
			result[tiles %over% snpBuffer[which(snpBuffer$score == 1 & snpBuffer$geno == 'homo')]] <- rgb(178, 90, 39, maxColorValue=255)
			result[tiles %over% snpBuffer[which(snpBuffer$score == 1 & snpBuffer$geno == 'hetero')]] <- rgb(84, 45, 135, maxColorValue=255)
		}

		results[i, ] <- result
	}

	if (cpg | gpc){
		# Rescale coordinates of points
		x2 <- base::do.call(c, x2.points) / length(tiles)
		y2 <- base::do.call(c, y2.points) / (n * (rep+1)-1)

		x2.5 <- base::do.call(c, x2.5.points) / length(tiles)
		y2.5 <- base::do.call(c, y2.5.points) / (n * (rep+1)-1)
	}

	# Expand results matrix
	results2 <- matrix(data=NA, nrow=nrow(gcReads) * (rep + 1), ncol=length(tiles))
	for(i in seq(rep)){
		results2[seq(i, nrow(results2), by=rep + 1), ] <- results
	}
	results2[is.na(results2)] <- "#FFFFFF"
	
	if(!plot){
		return(results2)
	}

	# if (runValue(strand(region)) == '-'){
	#     x <- -x + length(tiles) + 1
	#     x2 <- -x2 + length(tiles) + 1
	#     x2.5 <- -x2.5 + length(tiles) + 1
	#     results2 <- results2[,ncol(results2):1]
	# }

	results2 <- results2[-nrow(results2),]
	grid.raster(results2, interpolate=F, width=unit(1,'npc'), height=unit(1,'npc'))

	if(cpg | gpc){
		if(length(x2) > 0){
		   grid.points(x2, y2, pch=21, size=unit(rep*0.66/(n*(rep+1)), 'npc'), gp=gpar(fill='#ed2127', col='#ed2127', alpha=0.5)) 
		}
		if(length(x2.5) > 0){
			grid.points(x2.5, y2.5, pch=21, size=unit(rep*0.66/(n*(rep+1)), 'npc'), gp=gpar(fill='#0b529d', col='#0b529d', alpha=0.5))        
		}
	}

	# Highlight peaks if wanted
	if(!is.null(peaks)){

		.highlight <- function(peaks, tiles, n, rep){
			hits <- findOverlaps(peaks, tiles)
			
			.findBounds <- function(idx){
				start <- vector('integer')
				end <- vector('integer')
				
				i <- 1
				while(TRUE){
					start <- c(start, idx[i])
					search <- TRUE
					while(search & i < length(idx)){
						if (idx[i] != idx[i+1] - 1){
							end <- c(end, idx[i])
							search <- FALSE
						}
						i = i+1
					}
					if (search){
						end <- c(end, idx[i])
						break
					}
				}
				
				return(cbind(start, end))
			}
			
			bounds <- .findBounds(unique(subjectHits(hits)))
			
			for (i in seq(nrow(bounds))){
				s <- bounds[i,1]
				e <- bounds[i, 2]
				
				x0 <- s-1 
				x1 <- e+1
				y0 <- n * (rep+1)
				
				x <- c(x0, x1, x1, x0)
				y <- c(y0, y0, 0, 0)
				
				grid.polygon(x=x / length(tiles), y=y / (n * (rep + 1)), gp=gpar(col='red', fill=NA, lwd=2, lty='dashed'))
				
			}
		}

		.highlight(peaks, tiles, n, rep)
	}
}

sideColorBar <- function(colors, rep=3, width=0.5){
	# Plots colored rectangles to be paired with single reads
	#
	# Args:
	#	colors: the (ordered) colors to be plotted
	#	rep: the number lines (height) of the rectangle
	#	width: number indicating the fraction of viewport to fill in
	#
	# Returns:
	# 	Plots the provided colors with the specified number of lines
	
	# Repeat colors according to rep
	result <- rep(colors, rep(rep+1, length(colors)))
	
	# Set every rep+1 line to white
	result[seq(0,length(colors)*(rep+1),by=rep+1)[-1]] <- "#FFFFFF"

	# Drop last line
	result <- result[-length(result)]
	
	# Plot
	grid.raster(result, x=unit(1,'npc'), width=unit(width, 'npc'), height=unit(1,'npc'), just='right', interpolate=F)
}

context <- function(locus, peaks, resize.width, reverse=F){
	# Draws call-outs to provide context
	#
	# Args:
	#	locus: GRanges object of the entire region being plotted
	#	peaks: GRanges object of the peaks being called out
	#	resize.width: size to resize the peak widths
	#	reverse: boolean indicating whether to flip the triangles
	#
	# Returns:
	# 	Plots polygons to connect reads plot to Gviz peaks
	
	# Identify ranges
	range1 <- length(peaks)
	range2 <- width(locus)

	# Set y coordinates
	if(reverse){
		y0 <- 0.1
		y1 <- 0.9
	} else {
		y0 <- 0.9
		y1 <- 0.1
	}
		
	widths <- width(peaks) / (resize.width * 2)
	for (i in seq(length(peaks))){
		# Find centers
		c <- (2 * i - 1) / (2 * range1)

		x1 <- c - (widths[i] / range1)
		x2 <- c + (widths[i] / range1)
		x3 <- (end(peaks[i]) - start(locus)) / range2
		x4 <- (start(peaks[i]) - start(locus)) / range2
				
		grid.polygon(x=c(x1, x2, x3, x4), y=c(y0, y0, y1, y1), gp=gpar(fill='gray'))
		
	}    
}

plotTriMatrix <- function(mat, col_pal) {
	# Plots upper triangle of square matrix
	#
	# Args:
	#	mat: square, symmetric matrix with values between 0 and 1
	#	col_pal: array of colors corresponding to the scale
	#
	# Returns:
	# 	Plots the upper triangle in a (presumed) rotated viewport
	#	Returns 11-step color scale for further plotting

	# Bake in some warnings
	if (max(mat, na.rm=T) > 1){
		warning('Maximum value exceeds 1')
	}
	if (min(mat, na.rm=T) < 0){
		warning('Minimum value less than 0')
	}

	# Initialize result
	x <- matrix(data=NA, nrow=nrow(mat), ncol=ncol(mat))

	# Delete lower triangle
	mat[lower.tri(mat)] <- NA
	idx <- !is.na(mat)

	# Convert values to colors
	col_fun <- colorRamp(col_pal)
	# col_fun <- colorRamp(rev(brewer.pal(11, col_pal)))
	# col_fun <- colorRamp(brewer.pal(11, col_pal))
	x[idx] <- rgb(col_fun(mat[idx]), maxColorValue=255)
	x[!idx] <- rgb(255,255,255, maxColorValue=255)

	grid.raster(x, interpolate=F)
	return(rgb(col_fun(rev(0:11)/11), maxColorValue=255))
}

# Code examining a potential bug with rotated viewports
# pdf('test.pdf', width=12, height=12)
# # jpeg('test.jpg', width=2000, height=2000)
# grid.newpage()
# pushViewport(viewport(width=unit(1 / sqrt(2),'snpc'), height=unit(1 / sqrt(2),'snpc')))
# grid.raster(test2, width=unit(1,'npc'), height=unit(1,'npc'), interpolate=F)
# n <- nrow(test2)
# grid.polyline(unit(rep(0:n, rep(2,n+1))/n, 'npc'), unit(c(rep(c(0, 1,1, 0), n/2),0,1), 'npc'), gp=gpar(col='red', lwd=2))
# grid.polyline(unit(c(rep(c(0, 1,1, 0), n/2),0,1), 'npc'), unit(rep(0:n, rep(2,n+1))/n, 'npc'), gp=gpar(col='red', lwd=2))

# grid.newpage()
# pushViewport(viewport(width=unit(1 / sqrt(2),'snpc'), height=unit(1 / sqrt(2),'snpc'), angle=45))
# grid.raster(test2, width=unit(1,'npc'), height=unit(1,'npc'), interpolate=F)
# n <- nrow(test2)
# grid.polyline(unit(rep(0:n, rep(2,n+1))/n, 'npc'), unit(c(rep(c(0, 1,1, 0), n/2),0,1), 'npc'), gp=gpar(col='red', lwd=2))
# grid.polyline(unit(c(rep(c(0, 1,1, 0), n/2),0,1), 'npc'), unit(rep(0:n, rep(2,n+1))/n, 'npc'), gp=gpar(col='red', lwd=2))
# dev.off()


arcs <- function(corMat, n=5, peak.ind=NULL, min.height=0.5, hl.ratio){
	# Draws arcs underneath triangular heatmap
	#
	# Args:
	#	corMat: square, symmetric matrix with values between 0 and 1
	#	n: number of top interactions to plot
	#	peak.ind: array of indices that contain peaks for limiting possible interactions
	#	hl.ratio: ratio for scaling curvature
	#
	# Returns:
	# 	Plots n arcs in current viewport

	# Strip redundant info
	corMat[lower.tri(corMat)] <- NA
	diag(corMat) <- NA
	range <- nrow(corMat)

	# Get all possible pairwise interactions and order them
	if(is.null(peak.ind)){
		ranked <- sort(as.vector(corMat), na.last=T, decreasing=T)
		connections <- (which(corMat >= ranked[n], arr.ind=T) - 0.5) / range
	} else {
		combos <- combn(peak.ind, 2)
		order <- order(as.vector(corMat[t(combos)]), na.last=T, decreasing=T)
		connections <- (t(combos)[order[1:n],] - 0.5) / range
		ranked <- sort(as.vector(corMat[t(combos)]), na.last=T, decreasing=T)
	}

	# Draw arcs
	for(i in seq(n)){
		x1 <- connections[i,1]
		x2 <- connections[i,2]
		c <- 2 * hl.ratio
		grid.curve(x1, 0, x2, 0, ncp=10, curvature=-c, square=F, gp=gpar(lwd=3*ranked[i]))
	}
}

highlight <- function(diNT, peaks, lwd=1.5){
	# Highlights the intersection of peaks in the triangle plot
	#
	# Args:
	#	diNT: GRanges object corresponding rows/cols of square matrix
	#	peaks: GRanges object corresponding to regions of interest
	#	lwd: the width of the highlighting line
	#
	# Returns:
	# 	Draws green rectangles on existing raster plot indicating peaks
	
	# Identify ranges
	range <- length(diNT)

	.findBounds <- function(idx){
		start <- vector('integer')
		end <- vector('integer')
		
		i <- 1
		while(TRUE){
			start <- c(start, idx[i])
			search <- TRUE
			while(search & i < length(idx)){
				if (idx[i] != idx[i+1] - 1){
					end <- c(end, idx[i])
					search <- FALSE
				}
				i = i+1
			}
			if (search){
				end <- c(end, idx[i])
				break
			}
		}
		
		return(cbind(start, end))
	}

	# Identify peak overlaps
	hits <- findOverlaps(diNT, peaks)
	bounds <- .findBounds(unique(queryHits(hits)))
	
	l <- bounds[,2] - bounds[,1] + 1
	x <- (bounds[,1] + l/2) - 1
	coords <- expand.grid(x, x) / range
	dims <- expand.grid(l, l) / range

	nonSkip <- coords[,1] < coords[,2]

	for (i in seq(sum(nonSkip))){
		x1 <- coords[nonSkip,][i,2]
		y1 <- 1 - coords[nonSkip,][i,1]
		
		w <- dims[nonSkip,][i,2]
		h <- dims[nonSkip,][i,1]
		grid.rect(x=x1, y=y1, width=w, height=h, gp=gpar(fill=NA, col='green', lwd=lwd))    
	}    
}
