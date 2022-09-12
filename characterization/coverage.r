# This code examines the coverage of targeted NanoNOMe.
# It looks at each experiment that was run and plots the
# coverage for each loci as well as a control whole genome experiment.

# libraries and functions
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(ggplot2)
library(cowplot)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

genome <- BSgenome.Hsapiens.UCSC.hg38

# Paths to data and references
loci_gm.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
loci_h9.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/210603_loci.bed'
loci_tcell.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/TcellNov2020_loci.bed'
loci_stim.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_march2020/Tstim/donor73/200326_loci.bed'

bam_gm.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_GM12878/workspace/200913_GM12878.filtered.bam'
bam_k562.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/workspace/200913_K562.filtered.bam'
bam_h9.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9.filtered.bam'
bam_hsmm.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/HSMM/workspace/HSMM.filtered.bam'
bam_0h.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/t0h.filtered.bam'
bam_24h.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t24h/workspace/t24h.filtered.bam'
bam_48h.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/t48h.filtered.bam'
bam_stim.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_march2020/Tres/donor72/workspace/200304.filtered.bam'
bam_wg <- '/seq/epiprod02/Battaglia/NanoNOMe/190330_GpC/190330.filtered.bam'

# Functions
getRandomRegion <- function(width, avoid, hit){
	set.seed(1234)
	
	# Try until random region is obtained
	while(TRUE){
		# Choose random chromosome
		chr <- sample(paste0('chr', c(seq(22), 'X', 'Y')), 1)
		len <- seqlengths(seqinfo(genome)[chr])

		# Pick random region somewhere within bounds of chromosome
		seed <- sample((width/2):(len - width/2), 1)

		gr <- GRanges(chr, IRanges(seed-width/2, seed+width/2 - 1))

		if (gr %outside% avoid & gr %over% hit){
			# Success if random region is outside avoid and within hit
			return(gr)
		}
	}
}

getCoverage <- function(gcReads, GpC, cgReads, CpG){  

	chr <- unique(seqnames(GpC))

	# Determine extent of reads
	gcBounds <- apply(gcReads, 1, function(row){
		buffer <- which(!is.na(row))
		return(cbind(min(buffer), max(buffer)))
	})
	gcStarts <- start(GpC[gcBounds[1,]])
	gcEnds <- end(GpC[gcBounds[2,]])
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
	}

	return(reads)
}

expandText <- function(x){
	gsub('_', '\n', x)
}

## Coverage Plots ##

# Whole Genome coverage for reference
bam_wg <- import(bam_wg)
bam_wg.gr <- GRanges(seqnames(bam_wg), IRanges(start(bam_wg), end(bam_wg)))

# Total estimated genome length: 3,110,748,599
# Probably downsample this number
message(sprintf('%sx genome-wide coverage', sum(width(bam_wg)) / 3110748599))

# Examine each experiment
tab <- data.frame(loci=c(rep(loci_gm.path,2), rep(loci_h9.path,2), 
						 rep(loci_tcell.path, 3), loci_stim.path),
				  bam=c(bam_gm.path, bam_k562.path, bam_h9.path, bam_hsmm.path,
						bam_0h.path, bam_24h.path, bam_48h.path, bam_stim.path),
				  stringsAsFactors=F)

pdf('loci_cov.pdf', width=24, height=12)
experiment.list <- mapply(function(bam.path, loci.path){
	bam <- import(bam.path)
	bam.gr <- GRanges(seqnames(bam), IRanges(start(bam), end(bam)))
	
	loci <- import(loci.path)
	loci <- loci[order(width(loci), decreasing=F)]
	
	# Whole genome reference
	rand_size <- 130000
	rand <- getRandomRegion(rand_size, loci, bam_wg.gr)

	cov.df.wg.ds <- data.frame(width=sort(width(bam_wg.gr[bam_wg.gr %over% rand])), locus='WG')
	cov.df.wg.ds2 <- data.frame(width=sort(width(bam_wg.gr[bam_wg.gr %over% rand])) / rand_size, locus='WG')
	cov.df.wg.ds$x <- seq(nrow(cov.df.wg.ds))
	cov.tot.df.wg  <- data.frame(cov=sum(width(bam_wg)) / 3110748599, locus='WG')

	# Get reads over loci
	all_loci.list <- lapply(seq_along(loci), function(i){
		locus <- loci[i]
		return(bam.gr[bam.gr %over% locus])
	})

	# Convert reads to data frame
	cov.list <- lapply(seq_along(loci), function(i){
		locus <- loci[i]
		
		buffer <- all_loci.list[[i]]
		
		# Clip reads
		start(buffer[start(buffer) < start(locus)]) <- start(locus)
		end(buffer[end(buffer) > end(locus)]) <- end(locus)
		
		df <- data.frame(width=width(buffer), locus=locus$name)
		df <-  df[order(df$width, decreasing=F), ]
		df$x <- seq(nrow(df))
		
		return(df)
	})
	cov.df <- do.call(rbind, cov.list)
	cov.df <- rbind(cov.df.wg.ds, cov.df)
	
	# Percent full length
	cov.list2 <- lapply(seq_along(loci), function(i){
		locus <- loci[i]
		
		buffer <- all_loci.list[[i]]
		# Clip reads
		start(buffer[start(buffer) < start(locus)]) <- start(locus)
		end(buffer[end(buffer) > end(locus)]) <- end(locus)

		df <- data.frame(width=width(buffer) / width(locus), locus=locus$name)
		df <-  df[order(df$width, decreasing=F), ]
		df$x <- seq(nrow(df))

		return(df)
	})

	top <- max(sapply(cov.list2, function(el){
		return(mean(el$width == 1))
	}))
	message(sprintf('Max full coverage: %s', top))

	frac.full.length.list <- lapply(append(list(cov.df.wg.ds2), cov.list2) , function(el){
		return(data.frame(pct=mean(el$width == 1), locus=unique(el$locus)))
	})
	frac.full.length <- do.call(rbind, frac.full.length.list)

	cov.tot.list <- lapply(seq_along(all_loci.list), function(i){
		widths <- cov.list[[i]]$width
		data.frame(cov=sum(widths) / max(widths), locus=loci[i]$name)
	})
	cov.tot.df <- do.call(rbind, cov.tot.list)
	cov.tot.df <- rbind(cov.tot.df.wg, cov.tot.df)

	# Plot results
	p1 <- ggplot(cov.df[cov.df$width > 0,], aes(x=x, y=width)) + geom_area(fill='cornflowerblue', stat='identity', color='black') + 
		facet_grid(. ~ locus, scales='free_x') + scale_y_continuous(breaks=seq(0,200000, by=10000), labels=seq(0,200,by=10), limits=c(0,140000)) +
		ylab('Read Length (kB)') + theme_bw(base_size=24) + 
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=16), panel.grid.major.x = element_blank(),
			strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank(), panel.spacing.x=unit(0,'lines'))

	p2 <- ggplot(cov.tot.df, aes(x=locus, y=cov)) + geom_bar(stat='identity', fill='lightcoral') +
		scale_y_continuous(breaks=seq(0,1000, by=100), labels=paste0(seq(0,1000,by=100), 'X')) + scale_x_discrete(label=expandText) + 
		geom_errorbar(aes(x=locus, ymax = 100, ymin = 100), size=0.5, linetype = "longdash", inherit.aes = F, width = 1) + 
		theme_bw(base_size=24) + xlab('') + ylab('Depth of Coverage') +
		theme(axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank())

	p3 <- ggplot(frac.full.length, aes(x=locus, y=pct)) + geom_bar(stat='identity') + scale_x_discrete(label=expandText) +
		scale_y_continuous(limits=c(0,0.35), breaks=seq(0,0.35, by=0.1), labels=paste0(seq(0,35, by=10), '%')) + theme_bw(base_size=24) + xlab('Locus') + ylab('Percent Full Length') +
		theme(axis.text.y = element_text(size=16), axis.text.x = element_text(size=14), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank())

	print(plot_grid(p1, p2, p3, align = "v", axis='l', ncol=1, rel_heights = c(2, 1, 1)))

},  bam.path=tab$bam, loci.path=tab$loci)
dev.off()
