# This code tallies the the fold enrichment of each experiment.
# It sums the total lenghts of reads in each locus and
# divides it by the width of the locus. Outputs a text file.

# libraries and functions
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(ggplot2)
library(cowplot)

genome <- BSgenome.Hsapiens.UCSC.hg38

# Paths to data and references
loci_gm.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
loci_h9.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/210603_loci.bed'
loci_tcell.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/TcellNov2020_loci.bed'
loci_stim.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_march2020/Tstim/donor73/200326_loci.bed'

bam_gm.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_GM12878/workspace/200913_GM12878.bam'
bam_k562.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/workspace/200913_K562.bam'
bam_h9.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/H9/workspace_hg38/H9.bam'
bam_hsmm.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Imprinting_May2021/HSMM/workspace/HSMM.bam'
bam_0h.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/t0h.bam'
bam_24h.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t24h/workspace/t24h.bam'
bam_48h.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/t48h.bam'
bam_stim.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_march2020/Tres/donor72/workspace/200304.reduced.bam'
bam_wg <- '/seq/epiprod02/Battaglia/NanoNOMe/190330_GpC/190330.filtered.bam'

## Fold Enrichment ##
tab <- data.frame(loci=c(rep(loci_gm.path,2), rep(loci_h9.path,2), 
						 rep(loci_tcell.path, 3), loci_stim.path),
				  bam=c(bam_gm.path, bam_k562.path, bam_h9.path, bam_hsmm.path,
						bam_0h.path, bam_24h.path, bam_48h.path, bam_stim.path),
				  stringsAsFactors=F)

cov.list <- mapply(function(bam.path, loci.path){
	bam <- import(bam.path)
	loci <- import(loci.path)
	bam.gr <- GRanges(seqnames(bam), IRanges(start(bam), end(bam)))

	# Get reads over loci
	all_loci.list <- lapply(seq_along(loci), function(i){
		locus <- loci[i]
		return(bam.gr[bam.gr %over% locus])
	})

	# Sum up coverage as fraction of locus width
	cov.list <- sapply(seq_along(loci), function(i){
		locus <- loci[i]
		buffer <- all_loci.list[[i]]

		# Clip reads
		start(buffer[start(buffer) < start(locus)]) <- start(locus)
		end(buffer[end(buffer) > end(locus)]) <- end(locus)

		return(sum(width(buffer))/width(locus))
	})
	names(cov.list) <- loci$name

	return(cov.list)
	# # Target enrichments over rest of the genome
	# off_width <- sum(width(bam[bam %outside% loci]))
	# off_enrich <- off_width / (3110748599 - sum(width(loci)))
	# off_enrich
	# max(cov.list / off_enrich)
	# loci[which.max(cov.list)]
},  bam.path=tab$bam, loci.path=tab$loci)

# Write table of coverages
lapply(cov.list, write.table, 'fold_enrichments.txt', append=T, quote=F)