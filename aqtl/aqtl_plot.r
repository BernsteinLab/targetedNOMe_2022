# This code examines an aQTL in the IKZF3 locus.
# It looks at haplotyped GM12878 data and identifies
# a point with a large difference in accessibility.

# libraries and functions
library(rtracklayer)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
hap1.path <- '/seq/epiprod02/Battaglia/NanoNOMe/eQTL_GM12878/workspace/GM12878_all_hap1.txt'
hap2.path <- '/seq/epiprod02/Battaglia/NanoNOMe/eQTL_GM12878/workspace/GM12878_all_hap2.txt'
rdata.path <- '/seq/epiprod02/Battaglia/NanoNOMe/eQTL_GM12878/workspace/rdata'

# Load in data
locus <- GRanges('chr17', IRanges(39853503,39906106))
locus$name <- 'eQTL'

hap1.names <- read.table(hap1.path, stringsAsFactors=F)$V1
hap2.names <- read.table(hap2.path, stringsAsFactors=F)$V1

pdf('IKZF3_aqtl.pdf')
all_loci.list <- lapply(seq_along(locus), function(i){
	locus <- locus[i]
	message(locus$name)
	load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
	load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))

	# Filter reads
	message(sprintf('%s total reads', nrow(runMat.open)))
	readCov <- rowMeans(!is.na(runMat.open))
	select <- open_run.qc < 2
	message(sprintf('%s total reads after filtering', sum(select)))
	hap1.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap1.names
	hap2.idx <- chopStr(rownames(meMat.gc.f[select,]), 2) %in% hap2.names
	message(sprintf('%s total reads assigned to a haplotype', sum(hap1.idx) + sum(hap2.idx)))
	
	# Calculate beta values and coverage
	if(sum(hap1.idx) > 1){
		GpC_fil$hap1 <- colMeans(meMat.gc.f[select,][hap1.idx,], na.rm=T)
		GpC_fil$cov1 <- colSums(!is.na(meMat.gc.f[select,][hap1.idx,]), na.rm=T)
	} else if(sum(hap1.idx) == 1){
		GpC_fil$hap1 <- meMat.gc.f[select,][hap1.idx,]
	}
	if(sum(hap2.idx) > 1){
		GpC_fil$hap2 <- colMeans(meMat.gc.f[select,][hap2.idx,], na.rm=T)
		GpC_fil$cov2 <- colSums(!is.na(meMat.gc.f[select,][hap2.idx,]), na.rm=T)
	} else if(sum(hap2.idx) == 1){
		GpC_fil$hap2 <- meMat.gc.f[select,][hap2.idx,]
	}
	
	# Require coverage to be greater than 10
	GpC_fil <- GpC_fil[GpC_fil %over% locus & GpC_fil$cov1 > 9 & GpC_fil$cov2 > 9]

	# Plot results
	plot(GpC_fil$hap1, GpC_fil$hap2, xlab='Haplotype 1 GpC Methylation', ylab='Haplotype 2 GpC Methylation', xlim=c(0,1), ylim=c(0,1), main=locus$name, pch=19)
	abline(0.4, 1, col='red', lty='dashed', lwd=1.5)
	abline(-0.4, 1, col='red', lty='dashed', lwd=1.5)
	print(GpC_fil[abs(GpC_fil$hap1 - GpC_fil$hap2) > 0.4])
	return(GpC_fil)
})
dev.off()