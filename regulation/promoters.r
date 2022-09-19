# This code examines data around transcription start sites.
# It uses GM12878, K562, and activated T-cell data.
# Metaplots are produced aggregating open run and CpG signal
# around TSSs and stratifying by expression. In addition,
# heatmaps are produced to visualize all TSSs in T-cells.

# libraries and functions
library(ggplot2)
library(rtracklayer)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
loci_tcell.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/TcellNov2020_loci.bed'
rdata.gm.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/200913_GM12878/rdata'
rdata.k562.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/200913_K562/rdata'
rdata.0.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t0h/workspace/rdata/'
rdata.48.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/t48h/workspace/rdata/'

tss_gm_k562.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/tss_gm_k562.rds'
tss_tcell.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/tss_tcell.rds'

get_averages <- function(loci, rdata.path, tss, variable, var, cpg=F, pad=F){
	.get_scores <- function(peaks, tiles, runMat, variable, var){
		scores.list <- lapply(seq_along(peaks), function(i,  tiles, runMat, variable, var){
			peak <- peaks[i]
			idx <- tiles %over% peak
			x <- start(tiles[idx]) - start(resize(peak, width=1, fix='center'))
			if(sum(idx) > 0){
				goodReads <- !is.na(rowMeans(runMat[,idx], na.rm=T))
				
				vals <- colMeans(runMat[goodReads,idx], na.rm=T)

				if (runValue(strand(peak)) == '-'){
					vals <- rev(vals)
				}
				vals.df <- as.data.frame(t(vals))
			}

			if(sum(idx) < 801){
				if(pad){
					pad_n <- 801 - sum(idx)
					pad <- as.data.frame(matrix(NA, 1, pad_n))
					if(min(x) > -1995 & runValue(strand(peak)) == '+' | 
					   max(x) < 1995 & runValue(strand(peak)) == '-'){
						vals.df <- cbind(pad, vals.df)
					} else {
						vals.df <- cbind(vals.df, pad)
					}
					colnames(vals.df) <- paste0('V', seq(801))
					idx <- rep(TRUE, 801)
				}
			}

			if(sum(idx) == 801){
				vals.df$TPM <- peak@elementMetadata@listData[var][[1]]
				vals.df$cell <- variable
				vals.df$id <- peak$id
				return(vals.df)
			}
		}, tiles=tiles, runMat=runMat, variable=variable, var=var)
		scores <- do.call(rbind, scores.list)
		return(scores)
	}

	averages.list <- lapply(seq_along(loci), function(i, rdata.path, tss, variable, var, cpg){
		locus <- loci[i]
		message(locus$name)
		load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))
		
		buffer <- resize(tss[tss %over% locus], width=4005, fix='center')
		if(length(buffer) == 0){
			return(NULL)
		}

		if(!cpg){
			runMat.open <- runMat.open[open_run.qc < 2,]
			scores <- .get_scores(buffer, tiles, runMat.open, variable, var)
		} else {
			load(sprintf('%s/%s-run.RData', rdata.path, locus$name))
			col_idx <- nearest(tiles, CpG_fil)
			runMat.cpg <- meMat.cg.f[,col_idx]
			colnames(runMat.cpg) <- NULL
			scores <- .get_scores(buffer, tiles, runMat.cpg, variable, var)
		}
		return(scores)

	}, rdata.path=rdata.path, tss=tss, variable=variable, var=var, cpg=cpg)

	averages <- do.call(rbind, averages.list)
}

plot_tss <- function(matrix, df, cellType, matType, lines=F){
	x <- seq(-2000, 2000,by=5)
	colIdx <- 1:801
	rowIdx <- df$cell == cellType
	par(mfrow=c(1,2))
	par(fig=c(0,2/3,0,1))
	plot(x, colMeans(matrix[df$TPM < 1 & rowIdx, colIdx], na.rm=T),
		 ylim=c(0,1), col='red', type='l', lwd=1.5, xlab='Distance from TSS', 
		 ylab=sprintf('%s Signal', matType), 
		 main=sprintf('%s Meta-plot (%s)', matType, cellType))
	lines(x, colMeans(matrix[df$TPM >= 1 & df$TPM < 10 & rowIdx, colIdx], na.rm=T), 
		  ylim=c(0,1), col='orange', type='l', lwd=1.5)
	lines(x, colMeans(matrix[df$TPM >= 10 & rowIdx, colIdx], na.rm=T), 
		  ylim=c(0,1), col='green', type='l', lwd=1.5)
	legend("topright", legend=c('TPM < 1', 'TPM < 10', 'TPM > 10'), 
		   col=c('red', 'orange', 'green'), lwd=2)
	
	if (matType == 'Open Run'){
		if (lines){
			abline(v=-210,col='red',lty='dashed', lwd=1.5)
			abline(v=90,col='red',lty='dashed', lwd=1.5)
		}
		bpPos <- 388
		bpText <- '-60'
	} else if (matType == 'CpG'){
		if (lines){
			abline(v=-50,col='red',lty='dashed', lwd=1.5)
			abline(v=250,col='red',lty='dashed', lwd=1.5)
		}
		bpPos <- 420
		bpText <- '+100'
	}
	
	par(new=T)
	par(fig=c(2/3,1,0,1))
	boxplot(matrix[df$TPM < 1 & rowIdx, bpPos], 
			matrix[df$TPM < 10 & df$TPM >=1& rowIdx, bpPos], 
			matrix[df$TPM >= 10 & rowIdx, bpPos], 
			names=c('TPM < 1', '1 < TPM < 10', 'TPM > 10'), 
			col=c('red', 'orange', 'green'), ylim=c(0,1), 
			main=sprintf('300 bp Smoothed Distributions at %s Bp', bpText))
	par(mfrow=c(1,1))
}

# Load in data
tss.f <- readRDS(tss_gm_k562.path)
tss.l <- readRDS(tss_tcell.path)

loci <- import.bed(loci.path)
loci_tcell <- import.bed(loci_tcell.path)

## GM12878 & K562 ##
tss_open.gm <- get_averages(loci, rdata.gm.path, tss.f, 'GM12878', 'GM')
tss_open.k562 <- get_averages(loci, rdata.k562.path, tss.f, 'K562', 'K562')

tss_open <- rbind(tss_open.gm, tss_open.k562)

tss_open.s <- apply(tss_open[,1:801], 1, function(row){
	smooth(as.numeric(row), 59)
})

pdf('TSS_open_run.pdf', width=15, height=10)
plot_tss(t(tss_open.s), tss_open, 'GM12878', 'Open Run', lines=F)
plot_tss(t(tss_open.s), tss_open, 'K562', 'Open Run', lines=F)
dev.off()

# CpG
tss_cpg.gm <- get_averages(loci, rdata.gm.path, tss.f, 'GM12878', 'GM', cpg=T)
tss_cpg.k562 <- get_averages(loci, rdata.k562.path, tss.f, 'K562', 'K562', cpg=T)

tss_cpg <- rbind(tss_cpg.gm, tss_cpg.k562)
tss_cpg.s <- apply(tss_cpg[,1:801], 1, function(row){
	smooth(as.numeric(row), 59)
})

pdf('TSS_cpg.pdf', width=15, height=10)
plot_tss(t(tss_cpg.s), tss_cpg, 'GM12878', 'CpG', lines=F)
plot_tss(t(tss_cpg.s), tss_cpg, 'K562', 'CpG', lines=F)
dev.off()

# Scatter plot
rowIdx <- tss_open$cell == 'GM12878'
scatter.df <- data.frame(cpg=t(tss_cpg.s)[rowIdx, 420], gpc=t(tss_open.s)[rowIdx, 388], 
						 id=tss_open$id[rowIdx], active=tss_open$TPM[rowIdx] > 1, 
						 tpm = tss_open$TPM[rowIdx])
scatter.df$principal <- tss.f$principal[scatter.df$id]
scatter.df$symbol <- tss.f$symbol[scatter.df$id]

pdf('GM12878_TSS_Scatter.pdf', width=10, height=8)
ggplot(scatter.df[!is.na(scatter.df$active) & scatter.df$principal,], 
	   aes(x=cpg, y=gpc, col=active)) + 
	geom_point(size=4) + geom_text(aes(label=symbol), hjust= -0.5, vjust=-0.5) + 
	xlab('CpG Methylation') + ylab('GpC Methylation') + ylim(0,0.8) + xlim(0,0.8) +
	scale_color_discrete(name='TPM > 1') + theme_bw(base_size=20) + 
	guides(colour = guide_legend(reverse=T))
dev.off()

## T-cell ##
tss_open.0 <- get_averages(loci_tcell, rdata.0.path, tss.l, 't0', 'V0', pad=T)
tss_open.48 <- get_averages(loci_tcell, rdata.48.path, tss.l, 't48', 'V48', pad=T)

tss_open_tcell <- rbind(tss_open.0, tss_open.48)

tss_open_tcell.s <- apply(tss_open_tcell[,1:801], 1, function(row){
	smooth(as.numeric(row), 59)
})

pdf('TSS_open_run_tcell.pdf', width=15, height=10)
plot_tss(t(tss_open_tcell.s), tss_open_tcell, 't0', 'Open Run', lines=F)
plot_tss(t(tss_open_tcell.s), tss_open_tcell, 't48', 'Open Run', lines=F)
dev.off()

# Poised promoters
tss_open.0$active <- tss_open.0$TPM > 1
tss_open.48$active <- tss_open.48$TPM > 1

tss_open.0.s <- apply(tss_open.0[,1:801], 1, function(row){
	smooth(as.numeric(row), 59)
})
tss_open.48.s <- apply(tss_open.48[,1:801], 1, function(row){
	smooth(as.numeric(row), 59)
})

pdf('tcell_poised_promoter.pdf', width=8, height=8)
plot(seq(-2000, 2000,by=5), 
	 colMeans(t(tss_open.0.s)[which(!tss_open.0$active & !tss_open.48$active),], na.rm=T), 
	 ylim=c(0,1), col='red', type='l', lwd=1.5, xlab='Distance from TSS', 
	 ylab='Open Run Signal', main='Open Run Meta-plot (t0)')
lines(seq(-2000, 2000,by=5), 
	  colMeans(t(tss_open.0.s)[which(!tss_open.0$active & tss_open.48$active),], na.rm=T), 
	  ylim=c(0,1), col='orange', type='l', lwd=1.5)
legend("topright", legend=c('Remains Inactive', 'Becomes Active'), 
	   col=c('red', 'orange'), lwd=2)
dev.off()

# Heatmap 
tss_mat.0 <- t(tss_open.0.s)
tpm.0 <- tss_open.0$TPM
cols.0 <- vector(mode='character', length=length(tpm))
cols.0[tpm.0 < 1] <- 'red'
cols.0[tpm.0 >= 1] <- 'green'
rownames(tss_mat.0) <- paste0(tss.l$symbol[tss_open.0$id], ' (', round(tpm.0, 2), ')')
colnames(tss_mat.0) <- seq(-2000, 2000,by=5)

tss_mat.48 <- t(tss_open.48.s)
tpm.48 <- tss_open.48$TPM
cols.48 <- vector(mode='character', length=length(tpm))
cols.48[tpm.48 < 1] <- 'red'
cols.48[tpm.48 >= 1] <- 'green'
rownames(tss_mat.48) <- paste0(tss.l$symbol[tss_open.48$id], ' (', round(tpm.48, 2), ')')
colnames(tss_mat.48) <- seq(-2000, 2000,by=5)

# Reorder by moving low expression to the bottom
order <- order(tpm.0, decreasing=F)
split <- min(which(tpm.0[order] >= 1 ))
num <- length(tpm)
order <- c(order[1:(split-1)][tss_open.48$TPM[order][1:(split-1)] < 1], 
		   order[1:(split-1)][tss_open.48$TPM[order][1:(split-1)] >= 1], 
		   order[split:num][tss_open.48$TPM[order][split:num] < 1], 
		   order[split:num][tss_open.48$TPM[order][split:num] >= 1])

pdf('tcell_promoter_heatmap.pdf', width=9, height=9)
heatmap(tss_mat.0[order, ], Rowv=NA, Colv=NA, main='t0h', scale='none', 
		RowSideColors=cols.0[order], margins=c(3,7.5))
heatmap(tss_mat.48[order, ], Rowv=NA, Colv=NA, main='t48h', scale='none', 
		RowSideColors=cols.48[order], margins=c(3,7.5))
dev.off()
