# This code compares the open run signal with DNase.
# It does this for both GM12878 and K562.

# libraries and functions
library(cowplot)
library(dplyr)
library(ggplot2)
library(pROC)
library(rtracklayer)
source('../util/helper_functions.r')
source('../util/plotting_functions.r')

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
rdata.gm.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/200913_GM12878/rdata'
rdata.k562.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/200913_K562/rdata'

dnase.1.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/ENCODE/gm12878_dnase_hg38_rep1_ENCFF598KWZ.bed'
dnase.2.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/ENCODE/gm12878_dnase_hg38_rep2_ENCFF073ORT.bed'
dnase.k562.path <- '/seq/epiprod02/Battaglia/ENCODEdata/K562/hg38/ENCFF433TIR_DNase_k562_narrowPeak_JS_2.bed'

# Load in data
loci <- import.bed(loci.path)

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
						  qValue = "numeric", peak = "integer")
dnase.1 <- import(dnase.1.path, extraCols=extraCols_narrowPeak)
dnase.2 <- import(dnase.2.path, extraCols=extraCols_narrowPeak)
dnase.gm <- mergePeaks(dnase.1, dnase.2)
dnase.gm$idx <- seq(length(dnase.gm))

dnase.k562 <- import(dnase.k562.path, extraCols=extraCols_narrowPeak)
dnase.k562$idx <- seq(length(dnase.k562))

generate_controls <- function(loci, dnase){
	# Fixed random seed
	set.seed(1234)

	control.list <- lapply(seq_along(loci), function(i, dnase){
		locus <- loci[i]
		control <- GRanges()
		
		num <- sum(dnase %over% locus)
		if(num == 0){
			return(control)
		}

		# Avoid anything within 500 bp of peak center
		wide_dnase <- resize(dnase, width=1000, fix='center')
		for(j in seq(num)){
			while(TRUE){
				seed <- sample((start(locus)+75):(end(locus)-75), 1)
				try <- resize(GRanges(seqnames(locus), IRanges(start=seed, end=seed)), width=150, fix='center')
				if(try %outside% wide_dnase & try %outside% control){
					break
				}
			}
			control <- c(control, try)
		}

		return(control)
	}, dnase=dnase)
	control <- do.call(c, control.list)
	control$idx <- seq(length(control))
	return(control)
}

get_averages <- function(loci, rdata.path, peaks, controls){
	.get_scores <- function(peaks, tiles, runMat, variable, locus){
		scores.list <- lapply(seq_along(peaks), function(i, tiles, runMat, variable, locus){
			peak <- peaks[i]
			idx <- tiles %over% peak
			
			mergedMat <- mergeMat(runMat, tiles, peak)

			return(data.frame(variable=variable, value=mean(mergedMat, na.rm=T), 
							cov=sum(!is.na(mergedMat)), id=peak$idx, region=locus$name))
		
		}, tiles=tiles, runMat=runMat, variable=variable, locus=locus)

		scores <- do.call(rbind, scores.list)
		return(scores)
	}

	all_loci.list <- lapply(seq_along(loci), function(i, loci, rdata.path, peaks, controls){
		locus <- loci[i]
		message(locus$name)
		load(sprintf('%s/%s-runMat.RData', rdata.path, locus$name))
		runMat.open <- runMat.open[open_run.qc < 2,]
		dnase.buffer <- peaks[peaks %over% locus]
		control.buffer <- controls[controls %over% locus]

		if(length(dnase.buffer) == 0){
			return(NA)
		}

		dnase_scores <- .get_scores(dnase.buffer, tiles, runMat.open, 'dnase', locus)
		control_scores <- .get_scores(control.buffer, tiles, runMat.open, 'control', locus)
		
		return(rbind(dnase_scores, control_scores))
	}, loci=loci, rdata.path=rdata.path, peaks=peaks, controls=controls)

	averages <- do.call(rbind, all_loci.list)
	return(averages)
}

## DNase Comparisons ##
control.gm <- generate_controls(loci, dnase.gm)
averages.gm <- get_averages(loci, rdata.gm.path, dnase.gm, control.gm)

control.k562 <- generate_controls(loci, dnase.k562)
averages.k562 <- get_averages(loci, rdata.k562.path, dnase.k562, control.k562)

roc.gm <- roc(as.factor(averages.gm$variable), averages.gm$value)
roc.k562 <- roc(as.factor(averages.k562$variable), averages.k562$value)

pdf('GM12878_K562_ROC.pdf', width=8, height=8)
plot(roc.gm, main='ROC', col='#7671b4', lwd=2.5)
lines(roc.k562, col='#e7218c', lwd=2.5)
legend('bottomright', col=c('#7671b4', '#e7218c'), lwd=2,
	   legend=c(sprintf('GM12878 (AUC = %s)', round(roc.gm$auc,2)),
				sprintf('K562 (AUC = %s)', round(roc.k562$auc,2))))
message(sprintf('Optimal GM12878 threshold: %s', coords(roc.gm, x='best')$threshold))
message(sprintf('Optimal K562 threshold: %s', coords(roc.k562, x='best')$threshold))
dev.off()

# Combined
dnase.c <- mergePeaks(dnase.gm[dnase.gm %over% loci], dnase.k562[dnase.k562 %over% loci])
dnase.c$idx <- seq(length(dnase.c))

control.c <- generate_controls(loci, dnase.c)
averages.gm.c <- get_averages(loci, rdata.gm.path, dnase.c, control.c)
averages.gm.c$cell <- 'GM12878'

averages.k562.c <- get_averages(loci, rdata.k562.path, dnase.c, control.c)
averages.k562.c$cell <- 'K562'

averages.consensus <- inner_join(averages.gm.c, averages.k562.c, by=c('variable', 'id'))
averages.consensus$active.x <- mapply(function(var, id){
	if(var == 'control'){
		return('OFF')
	}
	if(dnase.c[id] %over% dnase.gm){
		return('ON')
	}
	return('OFF')
}, averages.consensus$variable, averages.consensus$id)
averages.consensus$active.y <- mapply(function(var, id){
	if(var == 'control'){
		return('OFF')
	}
	if(dnase.c[id] %over% dnase.k562){
		return('ON')
	}
	return('OFF')
}, averages.consensus$variable, averages.consensus$id)
averages.consensus$int <- interaction(averages.consensus$active.x, 
									  averages.consensus$active.y, sep='/')

pdf('GM12878_K562_Scatter.pdf', width=12, height=8)
pmain <- ggplot(averages.consensus, aes(x=value.x, y=value.y, col=int)) + 
		geom_point(size=3, alpha=0.8) + 
		scale_color_manual(values=c("#a9a9a9","#7671b4", "#e7218c", "#65a743"),
						   name='GM12878/K562 Activity') + 
		ylim(0,0.8) + xlim(0,0.8) + theme_bw(base_size=20) + 
		xlab('GM12878 Open Run Signal') + ylab('K562 Open Run Signal') +
		theme(panel.grid.minor = element_blank())
print(pmain)

# Exclude Off/Off
averages.consensus <- averages.consensus[!(averages.consensus$active.x == 'OFF' & averages.consensus$active.y == 'OFF'),]
pmain <- ggplot(averages.consensus, aes(x=value.x, y=value.y, col=int)) + 
		 geom_point(size=3, alpha=0.8) + 
		 scale_color_manual(values=c("#a9a9a9","#7671b4", "#e7218c", "#65a743"),
						   name='GM12878/K562 Activity') + 
		 ylim(0,0.8) + xlim(0,0.8) + theme_bw(base_size=20) + 
		 xlab('GM12878 Open Run Signal') + ylab('K562 Open Run Signal') +
		 theme(panel.grid.minor = element_blank())
xdens <- axis_canvas(pmain, axis='x') + 
		 geom_density(data=averages.consensus[averages.consensus$cell.x=='GM12878',], 
					  aes(x=value.x, fill=active.x), alpha=0.7, size=0.2) + 
		 scale_fill_manual(values=c("gray","black"), name='GM12878/K562 Activity')
ydens <- axis_canvas(pmain, axis='y', coord_flip=T) + 
		 geom_density(data=averages.consensus[averages.consensus$cell.y=='K562',], 
		 			  aes(x=value.y, fill=active.y), alpha=0.7, size=0.2) +
		 scale_fill_manual(values=c("gray","black"), name='GM12878/K562 Activity') + 
		 coord_flip()

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(0.2, 'null'), position='top')
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(0.2, 'null'), position='right')
ggdraw(p2)

dev.off()
