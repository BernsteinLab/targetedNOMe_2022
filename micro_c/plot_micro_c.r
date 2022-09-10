# This code plots microC matrices extracted from coolers.

# libraries and functions
library(RColorBrewer)
library(grid)
source('../util/plotting_functions.r')

plotMicroC <- function(path, percentile, col_pal, title){
	# Plots Micro-C matrix
	#
	# Args:
	#	path: path to Micro-C matrix extracted by cooltools
	#	percentile: percentile at which the scale saturates
	#	col_pal: array of colors corresponding to the scale
	#	title: name of plot
	#
	# Returns:
	# 	Plots one half of the Micro-C matrix 

	cool <- as.matrix(read.table(path))

	# Saturate at given percentile
	thresh <- quantile(as.vector(cool[upper.tri(cool, diag=T)]), probs=percentile, na.rm=T)
	cool[cool > thresh] <- thresh
	cool <- cool / thresh

	# Rough numbers to define viewports
	tot <- 20
	margin <- 1
	width <- tot - (4*margin)
	height <- width / 2

	# New plot
	grid.newpage()
	pushViewport(viewport(width=unit(1, 'snpc'), height=unit(1,'snpc'))) # square

	# Triangle
	pushViewport(viewport(x=(3*margin+width/2)/tot, y=(2.5*margin+height)/tot, width=unit(height/tot* 2/sqrt(2),'snpc'), height=unit(height/tot * 2/sqrt(2),'snpc'), angle=45))
	scale <- plotTriMatrix(as.matrix(cool), col_pal)
	popViewport()

	pushViewport(viewport(layout = grid.layout(5,3, heights=c(margin*1.5,height,1.5*margin,height,margin), widths=c(3*margin, width, margin))))

	# Title
	pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	grid.text(title, gp=gpar(fontsize=20))
	popViewport()

	# Scale bar
	pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
	pushViewport(viewport(layout = grid.layout(3,1, heights=c(1.5,2,1))))

	pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
	grid.raster(scale)
	grid.text(label=paste0('-  ', seq(0,percentile,l=4)), x=unit(0.5+1/3*3*margin/tot, 'npc'), y=seq(0,1,l=4), just='left', gp=gpar(fontsize=10))
	popViewport()

	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	grid.text('Scale', y=1/6, just='bottom',  gp=gpar(fontsize=12))
	popViewport()

	popViewport()
	popViewport()
}

# Paths to matrices extracted by cooltools
path1 <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/microC_merged/coolMatrix_chm13_merged_all.txt'
path2 <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/microC_merged/coolMatrix_repaired_merged_all.txt'
path3 <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/microC/coolMatrix_paper_all.txt'
path4 <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/assembly/microC_merged/coolMatrix_repaired_merged.txt'
paths <- c(path1, path2, path3, path4)
titles <-c('CHM13', 'Repaired', 'hg38', 'Repaired')

# Color scale
col_pal <- brewer.pal(9, 'Reds')
percentile <- 0.9

pdf('microc.pdf', width=12, height=12)
for (i in seq_along(paths)){
	plotMicroC(paths[i], percentile, col_pal, titles[i])
}
dev.off()