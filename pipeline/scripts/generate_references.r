library("BSgenome.Hsapiens.UCSC.hg38")
genome <- BSgenome.Hsapiens.UCSC.hg38

# Save file paths
gpc.out <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_gpc.rds'
gpc_tiles.out <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_gpc_tiles.rds'
cpg.out <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_cpg.rds'
cpg_tiles.out <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/all_cpg_tiles.rds'
gcg.out <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/gcg_strands.rds'

# Identify GpCs (only one strand necessary due to symmetry)
gpc.list <- lapply(paste0("chr", c(seq(22), "X", "Y", "M")), function(chr) {
	message(chr)
	result <- GRanges(chr, matchPattern("GC", genome[[chr]], fixed='subject')@ranges)
	result <- resize(result, width=1, fix='end')
	return(result)
})
names(gpc.list) <- paste0("chr", c(seq(22), "X", "Y", "M"))

gpc <- sort(Reduce(c, gpc.list))
seqinfo(gpc) <- BSgenome.Hsapiens.UCSC.hg38@seqinfo
gpc <- keepSeqlevels(gpc, names(gpc.list))
saveRDS(gpc, file=gpc.out)

# Extend ranges
gpc_tiles.list <- lapply(gpc.list, function(gr){
	mid <- floor((start(gr)[-1] - end(gr)[-length(gr)]) / 2 + end(gr)[-length(gr)])
	end(gr)[-length(gr)] <- mid
	start(gr)[-1] <- mid + 1
	return(gr)
})
names(gpc_tiles.list) <- paste0("chr", c(seq(22), "X", "Y", "M"))

gpc_tiles <- sort(Reduce(c, gpc_tiles.list))
seqinfo(gpc_tiles) <- BSgenome.Hsapiens.UCSC.hg38@seqinfo
gpc_tiles <- keepSeqlevels(gpc_tiles, names(gpc_tiles.list))
saveRDS(gpc_tiles, file=gpc_tiles.out)

# Identify CpGs (only one strand necessary due to symmetry)
cpg.list <- lapply(paste0("chr", c(seq(22), "X", "Y", "M")), function(chr) {
  	message(chr)
	result <- GRanges(chr, matchPattern("CG", genome[[chr]], fixed='subject')@ranges)
	result <- resize(result, width=1, fix='start')
	return(result)
})
names(cpg.list) <- paste0("chr", c(seq(22), "X", "Y", "M"))

cpg <- sort(Reduce(c, cpg.list))
seqinfo(cpg) <- BSgenome.Hsapiens.UCSC.hg38@seqinfo
cpg <- keepSeqlevels(cpg, names(cpg.list))
saveRDS(cpg, file=cpg.out)

# Extend ranges
cpg_tiles.list <- lapply(cpg.list, function(gr){
	mid <- floor((start(gr)[-1] - end(gr)[-length(gr)]) / 2 + end(gr)[-length(gr)])
	end(gr)[-length(gr)] <- mid
	start(gr)[-1] <- mid + 1
	return(gr)
})
names(cpg_tiles.list) <- paste0("chr", c(seq(22), "X", "Y", "M"))

cpg_tiles <- sort(Reduce(c, cpg_tiles.list))
seqinfo(cpg_tiles) <- BSgenome.Hsapiens.UCSC.hg38@seqinfo
cpg_tiles <- keepSeqlevels(cpg_tiles, names(cpg_tiles.list))
saveRDS(cpg_tiles, file=cpg_tiles.out)

# Identify GCGs on plus and minus strand
gcg.list <- lapply(paste0("chr", c(seq(22), "X", "Y", "M")), function(chr) {
	message(chr)
	plus <- GRanges(chr, matchPattern("GCG", genome[[chr]], fixed='subject')@ranges, '+')
	minus <- GRanges(chr, matchPattern("CGC", genome[[chr]], fixed='subject')@ranges, '-')
	c(plus, minus)
})
names(gcg.list) <- paste0("chr", c(seq(22), "X", "Y", "M"))

gcg <- sort(Reduce(c, gcg.list))
seqinfo(gcg) <- BSgenome.Hsapiens.UCSC.hg38@seqinfo

gcg <- keepSeqlevels(gcg, names(gcg.list))

gcg <- resize(gcg, width=1, fix='center')

saveRDS(gcg, file=gcg.out)