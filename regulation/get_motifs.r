# This code is used to get motifs of interest and save them.

# libraries and functions
library(rtracklayer)

# Paths to data and references
loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/200913_K562/200913_Loci.bed'
# Folder with motifs split by chromosome
motifs.path <- '/seq/epiprod02/jingyi/nano_nome/210112_motif'
save.filename <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/gm_k562_motifs.RDS'

# Load in data
loci <- import.bed(loci.path)

chrs <- as.character(unique(seqnames(loci)))

# Keep only motifs within enriched loci
motifs.list <- lapply(chrs, function(chr){
	buffer <- import.bed(sprintf('%s/hg38_motif.%s.txt', motifs.path, chr))

	return(buffer[buffer %over% loci])
})

motifs <- do.call(c, motifs.list)
saveRDS(motifs, save.filename, version=2)
