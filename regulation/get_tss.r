# This code is used to get transcription start sites and save them.
# It gets transcripts from RefSeq and adds gene names and symbols
# from Ensembl. In addition, it adds APPRIS annotations as well as
# transcript quantifications for GM12878 and K562 (from ENCODE) as 
# well as t0 and t48 T-cells (from RNA-seq)

# libraries and functions
library(biomaRt)
library(data.table)
library(rtracklayer)

# Paths to data and references
refseq.path <- '/seq/epiprod02/kdong/ncbi-genomes-2020-10-13/GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bed'
chromosomes.path <- '/seq/epiprod02/kdong/ncbi-genomes-2020-10-13/GCF_000001405.39_GRCh38.p13_assembly_chromosomes.txt'
biomart.path <- '/seq/epiprod02/kdong/ncbi-genomes-2020-10-13/biomart_xref.tsv'

txQuant.gm.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/ENCODE/gm12878_rna_hg38_rep1_ENCFF688XYV.tsv'
txQuant.k562.path <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/ENCODE/k562_rna_hg38_rep1_ENCFF205KYZ.tsv'

loci.path <- '/seq/epiprod02/Battaglia/NanoNOMe/Tcells_nov2020/TcellNov2020_loci.bed'

txQuant.0.path <- '/seq/epiprod02/jingyi/nome/210316_RNAseq/hg38/0h/isoforms.fpkm_tracking'
txQuant.48.path <- '/seq/epiprod02/jingyi/nome/210316_RNAseq/hg38/48h/isoforms.fpkm_tracking'

tss.f.out <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/tss_gm_k562.rds'
tss.l.out <- '/seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/references/tss_tcell.rds'

# Transcript quantifications
txQuant.gm <- read.table(txQuant.gm.path, header=T, stringsAsFactors=F)
txQuant.gm$stable <- gsub('\\..*','', txQuant.gm$transcript_id)
txQuant.k562 <- read.table(txQuant.k562.path, header=T, stringsAsFactors=F)
txQuant.k562$stable <- gsub('\\..*','', txQuant.k562$transcript_id)

loci <- import.bed(loci.path)
txQuant.0 <- read.table(txQuant.0.path, header=T, stringsAsFactors=F)
txQuant.0.gr.list <- lapply(txQuant.0$locus, function(string){
	split1 <- strsplit(string, ':')
	chr <- split1[[1]][1]
	split2 <- strsplit(split1[[1]][2], '-')
	start <- as.numeric(split2[[1]][1])
	stop <- as.numeric(split2[[1]][2])
	return(GRanges(chr, IRanges(start, stop)))
})
txQuant.0.gr <- do.call(c, txQuant.0.gr.list)

txQuant.48 <- read.table(txQuant.48.path, header=T, stringsAsFactors=F)
txQuant.48.gr.list <- lapply(txQuant.48$locus, function(string){
	split1 <- strsplit(string, ':')
	chr <- split1[[1]][1]
	split2 <- strsplit(split1[[1]][2], '-')
	start <- as.numeric(split2[[1]][1])
	stop <- as.numeric(split2[[1]][2])
	return(GRanges(chr, IRanges(start, stop)))
})
txQuant.48.gr <- do.call(c, txQuant.48.gr.list)

txQuant.0.gr$tx_id <- txQuant.0$tracking_id
txQuant.0.gr$gene_id <- txQuant.0$gene_id
txQuant.0.gr$symbol <- txQuant.0$gene_short_name
txQuant.0.gr$cov <- txQuant.0$coverage
txQuant.0.gr$FPKM <- txQuant.0$FPKM

txQuant.48.gr$tx_id <- txQuant.48$tracking_id
txQuant.48.gr$gene_id <- txQuant.48$gene_id
txQuant.48.gr$symbol <- txQuant.48$gene_short_name
txQuant.48.gr$cov <- txQuant.48$coverage
txQuant.48.gr$FPKM <- txQuant.48$FPKM

txQuant.0.l.list <- lapply(unique(txQuant.0.gr$symbol[txQuant.0.gr %over% loci]), function(gene){
	buffer <- txQuant.0.gr[txQuant.0.gr$symbol == gene]
	return(buffer[which.max(buffer$FPKM)])
})
txQuant.0.l <- do.call(c, txQuant.0.l.list)

txQuant.48.l.list <- lapply(unique(txQuant.48.gr$symbol[txQuant.48.gr %over% loci]), function(gene){
	buffer <- txQuant.48.gr[txQuant.48.gr$symbol == gene]
	return(buffer[which.max(buffer$FPKM)])
})
txQuant.48.l <- do.call(c, txQuant.48.l.list)

txQuant.0.l$name <- txQuant.0.l$symbol
txQuant.48.l$name <- txQuant.48.l$symbol
txQuant.0.l$score <- txQuant.0.l$FPKM
txQuant.48.l$score <- txQuant.48.l$FPKM


# Search https://www.ncbi.nlm.nih.gov/assembly/ for GRCh38.p13
# Click top result for latest assembly then "Download Assembly"
# In the subsequent prompt, select "RefSeq transcript alignments" for the "File type"
# Download and untar (tar -xvf)
# Retrieve the relevant info with the following command:
# bedtools bamtobed -i GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bam > GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bed
refseq <- import.bed(refseq.path)

# Download sequence report for the same assembly (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt)
# Retrieve just the relevant chromosomes with the following command:
# grep 'assembled-molecule' GCF_000001405.39_GRCh38.p13_assembly_report.txt | cut -f7,9,10 > GCF_000001405.39_GRCh38.p13_assembly_chromosomes.txt
chromosomes <- read.table(chromosomes.path)
chromosomes <- chromosomes[-nrow(chromosomes),] # drop chrM
colnames(chromosomes) <- c('refseq', 'size', 'name')
si <- Seqinfo(as.character(chromosomes$name), chromosomes$size, 
			  isCircular=rep(F, nrow(chromosomes)), genome='hg38')

# Fix seqinfo objects
refseq.f <- keepSeqlevels(refseq, chromosomes$refseq, pruning.mode='coarse')
seqinfo(refseq.f, new2old=seq(length(si))) <- si

refseq.f$stable <- gsub('\\..*','', refseq.f$name)

tss <- c(resize(refseq.f, width=1, fix='start'), ignore.mcols=T)
tss.u <- unique(tss)
hits <- findOverlaps(tss.u, tss)

# Cross-reference downloaded from Biomart
# http://uswest.ensembl.org/biomart/martview/d2d026294f6c270d7d94cfd96f01881b
# Choose "Ensembl Genes 101" then "Human Genes (GRCh38.p13)"
# "Attributes" -> "GENES" -> "Gene Stable ID", "Transcript Stable ID", "Chromosome/scaffold name", "Transcription start site (TSS)", "Strand", "APPRIS annotation"
# "Attributes" -> "EXTERNAL" -> "HGNC symbol", "RefSeq mRNA ID", "RefSeq ncRNA ID"
# Export all results to TSV

biomart <- read.table(biomart.path, header=T, stringsAsFactor=F, sep='\t')
colnames(biomart) <- c('gene_id', 'tx_id', 'chr', 'tss', 'strand', 'appris', 
					   'symbol', 'refseq_m', 'refseq_r')

## GM and K562 ##
tss_meta.list <- sapply(seq(length(tss.u)), function(i){
	r_ids <- refseq.f$stable[subjectHits(hits)[queryHits(hits) == i]]
	biomart_ids <- biomart$refseq_m %in% r_ids | biomart$refseq_r %in% r_ids
	e_ids <- unique(biomart$tx_id[biomart_ids])
	row_idx <- txQuant.gm$stable %in% e_ids
	if(sum(row_idx) == 0){
		return(c(NA, NA, NA))
	}
	gene <- unique(biomart$symbol[biomart_ids])
	tpm1 <- sum(txQuant.gm$TPM[row_idx])
	tpm2 <- sum(txQuant.k562$TPM[row_idx])
	return(c(gene, tpm1, tpm2))
})

tss_meta_appris.list <- sapply(seq(length(tss.u)), function(i){
	r_ids <- refseq.f$stable[subjectHits(hits)[queryHits(hits) == i]]
	biomart_ids <- biomart$refseq_m %in% r_ids | biomart$refseq_r %in% r_ids
	return(unique(biomart$appris[biomart_ids]))
})

meta.length <- sapply(tss_meta.list, function(el){length(el)})
tss_meta <- do.call(rbind, tss_meta.list[meta.length==3])
tss.f <- tss.u[meta.length==3]
tss.f$symbol <- tss_meta[,1]
tss.f$GM <- as.numeric(tss_meta[,2])
tss.f$K562 <- as.numeric(tss_meta[,3])

principal <- sapply(tss_meta_appris.list[meta.length == 3], function(el){
	return(sum(grepl('principal', el, ignore.case=T)) > 0)
})
tss.f$principal <- principal

saveRDS(tss.f, tss.f.out, version=2)

## T-cell ##
hits <- findOverlaps(txQuant.0.l, txQuant.48.l, type='equal')

# Remove alternative transcripts and non-coding RNA
picked <- c('ENST00000398532.7', 'ENST00000392320.5', 'ENST00000379888.7', 'ENST00000337386.8', 'ENST00000440181.4', 'ENST00000305798.6', 'ENST00000342854.8', 'ENST00000420343.1', 'ENST00000526143.2', 'ENST00000394089.5')
tss.l <- c(txQuant.0.l[queryHits(hits)],  txQuant.0.gr[match(picked, txQuant.0.gr$tx_id)])
exclude <- grepl('RP', tss.l$symbol) | grepl('CT', tss.l$symbol) | grepl('AC', tss.l$symbol) | tss.l$symbol == 'SCARNA5' | tss.l$symbol == 'SNORA14' | tss.l$tx_id == 'ENST00000577454.4'
tss.l <- tss.l[!exclude]
tss.l.list <- lapply(seq_along(tss.l), function(i){
	buffer <- tss.l[i]
	e_id <- gsub('\\..*','', buffer$gene_id)
	s <- unique(biomart$strand[grepl(e_id, biomart$gene_id)])

	if(length(s) == 0){
		return(NULL)
	}
	if(s == 1){
		end(buffer) <- start(buffer)
		strand(buffer) <- '+'
	} else {
		start(buffer) <- end(buffer)
		strand(buffer) <- '-'
	}
	return(buffer)
})
tss.l2 <- do.call(c, tss.l.list)
tss.l2$V0 <- tss.l2$FPKM
tss.l2$V48 <- txQuant.48.gr$FPKM[match(tss.l2$tx_id, txQuant.48.gr$tx_id)]
tss.l2$id <- seq(length(tss.l2))

saveRDS(tss.l2, tss.l.out, version=2)