# Call this script with loci BED filename (args[1])
# and methylation calls file (args[2]) 

args <- commandArgs(trailingOnly=T)

suppressMessages(library("data.table"))
suppressMessages(library("rtracklayer"))

loci <- import.bed(args[1])
calls <- fread(args[2])

gr <- GRanges(calls$chromosome, IRanges(calls$start + 1, calls$start))

subset <- calls[gr %over% loci,]

write.table(subset, file=basename(args[2]), row.names=F, quote=F, sep='\t')