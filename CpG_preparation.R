genome <- commandArgs(TRUE)[1]

#Find all CpG sites
library(Biostrings)
library(GenomicRanges)
fasta <- readDNAStringSet(paste(genome, ".fa", sep=""))
CpGs <- vmatchPattern("CG", fasta)
CpGs <- GRanges(rep(names(CpGs), elementLengths(CpGs)), IRanges(unlist(CpGs@ends), width=1))

#Download CpG islands
library(rtracklayer)

session <- browserSession("UCSC")
genome(session) <- genome
CpG <- track(ucscTableQuery(session, "cpgIslandExt"), asRangedData=FALSE)
CpGshores <- setdiff(resize(CpG, width(CpG)+4000, fix="center"), CpG)

CpG_annot <- split(CpGs, ifelse(overlapsAny(CpGs, CpG), "Island",
                         ifelse(overlapsAny(CpGs, CpGshores), "Shore", "Other")))
save(CpG_annot, file=paste(genome, "_CpG_annot.Rdata", sep=""))

