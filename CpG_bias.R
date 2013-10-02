library(Rsamtools)
library(GenomicRanges)

input <- commandArgs(TRUE)[1]
genome <- commandArgs(TRUE)[2]

stopifnot(file.exists(input))
stopifnot(file.exists(paste("/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/", genome, "/bismark_2/", genome, "_CpG_annot.Rdata", sep="")))

rs <- unlist(grglist(readBamGappedAlignments(input)))
load(paste("/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/", genome, "/bismark_2/", genome, "_CpG_annot.Rdata", sep=""))

cat(sapply(CpG_annot, function(x) mean(countOverlaps(x, rs)))[c("Island", "Shore", "Other")], "\n")

