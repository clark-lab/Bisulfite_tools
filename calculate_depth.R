input <- commandArgs(TRUE)[1]

stopifnot(file.exists(input))

library(Rsamtools)

rs <- unlist(grglist(readBamGappedAlignments(input)))
cat(sum(as.numeric(width(rs))), "\n")

