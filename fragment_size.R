input <- commandArgs(TRUE)[1]

stopifnot(file.exists(input))

library(Rsamtools)

rs <- readBamGappedAlignments(input, use.names=TRUE, param=ScanBamParam(scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)))

#Check we have paired end reads, and an even number of them
stopifnot(length(rs)>0 & length(rs)%%2==0)

#hand carve paired fragments
idx <- seq.int(2, length(rs), by=2)
stopifnot(all(names(rs)[idx-1]==names(rs)[idx]))

frags <- IRanges(pmin(start(rs)[idx-1], start(rs)[idx], end(rs)[idx-1], end(rs)[idx]),
                 pmax(start(rs)[idx-1], start(rs)[idx], end(rs)[idx-1], end(rs)[idx]))
frags.widths <- table(width(frags))

widths <- data.frame("width"=1:max(width(frags)), times=0)

widths$times[as.integer(names(frags.widths))] <- frags.widths

write.table(widths, row.names=FALSE, col.names=FALSE)

