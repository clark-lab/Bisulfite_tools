#supply bedGraph filename
stopifnot(length(commandArgs(TRUE))==1)
filename <- commandArgs(TRUE)[1]
stopifnot(file.exists(filename))

library(BSgenome.Hsapiens.UCSC.hg19)

#Find all CpGs
CpGs.hg19 <- vmatchPattern("CG", Hsapiens, fixed=TRUE, asRangedData=FALSE)
CpGs.hg19 <- CpGs.hg19[strand(CpGs.hg19)=="+"]

overlapSums <- function(x, y, z) {
    stopifnot(class(x)=="GRanges")
    stopifnot(class(y)=="GRanges")
    stopifnot(class(z)=="numeric" | class(z)=="integer")
    ov <- as.matrix(findOverlaps(y, x))
    ovSums <- rep(as.numeric(0), length(y))
    ovSums[unique(ov[,1])] <- viewSums(Views(z[ov[,2]], ranges(Rle(ov[,1]))))
    ovSums
}

tab <- read.table(filename, header=FALSE, stringsAsFactors=FALSE)
tab <- GRanges(tab$V1, IRanges(tab$V2+1, width=1), C=tab$V5, T=tab$V6, cov=tab$V5+tab$V6)

CpGs <- CpGs.hg19
values(CpGs)$C <- overlapSums(tab, CpGs, values(tab)$C)
values(CpGs)$T <- overlapSums(tab, CpGs, values(tab)$T)
CpGs$cov <- CpGs$C+CpGs$T

save(CpGs, file=gsub(".bedGraph", ".Rdata", filename))

