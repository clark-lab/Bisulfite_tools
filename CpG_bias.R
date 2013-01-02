library(Rsamtools)
library(GenomicRanges)

genomes <- c("hg18"="BSgenome.Hsapiens.UCSC.hg18",
             "hg19"="BSgenome.Hsapiens.UCSC.hg19")

input <- commandArgs(TRUE)[1]

stopifnot(file.exists(input))

genome <- commandArgs(TRUE)[2]

stopifnot(genome %in% names(genomes))

rs <- unlist(grglist(readBamGappedAlignments(input)))

#Once off create CpG annotations
if (!file.exists(paste(genome, "CpG_annot.Rdata", sep="_"))) {
    library(rtracklayer)
    library(genomes[genome], character.only=TRUE)

    #Find CpG sites, restrict to chr1-22
    CpGsites <- resize(vmatchPattern("CG", Hsapiens, fixed=TRUE, asRangedData=FALSE), 1, fix="start")
    seqlevels(CpGsites, force=TRUE) <- seqlevels(CpGsites)[1:22]
    CpGsites <- CpGsites[strand(CpGsites)=="+"]
    strand(CpGsites) <- "*"

    #Split into Islands/Shores/Other
    session <- browserSession("UCSC")
    genome(session) <- genome
    CpG <- track(ucscTableQuery(session, "cpgIslandExt"), asRangedData=FALSE)
    seqlevels(CpG, force=TRUE) <- seqlevels(CpG)[1:22]
    CpGshores <- setdiff(resize(CpG, width(CpG)+4000, fix="center"), CpG)

    CpG_annot <- split(CpGsites, ifelse(CpGsites %in% CpG, "Island",
                                 ifelse(CpGsites %in% CpGshores, "Shore", "Other")))

    rm(CpG, CpGshores, CpGsites)
    save(CpG_annot, file=paste(genome, "CpG_annot.Rdata", sep="_"))
} else load(paste(genome, "CpG_annot.Rdata", sep="_"))

cat(sapply(CpG_annot, function(x) mean(countOverlaps(x, rs)))[c("Island", "Shore", "Other")], "\n")

