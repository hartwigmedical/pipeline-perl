### run_QDNAseq.R
### Wrapper to perform CNV analysis on a set of bam files using the QDNAseq R package.

library(GetoptLong)
library(devtools)

library(Biobase)

GetoptLong(c(
    "bams|b=s@", "Bam files",
    "sampleNames|s=s@", "Sample names, same order as bams",
    "qdnaseq_path|qdnaseq=s", "Path to QDNAseq source"
))

## to have easier control over versions
load_all(qdnaseq_path)

data("b5.cl80")
bins <- b5.cl80

if (!file.exists("readCountsFiltered.rds")) {
    readCounts <- binReadCounts(bins, bams, cache=F, chunkSize=TRUE)

    readCountsFiltered <- readCounts
    readCountsFiltered <- applyFilters(readCountsFiltered, residual=FALSE, blacklist=TRUE, mappability=FALSE, bases=FALSE)
    readCountsFiltered <- estimateCorrection(readCountsFiltered)
    readCountsFiltered <- applyFilters(readCountsFiltered, residual=FALSE, blacklist=TRUE, mappability=FALSE, bases=FALSE, chromosomes=NA)

    readCountsFiltered <- correctBins(readCountsFiltered)
    readCountsFiltered <- normalizeBins(readCountsFiltered)
    readCountsFiltered <- smoothOutlierBins(readCountsFiltered)

    #readCountsFiltered <- normalizeSegmentedBins(segmentBins(readCountsFiltered, alpha=0.01, undo.splits="none"))
    readCountsFiltered <- normalizeSegmentedBins(segmentBins(readCountsFiltered, segmentStatistic="seg.median", alpha=0.01, undo.splits="none"))
    saveRDS(readCountsFiltered, "readCountsFiltered.rds")
} else {
    message("using cached read counts")
    readCountsFiltered <- readRDS("readCountsFiltered.rds")
}

png("copyNumberSegmented.png", width=1024, height=1024)
plot(readCountsFiltered)
dev.off()

exportBins(readCountsFiltered, "copynumber.igv", type="copynumber", format="igv")
exportBins(readCountsFiltered, "segments.igv", type="segments", format="igv")

allCalls <- callBins(readCountsFiltered, method="cutoff")
exportBins(allCalls, "calls.igv", type="calls", format="igv", logTransform=FALSE)
exportVCF(allCalls)

png("calls.png", width=1024, height=1024)
plot(allCalls)
dev.off()
