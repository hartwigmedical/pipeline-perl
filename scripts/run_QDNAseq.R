### run_QDNAseq.R
### Wrapper to perform CNV analysis on a set of bam files using the QDNAseq R package.

library(GetoptLong)
library(devtools)

library(Biobase)

# Command line arguments
GetoptLong(c(
	"bams|b=s@", "Bam files",
	"sampleNames|s=s@", "Sample names, same order as bams",
	"qdnaseq_path|qdnaseq=s", "Path to QDNAseq source"
))

## Load qdnaseq, don't use normal package loading using library to have easier control over versions
load_all(qdnaseq_path)

## Load callable loci bins
data("b5.cl80")
bins <- b5.cl80

#Open bam files and get read counts
readCounts <- binReadCounts(bins, bams, cache=F)

# Filter optimisation
if (!file.exists("readCountsFiltered.rds")) {
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
	readCountsFiltered <- readRDS("readCountsFiltered.rds")
}

#Export bins
exportBins(readCountsFiltered, "copynumber.igv", type="copynumber", format="igv")
exportBins(readCountsFiltered, "segments.igv", type="segments", format="igv")

# Get cnv calls
allCalls <- callBins(readCountsFiltered, method="cutoff")
exportBins(allCalls, "calls.igv", type="calls", format="igv", logTransform=FALSE)
exportVCF(allCalls)
