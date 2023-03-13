##Load required packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(SummarizedExperiment)
library(GenomicRanges)
library(tidyr)

##Function definitions
StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

##Start main
counts_file <- "../data/ATAC_fragcounts_raw.txt"
counts.table <- as.matrix(read.table(counts_file)) #Read the fragment counts
counts.table <- counts.table[which(rowMins(counts.table)>0), ] #Filter to peaks with non-zero fragments in all samples
peaks <- rownames(counts.table)
rownames(counts.table) <- NULL
parent_samples <- paste('par',c('1':'10'), sep='_')
counts.table <- counts.table[ , parent_samples]
fragment_counts <- SummarizedExperiment(assays = list(counts = counts.table),
                                        rowRanges = StringToGRanges(peaks))
counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
#bg <- getBackgroundPeaks(counts)
bg_gc <- getBackgroundPeaks(counts, bias = rowRanges(counts)$bias)
write.table(bg_gc, "bg_gc.dat", quote = FALSE, sep = "\t")
#perm_counts <- getPermutedData(counts, niterations = 2)