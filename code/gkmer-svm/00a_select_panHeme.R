library(BuenColors)
library(dplyr)
library(DESeq2)
library(BiocParallel)
library(data.table)
library(annotables)

register(MulticoreParam(2))

# Import raw expression values
ATAC.counts <- read.table("../../data/atac-svm/26August2017_EJCsamples_allReads_250bp.counts.txt.gz", header = TRUE)
peaks <- read.table("../../data/atac-svm/26August2017_EJCsamples_allReads_250bp.bed.gz")

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}

# get quantiles for Z scores
ATAC.counts.cpm <- sweep(ATAC.counts, 2, colSums(ATAC.counts), FUN="/") * 1000000

ATAC.counts.filt.Z <- scaleRows(ATAC.counts.cpm)
Qs <- apply(ATAC.counts.filt.Z,2,function(x) {quantile(x,0.97)})

pops <- colnames(ATAC.counts.filt.Z)

lapply(1:ncol(ATAC.counts.filt.Z), function(i){
  write.table(peaks[which(ATAC.counts.filt.Z[,i] > Qs[i]),], 
              file = paste0("../../output/processed_peaks_svm/", pops[i], ".bed"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t") 
  sum(ATAC.counts.filt.Z[,i] > Qs[i])
})

