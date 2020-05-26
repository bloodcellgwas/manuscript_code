library(BuenColors)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(cowplot)
library(ggrastr)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
"%ni%" <- Negate("%in%")

# Load deltaSVM data (takes ~20 seconds)
source("99_merge_deltaSVM_results.R")

# Load motifbreakr data
CS_mbreaker <- readRDS("../../data/motifs/ukbb_CS_all_mbreaker.rds")
# Convert ID back to UKID format
names(CS_mbreaker) <- sapply(CS_mbreaker$var, function(var){
  temp <- str_split_fixed(gsub("chr","",var),":",4)
  paste0(temp[1],":",temp[2],"_",temp[3],"_",temp[4])
})

# Load original CS
CS.df <- fread(cmd="zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz")

# Check deltaSVM scores
ukid <- "6:90976768_G_A"
rsid <- CS.df %>% filter(UKID == ukid) %>% pull(rsid) %>% unique

CS.df %>% filter(UKID == ukid)
all %>% filter(var == ukid)

# Plot motifbreakr plot for specific variant
pdf(file=paste0("../../output/processed_peaks_svm/motifbreakr_pwms/",rsid,".MBplot.pdf"), width = 4, height = 6)  
par(cex.main=0.8,mar=c(1,1,1,1))
plotMB(results = CS_mbreaker, rsid = ukid)
dev.off()
