library(data.table)
library(tidyverse)
library(GenomicRanges)

# Read in original CS file
CS.df <- fread(cmd="zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz")
all <- left_join(all,CS.df[,c("UKID","trait","PP","AF_Allele2")],by=c("var"="UKID"))

# Read in motifbreakr results
CS_mbreaker <- readRDS("../../data/motifs/ukbb_CS_all_mbreaker.rds")

# Convert ID back to UKID format
CS_mbreaker$var <- sapply(CS_mbreaker$var, function(var){
  temp <- str_split_fixed(gsub("chr","",var),":",4)
  paste0(temp[1],":",temp[2],"_",temp[3],"_",temp[4])
})
names(CS_mbreaker) <- NULL

# Merge with finemap PPs
mbreaker.df <- merge(CS.df[,c("UKID","PP","trait")],as.data.frame(CS_mbreaker),by.x="UKID",by.y="var") 
mbreaker.df$geneSymbol <- as.factor(mbreaker.df$geneSymbol)

write.table(mbreaker.df, file="../../output/ukbb_mbreaker_finemap_merged.txt",quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE) 
