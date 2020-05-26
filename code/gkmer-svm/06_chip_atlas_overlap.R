library(data.table)
library(tidyverse)

#----
# This script won't run due to large files that aren't sourced in this repository
# but the output file should be present
#---


#' Read in and munge ChIP-Atlas
# Read in meta data
chip_atlas.meta <- fread("../../../mpn-GWAS/data/transcription_factors/experimentList.clean.tab", sep = "\t", header = F)
names(chip_atlas.meta) <- c("id", "genome", "class1", "antigen", "class2", "celltype")
chip_atlas.meta <- chip_atlas.meta %>%
  dplyr::filter(genome == "hg19", class1 == "TFs and others", class2 == "Blood", 
                antigen %ni% c("-", "5-hmC", "5-mC", "Cyclobutane pyrimidine dimers", "Epitope tags", "HIV Tat", "MethylCap", "pFM2")) 

# Read in bed file (large file!)
chip_atlas.bed <- fread("zcat < /Volumes/broad_sankaranlab/ChIPAtlas/Oth.Bld.05.AllAg.AllCell.clean.bed.gz", sep = "\t", header = F)
names(chip_atlas.bed) <- c("seqnames", "start", "end", "id")
chip_atlas.bed <- chip_atlas.bed %>%
  dplyr::filter(id %in% chip_atlas.meta$id)

# Subset to 22 autosomes
setkey(setDT(chip_atlas.bed), id)
setkey(setDT(chip_atlas.meta), id)
chip_atlas.bed <- merge(chip_atlas.bed, chip_atlas.meta)
chip_atlas.bed.sub <- chip_atlas.bed %>% filter(seqnames %in% paste0("chr",seq(22)))

# Convert to GRanges since these are so much faster for overlaps
chip_atlas.bed.sub.gr <- GRanges(chip_atlas.bed.sub)

# Read in motifbreakR
CS_mbreaker <- readRDS("../../data/motifs/ukbb_CS_all_mbreaker.rds")
# Convert ID back to UKID format
CS_mbreaker$var <- sapply(CS_mbreaker$var, function(var){
  temp <- str_split_fixed(gsub("chr","",var),":",4)
  paste0(temp[1],":",temp[2],"_",temp[3],"_",temp[4])
})
names(CS_mbreaker) <- NULL

# Merge with fine-mapped variants
idx <- findOverlaps(CS_mbreaker, chip_atlas.bed.sub.gr)
mbreaker_chip.df <- data.frame(cbind(as.data.frame(CS_mbreaker[idx@from]), 
                               as.data.frame(chip_atlas.bed.sub.gr[idx@to])))

mbreaker_chip_match <- mbreaker_chip.df %>% filter(geneSymbol == antigen) %>%
  dplyr::select(-c(width,strand,seqnames.1,start.1,end.1,width.1,strand.1,genome,class2))

write.table(mbreaker_chip_match,file="../../output/chip_atlas_mbreaker_overlap.txt",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
