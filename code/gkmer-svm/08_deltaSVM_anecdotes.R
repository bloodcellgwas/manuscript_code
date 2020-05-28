library(BuenColors)
library(tidyverse)
library(data.table)
library(GenomicRanges)

# Load deltaSVM data
source("99_merge_deltaSVM_results.R")

# Choose Color Palette
accessibility_pallet <- jdb_palette("solar_rojos")

# Import ATAC data
peaksdf <- fread("../../data/atac-svm/26August2017_EJCsamples_allReads_250bp.bed.gz")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")

# Import counts and normalize
counts <-  data.matrix(fread("../../data/atac-svm/26August2017_EJCsamples_allReads_250bp.counts.txt.gz"))
ATAC.cpm <- t(t(counts)/colSums(counts))*1000000
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

# Read in original CS file
CS.df <- fread(cmd="zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz")
motif_breaker <- read.table("../../output/ukbb_mbreaker_finemap_merged.txt.gz", header = TRUE)
all.CS_merged <- left_join(all,CS.df[,c("UKID","rsid","trait","PP","AF_Allele2")],by=c("var"="UKID"))
all.CS_merged2 <- all.CS_merged %>% filter(PP > 0.5 & abs(deltaSVM) > quantile(abs(all$deltaSVM),0.99) )

all.CS_merged3 <- left_join(all.CS_merged2, motif_breaker[,c("UKID","geneSymbol","pctRef","pctAlt","effect")], by = c("var" = "UKID")) 

# Check for lineage specific
# Define lineage specific attributes
Ery <- c("HSC", "MPP", "CMP", "MEP", "Ery")
Meg <- c("HSC", "MPP", "CMP", "MEP", "Mega")
Mye <- c("HSC", "MPP", "CMP", "LMPP", "GMP-A", "GMP-B", "GMP-C", "Mono", "mDC")
Lymph <- c("HSC", "MPP", "LMPP", "CLP", "NK", "pDC", "CD4", "CD8", "B")

make2df <- function(trait, celltypes){
  return(data.frame(Celltype = celltypes, Trait = rep(trait, length(celltypes))))
}

lineageSpecificDF <- rbind(
  make2df("BASO", Mye),make2df("BASO_P", Mye),
  make2df("EOS", Mye),
  make2df("EO", Mye),make2df("EO_P", Mye),
  make2df("NEU", Mye),
  make2df("NEUT", Mye),make2df("NEUT_P", Mye),
  make2df("MONO", Mye),make2df("MONO_P", Mye),
  
  make2df("WBC", unique(c(Mye, Lymph))),
  make2df("LYM", Lymph),
  make2df("LYMPH", Lymph),make2df("LYMPH_P", Lymph),
  
  make2df("PLT", Meg),
  make2df("PCT", Meg),
  make2df("PDW", Meg),
  make2df("MPV", Meg),
  
  make2df("RDW", Ery),
  make2df("RDW_CV", Ery),
  make2df("HCT", Ery),
  make2df("RET", Ery),make2df("RET_P", Ery),
  make2df("MRV", Ery),
  make2df("HGB", Ery),
  make2df("HLR", Ery),make2df("HLR_P", Ery),
  make2df("IRF", Ery),
  make2df("MCH", Ery),
  make2df("MCV", Ery),
  make2df("MSCV", Ery),
  make2df("MCHC", Ery),
  make2df("RBC", Ery)
)

# Return T/F whether rows in df1 are in df2
rowCheck <- function(df1, df2){
  xx <- apply(df1, 1, paste, collapse = "_")
  yy <- apply(df2, 1, paste, collapse = "_")
  return(xx %in% yy)
}

all.CS_merged3$trait <- toupper(all.CS_merged3$trait)
all.CS_merged3$lineageSpecific <- rowCheck(all.CS_merged3[,c("celltype", "trait")], lineageSpecificDF)
all.CS_merged3 <- all.CS_merged3 %>% filter(lineageSpecific == TRUE)

# Require that motif score direction matches deltaSVM
all.CS_merged3 <- all.CS_merged3 %>% mutate(motif_delta_consistent = ifelse((pctAlt-pctRef)*deltaSVM > 0, "consistent","no")) %>%
  filter(motif_delta_consistent == "consistent")

df <- all.CS_merged3[,c("var", "rsid","celltype", "deltaSVM", "trait", "PP", "AF_Allele2","geneSymbol", "pctRef", "pctAlt","effect")]

# Overlap with heme peak set
# Add chromosome, start, end columns
df$chr <- paste0("chr",str_split_fixed(df$var,":",2)[,1])
df$end <- df$start<- as.integer(gsub("_.*","",str_split_fixed(df$var,":",2)[,2]))
idx <- findOverlaps(peaks,GRanges(df))
df_atac <- df[idx@to,]

# Find nearest gene
gencode_gr <- readRDS("../../data/annotations/gencode_filtered.rds")
df_atac$nearest_gene <- gencode_gr[nearest(GRanges(df_atac),gencode_gr,ignore.strand=TRUE),]$gene_name
df_atac <- df_atac %>% dplyr::select(var,rsid,celltype,deltaSVM,trait,PP,AF_Allele2,geneSymbol,nearest_gene,everything())

# Look for good anecdotes
df_atac %>% filter(PP>0.5,effect=="strong") %>% distinct(var,.keep_all = T)

df_atac %>% filter(var == "6:90976768_G_A") %>% distinct(geneSymbol,.keep_all = T)

make_pair_plot <- function(variant){
  
  cellCoordsDF <- data.frame(
    CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "mDC", "Ery", "Mega"),
    x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     4,      8,     10), 
    y = c(10,     8.3,      7,     5,       6.5,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,     2,      2,     2)
  )
  
  # Extract the variant position
  two <- as.numeric(strsplit(strsplit(variant, "_")[[1]][1], ":", 2)[[1]])
  variant_gr <- data.frame(chr = paste0("chr", as.character(two[1])), start = two[2], end = two[2]) %>% makeGRangesFromDataFrame()
  
  # Pull the ATAC se and delta SVM data; make sure that they are ordered properly
  ov <- subjectHits(findOverlaps(variant_gr, peaks))
  if(length(ov) == 1){
    vec <- ATAC.cpm[ov,]
  } else {
    vec <- 0
  }
  dSVM <- all.CS_merged[,c(1,2,3,4)] %>% filter(var == variant) %>% unique() %>% pull(deltaSVM)
  names(dSVM) <- sort(as.character(cellCoordsDF$CellLabel))
  
  # Make one dataframe for plotting
  plotdf <- data.frame(cellCoordsDF,
                       ATAC = vec[as.character(cellCoordsDF$CellLabel)] - min(vec),
                       deltaSVM = dSVM[as.character(cellCoordsDF$CellLabel)]
  )
  
  blank_theme <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  p1 <- ggplot(plotdf, aes(x = x, y = y, fill = ATAC)) + 
    geom_point(size = 6, color = "black", pch=21) + pretty_plot() + L_border() +
    geom_text(aes(label=CellLabel),hjust=0.5, vjust=3, color = "black", size = 2) + 
    scale_fill_gradientn(colors = accessibility_pallet) +
    scale_y_continuous(limits = c(1, 11)) + ggtitle(paste0( variant, " - ATAC")) + blank_theme
  
  # Make even color distribution
  biggest <- max(abs(dSVM))
  range_values <- (dSVM + biggest) / (max(dSVM) - min(dSVM)) # min / max normalize but for symmetry
  color_values <- scales::rescale(c(-1*biggest, range_values, 2*biggest))
  
  p2 <- ggplot(plotdf, aes(x = x, y = y, fill = deltaSVM)) + 
    geom_point(size = 6, color = "black", pch=21) + pretty_plot() + L_border() +
    geom_text(aes(label=CellLabel),hjust=0.5, vjust=3, color = "black", size = 2) + 
    scale_fill_gradient2(low = jdb_palette("solar_glare")[1:4], high = jdb_palette("solar_glare")[6:9]) +
    scale_y_continuous(limits = c(1, 11)) + ggtitle(paste0( variant, " - deltaSVM"))+ blank_theme
  
  cowplot::ggsave(cowplot::plot_grid(p1, p2), filename = paste0("../../output/plots_svm/accessibility_deltaSVM/", gsub(":","-",variant), ".deltaSVM.pdf"), 
                  width = 7, height = 3)
  
}


df %>% filter(PP.x > 0.7 & abs(deltaSVM) > 2 & AF_Allele2 < 0.05  ) %>%
  pull(var) %>% unique() -> allvars
lapply(allvars, make_pair_plot)

df %>% filter((deltaSVM) > 2 & AF_Allele2 < 0.01)
df %>% filter(var == "1:159174683_T_C")


make_pair_plot("6:90976768_G_A")

