library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(diffloop)
library(Matrix)
library(cowplot)
library(BuenColors)
library(tidyverse)

set.seed(1026)

# Read in ATAC data
peaksdf <- fread("../../data/panHeme/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/panHeme/29August2017_EJCsamples_allReads_500bp.counts.txt"))

# Create bulk Summarized Experiment from the ATAC data
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

# Import fine-mapped GWAS bedfile 
bcx <- importBedScore(rowRanges(SE), list.files("../../data/finemap_bedfiles/ukbb_v2/", full.names = TRUE, pattern = "*PP001.bed$"),
                      colidx=5)

# Compute weighted deviation scores using gchromVAR
bg <- getBackgroundPeaks(SE)
ukbb_wDEV <- computeWeightedDeviations(SE, bcx, background_peaks = bg)

# Reformat results
zscoreWeighted <- melt(t(assays(ukbb_wDEV)[["z"]]))
zscoreWeighted[,2] <- gsub(".PP001", "", zscoreWeighted[,2])
colnames(zscoreWeighted) <- c("Celltype","Trait","zscore")

zscoreWeighted$pval <- -log10(pnorm(zscoreWeighted$zscore, lower.tail = FALSE))
bonferroni_p <- -log10(0.05 / (length(unique(zscoreWeighted$Celltype))*length(unique(zscoreWeighted$Trait))))

# Export table
write.table(zscoreWeighted, "../..//output/gchromvar/gchromVAR_zscores_ukbb_v2.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Create plots ---------------------------------------------------
if (FALSE){
  # Read in g-chromVAR zscores
  zscores <- fread("../../output/gchromvar/gchromVAR_zscores_ukbb_v2.txt")
  bonferroni_p <- -log10(0.05 / (length(unique(zscores$Celltype))*length(unique(zscores$Trait))))
  
  # Facet bar plot
  p1 <- ggplot(zscoreWeighted, aes(x = Celltype, y = -log10(pval))) +
    geom_bar(width = 1, aes(fill = Celltype), colour="black",
             stat = "identity", position = position_dodge(width=1)) +
    pretty_plot(fontsize = 10) + 
    labs(x = "", y = "g-chromVAR Enrichment (-log10 p)", fill = "") +
    scale_fill_manual(values = ejc_color_maps) + 
    facet_wrap(~Trait, scales = "free_y") +
    theme(legend.position="bottom") +
    geom_hline(yintercept = -log10(0.05 / (18*length(unique(zscoreWeighted$Trait)))), linetype = 2) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(size=7)) + 
    guides(shape = guide_legend(override.aes = list(size = 2)))
  p1
  
  cowplot::ggsave2(p1, file="../../output/gchromvar/g-chromVAR_ukbb_v2.pdf",height=10,width=10)
  
  
  # Radial 
  # Cap logp at ceiling
  ylimit <- 10
  zscores <- zscores %>% mutate(pval = ifelse(pval > ylimit,ylimit,pval))
  
  # Add empty lines to the initial dataset
  if (F){
    to_add = zscores %>% distinct(Trait,.keep_all = TRUE) %>% mutate(Celltype = "", zscore=NA,pval=0)
    zscores=rbind(to_add,zscores)
    zscores=rbind(zscores, to_add)
  }
  
  # Reorder traits and cell types
  ct_order <- c("HSC","MPP","CMP","MEP","Ery","Mega",
                "GMP-A","GMP-B","GMP-C","Mono","mDC",
                "LMPP","CLP","CD4","CD8","B","NK","pDC")
  zscores$Celltype <- factor(zscores$Celltype,levels=ct_order)
  
  trait_order <- c("rbc","mcv","hct","hgb","hlr","irf","mch","mchc","mrv","mscv","rdw_cv","ret",
                   "plt","mpv","pct","pdw",
                   "wbc","baso","eo","mono","neut","lymph")
  zscores$Trait <- factor(zscores$Trait,levels=trait_order)
  
  p2 <- ggplot(zscores,aes(x=Trait,y=pval,group=Celltype)) + 
    geom_bar(aes(fill = Celltype),stat = "identity",position="dodge") +
    scale_fill_manual(values = ejc_color_maps) + 
    geom_hline(yintercept = bonferroni_p, linetype = 2) +
    ylim(-7,ylimit+0.1) +
    coord_polar(start = 0) + theme_minimal()+
    theme(axis.title.x =element_blank(), legend.position="right")
  
  cowplot::ggsave2(p2,file="../../output/gchromvar/g-chromVAR_bcx_ukbb_radial.pdf",height=4,width=5)
  
}
