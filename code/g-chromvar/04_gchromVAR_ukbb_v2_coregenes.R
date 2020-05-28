library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(diffloop)
library(Matrix)
library(cowplot)
library(BuenColors)
library(data.table)
library(tidyverse)

set.seed(1026)

run_gchromVAR <- function(SE,name,atac_data,df_path){
  bcx <- importBedScore(rowRanges(SE), list.files(paste0("../../data/finemap_bedfiles/ukbb_v2_gene_sets/",name,"/"), 
                                                  full.names = TRUE, pattern = "*.bed$"),colidx = 5)
  
  # Compute weighted deviation scores using gchromVAR
  bg <- getBackgroundPeaks(SE,niterations=100)
  ukbb_wDEV <- computeWeightedDeviations(SE, bcx, background_peaks = bg)
  
  zscoreWeighted <- melt(t(assays(ukbb_wDEV)[["z"]])) %>% filter(.,!grepl("dup",Var2))
  colnames(zscoreWeighted) <- c("Celltype","Trait","zscore")
  zscoreWeighted$logp <- -log10(pnorm(zscoreWeighted$zscore, lower.tail = FALSE))
  zscoreWeighted$Trait <- gsub(".PP001","",gsub(paste0(name,"_"),"",zscoreWeighted$Trait))
  
  bonferroni_p <- -log10(0.05 / (length(unique(zscoreWeighted$Celltype))*length(unique(zscoreWeighted$Trait))))
  
  p1 <- ggplot(zscoreWeighted, aes(x = Celltype, y = logp)) +
    geom_bar(width = 1, aes(fill = Celltype), colour="black",
             stat = "identity", position = position_dodge(width=1)) +
    pretty_plot(fontsize = 10) + 
    labs(x = "", y = "g-chromVAR Enrichment (-log10 p)", fill = "") +
    scale_fill_manual(values = colormap) + 
    facet_wrap(~Trait, scales = "free_y") +
    theme(legend.position="bottom") +
    geom_hline(yintercept = -log10(0.05 / (18*length(unique(zscoreWeighted$Trait)))), linetype = 2) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(size=7)) + 
    guides(shape = guide_legend(override.aes = list(size = 2)))
  
  cowplot::ggsave2(p1, file=paste0(df_path,"g-chromVAR_",name,"-",atac_data,".pdf"),height=6,width=6)
  
  # Make table and export
  write.table(zscoreWeighted, file=paste0(df_path,"gchromVAR_zscores_",name,"-",atac_data,".tsv"), 
              sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
  
  if (TRUE){
    # Radial barplot
    # Cap logp at ceiling
    ylimit <- 10
    zscoreWeighted <- zscoreWeighted %>% mutate(logp = ifelse(logp > ylimit,ylimit,logp))
    
    # Reorder traits and cell types
    ct_order <- c("HSC","MPP","CMP","MEP","Ery","Mega",
                  "GMP-A","GMP-B","GMP-C","Mono","mDC",
                  "LMPP","CLP","CD4","CD8","B","NK","pDC")
    zscoreWeighted$Celltype <- factor(zscoreWeighted$Celltype,levels=ct_order)
    
    p2 <- ggplot(zscoreWeighted,aes(x=Trait,y=logp,group=Celltype)) + 
      geom_bar(aes(fill = Celltype),stat = "identity",position="dodge") +
      scale_fill_manual(values = ejc_color_maps) + 
      geom_hline(yintercept = bonferroni_p, linetype = 2) +
      ylim(-7,ylimit+0.1) +
      coord_polar(start = 0) + theme_minimal()+
      theme(axis.title.x =element_blank(), legend.position="right")
    
    cowplot::ggsave2(p2,file=paste0(df_path,"g-chromVAR_radial_",name,"-",atac_data,".pdf"),
                     height=4,width=5)
  }
} 


# Read in ATAC data
peaksdf <- fread("../../data/panHeme/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/panHeme/29August2017_EJCsamples_allReads_500bp.counts.txt"))

# Create bulk Summarized Experiment from the ATAC data
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

core_gene_platforms <- c("BMF","BPD","SMD","any_genes")
colormap <- ejc_color_maps

# Run 
out_path <- "../../output/gchromvar/ukbb_gene_sets/"


for (i in 1:length(core_gene_platforms)){
  print(core_gene_platforms[i])
  run_gchromVAR(SE=SE,name=core_gene_platforms[i],atac_data=atac_data[1],df_path = out_path)
}

# compare enrichments for BPD before and after subsetting -----------------
colormap <- ejc_color_maps
traits_of_interest <- c("mpv","plt","pdw","pct")

all <- fread("../../output/gchromvar/gchromVAR_zscores_ukbb_v2.txt") %>% mutate(type = "all") %>%
  dplyr::rename(logp = "pval") %>% filter(Trait %in% traits_of_interest) %>% arrange(Celltype, Trait)

any <- fread("../../output/gchromvar/ukbb_gene_sets/gchromVAR_zscores_any_genes-heme_pops.tsv") %>% 
  mutate(type = "any_genes") %>% filter(Trait %in% traits_of_interest)%>% arrange(Celltype, Trait)

# BPD
BPD <- fread("../../output/gchromvar/ukbb_gene_sets/gchromVAR_zscores_BPD-heme_pops.tsv") %>% 
  mutate(type = "BPD")%>% filter(Trait %in% traits_of_interest)%>% arrange(Celltype, Trait)

merged <- bind_rows(all,BPD,any) %>% filter(Trait %in% traits_of_interest)%>% arrange(Celltype, Trait)

bonferroni_p <- -log10(0.05 / (length(unique(merged$Celltype))*length(unique(traits_of_interest))))

p1 <- ggplot(merged, aes(x = Celltype, y = logp)) +
  geom_bar(width = 1, aes(fill = Celltype), colour="black",
           stat = "identity", position = position_dodge(width=1)) +
  pretty_plot(fontsize = 10) + 
  labs(x = "", y = "g-chromVAR Enrichment (-log10 p)", fill = "") +
  scale_fill_manual(values = colormap) + 
  facet_grid(type ~ Trait)+
  theme(legend.position="bottom") +
  geom_hline(yintercept = bonferroni_p, linetype = 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=7)) + 
  guides(shape = guide_legend(override.aes = list(size = 2)))

p1

cowplot::ggsave2(p1, file="../../output/gchromvar/ukbb_gene_sets/platelet_trait_BPD_vs_all.pdf",
                 height = 5,width=4.5)

 
merged <- bind_rows(all,BPD,any) %>% filter(Trait %in% traits_of_interest)
merged_zscores <- pivot_wider(merged,id_cols = c(Celltype,Trait),
                              names_from = type, values_from = c(zscore,logp)) %>% as.data.frame()

# Plot z-scores
bonferroni_z <- 3.2

p1 <- ggplot(merged_zscores, aes(x = zscore_all, y = zscore_BPD, color = Celltype)) +
  labs(x = "g-chromVAR z-score (all variants)", y = "g-chromVAR z-score (BPD variants)", color = "") +
  geom_point() +
  scale_color_manual(values = ejc_color_maps) +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "bottom") + 
  geom_abline(intercept = 0, slope = 1, linetype = 1) +
  geom_hline(yintercept = bonferroni_z, linetype = 3) +
  geom_vline(xintercept = bonferroni_z, linetype = 3) 
  
p1

cowplot::ggsave2(p1, file="../output/gchromvar/ukbb_gene_sets/platelet_trait_BPD_vs_all_zscores.pdf",
                 height = 4,width=4.5)

