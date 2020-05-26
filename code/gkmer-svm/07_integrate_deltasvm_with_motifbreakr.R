library(BuenColors)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(cowplot)
library(ggrastr)
library(motifbreakR)
"%ni%" <- Negate("%in%")

# Load deltaSVM data (takes ~20 seconds)
source("99_merge_deltaSVM_results.R")

# Load motifbreakr data
mbreaker.df <- fread("../../output//ukbb_mbreaker_finemap_merged.txt.gz")

# Load chip-atlas/motifbreakr overlap
chipatlas_mbreaker_match <- fread("../../output/chip_atlas_mbreaker_overlap.txt") %>% dplyr::rename(UKID="var")

# See how motif scores track with deltaSVM scores
plot_motif_scores <- function(vars,title){
  motifs <- mbreaker.df %>% filter(UKID %in% vars) %>% distinct(UKID,providerName,.keep_all = T)
  toplot <- melt(motifs[,c("pctRef","pctAlt")],measure.vars=c("pctRef","pctAlt"))
  ggplot(toplot,aes(x=variable,y=value))+
    geom_violin(aes(fill=variable)) + 
    geom_boxplot(width=0.1) + 
    scale_fill_manual(values=jdb_palette("brewer_spectra")[c(7,1)]) +
    ggtitle(title)+
    labs(x="",y="") +
    pretty_plot(fontsize=6) + L_border() + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
}

threshold <- quantile(abs(all$deltaSVM),0.99)

gains <- all %>% filter(deltaSVM > threshold) %>% distinct(var) %>% .$var
lost <- all %>% filter(deltaSVM < -1* threshold) %>% distinct(var) %>% .$var

p1 <- plot_motif_scores(gains,"high deltaSVM") + labs(y="normalized motif score")
p2 <- plot_motif_scores(lost,"low deltaSVM")
combined <- cowplot::plot_grid(p1,p2)
if (FALSE){
  cowplot::ggsave(combined, file="../../output/plots_svm/gainers_motifbreakr_scores.pdf",width=2.5,height=1.5)
}

# Master TF motif creations --------------------------------------------------
plot_TF <- function(TF,ct,noncelltype,tosave = TRUE){
  # creations <- chipatlas_mbreaker_match %>% filter(geneSymbol == "GATA1", pctAlt > pctRef,effect=="strong") %>% distinct(UKID) %>%.$UKID
  gain_motifs <- mbreaker.df %>% filter(geneSymbol == TF, pctAlt > pctRef,effect=="strong") %>% 
    distinct(UKID) %>%.$UKID 
  lost_motifs <- mbreaker.df %>% filter(geneSymbol == TF, pctAlt < pctRef,effect=="strong") %>% 
    distinct(UKID) %>%.$UKID 
  
  # chip_motif <- all %>% filter(var %in% creations,celltype==celltype) %>% mutate(type= "motif + chip")
  gained <- all %>% filter(var %in% gain_motifs,celltype %in% ct) %>% mutate(type= "gain_motif")
  lost <- all %>% filter(var %in% lost_motifs,celltype %in% ct) %>% mutate(type= "lost_motif")
  nonmotif <- all %>% filter(celltype %in% ct,var %ni% c(gain_motifs,lost_motifs)) %>% mutate(type="non-motif")
  non_celltypespecific_gained <- all %>% filter(celltype %in% noncelltype,var %in% c(gain_motifs)) %>% mutate(type="non-lineage_gain")
  non_celltypespecific_lost <- all %>% filter(celltype %in% noncelltype,var %in% c(lost_motifs)) %>% mutate(type="non-lineage_lost")
  #ery_gained <- all %>% filter(celltype %in% c("Ery"),var %in% c(gain_motifs)) %>% mutate(type="ery_gain")
  #ery_lost <- all %>% filter(celltype %in% c("Ery"),var %in% c(lost_motifs)) %>% mutate(type="ery_lost")
  
  toplot <- NULL
  toplot <- do.call("rbind",list(gained,lost,nonmotif,non_celltypespecific_gained,non_celltypespecific_lost))
  toplot <- toplot %>% mutate(type = gsub("_"," ", type))
  
  toplot$type <- factor(toplot$type,levels=c("non-motif","gain motif","lost motif","non-lineage gain","non-lineage lost"))
  p1 <- ggplot(toplot,aes(x=type,y=deltaSVM))+
    geom_violin(aes(fill=type)) + 
    geom_boxplot_jitter(width=0.2,outlier.size=0.05, outlier.jitter.width = 0.08, outlier.alpha=1,
                        raster=T, raster.dpi = 400) +
    scale_fill_manual(values=jdb_palette("brewer_spectra")[c(3,7,9,2,1)]) +
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    ggtitle(paste(TF,ct,sep=" - "))+
    pretty_plot(fontsize=8) + L_border() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none")+
    labs(x="")
  
  if (tosave){
    ct <- ct[1]
    cowplot::ggsave(p1, file=paste0("../../output/plots_svm/TF_deltaSVM_boxplots/",TF,"-",ct,"_deltaSVM.pdf"),width=2.5,height=2.5)
  }
  
  return(p1)
}

plot_TF("GATA1",ct="P5",noncelltype=c("CD8","CD8","B"))
plot_TF("CEBPA",ct=c("GMP.A","GMP.B","GMP.C"),noncelltype=c("CD8","CD8","B"))
plot_TF("GABPA",ct=c("CD8"),noncelltype=c("P5","Mega"))

plot_TF("FLI1",ct=c("Mega"),noncelltype=c("Mono","mDC"))

plot_TF("TFAP4",ct=c("B"),noncelltype=c("P5","Mega"))

# Old code incorporating Chip-atlas data.. may not use
if (FALSE){
  gata1_creations <- chipatlas_mbreaker_match %>% filter(geneSymbol == "GATA1", pctAlt > pctRef,effect=="strong") %>% distinct(UKID) %>%.$UKID 
  gata1_motifs <- mbreaker.df %>% filter(geneSymbol == "GATA1", pctAlt > pctRef,effect=="strong") %>% distinct(UKID) %>%.$UKID 
  
  gata1 <- all %>% filter(var %in% gata1_creations,celltype=="Ery") %>% mutate(type= "motif + chip")
  gata1_motif <- all %>% filter(var %in% gata1_motifs,var %ni% gata1_creations,celltype=="Ery") %>% mutate(type= "motif")
  nongata1 <- all %>% filter(celltype=="Ery",var %ni% gata1$var) %>% mutate(type="non-motif")
  nonery <- all %>% filter(var %in% gata1_motifs,celltype==c("CD8","CD8","B")) %>% mutate(type= "non-ery motif")
  
  toplot <- NULL
  toplot <- do.call("rbind",list(gata1,gata1_motif,nongata1,nonery))
  toplot$type <- factor(toplot$type,levels=c("non-motif","motif","motif + chip","non-ery motif"))
  p1 <- ggplot(toplot,aes(x=type,y=deltaSVM))+
    geom_boxplot() + 
    ggtitle("GATA1")+
    pretty_plot() + L_border() +
    labs(x="")
  
  cowplot::ggsave(p1, file="../../output/plots_svm/TF_deltaSVM_boxplots/gata1_creations_deltaSVM.pdf",width=5,height=3)
  
}


# Which TF disruptions have highest deltaSVM changes per cell type? --------

# Define lineage specific attributes
# Ery <- c("HSC", "MPP", "CMP", "MEP", "Ery")
# Meg <- c("HSC", "MPP", "CMP", "MEP", "Mega")
# Mye <- c("HSC", "MPP", "CMP", "LMPP", "GMP-A", "GMP-B", "GMP-C", "Mono", "mDC")
# Lymph <- c("HSC", "MPP", "LMPP", "CLP", "NK", "pDC", "CD4", "CD8", "B")
Ery <- c("P5")
Meg <- c("Mega")
Mye <- c("Mono", "mDC")
Lymph <- c( "NK", "pDC", "CD4", "CD8", "B")

all_lineages <- list(Ery,Meg,Mye,Lymph)
names(all_lineages) <- c("Ery","Meg","Mye","Lymph")

top_motifs <- lapply(seq(1,length(all_lineages)),function(i){
  print(i)
  lineage <- names(all_lineages)[i]
  ct <- all_lineages[[i]]
  ct_specific <- all %>% filter(celltype %in% ct) %>% group_by(var) %>% summarise(maxSVM = max(abs(deltaSVM)))
  motif_merged <- mbreaker.df %>% filter(UKID %in% ct_specific$var) %>% left_join(., ct_specific,by=c("UKID"="var"))
  top_motifs <- motif_merged %>% filter(pctAlt > pctRef,effect=="strong") %>%
    group_by(geneSymbol) %>% summarise(medSVM = median(abs(maxSVM))) %>% arrange(desc(medSVM))
  top_motifs$lineage <- lineage
  return(top_motifs)
})

# Identify lineage-specific TFs
combined <- bind_rows(top_motifs) 
control <- combined %>% filter(lineage != "Lymph") %>% group_by(geneSymbol) %>% summarise(control = median(medSVM))  
left_join(top_motifs[[4]],control,by="geneSymbol") %>% mutate(diff = medSVM - control) %>% arrange(desc(diff))
  
