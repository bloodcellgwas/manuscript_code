library(BuenColors)
library(tidyverse)
library(data.table)

# Load deltaSVM data
source("99_merge_deltaSVM_results.R")

# Read in original CS file
CS.df <- fread(cmd="zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz")
all.CS_merged <- left_join(all,CS.df[,c("UKID","trait","PP","AF_Allele2")],by=c("var"="UKID"))

# Look at general distribution of deltaSVM scores -------------------------
unique_var <- all %>% distinct(var,.keep_all = T)

p1 <- ggplot(unique_var, aes(deltaSVM)) +
  geom_density() +
  pretty_plot(fontsize = 10)+
  L_border() +
  labs(x="deltaSVM",y="Density") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0, linetype = 2)

cowplot::ggsave(p1, file="../../output/plots_svm/deltaSVM_density.pdf",width=4.5,height=3.5)

# Distribution of significant delta scores across cell types
threshold <- quantile(abs(all$deltaSVM),0.99)
toplot <- all %>% filter(deltaSVM > threshold) %>% .$celltype %>% table() %>% as.data.frame() %>% mutate(category = "positive")  %>% dplyr::rename(celltype =".")
toplot <- rbind(toplot,all %>% filter(deltaSVM < -1* threshold) %>% .$celltype %>% table() %>% as.data.frame() %>% mutate(category = "negative",Freq=-1*Freq)%>% dplyr::rename(celltype="."))

p1 <- ggplot(toplot,aes(x=celltype,y=Freq))+
  geom_bar(aes(fill=category),stat="identity",position="identity") +
  pretty_plot() + L_border() + 
  labs(x="")+
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(7,2)])+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

cowplot::ggsave(p1, file="../../output/plots_svm/significant_scores_celltypes.pdf",width=5,height=3)

# Distribution across blood traits
toplot <- all.CS_merged %>% filter(deltaSVM > threshold,PP>0.10) %>% .$trait %>% table() %>% as.data.frame() %>% mutate(category = "positive")  %>% dplyr::rename(trait =".")
toplot <- rbind(toplot,all.CS_merged %>% filter(deltaSVM < -1* threshold,PP>0.10) %>% .$trait %>% table() %>% as.data.frame() %>% mutate(category = "negative",Freq=-1*Freq)%>% dplyr::rename(trait="."))

p1 <- ggplot(toplot,aes(x=trait,y=Freq))+
  geom_bar(aes(fill=category),stat="identity",position="identity") +
  pretty_plot() + L_border() + 
  labs(x="")+
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(7,2)])+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

cowplot::ggsave(p1, file="../../output/plots_svm/significant_scores_traits_PP10.pdf",width=5,height=3)

# Plot deltaSVM vs. MAF ---------------------------------------------------
median_vars <- all.CS_merged %>% mutate(medSVM = ave(abs(deltaSVM),var,FUN=median)) %>%
  distinct(var,.keep_all = T) %>%
  mutate(maf = ifelse(AF_Allele2 > 0.5, 1-AF_Allele2,AF_Allele2))

median_vars <- all.CS_merged %>% 
  group_by(var) %>%
  summarise(AF_Allele2 = median(AF_Allele2),maxSVM = max(abs(deltaSVM)))  %>%
  mutate(maf = ifelse(AF_Allele2 > 0.5, 1-AF_Allele2,AF_Allele2))

# Bin by MAF
bins = c(0,0.01,0.05,0.1,0.3,0.5)
median_vars$maf <- as.numeric(as.character(median_vars$maf))
median_vars$bin <- cut(median_vars$maf, bins,include.lowest=F)


summary(lm(maxSVM ~ maf, median_vars))

limits <- c(as.numeric(quantile(median_vars$maxSVM,0.1)),
            as.numeric(quantile(median_vars$maxSVM,0.9)))
medianline <- median_vars %>% group_by(bin) %>% summarise(median(maxSVM)) %>% .$`median(maxSVM)` %>% min()

p1 <- ggplot(median_vars,aes(x=bin,y=maxSVM)) +
  geom_violin(aes(fill=bin))+
  geom_boxplot(width =0.3)+
  scale_fill_manual(values =  jdb_palette("brewer_spectra")[-5]) +
  pretty_plot(fontsize = 8) + L_border() + 
  labs(x="MAF")+
  coord_cartesian(ylim = limits) +
  geom_hline(yintercept =medianline,linetype="dashed")+
  theme(legend.position="none")
p1

cowplot::ggsave(p1, file="../../output/plots_svm/deltaSVM_vs_MAF.pdf",width=2.5,height=2)

# Plot deltaSVM vs. PP ----------------------------------------------------
# Take the top PP for each variant and the max log2FC (between the ref and mut constructs)
max_vars <- all.CS_merged %>% group_by(var) %>%
  summarise(maxPP=max(PP),maxSVM = max(abs(deltaSVM)),medSVM = median(abs(deltaSVM))) 

summary(lm(maxSVM ~ maxPP, max_vars))
summary(lm(medSVM ~ maxPP, max_vars))

# Bin by PP
bins = c(0,0.01,0.1,0.50,0.75,1.0)
max_vars$PPbin <- cut(max_vars$maxPP, bins,include.lowest=T)

limits <- c(as.numeric(quantile(max_vars$maxSVM,0.2)),
            as.numeric(quantile(max_vars$maxSVM,0.8)))
medianline <- max_vars %>% group_by(PPbin) %>% summarise(median(maxSVM)) %>% .$`median(maxSVM)` %>% min()

p1 <- ggplot(max_vars,aes(x=PPbin,y=maxSVM)) +
  geom_boxplot(aes(color=PPbin)) + 
  scale_color_manual(values =  jdb_palette("brewer_spectra")[-5]) +
  pretty_plot(fontsize = 8) + L_border() + 
  labs(x="Fine-mapped PP")+
  coord_cartesian(ylim = limits) +
  geom_hline(yintercept =medianline,linetype="dashed")  +
  theme(legend.position="none")
p1

cowplot::ggsave(p1, file="../../output/plots_svm/deltaSVM_vs_PP.pdf",width=2.5,height=2)

# Look at erythroid variants ---------------------------------------------
threshold <- quantile(abs(all$deltaSVM),0.99)
high_ery <- all %>% filter(celltype=="Ery",abs(deltaSVM) >threshold)
table(high_ery$trait)

# Gain vs. loss of accessibility ---------------------------------------------------
all %>% filter(deltaSVM > threshold) %>% distinct(var,.keep_all = T)%>% nrow()
all %>% filter(deltaSVM < -1*threshold) %>% distinct(var,.keep_all = T) %>% nrow()
\