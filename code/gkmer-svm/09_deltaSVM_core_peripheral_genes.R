library(BuenColors)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(annotables)
"%ni%" <- Negate("%in%")

# Load deltaSVM data
source("99_merge_deltaSVM_results.R")

# Load finemap data and filter for PP>0.5
CS.df <- fread(cmd="zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz") %>% filter(PP>0.5)
# Merge finemap with deltaSVM
all.CS_merged <- left_join(CS.df[,c("UKID","trait","PP","AF_Allele2")],all,by=c("UKID"="var")) %>% unique()

# Exclude variants with coding consequences
coding_consequences <- c("missense_variant","synonymous_variant","frameshift_variant",
                         "splice_acceptor_variant","splice_donor_variant","splice_region_variant",
                         "inframe_insertion","stop_gained","stop_retained_variant",
                         "start_lost","stop_lost","coding_sequence_variant","incomplete_terminal_codon_variant")

# Read in Dragana's VEP worst consequence gene annotations for condind and/or PP>0.5 variants
v_to_gene <- fread("../../data/gene_annotations/Condind_vars_and_pp_gt_0.5_with_clinical_annotation_20_11_19.txt") %>%
  dplyr::rename(UKID="VAR",trait="pheno")  %>% filter(PP %ni% c("<0.5"),VEP_cons_most_serious %ni% coding_consequences)
v_to_gene[is.na(v_to_gene)] <- ""
v_to_gene <- v_to_gene %>% mutate(intergenic_or_not = ifelse(gene == "","intergenic","genic")) %>% 
  mutate(VEP_cons_most_serious = ifelse(VEP_cons_most_serious=="","intergenic_variant",VEP_cons_most_serious)) %>% unique()

all.CS_merged.gene <- merge(all.CS_merged,v_to_gene[,c("UKID","trait","beta","SD","VEP_cons_most_serious","gene","BRIDGE_gene","condind","intergenic_or_not")],by=c("UKID","trait"))

# Remove duplicate rows for the same variant saying it's both intergenic and genic
genic_vars <- all.CS_merged.gene %>% filter(intergenic_or_not == "genic") %>% .$UKID %>% unique()
intergenic_vars <- all.CS_merged.gene %>% filter(intergenic_or_not == "intergenic") %>% .$UKID %>% unique()
all.CS_merged.gene <-bind_rows(all.CS_merged.gene %>% filter(UKID %ni% intersect(genic_vars,intergenic_vars)),
                           all.CS_merged.gene %>% filter(UKID %in% intersect(genic_vars,intergenic_vars), intergenic_or_not == "genic"))

# Combine deltaSVM scores with gene assignments
max_vars <- all.CS_merged.gene %>% group_by(UKID) %>%
  summarise(maxPP=max(PP),maxSVM = max(abs(deltaSVM)),medSVM = median(abs(deltaSVM)),
            maxbeta = max(abs(beta)))
max_vars_nodups <- merge(max_vars, all.CS_merged.gene[,c("UKID","AF_Allele2","VEP_cons_most_serious","gene","BRIDGE_gene","intergenic_or_not")],by="UKID") %>% unique()

# Compare intergenic to genic variants
gene_deltaSVM <- maxvars_nodups %>% filter(intergenic_or_not == "genic") %>% .$maxSVM 
nogene_deltaSVM <- maxvars_nodups %>% filter(intergenic_or_not != "genic") %>% .$maxSVM 

t.test(x=gene_deltaSVM,y=nogene_deltaSVM,alternative="two.sided")

ggplot(maxvars_nodups, aes(x=maxSVM)) +
  geom_density(aes(fill=intergenic_or_not),alpha=0.4) +
  scale_fill_manual(values=jdb_palette("brewer_spectra")[c(1,4)]) +
  pretty_plot(fontsize = 8)+ L_border() +
  labs(x="deltaSVM",y="Density") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

# Compare core vs. peripheral variants
maxvars_nodups %>% filter(intergenic_or_not == "genic") %>%
  group_by(BRIDGE_gene) %>% summarise(mean(maxSVM))

# Stratify on type of disease caused by core genes
core_genes <- fread("../../data/gene_annotations/bridge_genes.txt")
BMF_genes <- core_genes %>% filter(Platform == "BMF") %>% .$Gene_Symbol_HGNC
BPD_genes <- core_genes %>% filter(Platform == "BPD") %>% .$Gene_Symbol_HGNC
SMD_genes <- core_genes %>% filter(Platform == "SMD") %>% .$Gene_Symbol_HGNC

all.CS_merged.gene %>% filter(gene %in% BPD_genes) %>% group_by(celltype) %>% 
  summarise(median(abs(deltaSVM))) %>% as.data.frame()

all.CS_merged.gene %>% filter(gene %in% SMD_genes) %>% group_by(celltype) %>% 
  summarise(median(abs(deltaSVM))) %>% as.data.frame()

all.CS_merged.gene %>% filter(gene %in% BMF_genes) %>% group_by(celltype) %>% 
  summarise(median(abs(deltaSVM))) %>% as.data.frame()

