library(BuenColors)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(annotables)
"%ni%" <- Negate("%in%")

# Load motifbreakr data
mbreaker.df <- fread("zcat < ../gwas_data/ukbb_mbreaker_finemap_merged.txt.gz")

# Load core BRIDGE flagship genes
core_genes <- read.table("../gwas_data/BRIDGE_Flagship_gene_list.txt")[,1]
peripheral_genes <- read.table("../gwas_data/List_gwas_genes_by_VEP_with_intergenic.txt")[,1] 
peripheral_genes <- peripheral_genes[complete.cases(peripheral_genes)]
peripheral_genes <- setdiff(peripheral_genes,core_genes)
intersect(core_genes,peripheral_genes)

# Read in Dragana's VEP worst consequence gene annotations for condind and/or PP>0.5 variants
v_to_gene <- fread("../gwas_data/Condind_vars_and_pp_gt_0.5_with_clinical_annotation_20_11_19.txt") %>%
  dplyr::rename(UKID="VAR",trait="pheno") 
v50 <- v_to_gene %>% filter(PP %ni% c("<0.5"))
v50_core <- v50 %>% filter(BRIDGE_gene == "yes")

# See if motif disrupting variants are enriched for being near core genes
CS.df <- fread("../../bcx-finemap/data/finemap_bedfiles/old/ukbb/bcx.CS.df") 
names(CS.df) <- c("seqnames","start","end","UKID","PP","region","trait")

bridge_motif <- mbreaker.df %>% filter(PP > 0.5, UKID %in% v50_core$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_motif <- mbreaker.df %>% filter(PP > 0.5, UKID %ni% v50_core$UKID) %>% .$UKID %>% unique() %>% length()
bridge_nomotif <- v50_core %>% filter(UKID %ni% mbreaker.df$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_nomotif <- CS.df %>% filter(PP > 0.5,UKID %ni% c(mbreaker.df$UKID,v50_core$UKID)) %>% .$UKID %>% unique() %>% length()

counts <- c("bm"=bridge_motif,"no_b_m"=nobridge_motif,"b_no_m"=bridge_nomotif,"no_b_no_m"=nobridge_nomotif)
counts
fisher.test(matrix(counts, nrow = 2))

# Compare new CS to old CS (motifbreaker was initially run on old CS)
CS.df.new <- fread(cmd="zcat < ../gwas_data/CS.ukid.ukbb_v2.PP0.001.txt.gz")
nobridge_nomotif <- CS.df %>% filter(PP > 0.5,UKID %ni% c(mbreaker.df$UKID,v50_core$UKID)) %>% .$UKID %>% unique()
nobridge_nomotif.2 <- CS.df.new %>% filter(PP > 0.5,UKID %ni% c(mbreaker.df$UKID,v50_core$UKID)) %>% .$UKID %>% unique()


# See if disrupted motifs are enriched for being core genes
bridge_motif <- mbreaker.df %>% filter(PP > 0.5, geneSymbol %in% v50_core$gene) %>% .$geneSymbol %>% unique() %>% length()
nobridge_motif <- mbreaker.df %>% filter(PP > 0.5, geneSymbol %ni% v50_core$gene) %>% .$geneSymbol %>% unique() %>% length()
bridge_nomotif <- v50_core %>% filter(gene %ni% mbreaker.df$geneSymbol) %>% .$UKID %>% unique() %>% length()
nobridge_nomotif <- v50 %>% filter(gene %ni% c(mbreaker.df$geneSymbol,BRIDGE_gene == "no")) %>% .$gene %>% unique() %>% length()

counts <- c("bm"=bridge_motif,"no_b_m"=nobridge_motif,"b_no_n"=bridge_nomotif,"no_b_no_m"=nobridge_nomotif)
counts
fisher.test(matrix(counts, nrow = 2))

