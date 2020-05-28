library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(diffloop)
library(Matrix)
library(cowplot)
library(BuenColors)
"%ni%" <- Negate("%in%")

# Read in variant -> gene annotations
v_to_gene <- fread("../../data/gene_annotations/Condind_vars_and_pp_gt_0.5_with_clinical_annotation_20_11_19.txt") %>%
  dplyr::rename(UKID="VAR",trait="pheno")  

v50 <- v_to_gene %>% filter(gene != "") # Exclude intergenic variants
v50_core <- v50 %>% filter(BRIDGE_gene == "yes")
v50_peripheral <- v50 %>% filter(BRIDGE_gene == "no")

# Read in rare variants
rare_vars_v2 <- readxl::read_xlsx("../../data/ukbb-rare/annotations_rare_vars_updated.xlsx") %>% as_tibble()

# Read in full CS
CS.df <- fread("zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz") %>% filter(PP>0.5)

# See if rare variants are enriched for being near core genes
gene_set <- v50_core
bridge_rare <- v_to_gene %>% filter(UKID %in% rare_vars_v2$var, UKID %in% gene_set$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_rare <- v_to_gene %>% filter(UKID %in% rare_vars_v2$var,UKID %ni% gene_set$UKID) %>% .$UKID %>% unique() %>% length()
bridge_norare <- gene_set %>% filter(UKID %ni% rare_vars_v2$var) %>% .$UKID %>% unique() %>% length()
nobridge_norare <- v_to_gene %>% filter(UKID %ni% c(rare_vars_v2$var,gene_set$UKID)) %>% .$UKID %>% unique() %>% length()

counts <- c("br"=bridge_rare,"no_b_r"=nobridge_rare,"b_no_r"=bridge_norare,"no_b_no_r"=nobridge_norare)
counts
fisher.test(matrix(counts, nrow = 2))

# See if rare variants are enriched for being near peripheral genes
gene_set <- v50_peripheral
bridge_rare <- v_to_gene %>% filter(UKID %in% rare_vars_v2$var, UKID %in% gene_set$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_rare <- v_to_gene %>% filter(UKID %in% rare_vars_v2$var,UKID %ni% gene_set$UKID) %>% .$UKID %>% unique() %>% length()
bridge_norare <- gene_set %>% filter(UKID %ni% rare_vars_v2$var) %>% .$UKID %>% unique() %>% length()
nobridge_norare <- v_to_gene %>% filter(UKID %ni% c(rare_vars_v2$var,gene_set$UKID)) %>% .$UKID %>% unique() %>% length()

counts <- c("br"=bridge_rare,"no_b_r"=nobridge_rare,"b_no_r"=bridge_norare,"no_b_no_r"=nobridge_norare)
counts
fisher.test(matrix(counts, nrow = 2))

# Exclude known pathogenic variants, re-assess enrichment
known_path <- readxl::read_xlsx("../../data/gene_annotations/T1_Known_variants_associated_penetrance_in_ukbb.xlsx") %>% .$Variant %>% unique

intersect(v_to_gene$rsID,known_path) %>% length()

v_to_gene_no_known <- v_to_gene %>% filter(rsID %ni% known_path)  

v50 <- v_to_gene_no_known %>% filter(gene != "") # Exclude intergenic variants
v50_core <- v50 %>% filter(BRIDGE_gene == "yes")
v50_peripheral <- v50 %>% filter(BRIDGE_gene == "no")

gene_set <- v50_core
bridge_rare <- v_to_gene_no_known %>% filter(UKID %in% rare_vars_v2$var, UKID %in% gene_set$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_rare <- v_to_gene_no_known %>% filter(UKID %in% rare_vars_v2$var,UKID %ni% gene_set$UKID) %>% .$UKID %>% unique() %>% length()
bridge_norare <- gene_set %>% filter(UKID %ni% rare_vars_v2$var) %>% .$UKID %>% unique() %>% length()
nobridge_norare <- v_to_gene_no_known %>% filter(UKID %ni% c(rare_vars_v2$var,gene_set$UKID)) %>% .$UKID %>% unique() %>% length()

counts <- c("br"=bridge_rare,"no_b_r"=nobridge_rare,"b_no_r"=bridge_norare,"no_b_no_r"=nobridge_norare)
counts
fisher.test(matrix(counts, nrow = 2))

