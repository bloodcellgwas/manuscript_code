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

# Read in ATAC data
peaksdf <- fread("../../data/panHeme/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../../data/panHeme/29August2017_EJCsamples_allReads_500bp.counts.txt"))
cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
log2cpm <- log2(cpm+1)
# Min / max scale
log2cpm.minmax <- log2cpm / rowMax(log2cpm)

# Exclude variants with coding consequences
coding_consequences <- c("missense_variant","synonymous_variant","frameshift_variant",
                         "splice_acceptor_variant","splice_donor_variant","splice_region_variant",
                         "inframe_insertion","stop_gained","stop_retained_variant",
                         "start_lost","stop_lost","coding_sequence_variant","incomplete_terminal_codon_variant")

# Read in merged gene / CS data
CS.df.merged <- fread("../../data/VEP/CS.ukid.ukbb_v2.VEP_merged.bed")%>%  mutate(chr = paste0("chr",chr)) %>%
  filter(Consequence %ni% coding_consequences) 
CS.gr <- GRanges(CS.df.merged)

# Find overlaps between peaks and all rare variants
idx <- findOverlaps(peaks,CS.gr)

# Construct overlap data.frame and granges with peaks, cpm, and finemap info
atac_overlap <- data.frame(
  CS.df.merged[idx@to,c("UKID","trait","PP","gene","BRIDGE")],
  log2cpm[idx@from,]
) 
atac_overlap_long <- atac_overlap %>% pivot_longer(.,cols = -c("UKID","trait","PP","gene","BRIDGE"),names_to="celltype",values_to = "log2cpm")

v_to_gene <- fread("../../data/gene_annotations/Condind_vars_and_pp_gt_0.5_with_clinical_annotation_20_11_19.txt") %>%
  dplyr::rename(UKID="VAR",trait="pheno")  %>% filter(PP %ni% c("<0.5"), VEP_cons_most_serious %ni% coding_consequences)
gene_atac_merged <- merge(atac_overlap_long,v_to_gene[,c("UKID","VEP_cons_most_serious")],by="UKID") 

# Gene vs. no gene
bridge_rare <- CS.df.merged %>% filter(UKID %in% atac_overlap$UKID, gene != "") %>% .$UKID %>% unique() %>% length()
nobridge_rare <- CS.df.merged %>% filter(UKID %in% atac_overlap$UKID,gene == "") %>% .$UKID %>% unique() %>% length()
bridge_norare <- CS.df.merged %>% filter(gene != "",UKID %ni% atac_overlap$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_norare <- CS.df.merged %>% filter(gene == "",UKID %ni% atac_overlap$UKID)  %>% .$UKID %>% unique() %>% length()

counts <- c("br"=bridge_rare,"no_b_r"=nobridge_rare,"b_no_r"=bridge_norare,"no_b_no_r"=nobridge_norare)
counts
fisher.test(matrix(counts, nrow = 2))

# Cell type specific analysis
bridge <- fread("../../data/gene_annotations//bridge_genes.txt")
table(bridge$Platform)

BMF_genes <- bridge %>% filter(Platform == "BMF") %>% .$Gene_Symbol_HGNC
BPD_genes <- bridge %>% filter(Platform == "BPD") %>% .$Gene_Symbol_HGNC
SMD_genes <- bridge %>% filter(Platform == "SMD") %>% .$Gene_Symbol_HGNC


gene_atac_merged %>% filter(gene %in% BPD_genes) %>% group_by(celltype) %>% summarise(mean(log2cpm))
gene_atac_merged %>% filter(gene %in% BMF_genes) %>% group_by(celltype) %>% summarise(mean(log2cpm))
gene_atac_merged %>% filter(gene %in% SMD_genes) %>% group_by(celltype) %>% summarise(mean(log2cpm))

# Read in variant -> gene annotations
v50 <- v_to_gene %>% filter(gene != "") # Exclude intergenic variants
v50_core <- v50 %>% filter(BRIDGE_gene == "yes")
v50_peripheral <- v50 %>% filter(BRIDGE_gene == "no")

# Read in Finemap data
CS.df <- fread("zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz") %>%  mutate(chr = paste0("chr",chr)) %>%
  filter(PP>0.5) 
CS.gr <- GRanges(CS.df)

# Merge atac overlap with gene annotations
gene_atac_merged <- merge(atac_overlap_long,v_to_gene[,c("UKID","VEP_cons_most_serious","gene")],by="UKID") 

# See if atac variants are enriched for being near core genes
gene_set <- v50
bridge_atac <- gene_set %>% filter(UKID %in% atac_overlap$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_atac <- v_to_gene %>% filter(UKID %in% atac_overlap$UKID,UKID %ni% gene_set$UKID) %>% .$UKID %>% unique() %>% length()
bridge_noatac <- gene_set %>% filter(UKID %ni% atac_overlap$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_noatac <- v_to_gene %>% filter(UKID %ni% c(atac_overlap$UKID,gene_set$UKID)) %>% .$UKID %>% unique() %>% length()

counts <- c("bm"=bridge_atac,"no_b_m"=nobridge_atac,"b_no_m"=bridge_noatac,"no_b_no_m"=nobridge_noatac)
counts
fisher.test(matrix(counts, nrow = 2))

# Core vs. peripheral enrichment
bridge_rare <- CS.df.merged %>% filter(UKID %in% atac_overlap$UKID, BRIDGE == "yes") %>% .$UKID %>% unique() %>% length()
nobridge_rare <- CS.df.merged %>% filter(UKID %in% atac_overlap$UKID,BRIDGE == "no") %>% .$UKID %>% unique() %>% length()
bridge_norare <- CS.df.merged %>% filter(BRIDGE == "yes",UKID %ni% atac_overlap$UKID) %>% .$UKID %>% unique() %>% length()
nobridge_norare <- CS.df.merged %>% filter(BRIDGE == "no",UKID %ni% atac_overlap$UKID)  %>% .$UKID %>% unique() %>% length()

counts <- c("br"=bridge_rare,"no_b_r"=nobridge_rare,"b_no_r"=bridge_norare,"no_b_no_r"=nobridge_norare)
counts
chisq.test(matrix(counts, nrow = 2))

