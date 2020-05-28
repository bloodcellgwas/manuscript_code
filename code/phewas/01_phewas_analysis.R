library(data.table)
library(tidyverse)
library(qvalue)
library(Matrix)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(BuenColors)
library(annotables)
library(ggrastr)
library(ggrepel)
"%ni%" <- Negate("%in%")

# Generate PheWAS input ---------------------------------------------------
if (FALSE){
  # Generate PheWAS input using rare variants supplied by Dragana
  rare_vars_v2 <- readxl::read_xlsx("../../data/ukbb-rare/annotations_rare_vars_updated.xlsx") 
  rare_vars_v2_phewas_input <- rare_vars_v2 %>% dplyr::select(var,chr,pos) %>% unique()
  range(rare_vars_v2$maf)
  fwrite(rare_vars_v2_phewas_input,file="../../data/ukbb-rare/rare_vars_updated_phewas_input.tsv",sep="\t")
  
  rare_vars_v2_phewas_input_leading_zero <- rare_vars_v2_phewas_input %>% mutate(var = ifelse(chr < 10,paste0("0",var),var)) 
  fwrite(rare_vars_v2_phewas_input_leading_zero,file="../../data/ukbb-rare/rare_vars_updated_phewas_input_leading_zeros.tsv",sep="\t")
}


# SAIGE pheWAS ------------------------------------------------------------
# Read phewas results for rare variants
icd <- fread("zcat < ../../output/rare_variants/phewas/SAIGE.allphenos.bcx_rarevars_updated.txt.gz") %>%
  mutate(var = paste(CHROM,POS,REF,ALT,sep=":"), maf = ifelse(af < 0.5, af, 1-af),
         expected_case_minor_AC =  2 * maf * num_cases, logp = -log10(pval)) %>% 
  dplyr::filter(expected_case_minor_AC > 25)

# Read and merge pheno codes 
phecodes <- readxl::read_xlsx("../../data/phewas/saige-phenotype-information.xlsx") %>% 
  dplyr::rename(pheno = `Phenotype Description`)
phecodes$PheCode <- as.double(phecodes$PheCode)
icd <- merge(icd,phecodes[,c("PheCode","pheno")],by.x="phecode",by.y="PheCode")

# Remove redundant phenotypes
to_remove <- c("Breast cancer [female]","Malignant neoplasm of female breast")
icd <- icd %>% filter(pheno %ni% to_remove)
unique(icd$pheno) %>% length()
icd$num_cases %>% range()

# Load sum stats
rare_vars <- readxl::read_xlsx("../../data/ukbb-rare/annotations_rare_vars_updated.xlsx")  %>% 
  mutate(var = gsub("_",":",var)) %>% dplyr::rename(trait = "pheno",blood_beta="beta")

# Assign blood traits to lineages
traits <- names(table(rare_vars$trait))
map <- c("BASO" = "BASO","BASO_P"="BASO", 
         "EO" = "EO","EO_P" = "EO", 
         "NEUT" = "NEUT", "NEUT_P"="NEUT",
         "RBC" = "RBC", "HCT" = "RBC", "HLR" = "RBC","HLR P" = "RBC", "HGB" = "RBC", "IRF" ="RBC","MCH" = "RBC", "MCHC" = "RBC", "MCV" = "RBC","MRV" = "RBC","MSCV" = "RBC","RDW_CV" = "RBC","RET" = "RBC","RET_P" = "RBC",
         "PLT" = "PLT", "MPV" = "PLT", "PCT" = "PLT", "PDW" = "PLT",
         "MONO" = "MONO", "MONO_P" = "MONO",
         "LYMPH" = "LYMPH","LYMPH_P" = "LYMPH",
         "WBC"="WBC") 
names(map) <- gsub(" ","_",names(map)) %>% tolower()
setdiff(traits,names(map))
lapply(traits, function(name){
  
  # Import and filter for PP > 0.1
  t <- rare_vars %>% filter(trait == name)
  
  # Annotate with trait and lineage
  lineage <- map[name]
  t$trait <- name
  t$lineage <- lineage
  t[complete.cases(t),]
}) %>% rbindlist () %>% as.data.frame() -> rare_vars_lineage

# Combine with sum stats 
rare_vars_traits_collapsed <- aggregate(trait ~ var, FUN=paste,data=rare_vars) %>%
  dplyr::select(var,trait) %>% dplyr::rename("all_traits"="trait")
rare_vars_all <- left_join(rare_vars_lineage,rare_vars_traits_collapsed,by="var")
merged <- merge(icd,rare_vars_all[,c("var","trait","effect_allele","blood_beta","all_traits","lineage","PP","VEP_cons_most_serious")],by="var") %>%
  mutate(CHROM=paste0("chr",CHROM))
merged.gr <- makeGRangesFromDataFrame(merged,keep.extra.columns = T,seqnames.field="CHROM",
                                      start.field = "POS",end.field = "POS")
# Add nearest gene
# Load protein coding annotations
grch38.pc <- grch38 %>% filter(biotype == "protein_coding")

# Read in gene body annotations
gencode_gr <- readRDS("../../data/annotations/gencode_filtered.rds")

# Find nearest gene
merged$nearest_gene <- gencode_gr[nearest(merged.gr,gencode_gr,ignore.strand=TRUE),]$gene_name

# PheWAS pvalue threshold
p_thresh <- .05/length(unique(merged$pheno))
strict_p_thresh <- .05/length(unique(merged$pheno))/ length(unique(merged$var))

# Make dot plot grouped by lineage
# Create data frame that only labels the top significant pheno/variant pair in each lineage
toplot <- merged %>% mutate(uniqueID = paste(ID,pheno,lineage,sep="-")) %>% distinct(uniqueID,.keep_all=TRUE)
tolabel <- toplot %>% group_by(lineage,pheno) %>% 
  dplyr::slice(which.max(logp)) %>%
  mutate(tolabel=ifelse(logp > -log10(p_thresh),pheno,"")) %>%
  arrange(lineage,pheno,desc(logp)) %>% as.data.frame() %>% dplyr::select(uniqueID,tolabel) %>%
  right_join(.,toplot,by="uniqueID") %>%
  mutate(tolabel = ifelse(is.na(tolabel),"",tolabel)) %>% arrange(pval) %>% 
  group_by(lineage) %>% mutate(lineage_rank = row_number()) %>% ungroup() %>%
  group_by(lineage,pheno) %>% mutate(mx = max(logp)) %>%    mutate(tolabel = ifelse(is.na(tolabel),"",tolabel)) %>% 
  arrange(lineage,desc(mx)) %>% ungroup() %>%   mutate(order = row_number()) %>% 
  mutate(tolabel = ifelse(lineage_rank <= 5 & logp > -log10(p_thresh),tolabel,"")) %>% as.data.frame() 


p2 <- ggplot(tolabel,aes(x=order,y=logp,label=tolabel))+
  geom_point_rast(aes(color=lineage,fill=lineage),alpha=1,size=0.25) +
  scale_color_manual(values=jdb_palette("lawhoops"))+
  geom_hline(yintercept = -log10(p_thresh), linetype = 2) +
  labs(y="pheWAS -log10(p)",x="Phenotype") +
  geom_text_repel(angle = 0,size=2) +  
  pretty_plot(fontsize=7) + L_border() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position ="none") +
  scale_y_continuous(expand = c(0.05, 0)) +scale_x_continuous(expand = c(0, 0))
# +  facet_grid(.~lineage,scales="free",switch="both")

cowplot::ggsave2(p2,file="../../output/rare_variants/phewas/phewas_dotplot_horizontal.pdf",
                width=4,height=2.5)

# Make table
phewas_table <- merged %>% filter(pval < p_thresh) %>% 
  dplyr::select(var,ID,pheno,trait,effect_allele,blood_beta,PP,VEP_cons_most_serious,num_cases,maf,beta,sebeta,pval,nearest_gene) %>% unique()
write.table(phewas_table,file="../../output/rare_variants/phewas/phewas_sig_phenos.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)

phewas_table <- merged %>% filter(pval < strict_p_thresh) %>% 
  dplyr::select(var,ID,pheno,trait,num_cases,maf,beta,sebeta,pval,nearest_gene)
phewas_table$trait <- vapply(phewas_table$trait, paste, collapse = ", ", character(1L))
write.table(phewas_table,file="../../output/rare_variants/phewas/phewas_sig_phenos_strict.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)

# Write table for VEP
if (FALSE){
  for_vep <- merged %>% filter(pval < p_thresh) %>%dplyr::select(CHROM,POS,REF,ALT) %>% 
    mutate(CHROM = gsub("chr","",CHROM),filler = ".",filler2=".",filler3=".",filler4=".") %>%
    dplyr::select(CHROM,POS,filler,everything())%>% unique()
  
  write.table(for_vep,file="../../output/rare_variants/phewas/phewas_sig_phenos.vcf",
              quote = FALSE, sep = "\t", col.names = F, row.names = F)
}
