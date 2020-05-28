library(data.table)
library(tidyverse)
"%ni%" <- Negate("%in%")

# UK Biobank --------------------------------------------------------------
# Standardize the positions of the ukbb finemap bed files
df_path <- "../../data/finemap_bedfiles/ukbb_v2/"
traits <- list.files(df_path,pattern = "*.bed$") %>% gsub("\\..*","",.)

ukbb <- lapply(traits,function(trait){
  df <- fread(paste0(df_path,trait,".PP001.bed"))
  bed <- df %>% mutate(V3 = V2) %>% filter(V5 > 0.001)
  write.table(bed,file=paste0(df_path,trait,".PP001.bed"),
              quote = FALSE, sep = "\t", col.names = F, row.names = FALSE)
  return(trait)
})

summary(lm(HSC_phenotype ~ rs1 + rs2 + rs1*rs2, data = ped.file))


# UK Biobank core gene variants -------------------------------------------
in_path <- "../../data/finemap_bedfiles/ukbb_v2/"
df_path <- "../../data/finemap_bedfiles/ukbb_v2_core_genes/"

# core genes 
bridge <- fread("../../data/gene_annotations/bridge_genes.txt")

# variant to gene annotations
CS.df.merged <- fread("../../data/VEP/CS.ukid.ukbb_v2.VEP_merged.bed")%>%  mutate(chr = paste0("chr",chr)) 

core_gene_platforms <- c("BMF","BPD","SMD")
platform_genes <- lapply(core_gene_platforms,function(i){
  print(i)
  bridge %>% filter(Platform == i) %>% .$Gene_Symbol_HGNC
})
names(platform_genes) <- core_gene_platforms

traits <- list.files(in_path,pattern = "*.bed$") %>% gsub("\\..*","",.)

# Write core gene bed files
ukbb <- lapply(traits,function(trait){
  
  temp <- lapply(1:3,function(i){
    df <- fread(paste0(in_path,trait,".PP001.bed"))
    bed <- df %>% mutate(V3 = V2) %>% filter(V5 > 0.001)
    
    # Filter for just core gene variants of a certain type of blood disorder
    vars <- CS.df.merged %>% filter(gene %in% platform_genes[[i]]) %>% .$rsid %>% unique
    bed_core <- bed %>% filter(V4 %in% vars)
    
    write.table(bed_core,file=paste0(df_path,names(platform_genes)[i],"_",trait,".PP001.bed"),
                quote = FALSE, sep = "\t", col.names = F, row.names = FALSE)
    
    names(bed_core) <- c("chr","start","end","var","PP","trait")
    bed_core <- bed_core %>% mutate(platform = names(platform_genes)[i])
    return(bed_core)
  }) %>% bind_rows()

  return(temp)
}) %>% bind_rows()

# Write any gene bed files
df_path <- "../../data/finemap_bedfiles/ukbb_v2_gene_sets/"

ukbb <- lapply(traits,function(trait){
  
    df <- fread(paste0(in_path,trait,".PP001.bed"))
    bed <- df %>% mutate(V3 = V2) %>% filter(V5 > 0.001)
    
    # Filter for any gene variants
    vars <- CS.df.merged %>% filter(gene != "") %>% .$rsid %>% unique
    bed_core <- bed %>% filter(V4 %in% vars)
    
    write.table(bed_core,file=paste0(df_path,"any_genes_",trait,".PP001.bed"),
                quote = FALSE, sep = "\t", col.names = F, row.names = FALSE)
    
    names(bed_core) <- c("chr","start","end","var","PP","trait")
    bed_core <- bed_core %>% mutate(platform = "any_genes")

  return(bed_core)
}) %>% bind_rows()


