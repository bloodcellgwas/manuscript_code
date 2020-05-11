library(data.table)
library(tidyverse)
"%ni%" <- Negate("%in%")

# UK Biobank --------------------------------------------------------------
# Standardize the positions of the ukbb finemap bed files
df_path <- "../data/finemap_bedfiles/ukbb_v2/"
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
in_path <- "../data/finemap_bedfiles/ukbb_v2/"
df_path <- "../data/finemap_bedfiles/ukbb_v2_core_genes/"

# core genes 
bridge <- fread("../data/gene_annotations/bridge_genes.txt")

# variant to gene annotations
CS.df.merged <- fread("../data/VEP/CS.ukid.ukbb_v2.VEP_merged.bed")%>%  mutate(chr = paste0("chr",chr)) 

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

# Non-core variants (OUTERSECT)
df_path <- "../data/finemap_bedfiles/ukbb_v2_noncore_genes/"

ukbb <- lapply(traits,function(trait){
  
  temp <- lapply(1:3,function(i){
    df <- fread(paste0(in_path,trait,".PP001.bed"))
    bed <- df %>% mutate(V3 = V2) %>% filter(V5 > 0.001)
    
    # Filter for just core gene variants of a certain type of blood disorder
    vars <- CS.df.merged %>% filter(gene %in% platform_genes[[i]]) %>% .$rsid %>% unique
    bed_core <- bed %>% filter(V4 %ni% vars)
    
    write.table(bed_core,file=paste0(df_path,names(platform_genes)[i],"_",trait,".PP001.bed"),
                quote = FALSE, sep = "\t", col.names = F, row.names = FALSE)
    
    names(bed_core) <- c("chr","start","end","var","PP","trait")
    bed_core <- bed_core %>% mutate(platform = names(platform_genes)[i])
    return(bed_core)
  }) %>% bind_rows()
  
  return(temp)
}) %>% bind_rows()

# Write any gene bed files
df_path <- "../data/finemap_bedfiles/ukbb_v2_any_genes/"

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

# Normalize by total number of variants in each trait
PP_thresh <- 0.001
ukbb_highPP <- ukbb %>% filter(PP > PP_thresh)
numvar <- CS.df.merged %>% filter(PP>PP_thresh) %>% dplyr::select(rsid, trait) %>% unique %>% group_by(trait) %>% summarise(numvar= n()) %>%
  filter(trait %in% traits)

# Summarize
breakdown <- data.frame(table(ukbb_highPP[,c("trait","platform")]))

normalized <- lapply(traits,function(t){
  denom <- numvar %>% filter(trait == t) %>% .$numvar
  out <- breakdown %>% filter(trait == t) %>% mutate(normalized = 100*Freq / denom)
  return(out)
}) %>% bind_rows()

ggplot(normalized, aes(x = trait, y = normalized)) +
  geom_bar(width = 1, aes(fill = trait), colour="black",
           stat = "identity", position = position_dodge(width=1)) +
  pretty_plot(fontsize = 10) + 
  labs(x = "", y = "Core gene variants (%)", fill = "") +
  scale_fill_manual(values = c(jdb_palette("lawhoops"),jdb_palette("Royal2"))) + 
  facet_wrap(~platform, scales = "free_y") +
  theme(legend.position="bottom") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=7)) + 
  guides(shape = guide_legend(override.aes = list(size = 2)))

  
# Trans-ethnic ------------------------------------------------------------

# Read in transethnic fine-mapped data
df_path <- "../data/finemap_bedfiles/transethnic_v2/trans/"
traits <- list.files(df_path,pattern = "*.txt") %>% gsub("_.*","",.)

trans <- lapply(traits,function(trait){
  df <- fread(paste0(df_path,trait,"_Trans_MRMEGA_20190214_0.95_CredibleSets.txt"))
  bed <- df %>% mutate(chr = paste0("chr",gsub(":.*","",MarkerName)), start = as.numeric(gsub("_.*","",gsub(".*:","",MarkerName))),end=start,trait=trait) %>%
    dplyr::select(chr,start,end,MarkerName,PP,trait)
  write.table(bed,file=paste0(df_path,trait,"_trans_95CS.bed"),
              quote = FALSE, sep = "\t", col.names = F, row.names = FALSE)
  return(bed)
}) %>% bind_rows()

# Transethnic ancestry-specific bed files
ancestries <- c("AA","EA","EAS","HA","SAS")

lapply(ancestries,function(a){
  print(a)
  df_path <- paste0("../data/finemap_bedfiles/transethnic_v2/",a,"/")
  traits <- list.files(df_path,pattern = "*.txt") %>% gsub("_.*","",.)
  remainder <-  list.files(df_path,pattern = "*.txt") %>% gsub(paste0(".*_",a),"",.) 
  
  trans <- lapply(seq(1,length(traits)),function(i){
    print(i)
    df <- fread(paste0(df_path,traits[i],"_",a,remainder[i]))
    if (nrow(df) == 0){return(NULL)}
    
    bed <- df %>% mutate(chr = paste0("chr",gsub(":.*","",MarkerName)), 
                         start = as.numeric(gsub("_.*","",gsub(".*:","",MarkerName))),
                         end=start,trait=traits[i],ancestry=a) %>%
      dplyr::rename("Sentinel_P"="Sentinel_P-   value_association") %>%
      dplyr::select(chr,start,end,MarkerName,PP,trait,ancestry,Sentinel_MarkerName, Sentinel_EAF,Sentinel_P)
    write.table(bed,file=paste0(df_path,traits[i],"_",a,"_95CS.bed"),
                quote = FALSE, sep = "\t", col.names = F, row.names = FALSE)
    return(bed)
  }) %>% bind_rows()
  return(trans)
}) %>% bind_rows() -> all_subancestries

all_trans_snps <- data.frame(union(all_subancestries$MarkerName,trans$MarkerName))
if (TRUE){
  fwrite(all_subancestries,file="../data/finemap_bedfiles/transethnic_v2/all_subancestry_CS.txt",sep="\t")
  fwrite(trans,file="../data/finemap_bedfiles/transethnic_v2/all_trans_CS.txt",sep="\t")
  fwrite(all_trans_snps,file="../data/finemap_bedfiles/transethnic_v2/all_trans_subancestry_SNPs.txt",sep="\t",col.names = FALSE)
}




# Sysmex batch 2 ----------------------------------------------------------

df_path <- "../data/finemap_bedfiles/sysmex/batch2/"
traits <- list.files(df_path,pattern = "*.bed$") %>% gsub("\\..*","",.)

sysmex_2 <- lapply(traits,function(trait){
  df <- fread(paste0(df_path,trait,".pp.bed"))
  bed <- df %>% mutate(V1 = paste0("chr",V1),pheno = trait) %>% filter(V4 > 0.001) %>%
    dplyr::select(V1,V2,V3,V5,V4,pheno)
  write.table(bed,file=paste0(df_path,trait,".PP001.bed"),
              quote = FALSE, sep = "\t", col.names = F, row.names = FALSE)
  return(trait)
})

# Combine into master CS
df_path <- "../data/finemap_bedfiles/sysmex/"
traits <- list.files(df_path,pattern = "*.PP001.bed$") %>% gsub("\\..*","",.)

sysmex_combined <- lapply(traits,function(trait){
  df <- fread(paste0(df_path,trait,".PP001.bed"))
}) %>% bind_rows() %>% setNames(.,c("seqnames","start","end","var","PP","trait"))

fwrite(sysmex_combined,file="../data/finemap_bedfiles/sysmex/sysmex_CS.bed",sep="\t")
