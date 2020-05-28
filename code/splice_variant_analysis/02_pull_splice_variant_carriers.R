library(tidyverse)
library(data.table)
"%ni%" <- Negate("%in%")

# Read in 1000G dosages for all splice variants
allvariants <- fread("../../data/splice_predictions/1000GP3_files/splice_variants_europeans.dosages.raw") %>% dplyr::select(IID,contains("rs")) 
allvariants_mat <- data.matrix(allvariants[,-1]); rownames(allvariants_mat) <- allvariants$IID
dim(allvariants_mat)
allvariants_mat[1:5,1:5]

if (FALSE){
  # Take out UBASH3A common variant
  splice_dosages <- fread("../../data/splice_predictions/1000GP3_files/splice_variants_europeans.dosages.raw") %>% dplyr::select(IID,contains("rs")) %>% dplyr::select(-c("rs1893592_C"))
  splice_dosage_mat <- data.matrix(splice_dosages[,-1]); rownames(splice_dosage_mat) <- splice_dosages$IID
}

# Check how many variants have genotype data 
splice_variants <- fread("../../output/splice_predictions/splice_variants.ukbb.PP001.delta0.2.tsv")
intersect(gsub("_.*","",colnames(allvariants_mat)),splice_variants$rsid) %>% length()
variants_to_test <- data.frame(var = intersect(gsub("_.*","",colnames(allvariants_mat)),splice_variants$rsid))


setdiff(gsub("_.*","",colnames(allvariants_mat)),splice_variants$rsid)  # One variant only from credible sets, not PP > 0.001

# Check how many variants have at least 1 carrier
table(colSums(allvariants_mat) > 0)

# Keep individuals who are a carrier for one or more variants
keep <- rowSums(allvariants_mat) > 0
table(keep)
allvariants_mat <- allvariants_mat[keep,]

# Pull fastq addresses and subset on only donors with RNA-seq data
fq <- fread("../../data/splice_predictions/1000GP3_files/file_addresses.txt")
fq %>% filter(Donor %in% rownames(allvariants_mat)) %>% dim()

splice_carriers <- allvariants_mat[rownames(allvariants_mat) %in% fq$Donor,]

# Write table for each variant
variants <- colnames(splice_carriers)
sapply(colnames(splice_carriers),function(y){
  subset <- splice_carriers[,y]
  subset <- subset[subset > 0] %>% as.data.frame() %>% setNames(., "allele_count")
  subset$ID <- rownames(subset) 
  
  fwrite(subset,paste0("../../data/splice_predictions/carriers_with_rnaseq/",y,"_carriers_with_rnaseq.txt"),sep="\t")
  return(nrow(subset))
})

# Pull fastq address
fq <- fread("../../data/splice_predictions/1000GP3_files/file_addresses.txt")
carrier_fq <- fq %>% filter(Donor %in% rownames(splice_carriers))
dim(carrier_fq)

fwrite(carrier_fq,"../../data/splice_predictions/1000GP3_files/carrier_fastq_addresses.txt",sep="\t")
