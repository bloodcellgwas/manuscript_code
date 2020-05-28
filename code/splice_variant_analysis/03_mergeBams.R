library(tidyverse)
library(data.table)
"%ni%" <- Negate("%in%")

# Read in all splice variants
splice_variants <- fread("../../output/splice_predictions/splice_variants.ukbb.PP001.delta0.2.tsv")

# Pick variants to merge bam files across all carriers and non-carriers (separately)
variants_to_merge <- c("rs8113779","rs12898397","rs139178017")

ids <- splice_variants %>% distinct(rsid,.keep_all = TRUE) %>% 
  dplyr::filter(rsid %in% variants_to_merge) %>% mutate(id = paste(rsid,alt,sep="_")) %>% .$id

all_individuals <- fread("../../data/splice_predictions/1000GP3_files/file_addresses.txt")$Donor

cluster_path <- "path to store bam files (too large for Github)"
# Merge carrier BAMs one variant at a time
for (variant in ids){
  name <- gsub("_.*","",variant)
  print(paste0("Starting to merge bams for ",name))
  
  # Pull out the region
  window <- 5000
  one_df <- splice_variants %>% distinct(rsid,.keep_all = TRUE) %>%
    dplyr::filter(rsid == name) %>% mutate(startbam = pos - window, endbam = pos + window)
  region <- paste0(one_df[1,"chr"], ":", one_df[1,"startbam"], "-",  one_df[1,"endbam"] )
  
  # Merge carrier BAMs
  df <- read.table(paste0("../../data/splice_predictions/carriers_with_rnaseq/",variant,"_carriers_with_rnaseq.txt"), header = TRUE)
  ids_to_merge <- df %>% pull(ID) 
  bams_to_merge <- paste0(cluster_path, ids_to_merge, ".bam")
  bams_to_merge2 <- paste0(cluster_path, ids_to_merge, ".", name,".bam")
  bams_to_merge_collapsed <- paste(bams_to_merge2, collapse = " ")
  
  new_bam <- paste0( cluster_path,name, ".merge.case.bam")
  
  # Make subset bams
  lapply(1:length(bams_to_merge), function(i){
    cmd_subset <- paste0("samtools view ", bams_to_merge[i], " ", region, " -b > ",  bams_to_merge2[i])
    system(cmd_subset)
    system(paste0("samtools index ", bams_to_merge2[i]))
  })
  
  samtools_call <- paste0("samtools merge -f ",  new_bam," ", bams_to_merge_collapsed)
  system(samtools_call)
  system(paste0("samtools index ", new_bam))
  
  lapply(bams_to_merge2, function(rmbam){
    system(paste0("rm ", rmbam))
    system(paste0("rm ", rmbam, ".bai"))
  })
  
  #-----------------------------------------------------------
  
  # Make another merged bam for all non-carriers
  control_ids <- all_individuals[all_individuals %ni% ids_to_merge]
  controls_to_merge <- paste0(cluster_path, control_ids, ".bam")
  controls_to_merge2 <- paste0(cluster_path, control_ids,".", name, ".bam")
  controls_to_merge_collapsed <- paste0(controls_to_merge2,collapse = " ")
  
  # Make subset bams
  lapply(1:length(controls_to_merge), function(i){
    cmd_subset <- paste0("samtools view ", controls_to_merge[i], " ", region, " -b > ",  controls_to_merge2[i])
    system(cmd_subset)
    system(paste0("samtools index ", controls_to_merge2[i]))
  })
  
  new_bam <- paste0(cluster_path, name, ".merge.control.bam")
  samtools_call <- paste0("samtools merge -f ",  new_bam," ", controls_to_merge_collapsed)
  system(samtools_call)
  system(paste0("samtools index ", new_bam))
  
  lapply(controls_to_merge2, function(rmbam){
    system(paste0("rm ", rmbam))
    system(paste0("rm ", rmbam, ".bai"))
  })
  
  print(paste0("Finished merging case and control bams for ",name))
}

