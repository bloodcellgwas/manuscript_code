library(data.table)
library(leafcutter)
library(tidyverse)
library(GenomicRanges)
"%ni%" <- Negate("%in%")

# Read in 1000G dosages for all splice variants
allvariants <- fread("../../data/splice_predictions/1000GP3_files/splice_variants_europeans.dosages.raw") %>% dplyr::select(IID,contains("rs")) 
allvariants_mat <- data.matrix(allvariants[,-1]); rownames(allvariants_mat) <- allvariants$IID
dim(allvariants_mat)

# Check how many variants have genotype data 
splice_variants <- fread("../../output/splice_predictions/splice_variants.ukbb.PP001.delta0.2.tsv")
intersect(gsub("_.*","",colnames(allvariants_mat)),splice_variants$rsid) %>% length()
variants_to_test <- data.frame(var = intersect(gsub("_.*","",colnames(allvariants_mat)),splice_variants$rsid))

splice_vars_gr <- splice_variants %>% filter(rsid %in% variants_to_test[[1]]) %>%
  group_by(rsid) %>% top_n(1, wt = trait) %>%
  mutate(seqnames = paste0("chr", chr)) %>% 
  makeGRangesFromDataFrame(seqnames.field = "seqnames",
                           start.field  = "pos", end.field = "pos", keep.extra.columns = TRUE)
length(splice_vars_gr)

# Import the leaf cutter data
lc_df <- fread("../../data/splice_predictions/all_GV_perind_numers.counts.gz")

# Need this to be a counts matrix for downstream use
lc_mat <- data.matrix(data.frame(lc_df[,2:dim(lc_df)[2]]))
rownames(lc_mat) <- lc_df[[1]]

# Make into a Granges object ot overlap
lc_gr <- data.frame(str_split_fixed(lc_df[["V1"]], ":", 4), stringsAsFactors = FALSE) %>%
  mutate(start = as.numeric(X2), end = as.numeric(X3)) %>%
  makeGRangesFromDataFrame(seqnames.field = "X1",keep.extra.columns = TRUE)
clusters_id <- mcols(lc_gr)[["X4"]]

# Overlap the leafcutter clusters with the splice variants to see when they overlap
ov <- findOverlaps(lc_gr, splice_vars_gr)
length(unique(subjectHits(ov))) # 60 out of the 108 could be tested

pairs_df <- data.frame(
  variant = mcols(splice_vars_gr)$rsid[subjectHits(ov)],
  lc_cluster = mcols(lc_gr)$X4[queryHits(ov)]
) 
pairs_df_unique <- pairs_df[!duplicated(pairs_df),]

# Make sure highlighed variants are still present
pairs_df_unique[pairs_df_unique$variant %in% c("rs139178017","rs141143931"), ]

# With all of that in line, now run leafcutter for the specific clusters 
# Takes a few minutes to test everything in a loop
lapply(1:dim(pairs_df_unique), function(idx){
  print(idx)
  # Extract the variant being tested
  variant <- pairs_df_unique[idx,1] %>% as.character()
  alt_allele <- (splice_variants %>% filter(rsid == variant) %>% pull(alt))[1]
  alt_name <- paste0(variant, "_", alt_allele)
  
  # Partition the donors into carriers or not
  carriers <- rownames(allvariants_mat)[allvariants_mat[,alt_name] > 0]
  has_variant <- (colnames(lc_mat) %in% carriers)*1
  
  # Subset out the cluster for fast testing
  cluster <- pairs_df_unique[idx,2] %>% as.character()
  cluster_mat <- lc_mat[clusters_id == cluster,]
  pvalue <- tryCatch({
    pvalue <- differential_splicing(cluster_mat, has_variant)[[1]]$lrtp
    closeAllConnections()
    pvalue
  }, # necessary afrer some trouble shooting,
  error = function(c) NA
  )
  
  data.frame(variant, cluster, pvalue, n_carriers = sum(has_variant), n_clusters = sum(cluster_mat), idx = idx)
}) %>% rbindlist() %>% data.frame() -> all_df_lc_results

# Examine pvalues
all_df_lc_results$highlight <- all_df_lc_results$variant %in% c("rs139178017","rs141143931")
out_df <- all_df_lc_results %>% arrange(desc(highlight), pvalue)
out_df$adjp <- p.adjust(out_df$pvalue)

sig_leafcutter <- out_df %>%
  filter(adjp < 0.01) %>% pull(variant) %>% unique() %>% as.character()

not_sig_leafcutter <- out_df %>% filter(variant %ni% sig_leafcutter) %>% 
  filter(adjp > 0.01) %>% pull(variant) %>% unique() %>% as.character()

detected_rare <- out_df %>% filter(variant %ni% c(sig_leafcutter,not_sig_leafcutter)) %>% 
  pull(variant) %>% unique() %>% as.character()


109 - length(sig_leafcutter) - length(not_sig_leafcutter) - length(detected_rare)

if (FALSE){
  write.table(out_df, file = "../../output/splice_predictions/leafcutter_analysis_results.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)  
}


# Pick anecdotes ----------------------------------------------------------

out_df <- fread("../output/splice_predictions/leafcutter_analysis_results.tsv")

# Look for good examples
CS <- fread("zcat < ../data/finemap_bedfiles/ukbb_v2/credible_sets_all_with_ranks.txt.gz")
CS.df.ranks <- fread("zcat < ../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.with_ranks.bed.gz")

# Add a column for fine-mapped rank of variant in each credible set
splice_variants_ranks <- left_join(splice_variants,CS.df.ranks[,c("rsid","trait","CS_rank")], by = c( "rsid"="rsid","trait"="trait"))
splice_variants_ranks %>% filter(rsid %in% sig_leafcutter) %>% arrange(AF_Allele2) %>% filter(PP > 0.1)

delta_cutoff <- 0.5
leafcutter_vars <- splice_variants_ranks %>% filter(rsid %in% sig_leafcutter) %>% arrange(AF_Allele2) %>% 
  filter(DS_AG > delta_cutoff | DS_AL > delta_cutoff | DS_DG > delta_cutoff | DS_DL > delta_cutoff) %>% .$rsid %>% unique() 
splice_variants_ranks %>% filter(rsid %in% leafcutter_vars)
CS %>% filter(var %in% leafcutter_vars)

# See what region a specific variant is in
CS.df.ranks %>% filter(rsid == "rs12898397")
CS %>% filter(CS == "rs12908814")

CS.df.ranks %>% filter(rsid == "rs8113779")
CS %>% filter(var == "rs8113779")
CS %>% filter(CS == "rs1005165")
