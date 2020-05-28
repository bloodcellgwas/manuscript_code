library(tidyverse)
library(GenomicRanges)
library(data.table)
library(BuenColors)

"%ni%" <- Negate("%in%")

# Add region_rank to cred_set
if (FALSE){
  cred_set <- fread("zcat < ../../data/finemap_bedfiles/ukbb_v2/credible_sets_all.txt.gz")
  cred_set <- cred_set %>% group_by(trait,CS) %>% arrange(desc(PP)) %>% mutate(CS_rank = row_number()) %>% ungroup() %>% as.data.frame() 
  fwrite(cred_set, file = "../../data/finemap_bedfiles/ukbb_v2/credible_sets_all_with_ranks.txt.gz",compress = "auto")
  
  CS.df <- fread("zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz") %>% filter(PP>0.001)
  CS.df.ranks <- left_join(CS.df,cred_set[,c("var","trait","CS_rank")], by=c("rsid"="var","trait"="trait"))
  fwrite(CS.df.ranks, file = "../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.with_ranks.bed.gz",compress = "auto")
}

# Read in variants
CS.df <- fread("zcat < ../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz") %>% filter(PP>0.001)
cred_set <- fread("zcat < ../../data/finemap_bedfiles/ukbb_v2/credible_sets_all_with_ranks.txt.gz")
ukid_conversions <- fread("../../data/finemap_bedfiles/cred_set_rsids_ukids_ALL.txt") %>% group_by(rsid) %>%
  filter(n() == 1) %>% ungroup()
cred_set_ukids <- merge(cred_set,ukid_conversions,by.x="var",by.y= "rsid") %>% dplyr::rename(rsid = "var")

ukid_conversions_combined <- bind_rows(CS.df[,c("rsid","UKID")],ukid_conversions) %>% unique()

format_splice_stats <- function(splice){
  colnames(splice) <- as.character(splice[1,]) %>% gsub("=.*","",.)
  colnames(splice)[1:5] <- c("var","chr","pos","ref","alt")
  splice_reformatted <- cbind(splice[,1:5],apply(splice[,6:ncol(splice)], 2, function(y) gsub(".*=","",y)))
  return(splice_reformatted)
}
# Subset just on delta score, without incorporating PP
make_splice_bed <- function(splice,delta=0.2){
  splice_variants <- splice %>% filter(DS_AG > cutoff | DS_AL > cutoff | DS_DG > cutoff | DS_DL > cutoff) 
  splice_variants %>% mutate(end = pos) %>% dplyr::select(chr,pos,end,var) %>% arrange(chr,pos) %>% unique()
}

ukbb_splice <- fread("../../data/splice_predictions/BCX_PP001.out.vcf")%>% format_splice_stats()
ukbb_cred_splice <- fread("../../data/splice_predictions/BCX_ukbb_cred_set_vars.out.vcf")%>% format_splice_stats()
all_splice <- bind_rows(ukbb_splice,ukbb_cred_splice) %>% unique() %>% arrange(chr,pos)

# Write european splice variants to pull from 1000GP3
if (FALSE){
  cutoff <- 0.2
  splice_bed <- all_splice %>% filter(DS_AG > cutoff | DS_AL > cutoff | DS_DG > cutoff | DS_DL > cutoff) %>%unique() 
  merged <- merge(splice_bed,ukid_conversions_combined[,c("UKID","rsid")],by.x="var",by.y="UKID") %>% dplyr::select(var,rsid,everything())
  
  fwrite(merged,file="../../output/splice_predictions/splice_variants.ukbb.CS.delta0.2.tsv",sep="\t")
}

# Merge with UKBB fine-mapping stats
all_splice$var %>% unique %>% length()
merged <- merge(all_splice,CS.df[,c("UKID","rsid","PP","AF_Allele2","trait")],by.x="var",by.y="UKID") 
# Merge with CS 
# merged <- merge(all_splice,ukid_conversions_combined[,c("UKID","rsid")],by.x="var",by.y="UKID")

# Delta score of a variant ranges from 0 to 1, and can be interpreted as the probability of the variant being splice-altering. In the paper, a detailed characterization is provided for 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs. Delta position conveys information about the location where splicing changes relative to the variant position (positive values are upstream of the variant, negative values are downstream).
cutoff <- 0.2
PP_threshold <- 0.5
splice_variants <- merged %>% filter(DS_AG > cutoff | DS_AL > cutoff | DS_DG > cutoff | DS_DL > cutoff) %>% arrange(desc(PP)) %>% dplyr::select(var,rsid,everything())
splice_variants$var %>% unique %>% length()
splice_variants$SYMBOL %>% unique %>% length()

# Annotate frame preserving vs. frame altering
# Acceptors: distance between canonical splice acceptor and newly created acceptor, multiple of 3?
# Donors: distance between canonical splice donor and newly created donor, multiple of 3?
frame <- splice_variants %>% group_by(var) %>% dplyr::slice(which.max(PP)) %>%
  pivot_longer(data=.,cols=starts_with("DS"),names_to = "type",values_to = "delta") %>%
  group_by(var) %>% dplyr::slice(which.max(delta)) %>% ungroup() %>% mutate(DIST = as.integer(DIST),frame = ifelse(DIST %% 3 == 0,"in-frame","frameshift"))
table(frame$frame)
frame %>% group_by(frame) %>% summarise(median(PP))

splice_variants <- splice_variants %>% left_join(.,frame[,c("var","frame")])
splice_variants %>% filter(frame == "in-frame")

highPP_splice_variants <- splice_variants %>% filter(PP> PP_threshold) 
highPP_splice_variants$var %>% unique %>% length()

if (FALSE){
  fwrite(splice_variants,file=paste0("../../output/splice_predictions/splice_variants.ukbb.PP001.delta",cutoff,".tsv"),sep="\t")
  fwrite(highPP_splice_variants,file=paste0("../../output/splice_predictions/splice_variants.ukbb.PP50.delta",cutoff,".tsv"),sep="\t")
}

# How many are non-canonical splice variants?
splice_variants[,c("DIST","DP_AG","DP_AL","DP_DG","DP_DL")] <- data.matrix(splice_variants[,c("DIST","DP_AG","DP_AL","DP_DG","DP_DL")])
splice_variants %>% filter(abs(DIST) > 2) %>% .$var %>% unique() %>% length()

splice_variants %>% distinct(var,.keep_all = TRUE) %>% .$DIST %>% summary()
splice_variants %>% distinct(var,.keep_all = TRUE) %>% .$TYPE %>% table()

# How many splice alterations are frame-preserving vs. not?
distances <- splice_variants %>% distinct(var,.keep_all = TRUE) %>% .$DIST 
distances %% 3 %>% table()

# Write table of all splice gene targets
if (FALSE){
  write.table(splice_variants$SYMBOL %>% unique(),file="../output/splice_predictions/splice_genes_ukbb.delta0.2.PP001.tsv",quote = FALSE, sep = "\t", col.names = F, row.names = FALSE) 
}

# Check for specific kinds of splice variants
splice_variants %>% filter(DS_AL > 0.8) 
splice_variants %>% filter(DS_AG > 0.5)
splice_variants %>% filter(DS_DG > 0.2,PP > 0.01) %>% distinct(var,.keep_all = TRUE)
splice_variants %>% filter(DS_DL > 0.5)

cred_set %>% filter(CS=="rs7180484")
cred_set %>% filter(CS  == "rs112463197")

# KS test to determine if splice variants have higher PP and lower MAF
CS.df.PPmax <- CS.df %>% group_by(UKID) %>% dplyr::slice(which.max(PP)) %>% 
  ungroup() %>% arrange(PP) %>% mutate(rank = row_number(),type="all finemap",
                                       MAF = ifelse(AF_Allele2 > 0.5, 1-AF_Allele2,AF_Allele2)) %>% dplyr::select(rank,PP,MAF,type)

splice_variants.PPmax <- splice_variants %>% group_by(var) %>% dplyr::slice(which.max(PP)) %>% 
  ungroup() %>% arrange(PP) %>% mutate(rank = row_number(), type="splice variants",
                                       MAF = ifelse(AF_Allele2 > 0.5, 1-AF_Allele2,AF_Allele2))%>% dplyr::select(rank,PP,MAF,type)

ks.test(splice_variants.PPmax$PP,CS.df.PPmax$PP,alternative="two.sided")
ks.test(splice_variants.PPmax$MAF,CS.df.PPmax$MAF,alternative="two.sided")
wilcox.test(splice_variants.PPmax$MAF,CS.df.PPmax$MAF,alternative="two.sided")
wilcox.test(splice_variants.PPmax$PP,CS.df.PPmax$PP,alternative="two.sided")

max_merged <- bind_rows(CS.df.PPmax,splice_variants.PPmax)

# Compare PP 
p1 <- ggplot(max_merged, aes(x=type,y=PP,fill=type)) +
  geom_violin()+
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,4)]) +
  labs(x="") + 
  pretty_plot(fontsize=8) + L_border() +
  theme(legend.position="none")
p1

# Compare MAF
p2 <- ggplot(max_merged, aes(x=type,y=MAF,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,4)]) +
  labs(x="") + 
  pretty_plot(fontsize=8) + L_border()+
  theme(legend.position="none")
p2

# Density plots
ggplot(max_merged, aes(x=PP)) +
  geom_density(aes(fill=type),alpha=0.4) +
  scale_fill_manual(values=jdb_palette("brewer_spectra")[c(1,4)]) +
  pretty_plot(fontsize = 8)+ L_border() +
  labs(x="PP",y="Density") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
p2 <- ggplot(max_merged, aes(x=MAF)) +
  geom_density(aes(fill=type),alpha=0.4) +
  scale_fill_manual(values=jdb_palette("brewer_spectra")[c(1,4)]) +
  pretty_plot(fontsize = 8)+ L_border() +
  labs(x="MAF",y="Density") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
p2

if (TRUE){
  cowplot::ggsave2(p1, file="../../output/splice_predictions/splice_variants_PP_distribution.pdf",width=2.5,height=2.5)
  cowplot::ggsave2(p2, file="../../output/splice_predictions/splice_variants_MAF_distribution.pdf",width=3,height=2.5)
}
