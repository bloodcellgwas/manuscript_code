library(seqinr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gkmSVM)
library(dplyr)

snps_df <- readr::read_tsv("../../data/finemap_bedfiles/ukbb_v2/CS.ukid.ukbb_v2.PP0.001.bed.gz", col_names = TRUE) %>% data.frame()
snps_df$chr <- paste0("chr", as.character(snps_df$chr))
snps_df$BP <- snps_df$start
pad <- 9

# Make a GenomicRanges object of the regions of interest
snps_df$start <- snps_df$BP - pad
snps_df$end <- snps_df$BP + pad
snps_gr <- makeGRangesFromDataFrame(snps_df, keep.extra.columns = TRUE)

# Extract sequences from a reference genome
sequence <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, snps_gr))

# Annotate with the reference and the alternate alleles
ref <-  paste0(substring(sequence, 1,9), snps_df$Allele1, substring(sequence, 11,19))
alt <-  paste0(substring(sequence, 1,9), snps_df$Allele2, substring(sequence, 11,19))
names <- paste0( snps_df$chr, ":", snps_df$BP, "_", snps_df$Allele1, "-", snps_df$Allele2, "-", snps_df$trait, "-", as.character(round(snps_df$PP, 4)))

boo <- duplicated(names)

write.fasta(as.list(ref[!boo]), names[!boo], "../../output/fastas_svm/FMall_v2_ref.fasta",
            open = "w", nbchar = 60)

write.fasta(as.list(alt[!boo]), names[!boo], "../../output/fastas_svm/FMall_v2_alt.fasta",
            open = "w", nbchar = 60, as.string = FALSE)



