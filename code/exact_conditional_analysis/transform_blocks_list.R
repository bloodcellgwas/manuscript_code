#!/software/R-3.3.0/bin/Rscript
library(data.table)
OUT_DIR <- Sys.getenv('OUT_DIR')


block_chrs <- read.table(paste0(OUT_DIR, "/blocks/block_chrs.tsv"), header=T)
new_order <- NULL
blocks_each_chrom <- NULL
for (chr in 1:23) {
  block_chrs_s <- block_chrs[block_chrs$CHR_num == chr,]
  new_order <- rbind(new_order, block_chrs_s)
  blocks_each_chrom <- rbind(blocks_each_chrom, c(chr, nrow(block_chrs_s)))
}
blocks_each_chrom <- as.data.frame(blocks_each_chrom)
blocks_each_chrom$blocks_cumulative <- cumsum(blocks_each_chrom$V2)
names(blocks_each_chrom) <- c("chr", "blocks", "blocks_cumulative")

write.table(new_order, paste0(OUT_DIR, "/blocks/block_chrs_order.tsv"), row.names=F, quote=F)
write.table(blocks_each_chrom, paste0(OUT_DIR, "/blocks/blocks_each_chrom.tsv"), row.names=F, quote=F)
