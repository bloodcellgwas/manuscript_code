#!/software/R-3.3.0/bin/Rscript
library(data.table)
OUT_DIR = Sys.getenv('OUT_DIR')
args <- commandArgs(TRUE)
chrom <- args[1]

pull_table = fread(paste0(OUT_DIR, "/genfiles/gwsig_snps_chr", chrom, "_pos.tsv"), header=F, colClasses="character")
pull_table$V1 = gsub("^0*", "", pull_table$V1)
ids = fread(paste0(OUT_DIR, "/genfiles/ids_", chrom, ".tsv"), header=F, colClasses="character")
if (chrom %in% c("XY", "X")) {
    ids$hg19 = paste0(gsub("(.+):.+_.+_.+","\\1", ids$V1), ":", gsub(".+:(.+)_.+_.+","\\1", ids$V1))
} else {
    ids$hg19 = paste0(as.numeric(gsub("(.+):.+_.+_.+","\\1", ids$V1)), ":", gsub(".+:(.+)_.+_.+","\\1", ids$V1))
}

if (sum(pull_table$V1 %in% ids$hg19) != nrow(pull_table)) {
  stop(paste0("not all pull id variants have been extracted for chromosome ", chrom))
}
