library(data.table)
OUT_DIR = Sys.getenv('OUT_DIR')
args <- commandArgs(TRUE)
chrom <- args[1]

pull_table = read.table(paste0(OUT_DIR, "/astle_compare/genfiles/ids_", chr, ".tsv"), header=F, stringsAsFactors=F)


