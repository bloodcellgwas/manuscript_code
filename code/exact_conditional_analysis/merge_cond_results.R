#!/usr/bin/Rscript
library(RcppEigen)
library(data.table)
library(rhdf5)
library(doMC)
options(show.error.locations = TRUE)
OUT_DIR <- Sys.getenv('OUT_DIR')
args = commandArgs(TRUE)
pheno = args[1]

merged = NULL
for (chr in 1:23) {
    this_chr = read.table(paste0(OUT_DIR, "/condout/results_chrom/condsig_", pheno, "_gwas_normalised_chr_", chr, ".tsv"), header=T, sep="\t")
    merged = rbind(merged, this_chr)
}

write.table(merged, paste0(OUT_DIR, "/condout/results/condsig_", pheno, "_gwas_normalised.tsv"), quote=F, sep="\t", row.names=F)
