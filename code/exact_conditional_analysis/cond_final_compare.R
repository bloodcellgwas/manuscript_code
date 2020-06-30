#!/usr/bin/Rscript
library(RcppEigen)
library(data.table)
library(rhdf5)
library(doMC)
options(show.error.locations = TRUE)
OUT_DIR <- Sys.getenv('OUT_DIR')
traits = read.table(Sys.getenv('PHENO_NAMES'), header=F, stringsAsFactors=F)$V1

signals_before = c()
signals_after =  c()

for (trait in traits) {

    recent_condsig = read.table(paste0(OUT_DIR, "/condout/results/condsig_", trait, "_gwas_normalised.tsv"), sep="\t", header=T, stringsAsFactors=F)
    previous_condsig = read.table(paste0("/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_3/condout/results/condsig_", trait, "_gwas_normalised.tsv"), sep="\t", header=T, stringsAsFactors=F)

    print("Total number of variants")
    print(paste("current analysis: ", nrow(recent_condsig)))
    print(paste("previous analysis: ", nrow(previous_condsig)))
    signals_before = c(signals_before, nrow(previous_condsig))
    signals_after = c(signals_after, nrow(recent_condsig))

    
    print("Variants present in this study but not previous")
    print(recent_condsig$VARIANT[!is.element(recent_condsig$VARIANT, previous_condsig$VARIANT)])

    print("Variants present in previous study but not current")
    print(previous_condsig$VARIANT[!is.element(previous_condsig$VARIANT, recent_condsig$VARIANT)])

}

result = data.frame(trait=traits, condsig_variants_before=signals_before, condsig_variants_after=signals_after)
write.csv(result, paste0(OUT_DIR, "/ukbb_compare_previous_number_signals.csv"), row.names=F, quote=F)
