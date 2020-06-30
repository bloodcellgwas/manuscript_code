#!Rscript
OUT_DIR <- Sys.getenv('OUT_DIR')
args <- commandArgs(TRUE)
chr <- args[1]

pulltable = read.table(paste0(OUT_DIR, "/condout/pull_ids/pull_ids_chr", chr, ".tsv"), stringsAsFactors=F, header=T)
ids = read.table(paste0(OUT_DIR, "/condout/subgenhd5/chr_ids", chr, ".tsv"), stringsAsFactors=F, header=F)
if (chr == 23) {
    ## order gets messed up in chromosome 23 because its actually two chromosomes but as long as all the snps we want are there it's fine but more stringent for the other chromosomes because the order should never be messed up so if it is something has happened
    if (sum(pulltable$SNPID %in% ids$V1) == nrow(pulltable) & nrow(pulltable) == nrow(ids)) {
        print("done")
    } else {stop("doesnt match for chr 23")}
} else if (sum(ids$V1 == pulltable$SNPID) != nrow(ids) | nrow(ids) != nrow(pulltable)) {
  stop(paste0("The rows don't match for chromosome: ", chr))
}  
