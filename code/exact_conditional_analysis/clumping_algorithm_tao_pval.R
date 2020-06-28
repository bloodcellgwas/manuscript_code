#!/software/R-3.3.0/bin/Rscript
library(data.table)
OUT_DIR <- Sys.getenv('OUT_DIR')
PHENO_NAMES <- Sys.getenv("PHENO_NAMES")
traits <- as.vector(read.table(PHENO_NAMES, header=F, stringsAsFactors=F)$V1)

final_output = fread(paste0(OUT_DIR, "/BCX_final_output_raw.csv"))
clump_table = fread(paste(OUT_DIR, "/condout/ldclump_dosage/hg19_clump_table_tao.tsv", sep=""))
stop("a")
