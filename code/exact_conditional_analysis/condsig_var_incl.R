#!/software/R-3.2.2/bin/Rscript
OUT_DIR <- Sys.getenv('OUT_DIR')
MAPPING_FILE <- Sys.getenv('MAPPING_FILE')
traits <- as.vector(read.table(Sys.getenv("PHENO_NAMES"), header=F, stringsAsFactors=F)$V1)
options(scipen=999)
library(data.table)

#mapping <- fread(MAPPING_FILE)

get_condsig_vars <- function() {
  condsig_vars <- c()
  for (trait in traits) {
    condsig_var <- read.table(paste(OUT_DIR, '/condout/results/condsig_', trait, '_gwas_normalised.tsv', sep=""), header=TRUE, stringsAsFactors=FALSE)
    condsig_vars <- c(condsig_vars, condsig_var$VARIANT)
  }
  condsig_vars <- unique(condsig_vars)
  return(condsig_vars)
}

position_pull <- function(condsig_vars) {
  chr <- gsub("(.+):.+_.+_.+","\\1", condsig_vars)
  bp <- gsub(".+:(.+)_.+_.+","\\1", condsig_vars)
  chr_num = chr
  chr_num[chr_num %in% c("X", "XY")] = "23"
  sort <- as.numeric(chr_num)*10^9 + as.numeric(bp)
  pull_table = data.frame(chr=chr, chr_num=chr_num, bp=bp, sort=sort, stringsAsFactors=F)
  pull_table = pull_table[order(pull_table$sort),]
  pull_table$chr[!is.element(pull_table$chr, c("X", "XY"))] = sprintf("%02d", as.numeric(pull_table$chr[!is.element(pull_table$chr, c("X", "XY"))]))
  pull_table = pull_table[, names(pull_table) != "sort"]
  pull_table$pos = paste0(pull_table$chr, ":", pull_table$bp)
  write.table(pull_table[, c("pos"), drop=FALSE], paste0(OUT_DIR, "/condout/ldclump/pull_table.tsv"), row.names=F, quote=F, sep="\t", col.names=F)

}

pull_from_genfiles <- function(condsig_vars) { 
  ref <- gsub(".+:.+_(.+)_.+","\\1", condsig_vars)
  alt <- gsub(".+:.+_.+_(.+)","\\1", condsig_vars)
  chr <- gsub("(.+):.+_.+_.+","\\1", condsig_vars)
  bp <- gsub(".+:(.+)_.+_.+","\\1", condsig_vars)
  pull_table <- data.frame(mapping[match(condsig_vars, mapping$COORDID), c(2,3), with=FALSE],
                            chromosome=chr, position=bp, alleleA=ref, alleleB=alt)

  names(pull_table) <- c("SNPID", "rsid", "chromosome", "position", "alleleA", "alleleB")

  ### SORT THE SNPS BY LOCATION IN GENOME ###
  pull_table$sort=as.numeric(pull_table$position)+10^9*as.numeric(pull_table$chromosome)
  pull_table <- pull_table[order(pull_table$sort),]
  pull_table <- pull_table[, !(names(pull_table) == "sort")]

  write.table(pull_table, paste(OUT_DIR, "/condout/ldclump/pull_table.tsv", sep=""), row.names=F, quote=F, sep="\t")
}

condsig_vars <- get_condsig_vars()
ids = data.frame(VARIANT=condsig_vars, stringsAsFactors=F)
write.table(ids, paste0(OUT_DIR, "/condout/ldclump/pull_table_ids.tsv"), row.names=F, col.names=F, quote=F, sep="\t")
stop("aa")
##pull_from_genfiles(condsig_vars)
position_pull(condsig_vars)
