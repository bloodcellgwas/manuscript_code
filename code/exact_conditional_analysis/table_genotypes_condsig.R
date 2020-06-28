#!/software/R-3.1.2/bin/Rscript
library(data.table)
options(scipen=999)
options(warn=2)
PHENO_NAMES=Sys.getenv('PHENO_NAMES')
OUT_DIR=Sys.getenv('OUT_DIR')

pheno_table=read.table(file=PHENO_NAMES,stringsAsFactors=FALSE)[[1]]

variant_table=NULL
for(trait in pheno_table)
  {
    print(trait)
    trait_variants=fread(sprintf("%s/condout/trait/var_list_%s.tsv", OUT_DIR, trait), head=TRUE, sep=",", colClasses="character")
    variant_table=rbind(variant_table, trait_variants)
  }

variant_table_unique=data.table(COORDID=unique(variant_table$VARIANT))

variant_table_unique$CHR=gsub("(.+):(.+)_(.+)_(.+)","\\1",variant_table_unique$COORDID)
variant_table_unique$BP=gsub("(.+):(.+)_(.+)_(.+)","\\2",variant_table_unique$COORDID)
variant_table_unique$REF=gsub("(.+):(.+)_(.+)_(.+)","\\3",variant_table_unique$COORDID)
variant_table_unique$ALT=gsub("(.+):(.+)_(.+)_(.+)","\\4",variant_table_unique$COORDID)
variant_table_unique$CHR_num = variant_table_unique$CHR
variant_table_unique$CHR_num[variant_table_unique$CHR_num %in% c("X", "XY")] = 23
for(chr in 1:23)
  {
    pull_table_int=cbind(
      variant_table_unique$COORDID[variant_table_unique$CHR_num==chr],
      variant_table_unique$COORDID[variant_table_unique$CHR_num==chr],
      variant_table_unique$CHR_num[variant_table_unique$CHR_num==chr],
      variant_table_unique$BP[variant_table_unique$CHR_num==chr],
      variant_table_unique$REF[variant_table_unique$CHR_num==chr],
      variant_table_unique$ALT[variant_table_unique$CHR_num==chr])
    pull_table_int = as.data.frame(pull_table_int, stringsAsFactors=F)
    colnames(pull_table_int)=c("SNPID","rsid","chromosome","position","alleleA","alleleB")
    pull_table_int$sort <- as.numeric(pull_table_int$chromosome)*10^9 + as.numeric(pull_table_int$position)
    pull_table_int = pull_table_int[order(as.numeric(pull_table_int$sort)),]
    pull_table_int = pull_table_int[, !is.element(names(pull_table_int), "sort")]
    write.table(pull_table_int,sprintf("%s/condout/pull_ids/pull_ids_chr%d.tsv", OUT_DIR, chr), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

    ## make position version
    pull_table_int <- data.frame(ID=paste0(sprintf("%02d", as.numeric(pull_table_int$chromosome[pull_table_int$chromosome==chr])), ":", pull_table_int$position[pull_table_int$chromosome==chr]))
    write.table(pull_table_int, paste0(OUT_DIR, "/condout/pull_ids/pull_ids_chr", chr, "_pos.tsv"), row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)

  }

