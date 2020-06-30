#!/software/R-3.1.2/bin/Rscript
options(warn=2)
library(data.table)
OUT_DIR=Sys.getenv("OUT_DIR")
PHENO_NAMES=Sys.getenv("PHENO_NAMES")

block_table=fread(sprintf("%s/blocks/blocks.tsv", OUT_DIR))

variant_table=NULL
for(block in 1:max(block_table$BLOCK))
{
    block_variants=fread(sprintf("%s/condout/blocks/ind_var_block%d.tsv", OUT_DIR, block), head=TRUE, sep="\t", colClasses=c("character"))
    block_variants$VARIANT <- gsub("^0:*", "", block_variants$VARIANT)
    
    pull_id <- read.table(paste0(OUT_DIR, "/blocks/pull_ids_block_", block, ".tsv"), header=T, colClasses="character")
    pull_id$VARIANT <- paste0(as.character(pull_id$chromosome), ":", as.character(pull_id$position), "_", pull_id$alleleA, "_", pull_id$alleleB)
    if (sum(is.element(block_variants$VARIANT, pull_id$VARIANT)) != nrow(block_variants)) {
        print(head(block_variants$VARIANT))
        print(head(pull_id$VARIANT))
        stop(paste(block, "contains variants which aren't in the pull file"))
    }
    
    variant_table=rbind(variant_table,cbind(block_variants, rep(unique(block_table$TRAIT[block_table$BLOCK==block]), dim(block_variants)[1])))
}
colnames(variant_table)[2]="TRAIT"
traits=unique(variant_table$TRAIT)
pheno_table=read.table(file=PHENO_NAMES,stringsAsFactors=FALSE)[[1]]
for(trait in traits)
{
    print(trait)
    trait_list=data.table(VARIANT=variant_table$VARIANT[variant_table$TRAIT==trait])
    chr = gsub("(.+):(.+)_(.+)_(.+)","\\1",trait_list$VARIANT)
    chr[chr %in% c("X", "XY")] = 23
    trait_list$sort=10^9*as.numeric(chr)+as.numeric(gsub("(.+):(.+)_(.+)_(.+)","\\2",trait_list$VARIANT))
    setkey(trait_list,"sort")
    trait_list[,sort:=NULL]
    write.table(trait_list,sprintf("%s/condout/trait/var_list_%s.tsv", OUT_DIR, trait), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
}
