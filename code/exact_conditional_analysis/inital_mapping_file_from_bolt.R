#!/usr/bin/Rscript
options(warn=2)
library(data.table)
OUT_DIR = Sys.getenv('OUT_DIR')
traits=as.vector(read.table(Sys.getenv('PHENO_NAMES'), header=F, stringsAsFactors=FALSE))$V1

variants = c()

for (trait in traits) {
    assoc = fread(paste0(OUT_DIR, "/assoc_files/", trait, "_gwsig.assoc"), colClasses="character", select="VARIANT")$VARIANT
    variants = unique(c(variants, assoc))
}

generate_vcf_from_marker_file<-function(coordid_list, outputfile) {
    chrom = gsub("(.+):.+_.+_.+","\\1", coordid_list)
    chrom[chrom == "XY"] = "X"
    chrom_num = chrom
    chrom_num[chrom %in% c("X", "XY")] = "23"
    print(unique(chrom))
  vcf=data.table(
    CHROM=chrom,
    POS=as.integer(gsub(".+:(.+)_.+_.+","\\1", coordid_list)),
    ID=coordid_list,
    REF=gsub(".+:.+_(.+)_.+","\\1", coordid_list),
    ALT=gsub(".+:.+_.+_(.+)","\\1", coordid_list),
    QUAL=rep(99,length(coordid_list)),
    FILTER=rep("PASS", length(coordid_list)),
    INFO=rep(".", length(coordid_list)),
    FORMAT=rep(".", length(coordid_list)),
    sort_order=as.numeric(chrom_num)*10^9 +
               as.numeric(gsub("(.+):(.+)_(.+)_(.+)","\\2",coordid_list))
  )
  setkey(vcf, "sort_order")
  vcf[,sort_order:=NULL]
  
  setnames(vcf,"CHROM","#CHROM")
  vcf_header=c("##fileformat=VCFv4.1")
  write.table(vcf_header, file=outputfile, quote=FALSE, col.names=FALSE, row.names=FALSE)
  fwrite(vcf, outputfile, sep="\t", row.names=F, col.names=T, quote=F, append=T)
  return(vcf)
}

vcf.file <- generate_vcf_from_marker_file(variants, paste0(OUT_DIR, "/tables/mapping_file.vcf"))

