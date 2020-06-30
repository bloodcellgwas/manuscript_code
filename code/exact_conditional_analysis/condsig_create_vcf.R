#!/software/R-3.3.0/bin/Rscript
library(data.table)
generate_vcf_from_marker_file<-function(coordid_list, name)
{
	vcf=data.table(
          CHROM=gsub("(.+):.+_.+_.+","\\1", coordid_list),
          POS=as.integer(gsub(".+:(.+)_.+_.+","\\1", coordid_list)),
          ID=coordid_list,
          REF=gsub(".+:.+_(.+)_.+","\\1", coordid_list),
          ALT=gsub(".+:.+_.+_(.+)","\\1", coordid_list),
          QUAL=rep(99,length(coordid_list)),
          FILTER=rep("PASS", length(coordid_list)),
          INFO=rep(".", length(coordid_list)),
          FORMAT=rep(".", length(coordid_list))
        )

	setnames(vcf,"CHROM","#CHROM")
	vcf_header=c("##fileformat=VCFv4.1")
	write.table(vcf, paste(OUT_DIR, "/condout/vep/", name, sep=""), sep="\t",
                    quote=FALSE, row.names=FALSE)        
}

# loop through all the conditionally significant results and create a VCF file based on these
PHENO_NAMES <- Sys.getenv('PHENO_NAMES')
OUT_DIR <- Sys.getenv('OUT_DIR')
traits=as.vector(read.table(PHENO_NAMES, header=F, stringsAsFactors=FALSE))$V1
condsig_list <- c()
for (trait in traits) {
  print(trait)
  if (!file.exists(paste(OUT_DIR, "/condout/results/condsig_", trait, "_gwas_normalised.tsv", sep=""))) {stop("ERROR"); next}
  variant_list <- read.table(paste(OUT_DIR, "/condout/results/condsig_", trait, "_gwas_normalised.tsv", sep=""), header=TRUE, stringsAsFactors=FALSE)$VARIANT
  variant_list = gsub("XY:", "X:", variant_list) # change XY to X
  condsig_list <- c(condsig_list, variant_list)
  
}
condsig_list <- unique(condsig_list)
print(length(condsig_list))
generate_vcf_from_marker_file(condsig_list, "all_condsig.vcf")
