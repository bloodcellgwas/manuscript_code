#!/software/R-3.1.2/bin/Rscript
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
	write.table(vcf, paste(OUT_DIR, "/vcf/", name, sep=""), sep="\t", quote=FALSE, row.names=FALSE)        
}

# loop through all the conditionally significant results and create a VCF file based on these
PHENO_NAMES <- Sys.getenv('PHENO_NAMES')
OUT_DIR <- Sys.getenv('OUT_DIR')
args = commandArgs(TRUE)
trait = args[1]
print(trait)
cond_sig_table <- read.csv(paste(OUT_DIR, "/assoc_files/", trait, "_gwsig.assoc", sep=""), header=TRUE, stringsAsFactors=FALSE)
generate_vcf_from_marker_file(cond_sig_table$VARIANT, paste(trait, "_gwsig.vcf", sep=""))

