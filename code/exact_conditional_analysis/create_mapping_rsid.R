#!/software/R-3.3.0/bin/Rscript
library(data.table)
OUT_DIR <- Sys.getenv('OUT_DIR')

mapping.vcf <- fread(paste0(OUT_DIR, "/tables/mapping_file_rsids.vcf"), skip=5, drop=c(6,7,8))
names(mapping.vcf)[3] <- "dbSNP"
with.rsid <- grepl("rs*", mapping.vcf$dbSNP)
mapping.vcf$dbSNP[!with.rsid] <- "."
print(paste(sum(with.rsid), "SNPs have an rsID which leaves", nrow(mapping.vcf) - sum(with.rsid), "snps without an rsid"))

mapping.vcf$COORDID <- paste0(mapping.vcf$"#CHROM", ":", mapping.vcf$POS, "_", mapping.vcf$REF, "_", mapping.vcf$ALT)

write.csv(mapping.vcf[, c("COORDID", "dbSNP")], paste0(OUT_DIR, "/tables/mapping_file_rsid.csv"), row.names=F, quote=F)
write.csv(mapping.vcf, paste0(OUT_DIR, "/tables/mapping_file_all.csv"), row.names=F, quote=F)
