#!Rscript
library(data.table)
OUT_DIR = Sys.getenv('OUT_DIR')

## load the first and second sample files
sample_x = fread("/rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chrX.sample", colClasses="character")
first_line_x = sample_x[1,]
sample_x = sample_x[-1,]
sample_xy = fread("/rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chrXY.sample", colClasses="character")
first_line_xy = sample_xy[1,]
sample_xy = sample_xy[-1,]

## find the intersection of them both
ids_both = sample_x$ID_1[is.element(sample_x$ID_1, sample_xy$ID_2)]

## get sample ids for the intersection making sure both sample files are in the same order
sample_x = sample_x[match(ids_both, sample_x$ID_1),]
sample_xy = sample_xy[match(ids_both, sample_xy$ID_1),]

if (nrow(sample_x) != nrow(sample_xy)) {
    stop("rows aren't same")
} else if (sum(sample_x$ID_1 == sample_xy$ID_1) != nrow(sample_x)) {
    stop("aren't in the same order")
}
sample_x = rbind(first_line_x, sample_x)
sample_xy = rbind(first_line_x, sample_xy)
## save the sample file
fwrite(sample_x, paste0(OUT_DIR, "/sex_chrom_x.sample"), row.names=F, sep=" ", na="NA", quote=F)
fwrite(sample_xy, paste0(OUT_DIR, "/sex_chrom_xy.sample"), row.names=F, sep=" ", na="NA", quote=F)
