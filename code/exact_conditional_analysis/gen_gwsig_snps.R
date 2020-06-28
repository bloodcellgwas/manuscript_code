#!/software/R-3.1.2/bin/Rscript
library(data.table)
OUT_DIR <- Sys.getenv('OUT_DIR')

generate_positions = TRUE

# table of all SNPs and their minimum p value across all traits
data <- fread(paste(OUT_DIR, "/all_min.tsv", sep=""))
print(OUT_DIR)
# include genomewide significant SNPs only
data <- data[data$GWSIG,]

# generate REF and ALT columns from the VARIANT ID
# VARIANT ID is CHR:BP_REF_ALT
data$REF <- gsub("(.+):(.+)_(.+)_(.+)","\\3",data$VARIANT)
data$ALT <- gsub("(.+):(.+)_(.+)_(.+)","\\4",data$VARIANT)

#library(doMC)
#registerDoMC(cores=LOCAL_CORES)

# helps map rsID to VARIANT id
##mapping <- fread(MAPPING_FILE)

for (chr in c(1:22, "X", "XY"))
{
    if (chr %in% 1:22) {
        chromosome = sprintf("%02d", as.numeric(data$CHR[data$CHR==chr]))
    } else {
        chromosome = chr
    }
    pull_table <- data.frame(ID=paste0(chromosome, ":", data$BP[data$CHR==chr]))
    write.table(pull_table, paste0(OUT_DIR, "/genfiles/gwsig_snps_chr",chr,"_pos.tsv"), row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
    print(chr)
}
