glibrary(plyr)
library(rhdf5)
library(data.table)

OUT_DIR = Sys.getenv('OUT_DIR')

## get vector of all blocks which have completed
completed = list.files(paste0(OUT_DIR, "/tmp_hd5/"))
completed = completed[grepl(".h5", completed)]
completed = ldply(strsplit(completed, split="_"))[[2]]
completed = ldply(strsplit(completed, split=".h5"))[[1]]


loadhdf5data <- function(h5File) {
    listing <- h5ls(h5File)
    # Find all data nodes, values are stored in *_values and corresponding column
    # titles in *_items
    data_nodes <- grep("_values", listing$name)
    name_nodes <- grep("_items", listing$name)
    data_paths = paste(listing$group[data_nodes], listing$name[data_nodes], sep = "/")
    name_paths = paste(listing$group[name_nodes], listing$name[name_nodes], sep = "/")
    rows = list()
    for (idx in seq(data_paths)) {
          # NOTE: matrices returned by h5read have to be transposed to obtain
          # required Fortran order!
          data <- data.frame(h5read(h5File, data_paths[idx]))
          names <- h5read(h5File, name_paths[idx])
          entry <- data.frame(data)
          rownames(entry) <- names
          rows <- append(rows, entry)
        }

    data <- data.frame(rows)
    rownames(data) = h5read(h5File, "/genotype_data/axis0")
    colnames(data) = h5read(h5File, "/genotype_data/axis1")
    return(data)
}
read_add_data<-function(block_id) {
    add_genotypes <- h5read(sprintf("%s/tmp_hd5/block_%s.h5", OUT_DIR, block), "/add")
    int_ids=fread(sprintf("%s/tmp_gen/ids_%s.tsv",OUT_DIR, block), head=FALSE, sep="\t")
    colnames(add_genotypes) = c("intercept", int_ids$V1)
    print(paste(block_id, "loaded"))
    return(add_genotypes)
}
chrom_prev = 0
for (block in completed) {
    print(block)
    if (file.exists(paste0(OUT_DIR, "/tmp_hd5_maf_check/", block, "_check_done.txt"))) {
        next
    }
    model_genotypes = read_add_data(block)
    chrom = strsplit(colnames(model_genotypes)[2], ":")[[1]][1]
    if (chrom %in% c("X", "XY")) {next}
    if (chrom != chrom_prev) {
        snpstats_file = fread(paste0("/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_meta_analysis", "/snpstats/snpstats_",chrom,".csv"))
        chrom_prev = chrom
    }
    
    differences = c()
    snps = colnames(model_genotypes)[2:ncol(model_genotypes)]
    for (snpi in 2:ncol(model_genotypes)) {
        snp = colnames(model_genotypes)[snpi]
        if (!is.element(snp,snpstats_file$alternate_ids)) next
        hd5_maf = mean(model_genotypes[, snpi]+1)/2
        snpstats_maf = snpstats_file[snpstats_file$alternate_ids == snp,]$alleleB_frequency
        diff_pct = abs(hd5_maf-snpstats_maf)
        ## if the difference is larger than 1% throw an error
        if (diff_pct > 1) {stop(paste("Theres a mismatch", snpi, snp, hd5_maf, snpstats_maf, diff_pct))}
        differences = c(differences, diff_pct)
    }
    print(summary(differences))
    sink(paste0(OUT_DIR, "/tmp_hd5_maf_check/", block, "_check_done.txt"))
    print(summary(differences))
    sink()
}
