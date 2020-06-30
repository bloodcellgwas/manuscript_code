OUT_DIR = Sys.getenv('OUT_DIR')
n = read.csv(paste0(OUT_DIR, "/BCX_final_output_raw.csv"), stringsAsFactors=F)

for (trait in unique(n$trait)) {
    if (trait %in% c("ret")) {next}
    ns = n[n$trait == trait,]
    threshold = table(ns$clump_id)[table(ns$clump_id) > 1]
    if (length(threshold) > 1) {
        stop("a")
    }
}
