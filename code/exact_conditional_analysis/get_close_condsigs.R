source(paste0(Sys.getenv('SCRIPT_DIR'), "/general_functions.R"))
OUT_DIR = Sys.getenv('OUT_DIR')
final = read.csv(paste0(OUT_DIR, "/BCX_final_output_raw.csv"), stringsAsFactors=F)
final$c = paste(final$trait, final$clump)
print(final[duplicated(final$c),])


final$rownum = final$CHR*10^9 + final$BP
final = final[order(final$rownum),]
window_start = c()
window_end = c()
window_chr = c()
window_trait = c()
threshold = 250e3
for (trait in unique(final$trait)) {
    finals = final[final$trait == trait,]
    for (i in 1:nrow(finals)) {
        
    }
}
