library(data.table)
OUT_DIR = Sys.getenv('OUT_DIR')
bcx_trait_names = c("BASO", "EOS", "HCT", "HGB", "LYM", "MCHC", "MCH", "MCV", "MONO", "MPV", "NEU", "PLT", "RBC", "RDW", "WBC")
trait_names = c("baso", "eo", "hct", "hgb", "lymph", "mchc", "mch", "mcv", "mono", "mpv", "neut", "plt", "rbc", "rdw", "wbc")

for (i in 1:length(trait_names)) {
    ukbb_variants = fread(paste0(OUT_DIR, "/assoc_files/", trait_names[i], "_gwsig.assoc"), header=T)
    bcx_variants = fread(paste0("/rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/ukbb500k/FINAL_EA_META/", bcx_trait_names[i], "_Eur_GWAMA_11Jul2018.meta"), header=T)
    bcx_variants = bcx_variants[bcx_variants$"_-log10_p-value" > 8.080399,]
    print(summary(as.factor(gsub("(.+):.+_.+_.+","\\1", bcx_variants$rs_number))))
}
