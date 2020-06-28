# not bothering to save XY and X chroms, to do this copy style from: gen_gwsig_snps.R
OUT_DIR = Sys.getenv('OUT_DIR')
astle_variants = read.csv(paste0(OUT_DIR, "/astle_compare/mmc4.csv"), stringsAsFactors=F)$Unique.Variant.ID
chr = gsub("(.+):(.+)_(.+)_(.+)", "\\1", astle_variants)
chr = sprintf("%02d", as.numeric(chr))
bp = gsub("(.+):(.+)_(.+)_(.+)", "\\2", astle_variants)
pos = paste0(chr, ":", bp)

traits = read.table(Sys.getenv('PHENO_NAMES'), header=F, stringsAsFactors=F)$V1
for (trait in traits) {
    print(trait)
    condsig_var = read.table(paste0(OUT_DIR, "/condout/results/condsig_", trait, "_gwas_normalised.tsv"), sep="\t", stringsAsFactors=F, header=T)
    chr = gsub("(.+):(.+)_(.+)_(.+)", "\\1", condsig_var$VARIANT)
    chr = sprintf("%02d", as.numeric(chr))
    bp = gsub("(.+):(.+)_(.+)_(.+)", "\\2", condsig_var$VARIANT)
    pos = unique(c(pos, paste0(chr, ":", bp)))
}
# assuming no sex chromosomes
order = data.frame(chr_pos = pos, chr = gsub("(.+):(.+)", "\\1", pos), pos = gsub("(.+):(.+)", "\\2", pos), stringsAsFactors=F)
order$chr_num = order$chr
order$chr_num[is.element(order$chr_num, c("X", "XY"))] = 23
order$chr_num = as.numeric(order$chr_num)
order$pos = as.numeric(order$pos)
order$sort = order$pos + order$chr_num*10^9
order = na.omit(order[order(order$sort),])
for (chr in 1:22) {
    subset_order = order[order$chr_num == chr,]
    write.table(subset_order$chr_pos, paste0(OUT_DIR, "/astle_compare/pull_tables/", chr, ".csv"), row.names=F, quote=F, col.names=F)
}
