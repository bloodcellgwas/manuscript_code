library(data.table)
OUT_DIR = Sys.getenv('OUT_DIR')

plink_ld = fread(paste0(OUT_DIR, "/sysmex/condout/ldclump_compare/plink.ld"))
plink_ld = plink_ld[plink_ld$R2 >= 0.8,]
fwrite(plink_ld, paste0(OUT_DIR, "/sysmex/condout/ldclump_compare/plink_slim08.ld"))
