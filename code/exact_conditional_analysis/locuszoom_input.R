#!/software/R-3.3.0/bin/Rscript
library(data.table)
OUT_DIR = Sys.getenv('OUT_DIR')
args <- commandArgs(TRUE)
trait <- args[1]

all_assoc = fread(paste0(OUT_DIR, "/all_variants/", trait, ".assoc"), select=c("ID", "P"))
names(all_assoc) = c("MarkerName", "P-value")

fwrite(all_assoc, paste0(OUT_DIR, "/locuszoom/input/", trait, ".tsv"), sep="\t")
