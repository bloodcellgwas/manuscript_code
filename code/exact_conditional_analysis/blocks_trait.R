OUT_DIR = Sys.getenv('OUT_DIR')
library(data.table)
block_table=fread(sprintf("%s/blocks/blocks.tsv", OUT_DIR)) 
