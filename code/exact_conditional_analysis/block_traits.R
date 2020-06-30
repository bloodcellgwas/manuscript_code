#!/software/R-3.1.2/bin/Rscript
# SPLITS THE GWSIG SNPS INTO BLOCKS
library(data.table)
library(doMC)
# remove exponention 1.36e8 etc
options(scipen=999)
PHENO_NAMES <- Sys.getenv('PHENO_NAMES')
OUT_DIR <- Sys.getenv('OUT_DIR')
CORES <- Sys.getenv('CORES')
CORES=10
pheno_table=read.table(PHENO_NAMES)$V1
block_table=NULL
num_blocks=0

convert_variant_id_double_remove_exp <- function(variant) {
  bp <- as.integer(gsub(".+:(.+)_.+_.+","\\1", variant))
  ref <- gsub(".+:.+_(.+)_.+","\\1", variant)
  alt <- gsub(".+:.+_.+_(.+)","\\1", variant)
  chr <- sprintf("%02d", as.numeric(gsub("(.+):.+_.+_.+","\\1", variant)))
  variant_output <- paste0(chr, ":", bp, "_", ref, "_", alt)
  return(variant_output)
}
variant_id_remove_exp <- function(variant) {
  bp <- as.integer(gsub(".+:(.+)_.+_.+","\\1", variant))
  ref <- gsub(".+:.+_(.+)_.+","\\1", variant)
  alt <- gsub(".+:.+_.+_(.+)","\\1", variant)
  chr <- as.numeric(gsub("(.+):.+_.+_.+","\\1", variant))
  variant_output <- paste0(chr, ":", bp, "_", ref, "_", alt)
  return(variant_output)
}

for(pheno in pheno_table)
{
	summary_stats=fread(paste(OUT_DIR, "/assoc_files/", pheno, "_gwsig.assoc", sep=""), colClasses=c("character"))
	summary_stats$dist=10^9*as.integer(summary_stats$CHR_num)+as.integer(summary_stats$BP)
	summary_stats$BLOCK=c(0,cumsum(diff(summary_stats$dist)>5*10^6))+1+num_blocks
	num_blocks=max(summary_stats$BLOCK)
	summary_stats$TRAIT=pheno
	block_table=rbind(block_table, summary_stats[,c("VARIANT","TRAIT","BLOCK", "CHR", "CHR_num", "BP", "REF","ALT"),with=FALSE])
}


## loop through all the blocks, if any are larger than the max size then split and break the for loop and start looping through
## all the blocks again
print(paste("number of blocks before splitting large:", num_blocks))
max_block_size = 2500
done = FALSE
cur_position = 1
while (done == FALSE) {
    for (block in cur_position:max(block_table$BLOCK)) {
        if (nrow(block_table[block_table$BLOCK == block,]) > max_block_size) {
            ## push up the number of all following blocks because the current block will be split into two
            block_table$BLOCK[block_table$BLOCK > block] = block_table$BLOCK[block_table$BLOCK > block] + 1
            ## split the current block
            ## variants in this block which are going to be put into the next block
            this_block_shift = which(block_table$BLOCK == block)[(max_block_size+1):length(which(block_table$BLOCK == block))] 
            block_table$BLOCK[this_block_shift] = block_table$BLOCK[this_block_shift] + 1
            done = FALSE
            break
        }
        cur_position = block - 1
        done = TRUE
    }
}
num_blocks = max(block_table$BLOCK)
print(paste("number of blocks before splitting large:", num_blocks))

write.table(block_table, file=sprintf("%s/blocks/blocks.tsv", OUT_DIR), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(unique(block_table[,c("BLOCK","CHR", "CHR_num"), with=FALSE]), file=sprintf("%s/blocks/block_chrs.tsv", OUT_DIR), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)	

allblocks <- unique(block_table$BLOCK)

registerDoMC(cores=CORES)
blockloop <- foreach (i = 1:length(allblocks)) %dopar% {
  block = allblocks[i]
  pull_table_int=cbind(
    block_table$VARIANT[block_table$BLOCK==block],
    block_table$VARIANT[block_table$BLOCK==block],
    as.character(block_table$CHR[block_table$BLOCK==block]),
    as.character(as.numeric(block_table$BP[block_table$BLOCK==block])),
    block_table$REF[block_table$BLOCK==block],
    block_table$ALT[block_table$BLOCK==block]
  )
  colnames(pull_table_int)=c("SNPID","rsid","chromosome","position","alleleA","alleleB")
  pull_table_int = as.data.frame(pull_table_int)
  print(block)
  #write.table(pull_table_int,paste(OUT_DIR, "/blocks/pull_ids_block_", block, ".tsv", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  fwrite(pull_table_int, paste(OUT_DIR, "/blocks/pull_ids_block_", block, ".tsv", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}
# 1770
