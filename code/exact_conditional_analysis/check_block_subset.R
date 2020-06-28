#!/software/R-3.3.0/bin/Rscript
OUT_DIR <- Sys.getenv('OUT_DIR')
args <- commandArgs(TRUE)
BLOCK = args[1]
#BLOCK <- Sys.getenv('BLOCK')
NODE_NAME <- Sys.getenv('NODE_NAME')

output_block_ids <- read.table(paste0(OUT_DIR, "/tmp_gen/ids_", BLOCK, ".tsv"), header=F, colClasses = "character")
input_pull_ids <- read.table(paste0(OUT_DIR, "/blocks/pull_ids_block_", BLOCK, ".tsv"), header=T, colClasses="character")
input_pull_ids$VARIANT <- paste0(input_pull_ids$chromosome, ":", as.character(input_pull_ids$position), "_", input_pull_ids$alleleA, "_", input_pull_ids$alleleB)


if (nrow(output_block_ids) != nrow(input_pull_ids)) {
  extra_ids_extracted <- output_block_ids$V1[!is.element(output_block_ids$V1, input_pull_ids$VARIANT)]
  if (length(extra_ids_extracted) > 0) {
    print(paste("The following variant(s) exists in output and shouldn't", paste(extra_ids_extracted, collapse=", ")))
    error.df <- data.frame(block=rep(BLOCK, length(extra_ids_extracted)), variant_id=extra_ids_extracted, node_name=rep(NODE_NAME, length(extra_ids_extracted)))
    write.table(error.df, paste0(OUT_DIR, "/tmp_gen_mistakes/extra_extraction.tsv"), append=T, row.names=F, col.names=F, quote=F)
  }
  missing_ids_extracted <- input_pull_ids$VARIANT[!is.element(input_pull_ids$VARIANT, output_block_ids$V1)]
  if (length(missing_ids_extracted) > 0) {
    error.df <- data.frame(block=rep(BLOCK, length(missing_ids_extracted)), variant_id=missing_ids_extracted, node_name=rep(NODE_NAME, length(missing_ids_extracted)))
    write.table(error.df, paste0(OUT_DIR, "/tmp_gen_mistakes/missing_extraction.tsv"), append=T, row.names=F, col.names=F, quote=F)
    print(paste("The following variant(s) don't exist in the output and they should do", paste(missing_ids_extracted, collapse=", ")))
  }
  file.rename(paste0(OUT_DIR, "/tmp_gen/ids_", BLOCK, ".tsv"), paste0(OUT_DIR, "/tmp_gen/ids_", BLOCK, ".tsv_WRONG"))
  stop(paste("doesn't match, output has", nrow(output_block_ids), "and input has", nrow(input_pull_ids)))
} else if (sum(output_block_ids$V1 == input_pull_ids$VARIANT) != nrow(input_pull_ids)) {
    # order doesnt match here but all variants might still be in output file
    # fine as long as we dont use the pull file to label the genotypes
    if (sum(input_pull_ids$VARIANT %in% output_block_ids$V1) == nrow(input_pull_ids)) {
        print("order doesn't match but all variants there")
    } else {
        file.rename(paste0(OUT_DIR, "/tmp_gen/ids_", BLOCK, ".tsv"), paste0(OUT_DIR, "/tmp_gen/ids_", BLOCK, ".tsv_WRONG"))
        stop("Doesn't match")
    }
} else if (sum(is.element(output_block_ids$V1, input_pull_ids$VARIANT)) == nrow(input_pull_ids)) {
                                        # rename id file because we use this as a measure that the script has completed properly
  print("variants match")
}

                                        # 1231 02:112359395_G_A 646
                                        # 1235 03:24339836_C_A  365
                                        # 2433 02:112359395_G_A 
                                        # block: 1751 lsb jobindex: 1585 chromosome 9 
