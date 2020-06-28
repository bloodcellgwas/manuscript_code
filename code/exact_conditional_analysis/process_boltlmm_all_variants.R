#!/usr/bin/Rscript
library(data.table)
library(plyr)
library(dplyr)
## remove exponention 1.36e6
options(scipen=999)
BOLT_OUT_DIR <- Sys.getenv('BOLT_OUT_DIR')
OUT_DIR <- Sys.getenv('OUT_DIR')
args <- commandArgs(TRUE)
trait = args[1]
print(paste("pheno", trait))
##exclude.snps <- fread(paste0(OUT_DIR, "/VarsInOld_notInNewSNPStats.txt"), header=F)$V1

total_snps = 0
## chr.sumstats <- foreach (traiti = 1:length(traits)) %dopar% {
boltlmm.out <- fread(sprintf("%s/%s_gwas_normalised_imputed_full_panel.out", BOLT_OUT_DIR, trait), colClasses=c("character", "character", rep("numeric", 2), rep("character", 2), rep("numeric", 10)), nThread=15)
boltlmm.out$CHR_num = boltlmm.out$CHR
boltlmm.xy = fread(sprintf("%s/%s_gwas_normalised_imputed_full_panel_chrXY.out", BOLT_OUT_DIR, trait), colClasses=c("character", "character", rep("numeric", 2), rep("character", 2), rep("numeric", 10)), nThread=15)
boltlmm.xy$CHR = "XY"
  boltlmm.xy$CHR_num = 23
  boltlmm.x = fread(sprintf("%s/%s_gwas_normalised_imputed_full_panel_chrX.out", BOLT_OUT_DIR, trait), colClasses=c("character", "character", rep("numeric", 2), rep("character", 2), rep("numeric", 10)), nThread=15)
  boltlmm.x$CHR = "X"
  boltlmm.x$CHR_num = 23
  boltlmm.out = rbind(boltlmm.out, boltlmm.xy, boltlmm.x)
  rm(boltlmm.xy)
  rm(boltlmm.x)
  total_snps = total_snps + nrow(boltlmm.out)
  print(total_snps)
### AND INFOSCORE ###
  ## we will filter MAF at the end
  boltlmm.out <- boltlmm.out[boltlmm.out$INFO > 0.4,]
### done filter ###
  ## add variant id
  boltlmm.out$VARIANT <- paste(boltlmm.out$CHR, ":", boltlmm.out$BP, "_", boltlmm.out$ALLELE1, "_", boltlmm.out$ALLELE0, sep="")
  
  if (sum(duplicated(boltlmm.out$VARIANT)) != 0) {stop("there are duplicated COORDIDs")}
  ## remove duplicates
  ##dups <- read.table(paste("/lustre/scratch114/projects/ukbiobank_t151/duplicate_variants/int_dups_chr", chr, sep=""), header=T)
  ## filter out duplicates
  ##  boltlmm.out <- boltlmm.out[!is.element(boltlmm.out$VARIANT, exclude.snps),]
  ##dups.b <- boltlmm.out[is.element(boltlmm.out$VARIANT, dups$COORDID),]

  
### SORT THE SNPS BY LOCATION IN GENOME ###
  boltlmm.out$sort_chrom=boltlmm.out$CHR_num
  boltlmm.out$sort_order=as.numeric(gsub("(.+):(.+)_(.+)_(.+)", "\\2", boltlmm.out$VARIANT))+10^9*as.numeric(boltlmm.out$sort_chrom)
  ## set sortorder col to key
  setkey(boltlmm.out, "sort_order")
  ## clean up
  boltlmm.out[,sort_order:=NULL]
  boltlmm.out[,sort_chrom:=NULL]
### DONE SORT ###

### FLIP ALLELE1 AND ALLELE2, * MINUS THE BETA AND FLIP THE FREQUENCY ####
  boltlmm.out$ALT_FREQ=1-boltlmm.out$A1FREQ
  boltlmm.out$EFFECT=-boltlmm.out$BETA
                                        # get rid of old columns
  boltlmm.out[,A1FREQ:=NULL]
  boltlmm.out[,BETA:=NULL]
### DONE FLIP ######
                                        # GET RID OF NON-INF P VALUE
  boltlmm.out[,P_BOLT_LMM:=NULL]
  boltlmm.out$R2=2.0*as.numeric(boltlmm.out$ALT_FREQ)*(1.0-as.numeric(boltlmm.out$ALT_FREQ))*boltlmm.out$EFFECT^2
  boltlmm.out$MLOG10P=-pchisq((boltlmm.out$EFFECT/boltlmm.out$SE)^2,df=1,lower.tail=FALSE,log.p=TRUE)/log(10)
  boltlmm.out$GWSIG=boltlmm.out$MLOG10P>9-log10(8.31)
  setnames(boltlmm.out, "SNP", "ID")
  setnames(boltlmm.out, "P_BOLT_LMM_INF", "P")
  setnames(boltlmm.out, "ALLELE1", "REF")
  setnames(boltlmm.out, "ALLELE0", "ALT")

  boltlmm.out$DIRECTION <- c("-", "+")[as.numeric(boltlmm.out$EFFECT>0)+1]
  boltlmm.out$ALT_MINOR <- as.numeric(boltlmm.out$ALT_FREQ)<0.5

                                        # find the frequency of the minor allele, this will be the same as
                                        # frequency of alternate allele, except in cases where alternate allele is not minor
  boltlmm.out$MA_FREQ <- as.numeric(boltlmm.out$ALT_FREQ)
  boltlmm.out$MA_FREQ[!boltlmm.out$ALT_MINOR] <- 1 - boltlmm.out$MA_FREQ[!boltlmm.out$ALT_MINOR]

### FILTER MAF ###
  boltlmm.out <- boltlmm.out[boltlmm.out$MA_FREQ > 0.00005,]
### done maf filter ###
  print(head(boltlmm.out))

## filter out all the other columns
  boltlmm.out <- boltlmm.out[, c("VARIANT", "ID", "CHR", "CHR_num", "BP", "GENPOS", "REF", "ALT", "EFFECT", "P", "MLOG10P", "SE", "ALT_FREQ", "ALT_MINOR", "MA_FREQ", "R2", "GWSIG", "INFO")]
  setcolorder(boltlmm.out, c("VARIANT", "ID", "CHR", "CHR_num", "BP", "GENPOS", "REF", "ALT", "EFFECT", "P", "MLOG10P", "SE", "ALT_FREQ", "ALT_MINOR", "MA_FREQ", "R2", "GWSIG", "INFO"))
fwrite(boltlmm.out, file=paste0(OUT_DIR, "/all_variants/", trait, ".assoc", sep=""), row.names=F, quote=F)
fwrite(boltlmm.out[boltlmm.out$GWSIG == TRUE,], file=paste0(OUT_DIR, "/assoc_files/", trait, "_gwsig.assoc", sep=""), row.names=F, quote=F)
## just used the SNP ids that were already in the file, didn't add from table
print(total_snps)
