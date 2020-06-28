#!/software/R-3.3.0/bin/Rscript
options(warn=2)
SUPPORT_FILES <- Sys.getenv('SUPPORT_FILES')
options(show.error.locations = TRUE)
options(show.warning.locations = TRUE)
library(data.table)
library(doMC)

OUT_DIR <- Sys.getenv("OUT_DIR")
PHENO_NAMES <- Sys.getenv("PHENO_NAMES")
args <- commandArgs(TRUE)
COND_OR_GWSIG <- args[1]
SUPPORT_FILES <- Sys.getenv('SUPPORT_FILES')
MAPPING_FILE_RSID <- Sys.getenv('MAPPING_FILE_RSID')
COND_OR_GWSIG="cond"
traits <- as.vector(read.table(PHENO_NAMES, header=F, stringsAsFactors=F)$V1)
print(traits)
# load VARIANT -> dbSNP49 table
rsid_table <- fread(MAPPING_FILE_RSID, colClasses=c("character"))
names(rsid_table) <- c("VARIANT", "ID")

# load user friendly file name function
source("general_functions.R")
#traits = c("eo", "wbc", "baso", "mch", "hct", "hgb", "mcv")
# 3:50035015_T_C mrv
for (traiti in 1:length(traits))  {
  trait = traits[traiti]
  print(trait)
  
  if (COND_OR_GWSIG == "cond") {
      if (!file.exists(paste(OUT_DIR, "/condout/results/condsig_", trait, "_gwas_normalised.tsv", sep=""))) {next}
      cond_snps <- read.table(paste(OUT_DIR, "/condout/results/condsig_", trait, "_gwas_normalised.tsv", sep=""), header=TRUE, stringsAsFactors=FALSE)
      gwsig_snps <- read.csv(paste0(OUT_DIR, "/assoc_files/", trait, "_gwsig.assoc"), header=TRUE, stringsAsFactors=FALSE)
      ##gwsig_snps <- read.table(paste0(OUT_DIR, "/assoc_files/", trait, "_gwsig.assoc"), header=TRUE, stringsAsFactors=FALSE)
      cond_igv <- merge(x = cond_snps, y=gwsig_snps, by="VARIANT")

      ## MERGE IN CLUMP AS WELL
      ## IVE PLANNED THIS PART OF THE ANALYSIS AND LEAVING THE CODE HERE BUT WILL IMPLEMENT LATER
      cond_igv$hg19 <- paste(cond_igv$CHR, ":", cond_igv$BP, sep="")
      hg19_clump_table <- read.table(paste(OUT_DIR, "/condout/ldclump_dosage/hg19_clump_table.tsv", sep=""), header=T, stringsAsFactors=F)
      cond_igv <- merge(cond_igv, hg19_clump_table, by="hg19")
      if (nrow(cond_igv) != nrow(cond_snps)) {print("there are rows being thrown out");print(cond_snps[!is.element(cond_snps$VARIANT, cond_igv$VARIANT),])}
      ## read vep annotation
      vep.table <- read.table(paste(OUT_DIR, "/condout/vep/tabular_all_condsig.txt", sep=""), header=T, stringsAsFactors=F)
    
      ## quit if the resultant table has lost rows, means merge hasn't worked
      ## rows will have been lost because the rare ones dont get clumped because of the way PLINK created hard aclled files form the BGENs
      if (nrow(cond_igv) != nrow(cond_snps)) {
          ## there are some SNPs in cond_snps which do not exist in gwsig_snps (not possible)
          stop(paste("Some VARIANT ids in condsig_", trait, "_gwas_normalised.tsv which do not exist in .qassoc file", sep=""))
      }
  
  } else if (COND_OR_GWSIG == "gwsig") {
      cond_igv <- read.table(paste(OUT_DIR, "/assoc_files/", trait, "_gwsig.assoc", sep=""), header=TRUE, stringsAsFactors=FALSE, sep=",")
      ## read vep annotation
      vep.table <- read.table(paste(OUT_DIR, "/vep/tabular_", trait, ".txt", sep=""), header=T, stringsAsFactors=F)
  }

  # add infoscore
  cond_igv$INFO <- get_info_score(cond_igv$VARIANT)
  
  # rename id to variant, use this to merge with previous
  names(vep.table)[names(vep.table) == "ID"] <- "VARIANT"
  vep.table <- vep.table[, names(vep.table) %in% c("VARIANT", "Consequence", "IMPACT", "SYMBOL", "Gene")]
  cond_igv <- merge(cond_igv, vep.table, by="VARIANT", all.x=T)
  
  # rename the titles to be more user friendly
  names(cond_igv)[names(cond_igv) == "Consequence"] = "VEP_CONSEQUENCE"
  names(cond_igv)[names(cond_igv) == "IMPACT"] = "VEP_IMPACT"
  names(cond_igv)[names(cond_igv) == "SYMBOL"] = "VEP_GENE_SYMBOL"
  names(cond_igv)[names(cond_igv) == "Gene"] = "VEP_GENE_ENSEMBL_ID"  
  
  cond_igv$P[log10(cond_igv$P)<(-300)]=10^(-300)

  # delete the id column which corresponds to the old dbSNP
  cond_igv <- cond_igv[, !(names(cond_igv) == "ID")]
  cond_igv$ID = get_rsid_from_variantid(cond_igv$VARIANT)
  #cond_igv <- merge(cond_igv, rsid_table, by="VARIANT", all.x=TRUE)

  cond_igv$sort=as.numeric(cond_igv$CHR_num)*10^9+as.numeric(cond_igv$BP)
  cond_igv <- cond_igv[order(cond_igv$sort),]
  cond_igv <- cond_igv[, !(names(cond_igv) == "sort")]
  cond_igv$VARIANT = paste0(cond_igv$CHR, ":", cond_igv$BP, "_", cond_igv$REF, "_", cond_igv$ALT)
  
  if (COND_OR_GWSIG == "cond") {
    names(cond_igv)[names(cond_igv) == "JOINT_BETA"] <- "JOINT_EFFECT"
    names(cond_igv)[names(cond_igv) == "clump_id"] <- "CLUMP"
    n = c("VARIANT", "ID", "INFO", "CHR", "BP", "GENPOS", "REF", "ALT", "ALT_MINOR", "JOINT_EFFECT", "JOINT_SE", "JOINT_MLOG10P", "P", "EFFECT", "SE", "MLOG10P", "ALT_FREQ", "MA_FREQ", "R2", "VEP_CONSEQUENCE", "VEP_IMPACT", "VEP_GENE_SYMBOL", "VEP_GENE_ENSEMBL_ID")

    cond_igv <- cond_igv[c("VARIANT", "ID", "INFO", "CHR", "BP", "GENPOS", "REF", "ALT", "ALT_MINOR", "JOINT_EFFECT", "JOINT_SE", "JOINT_MLOG10P", "P", "EFFECT", "SE", "MLOG10P", "ALT_FREQ", "MA_FREQ", "R2", "VEP_CONSEQUENCE", "VEP_IMPACT", "VEP_GENE_SYMBOL", "VEP_GENE_ENSEMBL_ID")]
    # removed clump above
    # cond_igv <- cond_igv[c("VARIANT", "ID", "INFO", "CHR", "BP", "GENPOS", "REF", "ALT", "ALT_MINOR", "CLUMP", "JOINT_EFFECT", "JOINT_SE", "JOINT_MLOG10P", "P", "EFFECT", "SE", "MLOG10P", "ALT_FREQ", "MA_FREQ", "R2", "VEP_CONSEQUENCE", "VEP_IMPACT", "VEP_GENE_SYMBOL", "VEP_GENE_ENSEMBL_ID")]
    write.table(cond_igv, file=paste(OUT_DIR, "/condout/igv/", user_friendly_file_name(trait), "_condsig.gwas", sep=""), quote=FALSE, row.names=FALSE, na = "NA", sep="\t")
    write.table(cond_igv, file=paste(OUT_DIR, "/condout/igv_orig_names/", trait, "_condsig.gwas", sep=""), quote=FALSE, row.names=FALSE, na = "NA", sep="\t")
  } else if (COND_OR_GWSIG == "gwsig") {
    cond_igv <- cond_igv[c("VARIANT", "ID", "INFO", "CHR", "BP", "GENPOS", "REF", "ALT", "ALT_MINOR", "EFFECT", "SE", "P", "MLOG10P", "ALT_FREQ", "MA_FREQ", "R2",  "VEP_CONSEQUENCE", "VEP_IMPACT", "VEP_GENE_SYMBOL", "VEP_GENE_ENSEMBL_ID")]
    write.table(cond_igv, file=paste(OUT_DIR, "/assoc_files/igv/", user_friendly_file_name(trait), "_gwsig.gwas", sep=""), quote=FALSE, row.names=FALSE, na = "NA", sep="\t")
  }
}
