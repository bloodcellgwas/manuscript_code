#!/software/R-3.3.0/bin/Rscript
library(data.table)
OUT_DIR <- Sys.getenv('OUT_DIR')
PHENO_NAMES <- Sys.getenv("PHENO_NAMES")
traits <- as.vector(read.table(PHENO_NAMES, header=F, stringsAsFactors=F)$V1)

# read plink.ld file and create hg19 ids for each pair
plink.ld <- fread(paste0(OUT_DIR, "/condout/ldclump_dosage/dosage.ld"), header=T, stringsAsFactors=F)
plink.ld$hg19A <- paste0(plink.ld$CHR_A, ":", plink.ld$BP_A)
plink.ld$hg19B <- paste0(plink.ld$CHR_B, ":", plink.ld$BP_B)
plink.ld$R2 <- as.numeric(plink.ld$R2)
all.vars <- unique(c(plink.ld$hg19A, plink.ld$hg19B))
plink.ld <- plink.ld[plink.ld$R2 > 0.8,]

## loop through all the vars
##  if var is in LD > 0.8 with any previous clump, merge into that clump
##  if var is in LD > 0.8 with more than one previous clump, merge both

hg19_table <- data.frame(hg19=all.vars, clump_id=rep("NA", length(all.vars)), stringsAsFactors=F)
hg19_table$clump_id <- as.character(hg19_table$clump_id)
hg19_table$hg19 <- as.character(hg19_table$hg19)
chr_num = gsub("(.+):(.+)", "\\1", hg19_table$hg19)
chr_num[chr_num %in% c("XY", "X")] = 23
hg19_table$sort = as.numeric(chr_num)*10^9 + as.numeric(gsub("(.+):(.+)", "\\2", hg19_table$hg19))
hg19_table = hg19_table[order(hg19_table$sort),c("hg19", "clump_id")]
print(length(all.vars))

get_partners <- function(var) {
  partnersA <- plink.ld[plink.ld$hg19A == var & plink.ld$R2 > 0.8,]$hg19B
  partnersB <- plink.ld[plink.ld$hg19B == var & plink.ld$R2 > 0.8,]$hg19A
  partners <- unique(c(partnersA, partnersB))
  return(partners)
}
compare_item_clumps <- function(partners, clump) {
  if (sum(clump %in% partners) > 0) {return(TRUE)}
  else {return(FALSE)}
}
clump_function <- function(partners, clumps, item) {
  if (length(clumps) == 0) {
    clumps[[length(clumps)+1]] <- unique(c(partners, item))
    return(clumps)
  }
  matching_clumps <- c()
  if (length(partners) != 0) {
    # loop through all the clumps to check if any of the partners already exist in a clump
    for (clumpi in 1:length(clumps)) {
      clump <- clumps[[clumpi]]
      # if this clump contains a partner add the index for this clump to the matching_clumps vector
      if (compare_item_clumps(partners, clump) == TRUE) {matching_clumps <- c(matching_clumps, clumpi)}
    }
  }
  # if there are no matching clumps add this item & its partners to a new clump
  if (length(matching_clumps) == 0) {
    clumps[[length(clumps)+1]] <- unique(c(partners, item))
  } else if (length(matching_clumps) == 1) {
    clumps[[matching_clumps]] <- unique(c(clumps[[matching_clumps]], item, partners))
  } else if (length(matching_clumps) > 1) {
    # create a new clumps merging both previous
    clumps[[length(clumps)+1]] <- unique(c(item, partners, unlist(clumps[matching_clumps])))
    # delete the previously matching clumps
    clumps[matching_clumps] <- NULL
  } else {
    stop("clump_function(): error with matching clumps variable")
  }
  return(clumps)
}

clumps <- list()
i = 1
for (var in hg19_table$hg19) {
  if (i %% 100 == 0) {print(i)}
  i = i + 1
  partners <- get_partners(var)
  clumps <- clump_function(partners, clumps, var)
}

# now we need to recreate hg19 table based on the new clumps list object
hg19_table_clump_col <- c()
for (hg19 in hg19_table$hg19) {
  clumpi = which(rapply(clumps, function(x) hg19 %in% x))
  hg19_table_clump_col <- c(hg19_table_clump_col, paste0("chr:", clumpi))
  if (length(clumpi) == 0) {
    stop("wasn't assigned to any clump")
  } else if (length(clumpi) > 1) {
    stop("a snp was assigned to two clumps") }
}
hg19_table$clump_id <- hg19_table_clump_col

write.table(hg19_table, paste(OUT_DIR, "/condout/ldclump_dosage/hg19_clump_table.tsv", sep=""), sep="\t", quote=FALSE, row.names=F)
print(paste0("Overall we have ", length(clumps), " clumps"))
