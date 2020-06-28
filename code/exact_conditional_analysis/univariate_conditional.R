#!/software/R-3.3.3/bin/Rscript
# conditional analysis
.libPaths("/home/pa354/R_libs_3_3_3")
options(show.error.locations = TRUE)
block <- as.integer(Sys.getenv("BLOCK"))
print(block)
library(doMC)
library(data.table)
library(RcppEigen)
args <- commandArgs(TRUE)
CORES <- args[1]
if (CORES > 30) {CORES=20}
OUT_DIR <- Sys.getenv('OUT_DIR')
PHENO_FILE <- Sys.getenv('PHENO_FILE')
PHENO_NAMES <- Sys.getenv('PHENO_NAMES')

block_chrs = read.csv(paste0(OUT_DIR, "/blocks/block_chrs.tsv"), header=T, sep="\t")
block_chrs = block_chrs[block_chrs$CHR_num == 23,]

if (as.character(block) %in% block_chrs$BLOCK) {
    sex = TRUE
} else {
    sex = FALSE
}
if (sex == FALSE) {
    int_phenotypes<-fread(PHENO_FILE,head=TRUE)
} else {
    int_phenotypes = fread(paste0(OUT_DIR, "/genfiles/chrXY_merged.sample"))
}
block_table=fread(sprintf("%s/blocks/blocks.tsv", OUT_DIR))
clear_files = FALSE

# load list of genotypes & individuals in this block
library(rhdf5)
print(block)
if (!file.exists(sprintf("%s/tmp_hd5/block_%d.h5", OUT_DIR, block))) {
  stop(paste0("hd5 file doesn't exist block_", block, ".hd5 you should delete the text file"))
} else {
  add_genotypes <- h5read(sprintf("%s/tmp_hd5/block_%d.h5", OUT_DIR, block), "/add")
  print("DONE LOADING HD5 FILE")
}

# check we have the right number of rows
if (nrow(int_phenotypes[-1,]) != nrow(add_genotypes)) {
    stop(paste("the rows aren't equal", nrow(add_genotypes), "versus", nrow(int_phenotypes[-1,])))
}

# list of IDs in this block only (COORDID)
int_ids=fread(sprintf("%s/tmp_gen/ids_%d.tsv",OUT_DIR, block), head=FALSE, sep="\t")
number_of_variants=dim(add_genotypes)[2]-1 # minus one because of the intercept

# just to make sure that the hd5 file included the right number of variants
# everything is will be crosschecked by the python file
output <- data.frame(vars=c(1:number_of_variants))
# THIS IS WRITTEN TO FILE AS LAST LINE IN THE SCRIPT
block_ids <- read.table(paste0(OUT_DIR, "/blocks/pull_ids_block_", block, ".tsv"), header=T, colClasses="character")
block_ids$VARIANT <- paste0(block_ids$chromosome, ":", block_ids$position, "_", block_ids$alleleA, "_", block_ids$alleleB)
if (number_of_variants != nrow(block_ids)) {
  stop("Conditional analysis number of variants loaded from hd5 file doesn't equal number in block")
} else if (sum(is.element(block_ids$VARIANT, int_ids$V1)) != nrow(block_ids)) {
  difference <- nrow(block_ids) - sum(is.element(block_ids$VARIANT, int_ids$V1))
  if (clear_files == TRUE) {
    # remove the blocks_done file and the first stage output file
    if (file.exists(paste0(OUT_DIR, "/condout/blocks/ind_var_block", block, ".tsv"))) {
      file.remove(paste0(OUT_DIR, "/condout/blocks/ind_var_block", block, ".tsv"))
    }
    if (file.exists(paste0(OUT_DIR, "/condout/blocks_done/block_", block, ".tsv"))) {
      file.remove(paste0(OUT_DIR, "/condout/blocks_done/block_", block, ".tsv"))
    }
    # remove the hd5 file and the hd5 done file
    if (file.exists(paste0(OUT_DIR, "/tmp_hd5/block_", block, ".txt"))) {
      file.remove(paste0(OUT_DIR, "/tmp_hd5/block_", block, ".txt"))
    }
    if (file.exists(paste0(OUT_DIR, "/tmp_hd5/block_", block, ".hd5"))) {
      file.remove(paste0(OUT_DIR, "/tmp_hd5/block_", block, ".hd5"))
    }
    # remove the .gen file and the .gen ids file
    if (file.exists(paste0(OUT_DIR, "/tmp_gen/ids_", block, ".tsv"))) {
      file.remove(paste0(OUT_DIR, "/tmp_gen/ids_", block, ".tsv"))
    }
    if (file.exists(paste0(OUT_DIR, "/tmp_gen/block_", block, ".gen"))) {
      file.remove(paste0(OUT_DIR, "/tmp_gen/block_", block, ".gen"))
    }
  }
  stop(paste("Not all variants in block ids exist in gen id file the difference is:", difference))
}

all_variants=unlist(int_ids[,1,with=FALSE])

# prepare phenotypes
pheno_table=read.table(file=PHENO_NAMES,stringsAsFactors=FALSE)[[1]]
pheno_table <- paste(pheno_table, "_gwas_normalised", sep="")

# merged_phenotypes: extracts phenotype values for the traits we are studying
# gets rid of all the other extra columns in the PHENO_FILE
merged_phenotypes=apply(int_phenotypes[2:dim(int_phenotypes)[1], match(pheno_table, colnames(int_phenotypes)), with=FALSE],2,as.numeric)

# extracts covariate values for interval
int_numerical_covariates=int_phenotypes[2:dim(int_phenotypes)[1],4:13,with=FALSE]

setnames(int_numerical_covariates,as.character(1:10))
int_size=dim(int_phenotypes)[1]-1
rm(int_phenotypes)
merged_numerical_covs=int_numerical_covariates

merged_numerical_covs=apply(merged_numerical_covs, 2, as.numeric)		

selected_pheno=colnames(merged_phenotypes)==paste(unique(block_table$TRAIT[block_table$BLOCK==block]), "_gwas_normalised", sep="")

non_missing=!is.na(merged_phenotypes[,selected_pheno])

# adjust for covariates
# interaction term for PGs
adjusted_lm<-lm(merged_phenotypes[non_missing,selected_pheno]~merged_numerical_covs[non_missing,])
merged_phenotypes[non_missing,selected_pheno]=resid(adjusted_lm,na.action=na.exclude)

registerDoMC(cores=CORES)
compute_add_pvalues<-function(add_var) {
  if (skip_variants[add_var] == TRUE) {
      return(-1)
  }
  proposed_variants=current_variants
  proposed_variants[add_var]=TRUE
  ## log p value
  logpvalue=0
  ## if we haven't added anything because add_var was already in the model skip this step and return p value of one
  if(!all(proposed_variants==current_variants))
    {
      ## make sure the variant we are adding isnt highly correlated iwth another variant in the model
      ## once we have final set in model we will go back and look for proxies in the data
      rsq=cor(model_genotypes[,c(FALSE, current_variants)], model_genotypes[,add_var+1])^2
      ## reduce to 0.9
      if(max(rsq)<0.9)
        {		
          test_index=sum(proposed_variants[1:add_var])
          proposed_model<-fastLm(model_genotypes[,c(TRUE, proposed_variants)], model_phenotype)
          ## t statistic, ratio of beta hat to the standard erorr (t statistic) with mean zero, df = number of data points - n
          logpvalue=log(2)+pt(abs(proposed_model$coefficients[1+test_index]/proposed_model$se[1+test_index]), df=proposed_model$df.residual, log=TRUE,lower=FALSE)	
        }
    }
  ## if the p value is less than -log10 2.5 then add this to the skip list
  if (logpvalue > -2) {
      skip_variants[add_var] = TRUE
  }
  return(logpvalue)
}



# the number of variants in the block (-1 because of intercept)
nvar=dim(add_genotypes)[2]-1

model_phenotype=as.matrix(merged_phenotypes[non_missing,selected_pheno, drop=FALSE])
model_genotypes=add_genotypes[non_missing,]
rm(add_genotypes)
gc()
looping=TRUE
adding=TRUE
current_variants=rep(FALSE, nvar)	
## we add variants to the skip variants list if they fall below -log10 2
## this means we don't keep testing variants which will probably never win
skip_variants = rep(FALSE, nvar)
while(looping)	
  {
    adding_failed=FALSE
    print(sprintf("number of variants = %d / %d", sum(current_variants), length(current_variants)))
    if(adding)
      {
          ## adding each variant to the model in turn (multithreaded)
        logadd_pvalues=unlist(foreach(n=1:nvar, .combine="c") %dopar% compute_add_pvalues(n))
          if (length(logadd_pvalues) != nvar) {
              stop("the lengths are not equal, add more memory")
          }
          skip_variants[logadd_pvalues > -2] = TRUE
          ## change to log p
        logminp=min(logadd_pvalues)
        print(logminp)
        print(log(8.31)-9*log(10))
        if(logminp<=log(8.31)-9*log(10))
          {
            print("adding")
            current_variants[which.min(logadd_pvalues)]=TRUE
          }else {print("nothing addable");adding_failed=TRUE}
        adding=FALSE	
      }
    print(sprintf("number of variants = %d / %d", sum(current_variants), length(current_variants)))
    
    if(sum(current_variants)>0)
      {
        logdrop_pvalues=rep(-Inf,nvar)
        
        current_model<-fastLm(model_genotypes[,c(TRUE, current_variants)], model_phenotype)
        logdrop_pvalues[current_variants]=log(2)+pt(abs(current_model$coefficients[1+1:sum(current_variants)]/current_model$se[1+1:sum(current_variants)]), df=current_model$df.residual, log=TRUE,lower=FALSE)
        logmaxp=max(logdrop_pvalues)
        if(logmaxp>log(8.31)-9*log(10))
          {
					# use which .max
            current_variants[which.max(logdrop_pvalues)]=FALSE
            print("dropping")
            if(sum(current_variants)==0)
              {
                adding=TRUE
              }  
          }
        else
          {
            print("nothing dropable")
            if(sum(current_variants)<length(current_variants)){adding=TRUE}
            if((adding_failed)||(sum(current_variants)==length(current_variants))){looping=FALSE}
          }
      }
    else
      {
          print("empty model: nothing to drop")
          adding=TRUE
        if(adding_failed) {
            ## therefore we are saving with an empty model
            pvalues = data.frame(VARIANT = all_variants, pvalues = logadd_pvalues)
            write.csv(pvalues, paste0(OUT_DIR, "/missing_blocks_pvalue/block_", block, ".tsv"), row.names=F)
            ## add the sentinel variant anyway
            current_variants[which.min(logadd_pvalues)]=TRUE
            looping=FALSE
        }
      }
  }
                                        # redo with gen derived names
surviving_variants=all_variants[current_variants]
surviving_table=data.frame(VARIANT=surviving_variants)
write.table(surviving_table,sprintf("%s/condout/blocks/ind_var_block%d.tsv", OUT_DIR, block), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(output, paste0(OUT_DIR, "/condout/blocks_done/block_", block, ".tsv"), row.names=F, col.names=F)
q <- function (save = "no", status = 0, runLast = TRUE)
.Internal(quit(save, status, runLast))
