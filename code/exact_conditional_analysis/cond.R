#!/usr/bin/Rscript
library(RcppEigen)
library(data.table)
library(rhdf5)
library(doMC)
options(show.error.locations = TRUE)
OUT_DIR <- Sys.getenv('OUT_DIR')
args = commandArgs(TRUE)
pheno = args[1]
chr = args[2]
pheno = "mrv"; chr = 3
PHENO_FILE <- Sys.getenv('PHENO_FILE')
if (chr == '23') {
    int_phenotypes  = fread(paste0(OUT_DIR, "/genfiles/chrXY_merged.sample"), head=TRUE)
} else {
    int_phenotypes <-fread(PHENO_FILE, head=TRUE)
}
PHENO_NAMES <- Sys.getenv('PHENO_NAMES')
## load list of phenotypes being analysed
pheno_table=as.vector(read.table(file=PHENO_NAMES,stringsAsFactors=FALSE)$V1)
## loop through each phoenotype and postprocess conditional analysis results
print(pheno_table)
##  pheno <- pheno_table[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
print(pheno)
print(dim(int_phenotypes))
options(show.error.locations = TRUE)
CORES <- 1
OUT_DIR <- Sys.getenv('OUT_DIR')
PHENO_NAMES <- Sys.getenv('PHENO_NAMES')
##  load list of phenotypes
pheno_table=as.vector(read.table(file=PHENO_NAMES,stringsAsFactors=FALSE)$V1)
## extract a list of phenotype values which match the phenotype we are currently studying, convert this to numeric
## also ignore the first row of the phenotype file which contains data type information
## convert names to have _gwas_normalised
pheno_table <- paste(pheno_table, "_gwas_normalised", sep="")
oldpheno <- pheno
pheno <- paste(pheno, "_gwas_normalised", sep="")
merged_phenotypes=int_phenotypes[2:dim(int_phenotypes)[1], match(pheno_table, colnames(int_phenotypes)),with=FALSE]
merged_phenotypes=apply(merged_phenotypes, 2, as.numeric)
## extract covariates from the .sample file
int_numerical_covariates=int_phenotypes[2:dim(int_phenotypes)[1],4:13,with=FALSE]	
setnames(int_numerical_covariates,as.character(1:10))
  
## total size of dataset
int_size=dim(int_phenotypes)[1]-1
  
## merge covariates with other cohorts (we only have one) and set them to numeric
merged_numerical_covs=int_numerical_covariates
merged_numerical_covs=apply(merged_numerical_covs, 2, as.numeric)		
## which phenotypes in .sample file are we studying
selected_pheno=colnames(merged_phenotypes)==pheno 
## vector of samples which aren't missing any of the selected phenotype values
non_missing=!is.na(merged_phenotypes[,selected_pheno])
## adjust for covariates
## interaction term for PGs
adjusted_lm<-lm(merged_phenotypes[non_missing,selected_pheno]~merged_numerical_covs[non_missing, 1:10])
merged_phenotypes[non_missing,selected_pheno]=resid(adjusted_lm,na.action=na.exclude)

	
trait_variants<-read.table(sprintf("%s/condout/trait/var_list_%s.tsv", OUT_DIR, oldpheno), sep="\t", header=T, stringsAsFactors=F)
## we only want trait_variants which are on the chromosome we are currently studying
trait_variants$chr = gsub("(.+):.+_.+_.+","\\1", trait_variants$VARIANT)
trait_variants$chr[trait_variants$chr == "X"] = "23"
trait_variants$chr[trait_variants$chr == "XY"] = "23"
trait_variants = trait_variants[trait_variants$chr == chr,c("VARIANT")]
if (length(trait_variants) == 0) {
    coef_table_labeled = data.frame(matrix(ncol = 4, nrow = 0))
    colnames(coef_table_labeled)=c("VARIANT", "JOINT_BETA", "JOINT_SE", "JOINT_MLOG10P")
    write.table(coef_table_labeled, sprintf("%s/condout/results_chrom/condsig_%s_chr_%s.tsv", OUT_DIR, pheno, chr), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    stop("done")
}
variant_ids=NULL
print(chr)
this_chr.variant_ids <- fread(sprintf("%s/condout/subgenhd5/ids_pos_%s.tsv",OUT_DIR,chr), head=FALSE, sep="\t", colClasses=c("character"))
this_chr.pull_ids <- fread(paste0(OUT_DIR, "/condout/pull_ids/pull_ids_chr", chr, ".tsv"), colClasses=c("character"))
this_chr.pull_ids$VARIANT <- paste0(as.numeric(this_chr.pull_ids$chromosome), ":", this_chr.pull_ids$position, "_", this_chr.pull_ids$alleleA, "_", this_chr.pull_ids$alleleB)
this_chr.variant_ids$V2 = gsub("XY", "23", this_chr.variant_ids$V1)
this_chr.variant_ids$V2 = gsub("X", "23", this_chr.variant_ids$V2)
#this_chr.pull_ids = this_chr.pull_ids[this_chr.pull_ids$chromosome == chr,]
if (sum(this_chr.pull_ids$VARIANT %in% this_chr.variant_ids$V2) != nrow(this_chr.pull_ids)) {
    stop(paste(chr, "pull id doesn't match with gen variant ids"))
}
variant_ids<-c(variant_ids,unlist(this_chr.variant_ids$V1))

names(variant_ids)=NULL
if (sum(!is.element(trait_variants, variant_ids)) > 0) {
#    save(list=c("trait_variants", "variant_ids"), file=paste0(pheno, "_pull_check.Rdata"))
    stop(paste(pheno, "there are variants to analyse which aren't in the pull table"))
}

registerDoMC(cores=CORES)
read_add_data<-function(chr, trait_variants)
{
    print(paste(chr, "loaded"))
    add_genotypes<-h5read(sprintf("%s/condout/subgenhd5/chr%s_pos_clean_id.h5", OUT_DIR ,chr), "/add")
    add_genotypes<-add_genotypes[,-1]
    this_chr.variant_ids <- fread(sprintf("%s/condout/subgenhd5/ids_pos_%s.tsv",OUT_DIR,chr), head=FALSE, sep="\t", colClasses=c("character"))$V1
    colnames(add_genotypes) <- this_chr.variant_ids
    add_genotypes <- add_genotypes[, is.element(this_chr.variant_ids, trait_variants), drop=FALSE]
    print(paste(chr, nrow(add_genotypes)))
    return(add_genotypes)
}
read_dom_data<-function(chr)
{
    dom_genotypes<-h5read(sprintf("%s/condout/subgenhd5/chr%s.h5", OUT_DIR ,chr), "/dom")
    return(dom_genotypes)
}

gc()
#add_genotypes=foreach(n=1:23, .combine="cbind") %dopar% read_add_data(n, trait_variants)
#if (length(add_genotypes) < 23) {stop("add genotypes length not right")}
print(chr)
add_genotypes = read_add_data(chr, trait_variants)
gc()
variant_order <- colnames(add_genotypes)
if (sum(is.element(trait_variants, variant_order)) != length(trait_variants) ||
    length(trait_variants) != length(variant_order)) {
    save(list=c("trait_variants", "variant_order"), file=paste0(pheno, "_pull_check.Rdata"))
    paste(stop("not all trait variants are included in the pull table"))
}
colnames(add_genotypes) <- NULL
model_genotypes=cbind(rep(1, sum(non_missing)),add_genotypes[non_missing,])
model_phenotype=as.matrix(merged_phenotypes[non_missing,selected_pheno, drop=FALSE])
nvar=dim(model_genotypes)[2]-1

compute_add_pvalues<-function(add_var)
{
    proposed_variants=current_variants
    proposed_variants[add_var]=TRUE
    ## log p value
    logpvalue=0
    if(!all(proposed_variants==current_variants))
    {
        test_index=sum(proposed_variants[1:add_var])
        proposed_model<-fastLm(model_genotypes[,c(TRUE, proposed_variants)], model_phenotype)
        logpvalue=log(2)+pt(abs(proposed_model$coefficients[1+test_index]/proposed_model$se[1+test_index]), df=proposed_model$df.residual, log=TRUE,lower=FALSE)	
        
    }
    return(logpvalue)
}

looping=TRUE
adding=FALSE
current_variants=rep(TRUE, nvar)	

while(looping)	
{
    print(looping)
    adding_failed=FALSE
    print(sprintf("number of variants = %d / %d", sum(current_variants), length(current_variants)))
    if(adding)
    {
        logadd_pvalues=unlist(foreach(n=1:nvar, .combine="c") %dopar% compute_add_pvalues(n))
        if (length(logadd_pvalues) != nvar) {
            stop("the lengths are not equal, add more memory")
        }
        ## change to log p
        logminp=min(logadd_pvalues)
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
        
        summary(logmaxp)
        summary(log(8.31)-9*log(1))
        if(logmaxp>(log(8.31)-9*log(10)))
        {
            ## use which .max
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
        if(adding_failed){looping=FALSE}
    }    
}


	
check_model<-fastLm(model_genotypes[,c(TRUE, current_variants)], model_phenotype)

coef_table=coef(summary(check_model))
coef_table[,4]=as.character(-(log(2)+pt(abs(check_model$coefficients/check_model$se), df=check_model$df.residual, log=TRUE,lower=FALSE))/log(10))
dddff
coef_table_labeled=cbind(variant_order[current_variants], coef_table[2:dim(coef_table)[1],c(1:2,4)])

if (sum(current_variants) == 1) {
    coef_table_labeled = data.frame(VARIANT=variant_order[current_variants], JOINT_BETA=coef_table[2,1], JOINT_SE=coef_table[2,2], JOINT_MLOG10P=coef_table[2,4])
    rownames(coef_table_labeled) = NULL
} else if (sum(current_variants) == 0) {
    coef_table_labeled = data.frame(matrix(ncol = 4, nrow = 0))
    colnames(coef_table_labeled)=c("VARIANT", "JOINT_BETA", "JOINT_SE", "JOINT_MLOG10P")
} else {
    colnames(coef_table_labeled)=c("VARIANT", "JOINT_BETA", "JOINT_SE", "JOINT_MLOG10P")
}
print(summary(coef_table_labeled))

write.table(coef_table_labeled, sprintf("%s/condout/results_chrom/condsig_%s_chr_%s.tsv", OUT_DIR, pheno, chr), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
print("JOB COMPLETE")

