###########################################
### CLINVAR annotation for GWAS results ###
###########################################

# overlap our findings with ClinVar db to see which pathogenic variants we detect 
# and if we can re-assign pathogenicity based on effect sizes distribution.


setwd("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/Results_final_fine_mapping/")
res=read.table("FM_vars_PP_gt_0.95.txt", he=T, strings=F)  
#res=read.table("FM_vars_PP_gt_0.50.txt", he=T, strings=F)
head(res)
dim(res)

# read and format clinvar downloaded dataset
clinvar=read.csv("../../Clinvar_variant_summary_051118.txt", he=T, strings=F, sep="\t", encoding="UTF-8", fill=T, comment.char = "")
head(clinvar)
clinvar=clinvar[which(clinvar$Assembly=="GRCh37"),]
dim(clinvar)  # 459080
clinvar$STAR=NA
clinvar$STAR[which(clinvar$ReviewStatus=="practice guideline")]=4
clinvar$STAR[which(clinvar$ReviewStatus=="reviewed by expert panel")]=3
clinvar$STAR[which(clinvar$ReviewStatus=="criteria provided, multiple submitters, no conflicts")]=2
clinvar$STAR[which(clinvar$ReviewStatus=="criteria provided, conflicting interpretations")]=1
clinvar$STAR[which(clinvar$ReviewStatus=="criteria provided, single submitter")]=1
clinvar$STAR[which(is.na(clinvar$STAR))]=0
table(clinvar$ClinSigSimple,clinvar$STAR)
clinvar$VAR=paste(clinvar$Chromosome,clinvar$Start,sep=":")
clinvar$VAR=paste(clinvar$VAR, clinvar$ReferenceAllele, clinvar$AlternateAllele, sep="_")
dim(clinvar)
clinv=clinvar[,c("ClinSigSimple","PhenotypeList","OriginSimple","VAR","STAR","ReviewStatus")]
dim(clinv)
head(clinv)
table(clinv$STAR)
tail(clinv)
length(unique(clinv$VAR))
clinv[which(duplicated(clinv$VAR)),] 

# merge with results file
head(res)
res$VAR=paste(res$chromosome,res$position,sep=":")
res$VAR=paste(res$VAR, res$allele1, res$allele2, sep="_")
fin=merge(res, clinv, by="VAR", all.x=T)
dim(fin)
head(fin)
length(unique(fin$VAR))
table(fin$ClinSigSimple, fin$STAR)

# plot effect sizes by clinvar STAR rating
library(ggplot2)
tmp=fin[which(fin$ClinSigSimple %in% c(0,1)),]
ggplot(data=tmp, aes(x=factor(STAR), y=abs(beta)))+
  geom_violin()+geom_boxplot(width = 0.2,aes(fill = "grey"))+scale_y_continuous(trans='log10')+
  theme(axis.text=element_text(size=14),axis.title.x=element_blank())+theme(legend.position="none")+facet_wrap(~ClinSigSimple)


### annotate 95% credible sets 
setwd("/Users/dv3/Desktop/finamapping_ukbb500k_final_release/Results_final_fine_mapping")
res=read.table("Credible_sets/All_traits.95credSet_summStat", he=F, strings=F)
hedd=scan("Credible_sets/header_snp_files.txt", "char")
head(res)
dim(res)
hedd=c("Trait", hedd)
names(res)=hedd
res$VAR=paste(res$chromosome,res$position,sep=":")
res$VAR=paste(res$VAR, res$allele1, res$allele2, sep="_")
length(unique(res$VAR))

fin=merge(res, clinv, by="VAR", all.x=T)
dim(fin)
head(fin)
length(unique(fin$VAR))
table(fin$ClinSigSimple, fin$STAR)

write.table(fin, "Credible_sets/All_traits.95credSet_summStat_clinvar_annotated.txt", row.names=F, quote=F, sep="\t")


