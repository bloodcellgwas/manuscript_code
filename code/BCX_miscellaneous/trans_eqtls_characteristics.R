### BCX - revision 
## trans eQTLs

setwd("/Users/dv3/Desktop/gene_expression_analyses/")
trans=read.table("../trans_eQTL_in_blood_GWAS_sign_2018_11_14_extracted_290419.txt", he=T, strings=F)
mer=read.table("../finemapping_ukbb500k_final_release/2019_01_14_ukbb500k_condsig_with_LD_clumping_and_comparison_to_Astle.txt", he=T, strings=F, sep="\t", comment.char="")
head(mer)

uniq=mer[!duplicated(mer$Unique.Variant.ID),]
allt=scan("../trans_eqtl_list_of_tested_snps.txt", "char")
uniq$transqtl=NA
uniq$transqtl[which(uniq$rsID..where.available. %in% allt)]="no"
uniq$transqtl[which(uniq$rsID..where.available. %in% trans$SNP)]="yes"

table(uniq$transqtl)

# add closest genes for intergenic variants to VEP annotation # check the nearest TSS 
tss=read.table("/Users/dv3/Desktop/TSS.Ensemble_genes_downloaded_28_10_2019.txt", he=T, strings=F, sep="\t")
#J=which(uniq$Gene.Symbol.s..for.Most.Serious.Consequence=="")
#length(J)
for(j in 1:dim(uniq)[1]){
  if(uniq$Gene.Symbol.s..for.Most.Serious.Consequence[j]!="")
    next
  subs=tss[which(tss$Chromosome.scaffold.name==uniq$Chr..GRCh37.[j]),]
  k=which.min(abs(subs$Transcription.start.site..TSS.-uniq$BP..GRCh37.[j]))
  if(length(k)>1){
    warning(paste("check ",j))
  }
  while(is.na(subs$HGNC.symbol[k]) | subs$HGNC.symbol[k]==""){
    subs=subs[-k,]
    k=which.min(abs(subs$Transcription.start.site..TSS.-uniq$BP..GRCh37.[j]))
  }
  uniq$Gene.Symbol.s..for.Most.Serious.Consequence[j]=subs$HGNC.symbol[k]
}

mendelian=scan("../Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")
uniq$mendelian=NA
uniq$mendelian[which(uniq$Gene.Symbol.s..for.Most.Serious.Consequence!="")]="no"
tmp=sapply(uniq$Gene.Symbol.s..for.Most.Serious.Consequence, function(x){unlist(strsplit(x,","))})
k=unlist(lapply(tmp, function(x){length(which(x %in% mendelian))}))
uniq$mendelian[which(k>=1)]="yes"
#table(uniq$mendelian)

table(uniq$transqtl,uniq$mendelian, exclude=NULL)

hist(table(uniq$mendelian, uniq$transqtl, uniq$Most.Serious.VEP.Consequence.of.Variant))

table(uniq$transqtl[which(uniq$mendelian=="yes")])
table(uniq$transqtl[which(uniq$mendelian=="no")])

# trans-eQTL in Mendelian genes 48 variants (12%)
# non trans-eQTL in Mendelian 408 variants 
# trans-eQTL in peripheral 529, vs non trans in non mend 9738 (5.4%)


#chisq.test(table(uniq$transqtl, uniq$mendelian))
unique(uniq$Most.Serious.VEP.Consequence.of.Variant)->pio
subs=uniq[which(uniq$Most.Serious.VEP.Consequence.of.Variant %in% pio[c(1,2,3,7,9,10,11)]),]
res=data.frame(conseq=pio[c(1,2,3,7,9,10,11)], perc_trans_other=NA, perc_trans_mend=NA)
for(i in res$conseq){
  subs.s=subs[which(subs$Most.Serious.VEP.Consequence.of.Variant==i),]
  A=table(subs.s$transqtl, subs.s$mendelian)
  res[which(res$conseq==i),"perc_trans_other"]=A[2,1]/(A[1,1]+A[2,1])*100
  res[which(res$conseq==i),"perc_trans_mend"]=A[2,2]/(A[1,2]+A[2,2])*100
}

meltData=melt(res)
unique(meltData$conseq)
meltData$conseq=as.character(meltData$conseq)
meltData$conseq=gsub(meltData$conseq, pattern="3_prime_UTR_variant", replacement="3' UTR")
meltData$conseq=gsub(meltData$conseq, pattern="5_prime_UTR_variant", replacement="5' UTR")
meltData$conseq=gsub(meltData$conseq, pattern="downstream_gene_variant", replacement="downstream")
meltData$conseq=gsub(meltData$conseq, pattern="upstream_gene_variant", replacement="upstream")
meltData$conseq=gsub(meltData$conseq, pattern="intron_variant", replacement="intron")
meltData$conseq=gsub(meltData$conseq, pattern="variant", replacement="")
meltData$conseq=gsub(meltData$conseq, pattern="_", replacement="")

head(meltData)
meltData$conseq=factor(meltData$conseq)
gg=ggplot(meltData, aes(x=conseq, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")+
  theme_classic()+theme(axis.text.x = element_text(angle = 15, size=11), axis.text=element_text(size=11))+labs(x="VEP Consequence", y="%")+
  scale_fill_discrete(name="trans-eQTLs", labels=c("in peripheral genes","in Mendelian genes"))

pdf("Trans_eQTLs_in_mendelian_and_peripheral_genes.pdf", height=4, width=7)
print(gg)
dev.off()

eqtl=read.csv("finemapping_ukbb500k_final_release/variant_to_gene_reassignment/2019_07_09_ld_signif_eqtl_coloc_results.csv", he=T, strings=F)
head(eqtl)
dim(eqtl)

uniq$ciseqtl="no"
uniq$ciseqtl[which(uniq$Unique.Variant.ID %in% eqtl$condsig_var)]="yes"
table(uniq$ciseqtl,  uniq$mendelian)
table(uniq$transqtl[which(uniq$ciseqtl=="yes")])

## target genes enrichment
trans$mendel="no"
trans$mendel[which(trans$GeneSymbol %in% mendelian)]="yes"
table(trans$mendel)
dim(trans)

length(unique(trans$GeneSymbol[which(trans$mendel=="yes")]))
# 87 out of 135! (314) mendelian genes that we detect in the GWAS

length(intersect(trans$GeneSymbol[which(trans$mendel=="no")], uniq$Gene.Symbol.s..for.Most.Serious.Consequence))
# 1132 out of 4546 genes

### 
targets=unique(trans$GeneSymbol)
mend.gwas=intersect(unlist(strsplit(uniq$Gene.Symbol.s..for.Most.Serious.Consequence,",")), mendelian)
intersect(targets, mendelian)->mend.targ
geni=unique(unlist(strsplit(uniq$Gene.Symbol.s..for.Most.Serious.Consequence,",")))

#### check expression levels 

### dataset 1: Nath et al.  
expr=read.table("whole_blood_gene_expression_finnish_cohort_normalised.txt", he=T, strings=F)

library(RNOmni)
expr$z=rankNorm(expr$x)
plot(density(expr$z), ylim=c(0,0.5), main="whole blood expression")
points(density(expr$z[which(expr$Group.1 %in% mendelian)]), type="l", col="red")
points(density(expr$z[which(expr$Group.1 %in% geni)]), type="l", col="green")
points(density(expr$z[which(expr$Group.1 %in% targets)]), type="l", col="blue")
legend("topleft", c("all genes","Mendelian","GWAS","trans-eQTL\ntargets"), lty=1, col=c("black","red","green","blue"))

write.table(expr, "whole_blood_gene_expression_finnish_cohort_normalised.txt", row.names = F, quote=F)

expr$mendelian=NA
expr$mendelian[which(expr$Group.1 %in% mendelian)]=expr$z[which(expr$Group.1 %in% mendelian)]
expr$target=NA
expr$target[which(expr$Group.1 %in% targets)]=expr$z[which(expr$Group.1 %in% targets)]
expr$GWAS=NA
expr$GWAS[which(expr$Group.1 %in% geni)]=expr$z[which(expr$Group.1 %in% geni)]
#expr$mendelian=as.factor(expr$mendelian)
head(expr)
dim(expr)
library(reshape2)
library(ggplot2)
#names(expr)
#ggplot(expr, aes(x=Group.1, y=x))+geom_boxplot()+facet_wrap(~mendelian+target+GWAS)

meltData=melt(expr[,c(1,2,4,5,6)])
meltData=meltData[which(!is.na(meltData$value)),]
head(meltData[which(meltData$variable=="mendelian"),])
head(meltData)

table(meltData$variable)
ggplot(meltData, aes(x=variable, y=value))+geom_boxplot()

meltData2=melt(expr[,c(1,3,4,5,6)])
meltData2=meltData2[which(!is.na(meltData2$value)),]
ggplot(meltData2, aes(x=variable, y=value))+geom_boxplot()

expr$mend.tgt=NA
expr[which(!is.na(expr$target) & !is.na(expr$mendelian)),"mend.tgt"]=expr$z[which(!is.na(expr$target) & !is.na(expr$mendelian))]
expr$gwas.tgt=NA
expr[which(!is.na(expr$target) & !is.na(expr$GWAS)),"gwas.tgt"]=expr$z[which(!is.na(expr$target) & !is.na(expr$GWAS))]

meltData3=melt(expr[,c(1,3,4,5,6,7,8)])
meltData3=meltData3[which(!is.na(meltData3$value)),]
meltData3$variable=factor(meltData3$variable, levels=c("z","target","mendelian","GWAS","mend.tgt","gwas.tgt"))

gg=ggplot(meltData3, aes(x=variable, y=value))+geom_boxplot()+
  ylab("Normalized median gene expression")+theme_classic()+theme(axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_text(size=14),legend.position="none")

gg

pdf("median_expression_targets_mend_gwas_comparison.pdf", width=7, height=4)
print(gg)
dev.off()
getwd()

wilcox.test(expr$z[which(!is.na(expr$mend.tgt))], expr$z[which(!is.na(expr$gwas.tgt))])
wilcox.test(expr$z[which(!is.na(expr$mendelian))], expr$z[which(!is.na(expr$GWAS))])

#Wilcoxon rank sum test with continuity correction

#data:  expr$z[which(!is.na(expr$mend.tgt))] and expr$z[which(!is.na(expr$gwas.tgt))]
#W = 33310, p-value = 0.7638
#alternative hypothesis: true location shift is not equal to 0

# in the full eQTL gen:
chisq.test(matrix(c(109,(314-109), 1504, (4546-1504)),2,2))

#Pearson's Chi-squared test with Yates' continuity correction
#data:  matrix(c(109, (314 - 109), 1504, (4546 - 1504)), 2, 2)
#X-squared = 0.282, df = 1, p-value = 0.5954


### dataset 2: Urmo Vosa, eQTL gen data from one cohort (BIOS)
expr=read.table("BIOS_expression_summary_20200528.txt", he=T, strings=F)
head(expr)
tss=read.table("Desktop/TSS.Ensemble_genes_downloaded_28_10_2019.txt", he=T, strings=F, sep="\t")
tss=tss[!duplicated(tss$HGNC.symbol),]
head(tss)
dim(tss)
expr=merge(expr, tss[,c("Gene.stable.ID","HGNC.symbol")], by.x="ENSG", by.y="Gene.stable.ID", all.x=T)
head(expr)
dim(expr)

expr$mendelian=NA
expr$mendelian[which(expr$HGNC.symbol %in% mendelian)]=expr$median_exp[which(expr$HGNC.symbol %in% mendelian)]
expr$target=NA
expr$target[which(expr$HGNC.symbol %in% targets)]=expr$median_exp[which(expr$HGNC.symbol %in% targets)]
expr$GWAS=NA
expr$GWAS[which(expr$HGNC.symbol %in% geni)]=expr$median_exp[which(expr$HGNC.symbol %in% geni)]

## however the overall gene expression is not normally distributed - rankinverse normalize

library(reshape2)
library(ggplot2)

meltData=melt(expr[,c(3,10,11,12,13)])
meltData=meltData[which(!is.na(meltData$value)),]
head(meltData[which(meltData$variable=="mendelian"),])
head(meltData)

table(meltData$variable)
ggplot(meltData, aes(x=variable, y=value))+geom_boxplot()

#meltData2=melt(expr[,c(1,3,4,5,6)])
#meltData2=meltData2[which(!is.na(meltData2$value)),]
#ggplot(meltData2, aes(x=variable, y=value))+geom_boxplot()

expr$mend.tgt=NA
expr[which(!is.na(expr$target) & !is.na(expr$mendelian)),"mend.tgt"]=expr$median_exp[which(!is.na(expr$target) & !is.na(expr$mendelian))]
expr$gwas.tgt=NA
expr[which(!is.na(expr$target) & !is.na(expr$GWAS)),"gwas.tgt"]=expr$median_exp[which(!is.na(expr$target) & !is.na(expr$GWAS))]

meltData3=melt(expr[,c(3,10,11,12,13,14,15)])
meltData3=meltData3[which(!is.na(meltData3$value)),]
table(meltData3$variable)
meltData3$variable=factor(meltData3$variable, levels=c("median_exp","target","mendelian","GWAS","mend.tgt","gwas.tgt"))

gg=ggplot(meltData3, aes(x=variable, y=value))+geom_boxplot()+
  ylab("Median gene expression")+theme_classic()+theme(axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_text(size=14),legend.position="none")

gg

pdf("median_expression_targets_mend_gwas_comparison_BIOS_dataset.pdf", width=7, height=4)
print(gg)
dev.off()
getwd()

wilcox.test(expr$median_exp[which(!is.na(expr$mend.tgt))], expr$median_exp[which(!is.na(expr$gwas.tgt))])
wilcox.test(expr$median_exp[which(!is.na(expr$mendelian))], expr$median_exp[which(!is.na(expr$GWAS))])


# update

expr$z=rankNorm(expr$median_exp)

coex.g=unique(dimnames(coex)[[1]])
expr$coex=NA
expr$coex[which(expr$HGNC.symbol %in% coex.g)]=expr$z[which(expr$HGNC.symbol %in% coex.g)]
expr$mendelian=NA
expr$mendelian[which(expr$HGNC.symbol %in% mendelian)]=expr$z[which(expr$HGNC.symbol %in% mendelian)]
expr$target=NA
expr$target[which(expr$HGNC.symbol %in% targets)]=expr$z[which(expr$HGNC.symbol %in% targets)]
expr$GWAS=NA
expr$GWAS[which(expr$HGNC.symbol %in% geni)]=expr$z[which(expr$HGNC.symbol %in% geni)]
expr$GWAS_w=NA
expr$GWAS_w[which(expr$HGNC.symbol %in% gwas.w)]=expr$z[which(expr$HGNC.symbol %in% gwas.w)]


meltData=melt(expr[,c(10,11,12,13,14,15,16)])
meltData=meltData[which(!is.na(meltData$value)),]
head(meltData[which(meltData$variable=="mendelian"),])
head(meltData)

table(meltData$variable)
ggplot(meltData, aes(x=variable, y=value))+geom_violin()+geom_boxplot()
pairwise.wilcox.test(meltData$value, meltData$variable)

length(coex.g)
length(geni)

expr$strat=0
expr$strat[which(!is.na(expr$coex))]=1
m.in=expr[,c("strat","z","HGNC.symbol")]
m.in=na.omit(m.in)
m.out=matchit(strat~z, data=m.in)
m.data<-match.data(m.out)
dim(m.data)

random_uno=m.data$HGNC.symbol[which(m.data$strat==0)]
random_due=m.data$HGNC.symbol[which(m.data$strat==0)]
length(intersect(random_uno,random_due))

########################################################
trans$z=NA
tmp=merge(trans, expr[,c("HGNC.symbol","z")], by.x="GeneSymbol",by.y="HGNC.symbol", all.x=T)
dim(tmp)
head(tmp)
trans=tmp
trans$mendel[which(trans$mendel=="yes")]=1
trans$mendel[which(trans$mendel=="no")]=0
head(uniq)
tmp=merge(trans, uniq[,c("rsID..where.available.","Minor.Allele.Frequency")], by.x="SNP", by.y="rsID..where.available.", all.x=T)
head(tmp)
dim(tmp)
trans=tmp
table(trans$mendel)
m.in=trans[,c("SNP","GeneSymbol","Zscore","mendel","z","Minor.Allele.Frequency")]
m.in=na.omit(m.in)
m.out=matchit(mendel~Zscore+z+Minor.Allele.Frequency, data=m.in, ratio=20)
m.data=match.data(m.out)
dim(m.data)
head(m.data)
table(m.data$mendel)
A=data.frame(table(m.data$GeneSymbol))
head(A)
dim(A)
A$mend="no"
A$mend[which(A$Var1 %in% mendelian)]="yes"
aggregate(A$Freq, list(A$mend), summary)


# repeat on the farm 
setwd("/lustre/scratch115/realdata/mdt3/projects/ukbb500k_t151/coexpression/trans-eTQLs")
trans=read.table("/nfs/team151_data02/eQTL_gen/trans-eQTL_significant_20181017_with_var_ID.txt", he=T, strings=F)
mendelian=scan("BRIDGE_Flagship_gene_list.txt", "char")
mer=read.table("../2019_01_14_ukbb500k_condsig_with_LD_clumping_and_comparison_to_Astle.txt", he=T, strings=F, sep="\t",comment.char="")
gwas=unique(unlist(strsplit(mer$Gene.Symbol.s..for.Most.Serious.Consequence,",")))
trans$mendel="no"
trans$mendel[which(trans$GeneSymbol %in% mendelian)]="yes"

expr=read.table("BIOS_expression_summary_20200528_with_normalisation.txt", he=T, strings=F)
head(expr)

tmp=merge(trans, expr[,c("HGNC.symbol","z")], by.x="GeneSymbol",by.y="HGNC.symbol", all.x=T)
dim(tmp)
head(tmp)
trans=tmp
trans$mendel[which(trans$mendel=="yes")]=1
trans$mendel[which(trans$mendel=="no")]=0
uniq=mer[!duplicated(mer$Unique.Variant.ID),]
tmp=merge(trans, uniq[,c("rsID..where.available.","Minor.Allele.Frequency")], by.x="SNP", by.y="rsID..where.available.", all.x=T)
head(tmp)
dim(tmp)
trans=tmp
table(trans$mendel)

m.in=subs[,c("SNP","GeneSymbol","Zscore","mendel","z")]
m.in=na.omit(m.in)
m.out=matchit(mendel~Zscore+z, data=m.in, ratio=20)
m.data=match.data(m.out)
dim(m.data)
head(m.data)
table(m.data$mendel)
m.data$gwas=0
m.data$gwas[which(m.data$GeneSymbol %in% gwas)]=1


B=data.frame(table(m.data$GeneSymbol))
head(A)
dim(A)
B$gwas="no"
B$gwas[which(B$Var1 %in% gwas)]="yes"
aggregate(B$Freq, list(B$gwas), summary)
### this datset is calibrated for mendelian genes but not for gwas
subs=trans[which(trans$GeneSymbol %in% mendelian | trans$GeneSymbol %in% gwas),]

