### revision 
## trans eQTLs

setwd("/Users/dv3/Desktop/gene_expression_analyses/")
trans=read.table("../trans_eQTL_in_blood_GWAS_sign_2018_11_14_extracted_290419.txt", he=T, strings=F)
mer=read.table("../finemapping_ukbb500k_final_release/2019_01_14_ukbb500k_condsig_with_LD_clumping_and_comparison_to_Astle.txt", he=T, strings=F, sep="\t", comment.char="")
head(mer)

uniq=mer[!duplicated(mer$Unique.Variant.ID),]
allt=scan("trans_eqtl_list_of_tested_snps.txt", "char")
uniq$transqtl=NA
uniq$transqtl[which(uniq$rsID..where.available. %in% allt)]="no"
uniq$transqtl[which(uniq$rsID..where.available. %in% trans$SNP)]="yes"

table(uniq$transqtl)

mendelian=scan("../Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")
uniq$mendelian="no"
tmp=sapply(uniq$Gene.Symbol.s..for.Most.Serious.Consequence, function(x){unlist(strsplit(x,","))})
k=unlist(lapply(tmp, function(x){length(which(x %in% mendelian))}))
uniq$mendelian[which(k>=1)]="yes"
table(uniq$mendelian)
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


#### check expression levels 
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