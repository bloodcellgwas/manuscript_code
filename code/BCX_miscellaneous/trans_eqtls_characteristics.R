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


