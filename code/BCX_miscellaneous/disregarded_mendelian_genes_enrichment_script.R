### Mendelian genes enrichment wrap-up script across multiple unrelated traits

setwd("/Users/dv3/Desktop/Mendelian_gene_enrichments_across_traits/")

filez=system("ls *_10_4_20.csv", intern=T)
traits=unlist(sapply(filez, function(x){unlist(strsplit(x,"_"))[1]}))
mendelian=scan("../Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")

# create a data frame to populate with observed characteristics
fin=data.frame(trait=traits, total_unique_SNPs=NA, mendelian_common_vars=NA, mendelian_rare_vars=NA, N_mendelian_common_genes=NA, N_mendelian_rare_genes=NA,mendelian_common_genes=NA,mendelian_rare_genes=NA) )

# for each trait overlap with mendelian genes and save the observed characteristics
for(i in 1:length(traits)){
res=read.csv(filez[i], he=T, strings=F)
res=res[!duplicated(res$Variant.and.risk.allele),]
res$mendelian="no"
tmp=sapply(res$Mapped.gene, function(x){unlist(strsplit(x,","))})
k=unlist(lapply(tmp, function(x){length(which(x %in% mendelian))}))
res$mendelian[which(k>=1)]="yes"
res$RAF=as.numeric(res$RAF)
res$maf=res$RAF
res$maf[which(res$RAF>0.5)]=1-res$RAF[which(res$RAF>0.5)]
res$freq_class="common"
res$freq_class[which(res$maf<0.01)]="rare"

fin[i,"total_unique_SNPs"]=dim(res)[1]
fin[i,"mendelian_common_vars"]=length(which(res$mendelian=="yes" & res$freq_class=="common"))
fin[i,"mendelian_rare_vars"]=length(which(res$mendelian=="yes" & res$freq_class=="rare"))
fin[i,"N_mendelian_common_genes"]=length(intersect(unique(unlist(strsplit(res$Mapped.gene[which(res$mendelian=="yes" & res$freq_class=="common")],","))), mendelian))
fin[i,"N_mendelian_rare_genes"]=length(intersect(unique(unlist(strsplit(res$Mapped.gene[which(res$mendelian=="yes" & res$freq_class=="rare")],","))), mendelian))
fin[i,"mendelian_common_genes"]=paste(intersect(unique(unlist(strsplit(res$Mapped.gene[which(res$mendelian=="yes" & res$freq_class=="common")],","))), mendelian), collapse=",")
fin[i,"mendelian_rare_genes"]=paste(intersect(unique(unlist(strsplit(res$Mapped.gene[which(res$mendelian=="yes" & res$freq_class=="rare")],","))), mendelian), collapse=",")

write.table(res[,c(8,10,13,14,15)], gsub(filez[i], pattern="_10_4_20.csv", replacement="_10_4_20_curated.txt"), row.names = F, quote=F, sep="\t")
}

write.table(fin, "Info_trait_overlaps_with_mendelian_genes.txt", row.names = F, quote=F, sep="\t")


## random permutation over all proteing coding genes
tss=read.table("../TSS.onlychr1.22.XY.MT_EnsembleGenes.v99_downloaded.13_04_20.txt", he=T, strings=F, sep="\t",fill=T)
head(tss)
tss=tss[!duplicated(tss$HGNC.symbol),]
dim(tss)
## 19,239 protein coding genes

# gene names in the files can have multiple genes separated by a comma
for(j in 1:10){
  res=read.csv(filez[j], he=T, strings=F,sep="\t", comment.char="")
  assign(paste("res",j,sep="_"),res)
  assign(paste("tmp",j,sep="_"), sapply(res$Mapped.gene, function(x){unlist(strsplit(x,","))}))
}

set.seed(123)

# generate data set to populate with permutation results
fin=data.frame(iteration=1:10000)
for(j in 1:10){
  x=rep(NA, 10000)
  fin=cbind(fin, x, x, x, x)
  names(fin)[dim(fin)[2]-3]=paste(traits[j],"_comm.var",sep="")
  names(fin)[dim(fin)[2]-2]=paste(traits[j],"_rare.var",sep="")
  names(fin)[dim(fin)[2]-1]=paste(traits[j],"_comm.gene",sep="")
  names(fin)[dim(fin)[2]]=paste(traits[j],"_rare.gene",sep="")
}

# permutation 
for(i in 1:10000){
  random=sample(tss$HGNC.symbol, 314)
    for(j in 1:10){
  res=get(paste("res",j,sep="_"))
  tmp=get(paste("tmp",j,sep="_"))
  k=unlist(lapply(tmp, function(x){length(which(x %in% random))}))
  fin[i,paste(traits[j],"_comm.var",sep="")]=dim(res[which(res$freq_class=="common" & k>0),])[1]
  fin[i,paste(traits[j],"_rare.var",sep="")]=dim(res[which(res$freq_class=="rare" & k>0),])[1]
  fin[i,paste(traits[j],"_comm.gene",sep="")]=length(intersect(unique(unlist(strsplit(res$Mapped.gene[which(res$freq_class=="common")],","))), random))
  fin[i,paste(traits[j],"_rare.gene",sep="")]=length(intersect(unique(unlist(strsplit(res$Mapped.gene[which(res$freq_class=="rare")],","))), random))
  }
}


### generate enrichment results and attach to the dataset generated above 
# these are statstics from randomly generated permutations

pio=read.table("Info_trait_overlaps_with_mendelian_genes.txt", he=T, strings=F, sep="\t")
head(pio)
names(pio)
pio$comm_vars_mean=NA
pio$comm_vars_sd=NA
pio$c.v.enrichm=NA
pio$c.v.p_value=NA

pio$rare_vars_mean=NA
pio$rare_vars_sd=NA
pio$r.v.enrichm=NA
pio$r.v.p_value=NA

pio$comm_gen_mean=NA
pio$comm_gen_sd=NA
pio$c.g.enrichm=NA
pio$c.g.p_value=NA

pio$rare_gen_mean=NA
pio$rare_gen_sd=NA
pio$r.g.enrichm=NA
pio$r.g.p_value=NA

for(j in 1:10){
  tmp=fin[,grep(traits[j], names(fin))]
  pio$comm_vars_mean[j]=mean(tmp[,1])
  pio$comm_vars_sd[j]=sd(tmp[,1])
  pio$c.v.enrichm[j]=pio$mendelian_common_vars[j]/mean(tmp[,1])
  pio$c.v.p_value[j]=length(which(tmp[,1]>pio$mendelian_common_vars[j]))/10000
  
  pio$rare_vars_mean[j]=mean(tmp[,2])
  pio$rare_vars_sd[j]=sd(tmp[,2])
  pio$r.v.enrichm[j]=pio$mendelian_rare_vars[j]/mean(tmp[,2])
  pio$r.v.p_value[j]=length(which(tmp[,2]>pio$mendelian_rare_vars[j]))/10000
  
  pio$comm_gen_mean[j]=mean(tmp[,3])
  pio$comm_gen_sd[j]=sd(tmp[,3])
  pio$c.g.enrichm[j]=pio$N_mendelian_common_genes[j]/mean(tmp[,3])
  pio$c.g.p_value[j]=length(which(tmp[,3]>pio$N_mendelian_common_genes[j]))/10000
  
  pio$rare_gen_mean[j]=mean(tmp[,4])
  pio$rare_gen_sd[j]=sd(tmp[,4])
  pio$r.g.enrichm[j]=pio$N_mendelian_rare_genes[j]/mean(tmp[,4])
  pio$r.g.p_value[j]=length(which(tmp[,4]>pio$N_mendelian_rare_genes[j]))/10000
  
}

write.table(pio,"Table_enrichment_results_all_gwas_traits_prot_coding_genes_13_04_20.txt", row.names = F, quote=F, sep="\t")

### enrichment for blood traits 
### original script for blood traits

res=read.table("../finemapping_ukbb500k_final_release/2019_01_14_ukbb500k_condsig_with_LD_clumping_and_comparison_to_Astle.txt", he=T, strings=F, sep="\t", comment.char="")
tmp=sapply(res$Gene.Symbol.s..for.Most.Serious.Consequence, function(x){unlist(strsplit(x,","))})

res$freq_class="common"
res$freq_class[which(res$Minor.Allele.Frequency<0.01)]="rare"

set.seed(123)
blood_comm.var=c()
blood_rare.var=c()
blood_comm.gene=c()
blood_rare.gene=c()
for(i in 1:10000){
  random=sample(tss$HGNC.symbol, 314)
    k=unlist(lapply(tmp, function(x){length(which(x %in% random))}))
    blood_comm.var=c(blood_comm.var,dim(res[which(res$freq_class=="common" & k>0),])[1])
    blood_rare.var=c(blood_rare.var,dim(res[which(res$freq_class=="rare" & k>0),])[1])
    blood_comm.gene=c(blood_comm.gene,length(intersect(unique(unlist(strsplit(res$Gene.Symbol.s..for.Most.Serious.Consequence[which(res$freq_class=="common")],","))), random)))
    blood_rare.gene=c(blood_rare.gene,length(intersect(unique(unlist(strsplit(res$Gene.Symbol.s..for.Most.Serious.Consequence[which(res$freq_class=="rare")],","))), random)))
}


k=unlist(lapply(tmp, function(x){length(which(x %in% mendelian))}))
mendelian_common_vars=dim(res[which(res$freq_class=="common" & k>0),])[1]
mendelian_rare_vars=dim(res[which(res$freq_class=="rare" & k>0),])[1]
N_mendelian_common_genes=length(intersect(unique(unlist(strsplit(res$Gene.Symbol.s..for.Most.Serious.Consequence[which(res$freq_class=="common")],","))), mendelian))
N_mendelian_rare_genes=length(intersect(unique(unlist(strsplit(res$Gene.Symbol.s..for.Most.Serious.Consequence[which(res$freq_class=="rare")],","))), mendelian))

comm_vars_med=median(blood_comm.var)
comm_vars_Ulim=quantile(blood_comm.var,0.975)
comm_vars_Llim=quantile(blood_comm.var,0.0275)
c.v.enrichm=mendelian_common_vars/median(blood_comm.var)
c.v.obs_perc=length(which(blood_comm.var>mendelian_common_vars))/10000

rare_vars_med=median(blood_rare.var)
rare_vars_Ulim=quantile(blood_rare.var,0.975)
rare_vars_Llim=quantile(blood_rare.var,0.0275)
r.v.enrichm=mendelian_rare_vars/median(blood_rare.var)
r.v.obs_perc=length(which(blood_rare.var>mendelian_rare_vars))/10000

comm_gen_med=median(blood_comm.gene)
comm_gen_Ulim=quantile(blood_comm.gene,0.975)
comm_gen_Llim=quantile(blood_comm.gene,0.0275)
c.g.enrichm=N_mendelian_common_genes/median(blood_comm.gene)
c.g.obs_perc=length(which(blood_comm.gene>N_mendelian_common_genes))/10000

rare_gen_med=median(blood_rare.gene)
rare_gen_Ulim=quantile(blood_rare.gene,0.975)
rare_gen_Llim=quantile(blood_rare.gene,0.0275)
r.g.enrichm=N_mendelian_rare_genes/median(blood_rare.gene)
r.g.obs_perc=length(which(blood_rare.gene>N_mendelian_rare_genes))/10000


### plot enrichment results 

pio=read.table("Table_enrichment_results_all_gwas_traits_prot_coding_genes_14_04_20.txt", he=T, strings=F, sep="\t")
pio$trait

## add confidence intervals to enrichment, consider the 95% percentile of the random distribution as upper limit
## if I want the CI to reflect when the enrichment is significant, then I should consider 95%-50% centile
## since this is only a one-sided test I am ignoring the left side (although I should consider it??)

fin=read.table("random_perms_results_all_traits.txt", he=T, strings=F)
prova=apply(fin[,grep("_comm.var", names(fin))], 2, function(x){quantile(x, 0.975)})
pio$c.v.Ulim=c(prova, comm_vars_Ulim)
prova=apply(fin[,grep("_comm.var", names(fin))], 2, function(x){quantile(x, 0.0275)})
pio$c.v.Llim=c(prova, comm_vars_Llim)
prova=apply(fin[,grep("_comm.var", names(fin))], 2, median)
pio$c.v.med=c(prova, comm_vars_med)

prova=apply(fin[,grep("_rare.var", names(fin))], 2, function(x){quantile(x, 0.975)})
pio$r.v.Ulim=c(prova, rare_vars_Ulim)
prova=apply(fin[,grep("_rare.var", names(fin))], 2, function(x){quantile(x, 0.0275)})
pio$r.v.Llim=c(prova, rare_vars_Llim)
prova=apply(fin[,grep("_rare.var", names(fin))], 2, median)
pio$r.v.med=c(prova, rare_vars_med)

prova=apply(fin[,grep("_comm.gen", names(fin))], 2, function(x){quantile(x, 0.975)})
pio$c.g.Ulim=c(prova, comm_gen_Ulim)
prova=apply(fin[,grep("_comm.gen", names(fin))], 2, function(x){quantile(x, 0.0275)})
pio$c.g.Llim=c(prova, comm_gen_Llim)
prova=apply(fin[,grep("_comm.gen", names(fin))], 2, median)
pio$c.g.med=c(prova, comm_gen_med)

prova=apply(fin[,grep("_rare.gene", names(fin))], 2, function(x){quantile(x, 0.975)})
pio$r.g.Ulim=c(prova, rare_gen_Ulim)
prova=apply(fin[,grep("_rare.gene", names(fin))], 2, function(x){quantile(x, 0.0275)})
pio$r.g.Llim=c(prova, rare_gen_Llim)
prova=apply(fin[,grep("_rare.gene", names(fin))], 2, median)
pio$r.g.med=c(prova, rare_gen_med)

# now that I have the basic stats I am calcualting the enrichment and the CI, separately for the two tails 
pio$c.v.enrichm=pio$mendelian_common_vars/pio$c.v.med
pio$c.v.FE_U[which(pio$c.v.enrichm>1)]=pio$c.v.enrichm[which(pio$c.v.enrichm>1)]+((pio$c.v.Ulim[which(pio$c.v.enrichm>1)]-pio$c.v.med[which(pio$c.v.enrichm>1)])/pio$c.v.med[which(pio$c.v.enrichm>1)])
pio$c.v.FE_L[which(pio$c.v.enrichm>1)]=pio$c.v.enrichm[which(pio$c.v.enrichm>1)]-((pio$c.v.Ulim[which(pio$c.v.enrichm>1)]-pio$c.v.med[which(pio$c.v.enrichm>1)])/pio$c.v.med[which(pio$c.v.enrichm>1)])
pio$c.v.FE_L[which(pio$c.v.enrichm<1)]=pio$c.v.enrichm[which(pio$c.v.enrichm<1)]+((pio$c.v.Llim[which(pio$c.v.enrichm<1)]-pio$c.v.med[which(pio$c.v.enrichm<1)])/pio$c.v.med[which(pio$c.v.enrichm<1)])
pio$c.v.FE_U[which(pio$c.v.enrichm<1)]=pio$c.v.enrichm[which(pio$c.v.enrichm<1)]-((pio$c.v.Llim[which(pio$c.v.enrichm<1)]-pio$c.v.med[which(pio$c.v.enrichm<1)])/pio$c.v.med[which(pio$c.v.enrichm<1)])
pio$c.v.p_value[which(pio$c.v.enrichm>1)]=as.numeric(pio$c.v.p_value[which(pio$c.v.enrichm>1)])*2
pio$c.v.p_value[which(pio$c.v.enrichm<1)]=(1-as.numeric(pio$c.v.p_value[which(pio$c.v.enrichm<1)]))*2

pio$r.v.enrichm=pio$mendelian_rare_vars/pio$r.v.med
pio$r.v.FE_U[which(pio$r.v.enrichm>1)]=pio$r.v.enrichm[which(pio$r.v.enrichm>1)]+((pio$r.v.Ulim[which(pio$r.v.enrichm>1)]-pio$r.v.med[which(pio$r.v.enrichm>1)])/pio$r.v.med[which(pio$r.v.enrichm>1)])
pio$r.v.FE_L[which(pio$r.v.enrichm>1)]=pio$r.v.enrichm[which(pio$r.v.enrichm>1)]-((pio$r.v.Ulim[which(pio$r.v.enrichm>1)]-pio$r.v.med[which(pio$r.v.enrichm>1)])/pio$r.v.med[which(pio$r.v.enrichm>1)])
pio$r.v.FE_L[which(pio$r.v.enrichm<1)]=pio$r.v.enrichm[which(pio$r.v.enrichm<1)]+((pio$r.v.Llim[which(pio$r.v.enrichm<1)]-pio$r.v.med[which(pio$r.v.enrichm<1)])/pio$r.v.med[which(pio$r.v.enrichm<1)])
pio$r.v.FE_U[which(pio$r.v.enrichm<1)]=pio$r.v.enrichm[which(pio$r.v.enrichm<1)]-((pio$r.v.Llim[which(pio$r.v.enrichm<1)]-pio$r.v.med[which(pio$r.v.enrichm<1)])/pio$r.v.med[which(pio$r.v.enrichm<1)])
pio$r.v.p_value[which(pio$r.v.enrichm>1)]=as.numeric(pio$r.v.p_value[which(pio$r.v.enrichm>1)])*2
pio$r.v.p_value[which(pio$r.v.enrichm<1)]=(1-as.numeric(pio$r.v.p_value[which(pio$r.v.enrichm<1)]))*2

pio$c.g.enrichm=pio$N_mendelian_common_genes/pio$c.g.med
pio$c.g.FE_U[which(pio$c.g.enrichm>1)]=pio$c.g.enrichm[which(pio$c.g.enrichm>1)]+((pio$c.g.Ulim[which(pio$c.g.enrichm>1)]-pio$c.g.med[which(pio$c.g.enrichm>1)])/pio$c.g.med[which(pio$c.g.enrichm>1)])
pio$c.g.FE_L[which(pio$c.g.enrichm>1)]=pio$c.g.enrichm[which(pio$c.g.enrichm>1)]-((pio$c.g.Ulim[which(pio$c.g.enrichm>1)]-pio$c.g.med[which(pio$c.g.enrichm>1)])/pio$c.g.med[which(pio$c.g.enrichm>1)])
pio$c.g.FE_L[which(pio$c.g.enrichm<1)]=pio$c.g.enrichm[which(pio$c.g.enrichm<1)]+((pio$c.g.Llim[which(pio$c.g.enrichm<1)]-pio$c.g.med[which(pio$c.g.enrichm<1)])/pio$c.g.med[which(pio$c.g.enrichm<1)])
pio$c.g.FE_U[which(pio$c.g.enrichm<1)]=pio$c.g.enrichm[which(pio$c.g.enrichm<1)]-((pio$c.g.Llim[which(pio$c.g.enrichm<1)]-pio$c.g.med[which(pio$c.g.enrichm<1)])/pio$c.g.med[which(pio$c.g.enrichm<1)])
pio$c.g.p_value[which(pio$c.g.enrichm>1)]=as.numeric(pio$c.g.p_value[which(pio$c.g.enrichm>1)])*2
pio$c.g.p_value[which(pio$c.g.enrichm<1)]=(1-as.numeric(pio$c.g.p_value[which(pio$c.g.enrichm<1)]))*2

pio$r.g.enrichm=pio$N_mendelian_rare_genes/pio$r.g.med
pio$r.g.FE_U[which(pio$r.g.enrichm>1)]=pio$r.g.enrichm[which(pio$r.g.enrichm>1)]+((pio$r.g.Ulim[which(pio$r.g.enrichm>1)]-pio$r.g.med[which(pio$r.g.enrichm>1)])/pio$r.g.med[which(pio$r.g.enrichm>1)])
pio$r.g.FE_L[which(pio$r.g.enrichm>1)]=pio$r.g.enrichm[which(pio$r.g.enrichm>1)]-((pio$r.g.Ulim[which(pio$r.g.enrichm>1)]-pio$r.g.med[which(pio$r.g.enrichm>1)])/pio$r.g.med[which(pio$r.g.enrichm>1)])
pio$r.g.FE_L[which(pio$r.g.enrichm<1)]=pio$r.g.enrichm[which(pio$r.g.enrichm<1)]+((pio$r.g.Llim[which(pio$r.g.enrichm<1)]-pio$r.g.med[which(pio$r.g.enrichm<1)])/pio$r.g.med[which(pio$r.g.enrichm<1)])
pio$r.g.FE_U[which(pio$r.g.enrichm<1)]=pio$r.g.enrichm[which(pio$r.g.enrichm<1)]-((pio$r.g.Llim[which(pio$r.g.enrichm<1)]-pio$r.g.med[which(pio$r.g.enrichm<1)])/pio$r.g.med[which(pio$r.g.enrichm<1)])
pio$r.g.p_value[which(pio$r.g.enrichm>1)]=as.numeric(pio$r.g.p_value[which(pio$r.g.enrichm>1)])*2
pio$r.g.p_value[which(pio$r.g.enrichm<1)]=(1-as.numeric(pio$r.g.p_value[which(pio$r.g.enrichm<1)]))*2

library(ggplot2)
require(gridExtra)

tmp=pio[,c(1, grep("c.v.", names(pio)))]
tmp=tmp[order(tmp$c.v.enrichm, decreasing = T),]
#head(tmp)
tmp$trait=factor(tmp$trait, levels=tmp$trait[order(tmp$c.v.enrichm)])
plot1=ggplot(tmp, aes(x=trait, y=c.v.enrichm))+geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes( ymin=c.v.FE_L, ymax=c.v.FE_U), position=position_dodge())+
  coord_flip() +theme_classic()+scale_fill_manual(values=c("blueviolet","#313fbd","#3182bd"))+
  geom_hline(yintercept = 1, col="black", linetype="dashed")+
  theme(axis.title=element_text(size=12),axis.text=element_text(size=12))+ylab("common variants enrichment")


tmp=pio[,c(1,7,grep("c.g.", names(pio)))]
tmp=tmp[order(tmp$c.g.enrichm, decreasing = T),]
#head(tmp)
tmp$trait=factor(tmp$trait, levels=tmp$trait[order(tmp$c.g.enrichm)])
plot2=ggplot(tmp, aes(x=trait, y=c.g.enrichm))+geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes( ymin=c.g.FE_L, ymax=c.g.FE_U), position=position_dodge())+
  coord_flip() +theme_classic()+scale_fill_manual(values=c("blueviolet","#313fbd","#3182bd"))+
  geom_hline(yintercept = 1, col="black", linetype="dashed") +
  theme(axis.text=element_text(size=14), )+ylab("common genes enrichment")


getwd()
pdf("Mend_genes_enrichment_in_other_GWAS.pdf", width=7,height=4)
grid.arrange(plot1, plot2, ncol=2)
dev.off()
ls()


## plot the blood-related enrichment for the other question about common and rare variants
pdf("Enrichment_by_common_and_rare_variants_and_genes.pdf", width=8, height=8)
par(mfrow=c(2,2))
print(plot(density(blood_comm.var), xlim=c(0,800), main="common",ylab="",xlab="", yaxt="n", bty="n"))
print(abline(v=mendelian_common_vars))
print(plot(density(blood_rare.var), xlim=c(0,187), main="rare",ylab="",xlab="", yaxt="n", bty="n"))
print(abline(v=mendelian_rare_vars))
print(plot(density(blood_comm.gene), xlim=c(20,110), main="common",ylab="",xlab="", yaxt="n", bty="n"))
print(abline(v=N_mendelian_common_genes))
print(hist(blood_rare.gene, xlim=c(0,40), main="rare",ylab="",xlab="", yaxt="n", bty="n"))
print(abline(v=N_mendelian_rare_genes))
dev.off()



## add total number of genes identified

fin$tot_genes=NA
for(j in 1:10){
  res=read.csv(filez[j], he=T, strings=F, comment.char="")
  tot=length(unique(unlist(strsplit(res$Mapped.gene,","))))
  fin[which(fin$trait==traits[j]),"tot_genes"]=tot
}

fin$chisq_p_genes=NA
for(j in 1:10){
 T=chisq.test(matrix(c(fin[j,"N_mendelian_common_genes"]+fin[j,"N_mendelian_rare_genes"], 314-fin[j,"N_mendelian_common_genes"]+fin[j,"N_mendelian_rare_genes"], fin[j,"tot_genes"]-fin[j,"N_mendelian_common_genes"]+fin[j,"N_mendelian_rare_genes"], 19239-fin[j,"tot_genes"]),2,2)) 
 fin[which(fin$trait==traits[j]),"chisq_p_genes"]=T$p.value
 }


