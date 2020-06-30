#DIR=
#/lustre/scratch115/realdata/mdt3/projects/ukbb500k_t151/risk_score/clin_annot

args=commandArgs(TRUE)
i=args[1]
print(i) 
i=as.numeric(i)

setwd("/lustre/scratch115/realdata/mdt3/projects/ukbb500k_t151/risk_score/clin_annot")

# create list of mendelian genes +- 100 kb
mend=read.table("ST7.txt",he=T, strings=F)
tss=read.table("TSS.Ensemble_genes_downloaded_28_10_2019.txt", he=T, strings=F,sep="\t")
tss=tss[!duplicated(tss$HGNC.symbol),]
mend=merge(mend, tss, by.x="Gene_Symbol_HGNC", by.y="HGNC.symbol")
mend$start=mend$Gene.start..bp.-100000
mend$end=mend$Gene.end..bp.+100000
pheno_list=scan("../results_ukbb_final/snp_lists/pheno_list.txt","char")

pheno=pheno_list[i]
print(pheno)

# create list of peripheral associations
gwas=read.table("2019_01_14_ukbb500k_condsig.txt", he=T, strings=F, sep="\t", comment.char="")
gwas$mend="no"
for(k in 1:dim(mend)[1]){
gwas$mend[which(gwas$Chr..GRCh37.==mend$Chromosome.scaffold.name[k] & gwas$BP..GRCh37.>mend$start[k] & gwas$BP..GRCh37.<mend$end[k])]="yes"
}

periph=gwas[which(gwas$mend=="no"),c("Chr..GRCh37.", "BP..GRCh37.","Unique.Variant.ID", "rsID..where.available."),]
periph=periph[!duplicated(periph),]


# read gwas results
  res=read.table(pipe(paste("awk \'{if($7>0.0005 && $7<(1-0.0005) && $8>0.4) print $2\" \"$3\" \"$16}\' Other_GWAS_traits_for_Mendelian_enrichment/",pheno,"_gwas_normalised_imputed_full_panel.out_pruned.txt",sep="")), he=F, strings=F)
  
# create expected to plot  
  p.vals=res$V3[order(res$V3)]
  ppoi=ppoints(p.vals)
  ppoi=sort(qchisq(1-ppoi, 1))
  ppoi=pchisq(ppoi, 1, lower.tail=F)
  ppoi=ppoi[order(ppoi)]
  
  png(paste("QQplot_with_mend_genes_enrichment_pruned_",pheno,".png",sep=""))
  plot(-log10(ppoi), -log10(p.vals), xlab="Expected",ylab="Observed",col="grey", pch=19, main=paste("QQ plot",pheno))
  abline(0,1,lty=2)
  
  ### plot SMD
  tmp=mend[which(mend$Platform=="SMD"),]
  fin=c()
  for(j in 1:dim(tmp)[1]){
    pio=res[which(res$V1==tmp$Chromosome.scaffold.name[j] & res$V2>tmp$start[j] & res$V2<tmp$end[j]),]
    fin=rbind(fin, pio)
  }
  fin=fin[!duplicated(fin),]
  p.vals=fin$V3[order(fin$V3)]
  ppoi=ppoints(p.vals)
  ppoi=sort(qchisq(1-ppoi, 1))
  ppoi=pchisq(ppoi, 1, lower.tail=F)
  ppoi=ppoi[order(ppoi)]
  
  points(-log10(ppoi), -log10(p.vals), col="blue", pch=19)
  
  ### plot BPD
  tmp=mend[which(mend$Platform=="BPD"),]
  fin=c()
  for(j in 1:dim(tmp)[1]){
    pio=res[which(res$V1==tmp$Chromosome.scaffold.name[j] & res$V2>tmp$start[j] & res$V2<tmp$end[j]),]
    fin=rbind(fin, pio)
  }
  fin=fin[!duplicated(fin),]
  p.vals=fin$V3[order(fin$V3)]
  ppoi=ppoints(p.vals)
  ppoi=sort(qchisq(1-ppoi, 1))
  ppoi=pchisq(ppoi, 1, lower.tail=F)
  ppoi=ppoi[order(ppoi)]
  
  points(-log10(ppoi), -log10(p.vals), col="violet", pch=19)
  
  
  ### plot BMF
  tmp=mend[which(mend$Platform=="BMF"),]
  fin=c()
  for(j in 1:dim(tmp)[1]){
    pio=res[which(res$V1==tmp$Chromosome.scaffold.name[j] & res$V2>tmp$start[j] & res$V2<tmp$end[j]),]
    fin=rbind(fin, pio)
  }
  fin=fin[!duplicated(fin),]
  p.vals=fin$V3[order(fin$V3)]
  ppoi=ppoints(p.vals)
  ppoi=sort(qchisq(1-ppoi, 1))
  ppoi=pchisq(ppoi, 1, lower.tail=F)
  ppoi=ppoi[order(ppoi)]
  
  points(-log10(ppoi), -log10(p.vals), col="orange", pch=19)
  
 #legend("topleft", c("SMD", "BPD", "BMF"), col=c("blue","violet","orange"), pch=19)

# plot peripheral 
res$var=paste(res$V1, res$V2, sep="_")
periph$var=paste(periph$Chr..GRCh37., periph$BP..GRCh37., sep="_")
fin=res[which(res$var %in% periph$var),]
fin=fin[!duplicated(fin),]
p.vals=fin$V3[order(fin$V3)]
ppoi=ppoints(p.vals)
ppoi=sort(qchisq(1-ppoi, 1))
ppoi=pchisq(ppoi, 1, lower.tail=F)
ppoi=ppoi[order(ppoi)]

points(-log10(ppoi), -log10(p.vals), col="black", pch=19)

legend("topleft", c("SMD", "BPD", "BMF", "Peripheral","GW"), col=c("blue","violet","orange","black","grey"), pch=19)


dev.off()


