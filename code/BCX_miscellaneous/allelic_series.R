### check allele series ### 
### define allelic series as loci with multiple independent signal and overlap this set with colocalisation results


# confront with colocalisation results
# read all files
setwd("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/")
coloc=read.csv("colocalizations/all_coloc_results.csv", he=T, strings=F)
head(coloc)
condsig=read.table("2019_01_14_ukbb500k_condsig_with_LD_clumping_and_comparison_to_Astle.txt", he=T, strings=F, sep="\t",comment.char="")

fm=read.table("Input_files/conditional_blocks/All_jma.txt", he=T, strings=F, sep="\t")
fm$Associated.Blood.Index=tolower(fm$Associated.Blood.Index)
fm$Associated.Blood.Index=gsub(fm$Associated.Blood.Index, pattern=" ", replacement="_")
fm$blockid=paste(fm$Associated.Blood.Index, fm$FINEMAP_Block, sep="_")

# select loci with multiple independent signals
fin=c()
for(i in 1:dim(coloc)[1]){
var=coloc$variant[i]
ph=coloc$trait[i]
fm[which(fm$Unique.Variant.ID==var & fm$Associated.Blood.Index==ph),]->tmp
k=dim(fm[which(fm$blockid==tmp$blockid[1]),])[1]
  if(k>1){
    fin=rbind(fin,tmp)
  }
}

write.table(fin, "/Users/dv3/Desktop/finemapping_ukbb500k_final_release/Colocalising_variants_falling_in_blocks_with_multiple_signals.txt", row.names=F, quote=F)

### manuscript plots
library(ggplot2)
setwd("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/")
coloc_as=read.table("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/Colocalising_variants_falling_in_blocks_with_multiple_signals.txt", he=T, strings=F)
tss=read.table("/Users/dv3/Desktop/TSS.Ensemble_genes_downloaded_28_10_2019.txt", he=T, strings=F, sep="\t")
head(coloc_as)
dim(coloc_as)
PLOTS=unique(coloc_as$variant)
coloc_as$disease[grep("IBS",coloc_as$disease)]="IBS_CD_UC"
coloc_as=coloc_as[!duplicated(coloc_as),]
coloc_as=coloc_as[order(abs(coloc_as$disease_sentinel_zscore), decreasing = T),] ## among multiple disease GWAS take the one with the largest Z score for the variant
coloc_as=coloc_as[which(!duplicated(coloc_as[,c("trait", "variant")])),]
coloc_as=coloc_as[order(coloc_as$PPA_3, decreasing = T),] #take the trait that colocaslises best
coloc_as=coloc_as[!duplicated(coloc_as$variant),]

# associate each region to the relevant fine-mapping block
blocks=read.table("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/Input_files/conditional_blocks/All_jma.txt", he=T, strings=F, sep="\t")
blocks$BP..GRCh37.=as.numeric(blocks$BP..GRCh37.)
blocks=read.table("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/Input_files/conditional_blocks/All_jma.txt", he=T, strings=F, sep="\t")
blocks$pheno=tolower(blocks$Associated.Blood.Index)
blocks$pheno=gsub(blocks$pheno, pattern=" ", replacement="_")

write.table(coloc_as,"Dataset_for_figure_5_FtoK.txt", row.names = F, quote=F, sep="\t")
PLOTS=unique(coloc_as$variant)

# produce plots
for(j in 1:5){
  var=PLOTS[j]
  ph=coloc_as$trait[which(coloc_as$variant==var)]
bl=blocks[which(blocks$Unique.Variant.ID==var & blocks$pheno %in% ph),]
 tmp=blocks[which(blocks$FINEMAP_Block==bl$FINEMAP_Block & blocks$pheno %in% ph),]
rrange=c(tmp$Chr..GRCh37.[1], tmp$BP..GRCh37.[1]-250000, tmp$BP..GRCh37.[dim(tmp)[1]]+250000)
rrange=as.numeric(rrange)
geni=tss[which(tss$Chromosome.scaffold.name==tmp$Chr..GRCh37.[1] & ((tss$Gene.end..bp.>rrange[2] & tss$Gene.end..bp.<rrange[3]) | (tss$Gene.start..bp.>rrange[2] & tss$Gene.start..bp.<rrange[3]))),]
geni=geni[which(geni$HGNC.symbol!=""),]
geni=geni[!duplicated(geni$HGNC.symbol),]
geni=geni[order(geni$Gene.start..bp.),]
if(min(geni$Gene.end..bp.<rrange[2]))
  geni$Gene.start..bp.[which.min(geni$Gene.start..bp.)]==rrange[2]
if(max(geni$Gene.end..bp.>rrange[3]))
  geni$Gene.end..bp.[which.max(geni$Gene.end..bp.)]=rrange[3]
geni

#bottom part of panel
gg=ggplot(data=tmp)+
  geom_point(data=tmp,aes(x=BP..GRCh37., y=0.8), shape=18,size=5)+
  scale_x_continuous(limits=c(min(min(geni$Gene.start..bp.),rrange[2]),max(max(geni$Gene.end..bp.),rrange[3])))+
  scale_y_continuous(limits=c(0.5,1.1))+
geom_rect(data=geni, mapping=aes(xmin=Gene.start..bp., xmax=Gene.end..bp., ymin=0.675, ymax=0.725), color="black", fill="white")+
 geom_text(data=geni, aes(x=Gene.start..bp.+(Gene.end..bp.-Gene.start..bp.)/2, y=c(rep(c(0.6,0.65),round(dim(geni)[1]/2)), rep(0.6, dim(geni)[1]%%2)), label=geni$HGNC.symbol), size=4)+
  theme(axis.text=element_text(size=12),axis.line.x =element_line(size=1.5), axis.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank())
  
pdf(paste(coloc_as$gene[j], "lower_panel.pdf", sep=""), width=7, height=4)
print(gg)
dev.off()
}

# upper panel
for(j in 1:5){
  var=PLOTS[j]
  ph=coloc_as$trait[which(coloc_as$variant==var)]
  bl=blocks[which(blocks$Unique.Variant.ID==var & blocks$pheno %in% ph),]
  tmp=blocks[which(blocks$FINEMAP_Block==bl$FINEMAP_Block & blocks$pheno %in% ph),]
  #print(tmp[which(tmp$Minor.Allele.Frequency<0.01),])
a=max(abs(tmp$EFFECT))
pdf(paste(coloc_as$gene[j], "upper_panel_new.pdf", sep=""), width=5, height=3, useDingbats=FALSE)
par(mar=c(5.1, 4.1, 4.1, 5.1))
print(plot(x=tmp$BP..GRCh37., y=tmp$EFFECT, pch=20,axes=FALSE,xlab= "", ylab = "Beta", col="red", xaxt="n", ylim=c(-a,a), col.lab="red"))
print(abline(h=0))
print(axis(side=2, col="red", at=round(seq(-a,a,(2*a/6)),2), col.ticks="red",  col.axis="red", xlab="Beta"))
#print(points(x=tmp$BP..GRCh37., y=tmp$EFFECT,pch=18, axes = FALSE, bty = "n", xlim=rrange[2:3],xlab= "", ylab = "", col="red"))
par(new = TRUE, mar=c(5.1, 4.1, 4.1, 5.1))
print(plot(x=tmp$BP..GRCh37., y=tmp$Minor.Allele.Frequency, axes = FALSE, pch=20, xaxt="n",ylab="",xlab="",ylim=c(-0.5, 0.5)))
#print(points(x=tmp$BP..GRCh37., y=tmp$Minor.Allele.Frequency,pch=18, xlim=rrange[2:3]))
print(axis(side=4, at=c(0,0.5)))
mtext("MAF",side=4, line=2)
dev.off()
}
j

tmp=read.table("Dataset_for_figure_5_K.txt", he=T, strings=F, sep="\t", comment.char = "")
a=max(abs(tmp$X.UNIVAR..Estimate.of.Additive.Allelic.Effect..REF.Baseline..ALT.Effect.))
pdf(paste("Asthma_", "upper_panel_new.pdf", sep=""), width=5, height=3, useDingbats=FALSE)
par(mar=c(5.1, 4.1, 4.1, 5.1))
print(plot(x=tmp$BP..GRCh37., y=tmp$X.UNIVAR..Estimate.of.Additive.Allelic.Effect..REF.Baseline..ALT.Effect., pch=20,axes=FALSE,xlab= "", ylab = "Beta", col="red", xaxt="n", ylim=c(-a,a), col.lab="red"))
print(abline(h=0))
print(axis(side=2, col="red", at=round(seq(-a,a,(2*a/6)),2), col.ticks="red",  col.axis="red", xlab="Beta"))
#print(points(x=tmp$BP..GRCh37., y=tmp$EFFECT,pch=18, axes = FALSE, bty = "n", xlim=rrange[2:3],xlab= "", ylab = "", col="red"))
par(new = TRUE, mar=c(5.1, 4.1, 4.1, 5.1))
print(plot(x=tmp$BP..GRCh37., y=tmp$Minor.Allele.Frequency, axes = FALSE, pch=20, xaxt="n",ylab="",xlab="",ylim=c(-0.5, 0.5)))
#print(points(x=tmp$BP..GRCh37., y=tmp$Minor.Allele.Frequency,pch=18, xlim=rrange[2:3]))
print(axis(side=4, at=c(0,0.5)))
mtext("MAF",side=4, line=2)
dev.off()
