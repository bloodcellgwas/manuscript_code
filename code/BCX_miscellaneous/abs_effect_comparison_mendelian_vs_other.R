### Mendelian variants effect sizes comparison
setwd("/Users/dv3/Desktop/")

head(mer)
mendelian=scan("Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")
mer$mendelian=NA
mer$mendelian[which(mer$Gene.Symbol.s..for.Most.Serious.Consequence!="")]="no"
tmp=sapply(mer$Gene.Symbol.s..for.Most.Serious.Consequence, function(x){unlist(strsplit(x,","))})
k=unlist(lapply(tmp, function(x){length(which(x %in% mendelian))}))
mer$mendelian[which(k>=1)]="yes"
table(mer$mendelian)

mer$abs_effect=abs(mer$X.UNIVAR..Estimate.of.Additive.Allelic.Effect..REF.Baseline..ALT.Effect.)

mer=mer[order(mer$abs_effect, decreasing = T),]
uniq=mer[!duplicated(mer$Unique.Variant.ID),]

# missense variants
# check maf first
library(MatchIt)

# create dataset of matched variants
m.in=uniq[which(uniq$Most.Serious.VEP.Consequence.of.Variant=="missense_variant"),c("Most.Serious.VEP.Consequence.of.Variant","Minor.Allele.Frequency","abs_effect","Unique.Variant.ID","mendelian")]
m.in$mendelian[which(m.in$mendelian=="no")]=0
m.in$mendelian[which(m.in$mendelian=="yes")]=1
m.out=matchit(mendelian~Minor.Allele.Frequency, data=m.in, ratio=2, method="optimal")
m.out
m.data<-match.data(m.out)
wilcox.test(m.data$Minor.Allele.Frequency~m.data$mendelian) # check that the selection worked well
wilcox.test(m.data$abs_effect~m.data$mendelian)
fin=(m.data)


# intron_variant 
m.in=uniq[which(uniq$Most.Serious.VEP.Consequence.of.Variant=="intron_variant"),c("Most.Serious.VEP.Consequence.of.Variant","Minor.Allele.Frequency","abs_effect","Unique.Variant.ID","mendelian")]
m.in$mendelian[which(m.in$mendelian=="no")]=0
m.in$mendelian[which(m.in$mendelian=="yes")]=1
m.out=matchit(mendelian~Minor.Allele.Frequency, data=m.in, ratio=2, method="optimal")
m.out
m.data<-match.data(m.out)
wilcox.test(m.data$Minor.Allele.Frequency~m.data$mendelian) # check that the selection worked well
wilcox.test(m.data$abs_effect~m.data$mendelian)
fin=rbind(fin, m.data)

# synonymous
m.in=uniq[which(uniq$Most.Serious.VEP.Consequence.of.Variant=="synonymous_variant"),c("Most.Serious.VEP.Consequence.of.Variant","Minor.Allele.Frequency","abs_effect","Unique.Variant.ID","mendelian")]
m.in$mendelian[which(m.in$mendelian=="no")]=0
m.in$mendelian[which(m.in$mendelian=="yes")]=1
m.out=matchit(mendelian~Minor.Allele.Frequency, data=m.in, ratio=2, method="optimal")
m.out
m.data<-match.data(m.out)
wilcox.test(m.data$Minor.Allele.Frequency~m.data$mendelian) # check that the selection worked well
wilcox.test(m.data$abs_effect~m.data$mendelian)
fin=rbind(fin, m.data)

# downstream
m.in=uniq[which(uniq$Most.Serious.VEP.Consequence.of.Variant=="downstream_gene_variant"),c("Most.Serious.VEP.Consequence.of.Variant","Minor.Allele.Frequency","abs_effect","Unique.Variant.ID","mendelian")]
m.in$mendelian[which(m.in$mendelian=="no")]=0
m.in$mendelian[which(m.in$mendelian=="yes")]=1
m.out=matchit(mendelian~Minor.Allele.Frequency, data=m.in, ratio=2, method="optimal")
m.out
m.data<-match.data(m.out)
wilcox.test(m.data$Minor.Allele.Frequency~m.data$mendelian) # check that the selection worked well
wilcox.test(m.data$abs_effect~m.data$mendelian)
fin=rbind(fin, m.data)

# upstream_gene_variant
m.in=uniq[which(uniq$Most.Serious.VEP.Consequence.of.Variant=="upstream_gene_variant"),c("Most.Serious.VEP.Consequence.of.Variant","Minor.Allele.Frequency","abs_effect","Unique.Variant.ID","mendelian")]
m.in$mendelian[which(m.in$mendelian=="no")]=0
m.in$mendelian[which(m.in$mendelian=="yes")]=1
m.out=matchit(mendelian~Minor.Allele.Frequency, data=m.in, ratio=2, method="optimal")
m.out
m.data<-match.data(m.out)
wilcox.test(m.data$Minor.Allele.Frequency~m.data$mendelian) # check that the selection worked well
wilcox.test(m.data$abs_effect~m.data$mendelian)
fin=rbind(fin, m.data)


# produce plots and summarize
annots=c("upstream_gene_variant","downstream_gene_variant", "synonymous_variant","intron_variant", "missense_variant" )

pdf("Abs_effect_by_mend_all_snps.pdf", width=10, height=2,useDingbats=FALSE)
tmp=uniq[which(!is.na(uniq$mendelian) & uniq$Most.Serious.VEP.Consequence.of.Variant %in% annots),]
gg=ggplot(tmp, aes(x=mendelian,y=abs_effect))+geom_violin(position=position_dodge(width=0.8))+
  geom_boxplot(position=position_dodge(width=0.5),width=0.2, aes(fill="grey") )+
  scale_y_continuous(trans = 'log10')+
  facet_wrap(~Most.Serious.VEP.Consequence.of.Variant, ncol=5)+
  theme(axis.line.y = element_blank(),axis.text=element_text(size=14),legend.title=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(), legend.position = "")
print(gg)  
dev.off()

pdf("Abs_effect_by_mend_matched_snps_by_maf.pdf", width=10, height=2,useDingbats=FALSE)
tmp=fin[which(!is.na(fin$mendelian) & fin$Most.Serious.VEP.Consequence.of.Variant %in% annots),]
gg=ggplot(tmp, aes(x=mendelian,y=abs_effect))+geom_violin(position=position_dodge(width=0.8))+
  geom_boxplot(position=position_dodge(width=0.5),width=0.2, aes(fill="grey") )+
  scale_y_continuous(trans = 'log10')+
  facet_wrap(~Most.Serious.VEP.Consequence.of.Variant, ncol=5)+
  theme(axis.line.y = element_blank(),axis.text=element_text(size=14),legend.title=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(), legend.position = "")
print(gg)
dev.off()

pdf("Check_MAF_by_mend_matched_snps_by_maf.pdf", width=10, height=2,useDingbats=FALSE)
gg=ggplot(tmp, aes(x=mendelian,y=Minor.Allele.Frequency))+geom_violin(position=position_dodge(width=0.8))+
  geom_boxplot(position=position_dodge(width=0.5),width=0.2, aes(fill="grey") )+
  facet_wrap(~Most.Serious.VEP.Consequence.of.Variant, ncol=5)+
  theme(axis.line.y = element_blank(),axis.text=element_text(size=14),legend.title=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(), legend.position = "")
print(gg)
dev.off()


fin=read.table("Desktop/matched_dataset_for_effect_size_comparison.txt", he=T, strings=F)
fin$abs_effect=as.numeric(fin$abs_effect)
for(a in annots){
  print(a)
 print((wilcox.test(fin$abs_effect[which(fin$Most.Serious.VEP.Consequence.of.Variant==a)]~fin$mendelian[which(fin$Most.Serious.VEP.Consequence.of.Variant==a)])))
}
