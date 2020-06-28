### VAR TO GENE

setwd("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/variant_to_gene_reassignment/")
fab=read.table("gene_reass_result.txt", he=T, strings=F)
head(fab)
dim(fab)

# read bridge genes list 
bridge=scan("../../Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")

fab$stand_bridge="no"
fab$stand_bridge[which(fab$stand_gene %in% bridge)]="yes"
fab$pc5_bridge="no"
fab$pc5_bridge[which(fab$pc5_gene %in% bridge)]="yes"
fab$pc10_bridge="no"
fab$pc10_bridge[which(fab$pc10_gene %in% bridge)]="yes"

table(fab$stand_bridge)
table(fab$pc5_bridge)
table(fab$pc10_bridge)

# read trans eQTL dataframe
trans=read.table("/Users/dv3/Desktop/trans_eQTL_in_blood_GWAS_sign_2018_11_14_extracted_290419.txt", he=T, strings=F)
trans=trans[,c("SNPChr.SNPPos_AssessedAllele_OtherAllele","SNPChr.SNPPos_OtherAllele_AssessedAllele","GeneSymbol")]
tmp1=merge(mer, trans, by.x="VAR", by.y="SNPChr.SNPPos_AssessedAllele_OtherAllele",all.x=T)
tmp2=merge(mer, trans, by.x="VAR", by.y="SNPChr.SNPPos_OtherAllele_AssessedAllele",all.x=T)
names(tmp1)[30]="transSNP"
names(tmp2)[30]="transSNP"
names(tmp1)[31]="trans_gene"
names(tmp2)[31]="trans_gene"
trans=rbind(tmp1[,c(1,31)], tmp2[,c(1,31)])
trans=trans[!duplicated(trans),]
trans=trans[which(!is.na(trans$trans_gene)),]
head(trans)
mmtr=merge(mer, trans, by.x="VAR", all.x=T)
fab$VAR=gsub(fab$VAR, pattern="chr", replacement="")
mmfab=merge(mmtr, fab, by="VAR", all.x=T)
mmfab$test_pc5=NA
mmfab$test_pc5[which(mmfab$trans_gene==mmfab$pc5_gene)]="correct"
mmfab$test_pc5[which(mmfab$trans_gene!=mmfab$pc5_gene)]="wrong"
table(mmfab$test_pc5, exclude=NULL)

###############################
### VEP vs network vs eQTLs ###
###############################

eqtl=read.csv("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/variant_to_gene_reassignment/2019_07_09_ld_signif_eqtl_coloc_results.csv", he=T, strings=F)
vep=read.table("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/Condind_and_PP_gt_0.5_vars_annotated_all_VEP_consequences.txt", he=T, strings=F)
head(vep)

# format variant ids
vep$VAR=sub(vep$ID,pattern="_",replacement=":")
vep$VAR=sub(vep$VAR,pattern="\\/",replacement="_")

#vep$VAR=paste(vep$CHROM,":",vep$POS,"_",vep$REF,"_",vep$ALT,sep="")
vep=vep[!duplicated(vep),]
vep[which(vep$VAR=="12:22610292_G_C"),]

# test how correct is VEP annotation
uniq=mer[!duplicated(mer$VAR),]
mm=merge(eqtl,uniq, by.x="condsig_var", by.y="VAR")
mm=mm[!duplicated(mm),]
mm$test="wrong"
for(i in 1:dim(mm)[1]){
if(mm$hgnc[i]==mm$gene[i])
  mm$test[i]="correct"
if(mm$hgnc[i] %in% unlist(strsplit(mm$gene[i],",")))
  mm$test[i]="correct"
}
table(mm$test) ## basic VEP annotation (worst consequence)

mm$test2="wrong"
for(i in 1:dim(mm)[1]){
  tmp=vep[which(vep$VAR==mm$condsig_var[i]),]
  if(mm$hgnc[i] %in% tmp$SYMBOL)
    mm$test2[i]="correct"
  if(mm$VEP_cons_most_serious[i]=="intergenic_variant" & mm$test[i]=="correct") #add this to account for closest TSS manual annotation
    mm$test2[i]="correct"
    
}
table(mm$test2)
mm[which(mm$test2=="wrong" & mm$test=="correct"),]

# summarize counts
A=data.frame(table(mm$test, mm$VEP_cons_most_serious))
A
A$perc_correct=NA
for(j in seq(1,22,2)){
  A$perc_correct[j]=A$Freq[j]/(A$Freq[j]+A$Freq[j+1])*100
}

B=data.frame(table(mm$test2, mm$VEP_cons_most_serious))
B
B$perc_correct=NA
for(j in seq(1,22,2)){
  B$perc_correct[j]=B$Freq[j]/(B$Freq[j]+B$Freq[j+1])*100
}

C=merge(A,B,by=c("Var1","Var2"))
C
names(C)=c("Var1","Var2","Freq_VEP_worst","P_correct_VEP_worst","Freq_VEP_all","P_correct_VEP_all")

# add closest genes for intergenic variants to VEP annotation # check the nearest TSS 
tss=read.table("/Users/dv3/Desktop/TSS.onlychr1.22.XY.MT_EnsembleGenes.v75_downloaded.07122017 copy.txt", he=T, strings=F)
J=which(mm$VEP_cons_most_serious=="intergenic_variant")
for(j in J){
  subs=tss[which(tss$chromosome_name==mm$chr[j]),]
  k=which.min(abs(subs$TSS-mm$pos[j]))
  if(length(k)>1){
  warning(paste("check ",j))
  }
  while(is.na(subs$hgnc_symbol[k])){
    subs=subs[-k,]
    k=which.min(abs(subs$TSS-mm$pos[j]))
  }
  mm$gene[j]=subs$hgnc_symbol[k]
}

# now repeat test assignments above

### add netwrok based re-assignment
coex=readRDS("/Users/dv3/Desktop/coexpression_networks/coexpression.rds")
edge_presence <- function(coexp, threshold) {
  yes_no <- matrix(0, ncol=ncol(coexp), nrow=nrow(coexp), dimnames=dimnames(coexp))
  yes_no[abs(coexp) > threshold] <- 1
  return(yes_no)
}
#mendelian=scan("/Users/dv3/Desktop/Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")
mat=edge_presence(coex, 0.2)
unlist(dimnames(mat)[1])->geni
length(geni)
length(unique(geni))
mat=mat[unique(geni),unique(geni)]
apply(mat, 1, sum)->k
j=which(k>2)
length(j)
mat=mat[j,j]
unlist(dimnames(mat)[1])->geni
#gene.pos=read.table("/Users/dv3/Desktop/iMK_Project/ChIP/TSS.onlychr1.22.XY.MT_EnsembleGenes.v75_downloaded.07122017.txt", he=T, strings=F)
gene.pos=tss[which(tss$hgnc_symbol %in% geni),]
gene.pos=gene.pos[,c(3,4,5,6)]
gene.pos=gene.pos[!duplicated(gene.pos),]
dim(gene.pos)

#re-assign genes as neighbours to mendelian genes 
mm$test_coex=NA
for(i in 1:dim(mm)[1]){
  chr=mm$chr[i]
  pos=mm$pos[i]
  subs=gene.pos[which(gene.pos$chromosome_name==chr & gene.pos$start_position-250000<pos & gene.pos$end_position+250000>pos),]
  if(dim(subs)[1]>0){
    subs$dist=pmin(abs(subs$start_position-pos), abs(subs$end_position-pos))
    mm$new_gene[i]=subs$hgnc_symbol[which.min(subs$dist)]
    if(mm$new_gene[i]==mm$hgnc[i])
      mm$test_coex[i]="correct"
    else
      mm$test_coex[i]="wrong"
  }
}
table(mm$test_coex)

D=data.frame(table(mm$test_coex, mm$VEP_cons_most_serious))
D
D$perc_correct=NA
for(j in seq(1,22,2)){
  D$perc_correct[j]=D$Freq[j]/(D$Freq[j]+D$Freq[j+1])*100
}
D
F=merge(C,D,by=c("Var1","Var2"))
F
names(F)[7]="Freq_coex"
names(F)[8]="P_correct_coex"
library(reshape2)
meltData=melt(F[which(F$Var1=="correct"),c("Var2", "P_correct_VEP_worst","P_correct_VEP_all","P_correct_coex")])
head(meltData)

ggplot(meltData, aes(x=Var2, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")+
  theme(axis.text.x = element_text(size=10,angle = 45),axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_blank())

F
# remove variants with less than 5 occurrences from plot 
aggregate(F$Freq_VEP_worst, list(F$Var2), sum)
meltData=meltData[which(meltData$Var2!="splice_region_variant"),]
meltData=meltData[which(meltData$Var2!="splice_acceptor_variant"),]
ggplot(meltData, aes(x=Var2, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")+
  theme(axis.text.x = element_text(size=10,angle = 45),axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_blank())


# if we were to choose a random consequence 
mm$random_vep=NA
for(i in 1:dim(mm)[1]){
  tmp=vep[which(vep$VAR==mm$condsig_var[i]),]
  if(dim(tmp)[1]>0){
  k=sample(1, c(1:dim(tmp)[1]))
  mm$random_vep[i]=tmp$SYMBOL[k]
  }else{
    if(mm$VEP_cons_most_serious[i]=="intergenic_variant")
      mm$random_vep[i]=mm$gene[i]
  }
}
mm$test_random="wrong"
mm$test_random[which(mm$random_vep==mm$hgnc)]="correct"
table(mm$test_random)
table(mm$test2)

E=data.frame(table(mm$test_random, mm$VEP_cons_most_serious))
E
E$perc_correct=NA
for(j in seq(1,22,2)){
  E$perc_correct[j]=E$Freq[j]/(E$Freq[j]+E$Freq[j+1])*100
}

G=F
G=merge(G,E,by=c("Var1","Var2"))
names(G)[11]="Freq_VEP_random"
names(G)[12]="P_correct_vep_random"

meltData=melt(G[which(G$Var1=="correct"),c("Var2", "P_correct_VEP_worst","P_correct_VEP_all","P_correct_coex","P_correct_vep_random")])
head(meltData)
meltData=meltData[which(meltData$Var2!="splice_region_variant"),]
meltData=meltData[which(meltData$Var2!="splice_acceptor_variant"),]


ggplot(meltData, aes(x=Var2, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")+scale_fill_brewer(palette="Blues")+
  theme(axis.text.x = element_text(size=10,angle = 45),axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_blank())

meltData=meltData[which(meltData$variable!="P_correct_coex"),]
ggplot(meltData, aes(x=Var2, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")+scale_fill_brewer(palette="Blues")+theme_bw()+
  theme(axis.text.x = element_text(size=10,angle = 45),axis.text=element_text(size=14),axis.title.x=element_blank(),axis.title.y=element_blank())

meltData$Var2=as.character(meltData$Var2)
meltData$Var2=gsub(meltData$Var2, pattern="_variant",replacement="")
meltData$Var2=gsub(meltData$Var2, pattern="_",replacement=" ")
meltData$Var2=gsub(meltData$Var2, pattern="3 prime",replacement="3'")
meltData$Var2=gsub(meltData$Var2, pattern="5 prime",replacement="5'")
unique(meltData$Var2)
meltData=meltData[which(meltData$Var2!="non coding transcript exon"),]
#meltData=meltData[which(meltData$variable!="P_correct_coex"),]
meltData$Var2=gsub(meltData$Var2, pattern=" gene",replacement="")

N=aggregate(G$Freq_VEP_all,list(G$Var2),sum )
N=N[c(1,2,3,4,5,6,10,11),]
N
#fix this manually
meltData[which(meltData$Var2=="intergenic"),]
meltData[which(meltData$Var2=="intergenic" & meltData$variable=="P_correct_vep_random"),"value"]=46.1538

pdf("VEP_annotation_summary_for_eQTLs.pdf",width=6, height=4)
ggplot(meltData, aes(x=Var2, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")+scale_fill_brewer(palette="Blues", labels=c("VEP worst conseq.","VEP all conseq.","VEP random conseq."))+theme_classic()+
  theme(axis.text.x = element_text(size=10,angle=90),axis.ticks.x = element_blank(),axis.text=element_text(size=14),axis.title.x=element_blank(),legend.title=element_blank())+scale_y_continuous(name="% correct")+
  annotate("text", label=paste("N=",N$x,sep=""),x=c(1:8),y=100)
dev.off()
getwd()


## add co-ex assignment
unique(meltData$variable)
meltData$variable=factor(meltData$variable, levels=c("P_correct_VEP_worst","P_correct_VEP_all","P_correct_vep_random","P_correct_coex"))
pdf("VEP_annotation_summary_for_eQTLs_with_network_based.pdf",width=6, height=4)
ggplot(meltData, aes(x=Var2, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")+scale_fill_brewer(palette="Reds",labels=c("VEP worst conseq.","VEP all conseq.","VEP random conseq.","Network based"))+theme_classic()+
  theme(axis.text.x = element_text(size=10,angle=90),axis.ticks.x = element_blank(),axis.text=element_text(size=14),axis.title.x=element_blank(),legend.title=element_blank())+scale_y_continuous(name="% correct")+
  annotate("text", label=paste("N=",N$x,sep=""),x=c(1:8),y=100)
dev.off()

## find how many ADDITIONAL variants we can assign to genes by using the co-expression network
head(mer)

# first fill VEP gaps with intergenic annotation to closest gene
tss=read.table("/Users/dv3/Desktop/TSS.onlychr1.22.XY.MT_EnsembleGenes.v75_downloaded.07122017 copy.txt", he=T, strings=F) # try to update this

J=which(mer$VEP_cons_most_serious=="intergenic_variant" | mer$VEP_cons_most_serious=="")
for(j in J){
  subs=tss[which(tss$chromosome_name==mer$chr[j]),]
  k=which.min(abs(subs$TSS-mer$pos[j]))
  if(length(k)>1){
    warning(paste("check ",j))
  }
  while(is.na(subs$hgnc_symbol[k])){
    subs=subs[-k,]
    k=which.min(abs(subs$TSS-mer$pos[j]))
  }
  mer$gene[j]=subs$hgnc_symbol[k]
}

#now repeat network based annotation for all variants
mer$gene_coex=NA
for(i in 1:dim(mer)[1]){
  chr=mer$chr[i]
  pos=mer$pos[i]
  subs=gene.pos[which(gene.pos$chromosome_name==chr & gene.pos$start_position-250000<pos & gene.pos$end_position+250000>pos),]
  if(dim(subs)[1]>0){
    subs$dist=pmin(abs(subs$start_position-pos), abs(subs$end_position-pos))
    mer$gene_coex[i]=subs$hgnc_symbol[which.min(subs$dist)]
  }
}

mer$coex_eq_vep=NA
mer$coex_eq_vep=apply(mer[,c("gene_coex","gene")],1,function(x){if(x[1] %in% unlist(strsplit(x[2],","))) return("yes") else return("no")})
table(mer$coex_eq_vep, mer$VEP_cons_most_serious)

getwd()
setwd("../finemapping_ukbb500k_final_release/")
rm(F)
write.table(mer, "Condind_and_PP_gt_0.5_vars_annotated_VEP_and_coexpression_plus_intergenic.txt", row.names = F, quote=F ,sep="\t")

gwas=unique(unlist(strsplit(mer$gene[which(mer$condind=="yes")],",")))
write.table(gwas, "List_gwas_genes_by_VEP_with_intergenic.txt",col.names = F, row.names = F, quote=F, sep="\n")

getwd()
gwas=read.table("../List_gwas_genes_by_VEP_with_intergenic.txt", he=F, strings=F, sep="\t")
mer=read.table("../Condind_and_PP_gt_0.5_vars_annotated_VEP_and_coexpression_plus_intergenic.txt", he=T, strings=F,sep="\t", comment.char="", fill=T)

## what happens to non-coding variants outside genes (intergenic, downstream and updtream) [?add intronic?]
nc=mer[which(mer$VEP_cons_most_serious %in% c("downstream_gene_variant","upstream_gene_variant","intergenic_variant")),]
nc=nc[which(nc$condind=="yes"),]
vep=unique(unlist(strsplit(nc$gene,",")))
netw=unique(nc$gene_coex)


