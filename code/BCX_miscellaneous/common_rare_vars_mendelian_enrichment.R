### PAPER REVISION ANALYSIS

setwd("/Users/dv3/Desktop/finemapping_ukbb500k_final_release/")
mer=read.table("2019_01_14_ukbb500k_condsig_with_LD_clumping_and_comparison_to_Astle.txt", he=T, strings=F, sep="\t", comment.char="")
coloc=read.csv("colocalizations/all_coloc_results.csv", he=T,strings=F)

mendelian=scan("../Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")
head(mer)
mer$mendelian="no"
tmp=sapply(mer$Gene.Symbol.s..for.Most.Serious.Consequence, function(x){unlist(strsplit(x,","))})
head(tmp)
length(tmp)
k=unlist(lapply(tmp, function(x){length(which(x %in% mendelian))}))
mer$mendelian[which(k>=1)]="yes"

mer$freq_class="common"
mer$freq_class[which(mer$Minor.Allele.Frequency<0.01)]="rare"
table(mer$mendelian)

## too few colocs for checking in this way

# check enrichement of Mendelian SNPs by MAF, compared to random genes 

tss=read.table("../TSS.Ensemble_genes_downloaded_28_10_2019.txt", he=T, strings=F, sep="\t")
comm.overlap=c()
rare.overlap=c()
set.seed(123)
res=mer[!duplicated(mer$Unique.Variant.ID),]

tmp=sapply(res$Gene.Symbol.s..for.Most.Serious.Consequence, function(x){unlist(strsplit(x,","))})

for(i in 1:10000){
  random=sample(unique(tss$HGNC.symbol[which(!tss$HGNC.symbol %in% mendelian)]), 314, replace = F)
  res$random="no"
  k=unlist(lapply(tmp, function(x){length(which(x %in% random))}))
  res$random[which(k>=1)]="yes"
  comm.overlap=c(comm.overlap, dim(res[which(res$random=="yes" & res$freq_class=="common"),])[1])
  rare.overlap=c(rare.overlap, dim(res[which(res$random=="yes" & res$freq_class=="rare"),])[1])
  }

comm.overlap.rate=comm.overlap/10213*100
rare.overlap.rate=rare.overlap/510*100

comm.mend.rate=length(which(res$mendelian=="yes" & res$freq_class=="common"))/10213*100
rare.mend.rate=length(which(res$mendelian=="yes" & res$freq_class=="rare"))/510*100

comm.p=length(which(comm.overlap.rate>comm.mend.rate))/10000
comm.enrich=comm.mend.rate/mean(comm.overlap.rate)
comm.enrich.CI=1.96*sd(comm.overlap.rate)/sqrt(10000)
rare.p=length(which(rare.overlap.rate>rare.mend.rate))/10000
rare.enrich=rare.mend.rate/mean(rare.overlap.rate)
rare.enrich.CI=1.96*sd(rare.overlap.rate)/sqrt(10000)

par(mfrow=c(1,1))
pdf("Common_and_rare_variants_enrichment_of_mendelian_genes_09_04_20.pdf", width=5,height=4)
plot(density(comm.overlap.rate), main="Variant based enrichment of mendelian genes",xlab="% of GWAS variants assigend to random sets of 314 genes",xlim=c(0,13))
points(density(rare.overlap.rate),type="l",col="red")
abline(v=comm.mend.rate)
abline(v=rare.mend.rate, col="red")
legend("topright", c("common", "rare"), lty=c(1,1), col=c("black", "red"))
dev.off()

## repeat this analysis with N. of genes overlapping in addition to the N. of variants which might be biased by the finding above
gwas.genes.comm=unique(unlist(strsplit(res$Gene.Symbol.s..for.Most.Serious.Consequence[which(res$freq_class=="common")],",")))
gwas.genes.rare=unique(unlist(strsplit(res$Gene.Symbol.s..for.Most.Serious.Consequence[which(res$freq_class=="rare")],",")))
overlap.mend.comm=length(intersect(gwas.genes.comm, mendelian))
overlap.mend.rare=length(intersect(gwas.genes.rare, mendelian))

comm.genes.overlap=c()
rare.genes.overlap=c()

for(i in 1:10000){
  random=sample(unique(tss$HGNC.symbol[which(!tss$HGNC.symbol %in% mendelian)]), 314, replace = F)
  comm.genes.overlap=c(comm.genes.overlap,length(intersect(random, gwas.genes.comm)))
  rare.genes.overlap=c(rare.genes.overlap,length(intersect(random, gwas.genes.rare)))
}

comm.p=length(which(comm.genes.overlap>overlap.mend.comm))/10000
rare.p=length(which(rare.overlap.rate>rare.mend.rate))/10000
comm.enrich=overlap.mend.comm/mean(comm.genes.overlap)
comm.enrich.CI=1.96*sd(comm.genes.overlap)/sqrt(10000)
rare.enrich=overlap.mend.rare/mean(rare.genes.overlap)
rare.enrich.CI=1.96*sd(rare.genes.overlap)/sqrt(10000)


pdf("Common_and_rare_gene_enrichment_of_mendelian_genes_09_04_20.pdf", width=10,height=4)
par(mfrow=c(1,2))
plot(density(comm.genes.overlap), main="Gene based enrichment of mendelian genes",xlab="GWAS genes overlap with random sets of 314 genes", xlim=c(0,115))
#points(density(rare.genes.overlap),type="h",col="red")
abline(v=overlap.mend.comm)
#abline(v=overlap.mend.rare, col="red")
#legend("topright", c("common", "rare"), lty=c(1,1), col=c("black", "red"))
plot(density(rare.genes.overlap), main="Gene based enrichment of mendelian genes",xlab="GWAS genes overlap with random sets of 314 genes", xlim=c(0,36))
abline(v=overlap.mend.rare)
#legend("topright", c("common", "rare"), lty=c(1,1), col=c("black", "red"))
dev.off()


