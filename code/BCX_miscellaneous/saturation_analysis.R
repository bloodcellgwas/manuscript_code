### saturation analysis 

#####################
### variant-based ###
#####################

sat=read.table("/Users/dv3/Desktop/Discovery_saturation/Table_N_variants_sign.txt", he=T, strings=F)
head(sat)
sat$N=NA
sat$N[which(sat$cohort=="INT")]=35000
sat$N[which(sat$cohort=="UKBiLEVE")]=43500
sat$N[which(sat$cohort=="UKBiobank")]=83000
sat$N[which(sat$cohort=="UKBB")]=400000
head(sat)

sat.t=sat[which(sat$pheno=="mono"),]
ggplot(sat, aes(x=N, y=N_sign_vars,group=pheno))+geom_line()+geom_point()

### interpolate based on INT, UKBiLEVE and UKBiobank 

head(sat)
mod=glm(sat.t$N_sign_vars~sat.t$N)
summary(mod)
mod$coefficients
y.p=mod$coefficients[2]*400000+mod$coefficients[1]
y.p

sat$N=as.numeric(sat$N)
sat$N_sign_vars=as.numeric(sat$N_sign_vars)
str(sat)
head(sat)
for(ph in unique(sat$pheno)){
 sat.t=sat[which(sat$pheno==ph),]
 sat.t$N=as.numeric(sat.t$N)
 sat.t$N_sign_vars=as.numeric(sat.t$N_sign_vars)
 sat$N=as.numeric(sat$N)
 sat$N_sign_vars=as.numeric(sat$N_sign_vars)
 mod=glm(sat.t$N_sign_vars[which(sat.t$cohort!="UKBB")]~sat.t$N[which(sat.t$cohort!="UKBB")])
 y.p=mod$coefficients[2]*400000+mod$coefficients[1]
 line=c("UKBB_pred", ph, y.p, 400000)
 sat=rbind(sat, line)
}

one=sat[which(sat$cohort!="UKBB_pred"),]
two=sat[which(sat$cohort!="UKBB"),]
one$N=as.numeric(one$N)
one$N_sign_vars=as.numeric(one$N_sign_vars)
ggplot()+geom_line(data=one, aes(x=N, y=N_sign_vars,color=pheno))+geom_point()+
  geom_line(data=two, aes(x=N, y=N_sign_vars,color=pheno),linetype="dashed")

two$N=as.numeric(two$N)
two$N_sign_vars=as.numeric(two$N_sign_vars)

ggplot(data=one, aes(x=N, y=N_sign_vars))+geom_line()+geom_point()+facet_wrap(~pheno)+geom_line(data=two, aes(x=N, y=N_sign_vars),linetype="dashed")

unique(one$pheno)
pls=c("pct","plt","pdw","mpv")
reds=c("mch","mchc","mcv","rbc","hct","hgb","hlr","hlr_p","irf","ret","ret_p")
phils=c("eo_p","eo","mono","mono_p","neut","neut_p","baso","baso_p")
lymph=c("lymph_p","lymph","wbc")

one$pheno=factor(one$pheno,levels=c(lymph,phils,pls,reds))
two$pheno=factor(two$pheno,levels=c(lymph,phils,pls,reds))

ggplot(data=one, aes(x=N/10000, y=N_sign_vars))+geom_line()+
  geom_point()+facet_wrap(~pheno,ncol=6)+
  geom_line(data=two, aes(x=N/10000, y=N_sign_vars),linetype="dashed")+
  theme(axis.text=element_text(size=15), strip.text.x=element_text(size=15))+xlab("")+ylab("")



setwd("/Users/dv3/Desktop/Discovery_saturation/")
dir()
options(scipen=999)

pdf("genome_wide_loci_overview.pdf")
for(ph in pheno){
res=read.table(paste(ph,"_LD_locus_based_res_summary.txt",sep=""), he=T, strings=F)
res$cumsum=cumsum(res$start/10000)
print(plot(x=res$cumsum, y=res$UKBB-1, pch=20, cex=3, col="grey", ylim=c(0,3.5), xlab="genome-wide loci", ylab="", xaxt="n",yaxt="n", main=ph))
print(points(x=res$cumsum[which(res$INT==1)], y=res$INT[which(res$INT==1)], pch=20, cex=3))
print(points(x=res$cumsum[which(res$UKBiLEVE==1)], y=res$UKBiLEVE[which(res$UKBiLEVE==1)]*2, pch=20, cex=3, col="blue"))
print(points(x=res$cumsum[which(res$UKBiobank==1)], y=res$UKBiobank[which(res$UKBiobank==1)]*3, pch=20, cex=3, col="orange"))
print(mtext(side=2, line=1, c("UKBB","INT","UKBiLEVE","UKBiobank"), at=c(0,1,2,3)))
}
dev.off()

###################
### locus-based ###
###################

N=c(30000,43500,83000,400000)
one=data.frame(cohort=rep(c("INT","UKBiLEVE","UKBiobank","UKBB"),26),N=rep(N,26), N_loci=NA, pheno=rep(pheno,each=4))
head(one)
two=data.frame(cohort=rep(c("INT","UKBiLEVE","UKBiobank","UKBB_pred"),25),N=rep(x,25), N_loci=NA,pheno=pheno)
head(two)

for(ph in pheno){
res=read.table(paste(ph,"_LD_locus_based_res_summary.txt",sep=""), he=T, strings=F)
y=unlist(apply(res[,c("INT","UKBiLEVE","UKBiobank","UKBB")],2,sum))
one[which(one$pheno==ph & one$cohort=="INT"), "N_loci"]=y["INT"]
one[which(one$pheno==ph & one$cohort=="UKBiLEVE"), "N_loci"]=y["UKBiLEVE"]
one[which(one$pheno==ph & one$cohort=="UKBiobank"), "N_loci"]=y["UKBiobank"]
one[which(one$pheno==ph & one$cohort=="UKBB"), "N_loci"]=y["UKBB"]
mod=lm(y[1:3]~c(30000,43500,83000))
y.p=mod$coefficients[2]*400000+mod$coefficients[1]
two[which(two$pheno==ph & two$cohort=="INT"), "N_loci"]=y["INT"]
two[which(two$pheno==ph & two$cohort=="UKBiLEVE"), "N_loci"]=y["UKBiLEVE"]
two[which(two$pheno==ph & two$cohort=="UKBiobank"), "N_loci"]=y["UKBiobank"]
two[which(two$pheno==ph & two$cohort=="UKBB_pred"), "N_loci"]=y.p
}

one$N=as.numeric(one$N)
one$N_loci=as.numeric(one$N_loci)
two$N=as.numeric(two$N)
two$N_loci=as.numeric(two$N_loci)
ggplot()+geom_line(data=one, aes(x=N, y=N_loci,color=pheno))+geom_point()+
  geom_line(data=two, aes(x=N, y=N_loci,color=pheno),linetype="dashed")

ggplot(data=one, aes(x=N, y=N_loci))+geom_line()+
  geom_point()+facet_wrap(~pheno)+geom_line(data=two, aes(x=N, y=N_loci),linetype="dashed")

pls=c("plt","pdw","mpv")
reds=c("mch","mchc","mcv","rbc","hct","hgb","hlr","hlr_p","irf","ret","ret_p")
phils=c("eo_p","eo","mono","mono_p","neut","neut_p","baso","baso_p")
lymph=c("lymph_p","lymph","wbc")

one$pheno=factor(one$pheno,levels=c(lymph,phils,pls,reds))
two$pheno=factor(two$pheno,levels=c(lymph,phils,pls,reds))

ggplot(data=one, aes(x=N/10000, y=N_loci))+geom_line()+
  geom_point()+facet_wrap(~pheno,ncol=6)+
  geom_line(data=two, aes(x=N/10000, y=N_loci),linetype="dashed")+
  theme(axis.text=element_text(size=15), strip.text.x=element_text(size=15))+xlab("")+ylab("")



##################
### gene-based ###
##################
getwd()
#how many genes are discovered? 
tab=read.table("Table_N_genes.txt", he=T, strings=F)
tab$N=NA
tab$N[which(tab$cohort=="INT")]=35000
tab$N[which(tab$cohort=="UKBiLEVE")]=43500
tab$N[which(tab$cohort=="UKBiobank")]=83000
tab$N[which(tab$cohort=="UKBB")]=400000

for(ph in unique(tab$pheno)){
  tab.t=tab[which(tab$pheno==ph),]
  tab.t$N_genes=as.numeric(tab.t$N_genes)
  tab.t$N=as.numeric(tab.t$N)
  mod=glm(tab.t$N_genes[which(tab.t$cohort!="UKBB")]~tab.t$N[which(tab.t$cohort!="UKBB")])
  y.p=mod$coefficients[2]*400000+mod$coefficients[1]
  line=c("UKBB_pred", ph, y.p,NA,NA,NA, 400000)
  tab=rbind(tab, line)
}

one=tab[which(tab$cohort!="UKBB_pred"),]
two=tab[which(tab$cohort!="UKBB"),]
one$N=as.numeric(one$N)
one$N_genes=as.numeric(one$N_genes)
two$N=as.numeric(two$N)
two$N_genes=as.numeric(two$N_genes)
ggplot(data=one, aes(x=N, y=N_genes))+geom_line()+
  geom_point()+facet_wrap(~pheno)+geom_line(data=two, aes(x=N, y=N_genes),linetype="dashed")

# sort phenotypes by cell type 
pls=c("plt","pdw","mpv")
reds=c("mch","mchc","mcv","rbc","hct","hgb","hlr","hlr_p","irf","ret","ret_p")
phils=c("eo_p","eo","mono","mono_p","neut","neut_p","baso","baso_p")
lymph=c("lymph_p","lymph","wbc")

one$pheno=factor(one$pheno,levels=c(lymph,phils,pls,reds))
two$pheno=factor(two$pheno,levels=c(lymph,phils,pls,reds))
ggplot(data=one, aes(x=N/10000, y=N_genes))+geom_line()+
  geom_point()+facet_wrap(~pheno,ncol=6)+
  geom_line(data=two, aes(x=N/10000, y=N_genes),linetype="dashed")+
  theme(axis.text=element_text(size=15), strip.text.x=element_text(size=15))+xlab("")+ylab("")

head(one)


### plot details for one phenotype ###
### add percentage plots instead of raw numbers so that curves are more comparable 

pheno="mpv"
#pheno="eo"
vars=read.table("Table_N_variants_sign.txt", he=T, strings=F)
vars=vars[which(vars$pheno==pheno),]
head(vars)
vars$N=c(33501, 43597, 83637, 391598)
vars$N_sign_vars_p=vars$N_sign_vars/vars$N_sign_vars[4]
vars
loci=read.table("Table_N_loci.txt", he=T, strings=F) #check this
loci=loci[which(loci$pheno==pheno),]
loci
loci$N=c(33501, 43597, 83637, 391598)
loci$N_loci_p=loci$N_loci/loci$N_loci[4]
geni=read.table("Table_N_genes.txt", he=T, strings=F)
geni=geni[which(geni$pheno==pheno),]
geni
geni$N=c(33501, 43597, 83637, 391598)
geni$N_genes_p=geni$N_genes/geni$N_genes[4]
geni
par(mfrow=c(1,3))
plot(vars$N, vars$N_sign_vars, pch=20, cex=2, cex.axis=1.5)
plot(loci$N, loci$N_loci, pch=20,cex=2, cex.axis=1.5)
plot(geni$N, geni$N_genes, pch=20,cex=2, cex.axis=1.5)
#dev.off()
plot(vars$N, vars$N_sign_vars_p, pch=20, cex=2, cex.axis=1.5,xlab="N",ylab="N variants")
plot(loci$N, loci$N_loci_p, pch=20,cex=2, cex.axis=1.5,xlab="N",ylab="N loci")
plot(geni$N, geni$N_genes_p, pch=20,cex=2, cex.axis=1.5,xlab="N",ylab="N genes")

mod1=lm(vars$N_sign_vars[1:3]~vars$N[1:3])
plot(vars$N, vars$N_sign_vars, pch=20, cex=2, cex.axis=1.5, xlab="Cohort size",ylab="N variants",cex.lab=1.5,main="MPV")
abline(mod1$coefficients[1:2], lty=2)
summary(mod1)
mod2=lm(vars$N_sign_vars~vars$N)
summary(mod2)
abline(mod2$coefficients, lty=2, col="red")
tmp=sqrt(vars$N)
mod3=lm(vars$N_sign_vars~tmp)
summary(mod3)
x=seq(min(vars$N), max(vars$N), 10000)
y=mod3$coefficients[2]*sqrt(x)+mod3$coefficients[1]
lines(x,y, col="red")

plot(loci$N, loci$N_loci, pch=20, cex=2, cex.axis=1.5, xlab="Cohort size",ylab="N loci",cex.lab=1.5,main="MPV")
mod1=lm(loci$N_loci[1:3]~loci$N[1:3])
abline(mod1$coefficients[1:2], lty=2)
mod2=lm(loci$N_loci~loci$N)
abline(mod2$coefficients, lty=2, col="red")
tmp=sqrt(loci$N)
mod3=lm(loci$N_loci~tmp)
x=seq(min(loci$N), max(loci$N), 10000)
y=mod3$coefficients[2]*sqrt(x)+mod3$coefficients[1]
lines(x,y, col="red")

geni$N=c(33501, 43597, 83637, 391598)
plot(geni$N, geni$N_genes, pch=20, cex=2, cex.axis=1.5, xlab="Cohort size",ylab="N genes",cex.lab=1.5,main="MPV")
mod1=lm(geni$N_genes[1:3]~geni$N[1:3])
abline(mod1$coefficients[1:2], lty=2)
mod2=lm(geni$N_genes~geni$N)
abline(mod2$coefficients, lty=2, col="red")
tmp=sqrt(geni$N)
mod3=lm(geni$N_genes~tmp)
x=seq(min(geni$N), max(geni$N), 10000)
y=mod3$coefficients[2]*sqrt(x)+mod3$coefficients[1]
lines(x,y, col="red")

# geni
x=seq(1,850000000,100000)
y=mod3$coefficients[2]*sqrt(x)+mod3$coefficients[1]
plot(x,y,lty=1, main="MPV Model projection up tp 400M participants",xlab="Cohort size",ylab="Genes discovered")
abline(h=20000,col="red")
tail(which(y<20000))
x[1066] #106500001

#loci
x=seq(1,850000000,100000)
y=mod3$coefficients[2]*sqrt(x)+mod3$coefficients[1]
plot(x,y,lty=1, main="Model projection up tp 400M participants",xlab="Cohort size",ylab="Loci discovered")
abline(h=1705, col="red")

##########################
### variants with meta ###
##########################

vars=read.table("Table_N_variants_sign_with_meta.txt", he=T, strings=F)
pheno="mpv"
vars=vars[which(vars$pheno==pheno),]
head(vars)
vars$N=c(33501, 43597, 83637, 160735,391598)
#loci=read.table("Table_N_loci.txt", he=T, strings=F) #check this
#loci=loci[which(loci$pheno==pheno),]
#loci
loci=loci[1:4,]
loci$N=c(33501, 43597, 83637, 391598)
geni=read.table("Table_N_genes.txt", he=T, strings=F)
geni=geni[which(geni$pheno==pheno),]
par(mfrow=c(1,3))
plot(vars$N, vars$N_sign_vars, pch=20, cex=2, cex.axis=1.5)
plot(loci$N, loci$N_loci, pch=20,cex=2, cex.axis=1.5)
plot(geni$N, geni$N_genes, pch=20,cex=2, cex.axis=1.5)
dev.off()

######################
### best model fit ###
######################

res=vars
head(res)
res$best_mod_vars=NA
res$best_mod_loci=NA
res$best_mod_geni=NA

#plot best fitted models forN of variants

vars.pred=data.frame(pheno=NA, x=NA, y=NA)
best.mod=data.frame(pheno=unique(res$pheno), model=NA)
pdf("Best_model_for_N_of_vars_per_pheno.pdf")
for(ph in unique(res$pheno)){
  vars.t=vars[which(vars$pheno==ph),]
  vars.t$N=c(30000,43500,83000,400000)
  mod1=lm(N_sign_vars~N, data=vars.t)
  mod2=lm(N_sign_vars~sqrt(N), data=vars.t)
  mod3=lm(N_sign_vars~sqrt(N)+N, data=vars.t)
  mod4=lm(N_sign_vars~log(N), data=vars.t)
  mod5=lm(N~N_sign_vars^2+N_sign_vars, data=vars.t)
  k=which.max(c(summary(mod1)$r.squared,summary(mod2)$r.squared,summary(mod3)$r.squared,summary(mod4)$r.squared,summary(mod5)$r.squared))
  
x=seq(1,400000,10000)
y=predict(get(paste("mod",k,sep="")), newdata=data.frame(N=x))
vars.pred.t=data.frame(pheno=ph, x=x, y=y)
mod.names=c("linear", "y~sqrt(x)", "y~sqrt(x)+x", "logarithmic", "quadratic in y")
best.mod[which(best.mod$pheno==ph),2]=mod.names[k]
print(plot(vars.t$N, vars.t$N_sign_vars, pch=20, cex=2, main=paste(ph, ", best model is ",mod.names[k],sep="")))
print(lines(x,y, col="red"))
vars.pred=rbind(vars.pred,vars.pred.t)
}
dev.off()

pls=c("pct","plt","pdw","mpv")
reds=c("mch","mchc","mcv","rbc","hct","hgb","hlr","hlr_p","irf","ret","ret_p")
phils=c("eo_p","eo","mono","mono_p","neut","neut_p","baso","baso_p")
lymph=c("lymph_p","lymph","wbc")

vars.pred=vars.pred[!is.na(vars.pred$pheno),]
vars.pred$pheno=factor(vars.pred$pheno,levels=c(lymph,phils,pls,reds))
vars$pheno=factor(vars$pheno,levels=c(lymph,phils,pls,reds))
ggplot(data=vars, aes(x=N/10000, y=N_sign_vars))+geom_line()+
  geom_point()+facet_wrap(~pheno,ncol=6)+
  geom_line(data=vars.pred, aes(x=x/10000, y=y),col="red")+
  theme(axis.text=element_text(size=15), strip.text.x=element_text(size=15))+xlab("")+ylab("")

# plot best fitted model for N.loci
loci=read.table("Table_N_loci.txt", he=T, strings=F)
loci.pred=data.frame(pheno=NA, x=NA, y=NA)
best.mod2=data.frame(pheno=pheno, model=NA)
pdf("Best_model_for_N_of_loci_per_pheno.pdf")
for(ph in pheno){
  loci.t=loci[which(loci$pheno==ph),]
  #loci.t$N=c(30000,43500,83000,400000)
  mod1=lm(N_loci~N, data=loci.t)
  mod2=lm(N_loci~sqrt(N), data=loci.t)
  mod3=lm(N_loci~sqrt(N)+N, data=loci.t)
  mod4=lm(N_loci~log(N), data=loci.t)
  mod5=lm(N~N_loci^2+N_loci, data=loci.t)
  k=which.max(c(summary(mod1)$r.squared,summary(mod2)$r.squared,summary(mod3)$r.squared,summary(mod4)$r.squared,summary(mod5)$r.squared))
  
  x=seq(1,400000,10000)
  y=predict(get(paste("mod",k,sep="")), newdata=data.frame(N=x))
  loci.pred.t=data.frame(pheno=ph, x=x, y=y)
  mod.names=c("linear", "y~sqrt(x)", "y~sqrt(x)+x", "logarithmic", "quadratic in y")
  best.mod2[which(best.mod$pheno==ph),2]=mod.names[k]
  print(plot(loci.t$N, loci.t$N_loci, pch=20, cex=2, main=paste(ph, ", best model is ",mod.names[k],sep="")))
  print(lines(x,y, col="red"))
  loci.pred=rbind(loci.pred,loci.pred.t)
}
dev.off()

loci.pred=loci.pred[which(!is.na(loci.pred$pheno)),]
loci.pred$pheno=factor(loci.pred$pheno,levels=c(lymph,phils,pls,reds))
loci$pheno=factor(loci$pheno,levels=c(lymph,phils,pls,reds))
ggplot(data=loci, aes(x=N/10000, y=N_loci))+geom_line()+
  geom_point()+facet_wrap(~pheno,ncol=6)+
  geom_line(data=loci.pred, aes(x=x/10000, y=y),col="red")+
  theme(axis.text=element_text(size=15), strip.text.x=element_text(size=15))+xlab("")+ylab("")

# plot best fit model for N. genes 

geni=read.table("Table_N_genes.txt", he=T, strings=F)
geni.pred=data.frame(pheno=NA, x=NA, y=NA)
best.mod3=data.frame(pheno=pheno, model=NA)
pdf("Best_model_for_N_of_geni_per_pheno.pdf")
for(ph in pheno){
  geni.t=geni[which(geni$pheno==ph),]
  geni.t$N=c(30000,43500,83000,400000)
  mod1=lm(N_genes~N, data=geni.t)
  mod2=lm(N_genes~sqrt(N), data=geni.t)
  mod3=lm(N_genes~sqrt(N)+N, data=geni.t)
  mod4=lm(N_genes~log(N), data=geni.t)
  mod5=lm(N~N_genes^2+N_genes, data=geni.t)
  k=which.max(c(summary(mod1)$r.squared,summary(mod2)$r.squared,summary(mod3)$r.squared,summary(mod4)$r.squared,summary(mod5)$r.squared))
  
  x=seq(1,400000,10000)
  y=predict(get(paste("mod",k,sep="")), newdata=data.frame(N=x))
  geni.pred.t=data.frame(pheno=ph, x=x, y=y)
  mod.names=c("linear", "y~sqrt(x)", "y~sqrt(x)+x", "logarithmic", "quadratic in y")
  best.mod3[which(best.mod$pheno==ph),2]=mod.names[k]
  print(plot(geni.t$N, geni.t$N_genes, pch=20, cex=2, main=paste(ph, ", best model is ",mod.names[k],sep="")))
  print(lines(x,y, col="red"))
  geni.pred=rbind(geni.pred,geni.pred.t)
}
dev.off()


geni.pred=geni.pred[which(!is.na(geni.pred$pheno)),]
geni.pred$pheno=factor(geni.pred$pheno,levels=c(lymph,phils,pls,reds))
geni$pheno=factor(geni$pheno,levels=c(lymph,phils,pls,reds))
ggplot(data=geni, aes(x=N/10000, y=N_genes))+geom_line()+
  geom_point()+facet_wrap(~pheno,ncol=6)+
  geom_line(data=geni.pred, aes(x=x/10000, y=y),col="red")+
  theme(axis.text=element_text(size=15), strip.text.x=element_text(size=15))+xlab("")+ylab("")


### PLOT VARS by freq 
vars=read.table("Table_N_variants_sign_by_MAF.txt", he=T, strings=F)
head(vars)
vars$N=rep(c(30000,43500,83000,400000),26)
meltData=melt(vars[,c("cohort","pheno","N_common","N_low_freq","N_rare")])
head(meltData)
meltData$N=rep(c(30000,43500,83000,400000),78)
meltData$cohort=as.factor(meltData$cohort)
meltData$pheno=factor(meltData$pheno,levels=c(lymph,phils,pls,reds))
meltData.t=meltData[which(meltData$pheno=="mpv"),]
meltData.t=meltData[which(meltData$variable!="N_common"),]

ggplot(meltData.t, aes(x=N, y=value, fill=variable))+geom_area()+facet_wrap(~pheno,ncol=6)

ggplot(meltData, aes(x=N, y=value, fill=variable))+geom_area()+facet_wrap(~pheno,ncol=6)

pheno="mpv"
meltData.t=meltData[which(meltData$pheno=="mpv"),]
ggplot(meltData.t, aes(x=N, y=value, fill=variable))+geom_area()+geom_vline(xintercept = 30000)+geom_vline(xintercept = 43500)+geom_vline(xintercept = 83000)+ggtitle("MPV")

head(meltData.t)
meltData.t$variable[which(is.na(meltData.t$variable))]="N_common"
meltData.t$value[which(meltData.t$variable=="N_common")]=meltData.t$value[which(meltData.t$variable=="N_common")]/1922
meltData.t$value[which(meltData.t$variable=="N_low_freq")]=meltData.t$value[which(meltData.t$variable=="N_low_freq")]/131
meltData.t$value[which(meltData.t$variable=="N_rare")]=meltData.t$value[which(meltData.t$variable=="N_rare")]/67

meltData.t
ggplot(meltData.t, aes(x=N, y=value, fill=variable))+geom_area()+geom_vline(xintercept = 30000)+geom_vline(xintercept = 43500)+geom_vline(xintercept = 83000)+ggtitle("MPV proportions")
head(vars)
vars$N_common_p=vars$N_common/vars$N_sign_vars
vars$N_low_freq_p=vars$N_low_freq/vars$N_sign_vars
vars$N_rare_p=vars$N_rare/vars$N_sign_vars

meltData.t=melt(vars[which(vars$pheno=="mpv"), c("pheno","cohort","N_common_p", "N_low_freq_p","N_rare_p")])
head(meltData.t)
meltData.t$N=rep(c(30000,43500,83000,400000),3)
ggplot(meltData.t, aes(x=N, y=value, fill=variable))+geom_area()+geom_vline(xintercept = 30000)+geom_vline(xintercept = 43500)+geom_vline(xintercept = 83000)+ggtitle("MPV proportions")


###
### further attempts with mendelian genes and coexpression ###
mpv$max_coex_mend=NA
mpv$closest_mend=NA
unlist(dimnames(mat)[1])->geni
for(i in 1:dim(mpv)[1]){
g=unlist(strsplit(mpv$gene[i],","))
if(length(intersect(g, geni))>0){
  g=intersect(g, geni)
  b=intersect(bridge, geni)
subs=coex[g, b]
mpv$max_coex_mend[i]=max(abs(subs))
if(length(g)>1)
mpv$closest_mend[i]=row.names(which(abs(subs)==max(abs(subs)), arr.ind = TRUE))
else
  mpv$closest_mend[i]=names(which.max(abs(subs)))
}
}


lymph$max_coex_mend=NA
lymph$closest_mend=NA
unlist(dimnames(mat)[1])->geni
for(i in 1:dim(lymph)[1]){
  g=unlist(strsplit(lymph$gene[i],","))
  if(length(intersect(g, geni))>0){
    g=intersect(g, geni)
    b=intersect(bridge, geni)
    subs=coex[g, b]
    lymph$max_coex_mend[i]=max(abs(subs))
    if(length(g)>1)
      lymph$closest_mend[i]=row.names(which(abs(subs)==max(abs(subs)), arr.ind = TRUE))
    else
      lymph$closest_mend[i]=names(which.max(abs(subs)))
  }
}

summary(mpv$max_coex_mend)
summary(mpv$max_coex_mend[which(mpv$INT_sign=="yes")])
summary(mpv$max_coex_mend[which(mpv$UKBiLEVE_sign=="yes")])
summary(mpv$max_coex_mend[which(mpv$UKBiobank_sign=="yes")])
summary(mpv$max_coex_mend[which(mpv$new=="yes")])
boxplot(mpv$max_coex_mend~mpv$new)
wilcox.test(mpv$max_coex_mend~mpv$new)

summary(lymph$max_coex_mend)
summary(lymph$max_coex_mend[which(lymph$INT_sign=="yes")])
summary(lymph$max_coex_mend[which(lymph$UKBiLEVE_sign=="yes")])
summary(lymph$max_coex_mend[which(lymph$UKBiobank_sign=="yes")])
summary(lymph$max_coex_mend[which(lymph$new=="yes")])
boxplot(lymph$max_coex_mend~lymph$new)
wilcox.test(lymph$max_coex_mend~lymph$new)

cut.offs=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)
for(k in cut.offs){
  print(k)
  print(table(mpv[which(mpv$max_coex_mend>=k),"new"])[1]/table(mpv[which(mpv$max_coex_mend>=k),"new"])[2])
  print(table(mpv[which(mpv$max_coex_mend<k),"new"])[1]/table(mpv[which(mpv$max_coex_mend<k),"new"])[2])
  print(chisq.test(c(table(mpv[which(mpv$max_coex_mend>=k),"new"]),table(mpv[which(mpv$max_coex_mend<k),"new"]))))
}


### heritability saturation ###
# Parsa's h2 estimates, on HPC
# DIR=/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/2018_11_28_heritability_prediction/rsquared_interval/raw

h2=data.frame(cohort=rep(c("INT","UKBiLEVE","UKBiobank","UKBB"),2), N=rep(c(30000,43500,83000,400000),2), pheno=rep(c("lymph","mpv"),each=4), h2=c(0.05148,0.0566,0.08494,0.143565,0.24733,0.25353,0.30484,0.3612))
plot(h2$N, h2$h2, pch=20, cex=2, cex.axis=1.5)

h2.t=h2[which(h2$pheno=="mpv"),]
mod1=lm(h2~N, data=h2.t)
mod2=lm(h2~sqrt(N), data=h2.t)
mod3=lm(h2~sqrt(N)+N, data=h2.t)
mod4=lm(h2~log(N), data=h2.t)
mod5=lm(N~h2^2+h2, data=h2.t)
k=which.max(c(summary(mod1)$r.squared,summary(mod2)$r.squared,summary(mod3)$r.squared,summary(mod4)$r.squared,summary(mod5)$r.squared))
k
plot(h2.t$N, h2.t$h2, pch=20, cex=2, cex.axis=1.5)


plot(h2.t$N, h2.t$h2, pch=20, cex=2, cex.axis=1.5, xlab="Cohort size",ylab="Heritability explained",cex.lab=1.5,main="MPV")
mod1=lm(h2.t$h2[1:3]~h2.t$N[1:3])
abline(mod1$coefficients[1:2], lty=2)
mod2=lm(h2.t$h2~h2.t$N)
abline(mod2$coefficients, lty=2, col="red")
tmp=sqrt(h2.t$N)
mod3=lm(h2.t$h2~tmp+h2.t$N)
mod3$coefficients
summary(mod3)
x=seq(min(h2.t$N), max(h2.t$N), 10000)
y=mod3$coefficients[2]*sqrt(x)+mod3$coefficients[1]+mod3$coefficients[3]*x
lines(x,y, col="red")


sat$N=NA
sat$N[which(sat$cohort=="INT")]=35000
sat$N[which(sat$cohort=="UKBiLEVE")]=43500
sat$N[which(sat$cohort=="UKBiobank")]=83000
sat$N[which(sat$cohort=="UKBB")]=400000
sat=sat[which(sat$pheno=="mpv"),]
tmp=sqrt(sat$N)
mod3=lm(sat$N_sign_vars~tmp+sat$N)

plot(sat$N, sat$N_sign_vars, pch=20, cex=2, cex.axis=1.5, xlab="Cohort size",ylab="N. of associated variants",cex.lab=1.5,main="MPV")
x=seq(min(sat$N), max(sat$N), 10000)
y=mod3$coefficients[2]*sqrt(x)+mod3$coefficients[1]+mod3$coefficients[3]*x
lines(x,y, col="red")

mod1=lm(sat$N_sign_vars[1:3]~sat$N[1:3])
mod2$coefficients
y=mod1$coefficients[2]*x+mod1$coefficients[1]
points(x,y, type="l", lty=2)

mod2=lm(sat$N_sign_vars~sat$N)
y=mod2$coefficients[2]*x+mod2$coefficients[1]
points(x,y, type="l", lty=2, col="red")

### plot MPV example for manuscript figure 
sat=read.table("/Users/dv3/Desktop/Discovery_saturation/Table_N_variants_sign.txt", he=T, strings=F)


### want to repeat the saturation analysis restricting to mendelian genes, to see if they are saturing (more wrt other genes)
## I need to get the number of sign vars, and genes for each cohort and across all phenotypes 

# saturation curves between core and peripheral genes
setwd("/Users/dv3/Desktop/Discovery_saturation/")

geni=read.table("Table_N_genes.txt", he=T, strings=F)
head(geni)
geni$N_mend_genes=NA  
mendelian=scan("../Haem_disease/BRIDGE_flagship_release/BRIDGE_Flagship_gene_list.txt","char")
vars=read.table("Table_N_variants_sign_by_MAF.txt", he=T, strings=F)
head(vars)
vars$N_mend_vats=NA

for(ph in unique(geni$pheno)){
  res=read.table(paste("discovery_saturation_tables/",ph,"_sign_genes_per_cohort.txt",sep=""),he=T, strings=F,sep="\ ")
  geni[which(geni$cohort=="UKBB" & geni$pheno==ph),"N_mend_genes"]=length(intersect(unlist(strsplit(res$gene[which(res$mendelian=="yes")],",")), mendelian))
  vars[which(vars$cohort=="UKBB" & vars$pheno==ph),"N_mend_vars"]=length(which(res$mendelian=="yes"))
  for(c in unique(geni$cohort)[1:3]){
    tmp=res[which(res[,paste(c,"_sign",sep="")]=="yes"),]
    geni[which(geni$cohort==c & geni$pheno==ph),"N_mend_genes"]=length(intersect(unlist(strsplit(tmp$gene[which(tmp$mendelian=="yes")],",")), mendelian))
    vars[which(vars$cohort==c & vars$pheno==ph),"N_mend_vars"]=length(which(tmp$mendelian=="yes"))
  }
}

library(ggplot2)
geni$N=rep(c(33501, 43597, 83637, 391598),26)
ggplot(data=geni, aes(x=N/10000, y=N_mend_genes))+geom_line()+
  geom_point()+facet_wrap(~pheno,ncol=6)+
  theme(axis.text=element_text(size=15), strip.text.x=element_text(size=15))+xlab("")+ylab("")

  geom_line(data=geni, aes(x=N/10000, y=N_genes),col="red")
  
dim(geni)
delta_geni=geni$N_genes[which(geni$cohort=="UKBB")]-geni$N_genes[which(geni$cohort=="UKBiobank")]
delta_geni=delta_geni/geni$N_genes[which(geni$cohort=="UKBB")]*100
delta_mend=geni$N_mend_genes[which(geni$cohort=="UKBB")]-geni$N_mend_genes[which(geni$cohort=="UKBiobank")]
delta_mend=delta_mend/geni$N_mend_genes[which(geni$cohort=="UKBB")]*100

# the above is not too convincing as the scale is different
# here I want to see % of mendelian vs non-mendelian genes detected per cohort and overall (keeping all the traits together)
setwd("/Users/dv3/Desktop/Discovery_saturation/")
res=read.table("discovery_saturation_tables/all_traits_genes_per_cohorts.txt", he=T, strings=F, sep=" ")
dim(res)
res=res[!duplicated(res$VAR),]

fin=data.frame(cohort=c("INT","UKBiLEVE","UKBiobank","UKBB"), N=c(33501, 43597, 83637, 391598), tot_genes=NA, mend_genes=NA, other_genes=NA, N_vars_in_mend=NA, N_vars_in_other=NA)
for(c in fin$cohort){
  if(c=="UKBB")
    tmp=res
  else
  tmp=res[which(res[,paste(c,"_sign",sep="")]=="yes"),]
  fin[which(fin$cohort==c), "N_vars_in_mend"]=table(tmp$mendelian)["yes"]
  fin[which(fin$cohort==c), "N_vars_in_other"]=table(tmp$mendelian)["no"]
  fin[which(fin$cohort==c), "other_genes"]=length(unique(tmp$gene[which(tmp$mendelian=="no")]))
  fin[which(fin$cohort==c), "mend_genes"]=length(intersect(unlist(strsplit(tmp$gene,",")), mendelian))        
    }

fin
fin$mend_genes_p=fin$mend_genes/(fin$mend_genes+fin$other_genes)*100
fin$N_vars_in_mend_p=fin$N_vars_in_mend/(fin$N_vars_in_mend+fin$N_vars_in_other)*100


library(reshape2)
meltData=melt(fin[,c(1,4,5)])
head(meltData)

ggplot(meltData, aes(x=cohort, y=value, fill=variable))+geom_bar(stat="identity", position="dodge")

pdf("Overall_saturation_of_mendelian_genes.pdf", height=7, width=7)
par(mar=c(4.1, 5.6, 5, 5.6))
plot(x=fin$N/10000, y=fin$other_genes, pch=20, xlim=c(0,50),ylim=c(0,5196), xlab="sample size x 10,000", ylab="N. of peripheral genes", cex.axis=1.4, cex.lab=1.4, main="Saturation of core vs peripheral genes")
tmp=sqrt(fin$N)
mod_o=lm(fin$other_genes~tmp+fin$N)
x=seq(0, 500000, 10000)
y=mod_o$coefficients[2]*sqrt(x)+mod_o$coefficients[1]+mod_o$coefficients[3]*x
lines(x/10000,y)

par(new=T)
plot(x=c(fin$N/10000,50), y=c(fin$mend_genes,135) , col=c("red", "red","red","red","white"),pch=1, ylab="", axes=F,ylim=c(0,135), xlim=c(0,50), xlab="")
axis(side = 4,cex.axis=1.4, col.axis="red")
mod_m=lm(fin$mend_genes~tmp+fin$N)
x=seq(0, 500000, 10000)
y=mod_m$coefficients[2]*sqrt(x)+mod_m$coefficients[1]+mod_m$coefficients[3]*x
lines(x/10000,y, col="red")
mtext(side = 4, line = 3, 'N. of Mendelian genes', cex=1.4, col="red")
dev.off()


