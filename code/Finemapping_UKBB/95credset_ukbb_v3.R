#credible sets
args=commandArgs(TRUE)
pheno=args[1]
N_blocks=args[2]
root_dir="/lustre/scratch115/projects/ukbb500k_t151/final_fine_mapping/"
for(j in c(1:N_blocks)){
  setwd(paste(root_dir,pheno,"/finemap_",pheno,"/",sep=""))
  res=try(read.table(paste(pheno,"_block",j,".config",sep=""), he=T, strings=F))
  if(class(res)!="try-error"){
    res$sum=cumsum(res$prob)
    if(max(res$sum, na.rm=T)<0.95) next
    snps=res[1:which(res$sum>=0.95)[1],"config"]
    if(length(snps)==0)
      snps=res[1,""]
    snps=unique(unlist(strsplit(snps, ",")))
    write.table(snps, paste(pheno,"_block",j,".95credSet",sep=""),col.names=F, row.names=F, quote=F)}
  else{
    warning(paste("check ", pheno," block ",j,sep="")) }
}



