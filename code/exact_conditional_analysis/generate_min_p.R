#!/software/R-3.1.2/bin/Rscript
library(data.table)
options(scipen=999)
options(warn=2)
OUT_DIR <- Sys.getenv('OUT_DIR')
PHENO_NAMES <- Sys.getenv('PHENO_NAMES')
traits=as.vector(read.table(PHENO_NAMES, header=F, stringsAsFactors=FALSE))$V1

data=NULL
for(trait in traits)
{
  	local_data=fread(paste(OUT_DIR, "/assoc_files/", trait, "_gwsig.assoc", sep=""), select=c("VARIANT", "ID","MLOG10P","MA_FREQ", "INFO", "CHR_num"), header=T)
        print(paste("loaded", trait))
	setkey(local_data,"VARIANT", "ID")
        # if this is the first trait we are iterating over:
	if(trait==traits[1])
	{
		data=local_data
		setnames(data, "MLOG10P","MLOG10P_MAX")      
	}
	else
	{
		data=merge(data, local_data[,c("VARIANT","ID","MLOG10P", "INFO", "CHR_num"),with=FALSE], by=c("VARIANT", "ID", "INFO", "CHR_num"), all.x=TRUE, all.y=TRUE)
                print(paste("merged", trait))
		data$MLOG10P_MAX=pmax(data$MLOG10P_MAX, data$MLOG10P,  na.rm=TRUE)
		data$MA_FREQ[is.na(data$MA_FREQ)]=local_data$MA_FREQ[match(data$VARIANT[is.na(data$MA_FREQ)],local_data$VARIANT)]
        	data[,"MLOG10P":=NULL]
	}
        print(paste("generate_min_p.R done", trait))


        data$sort_order=as.numeric(gsub("(.+):(.+)_(.+)_(.+)","\\2",data$VARIANT))+10^9*as.numeric(data$CHR_num)
	setkey(data, "sort_order")
        data[,sort_order:=NULL]
	data$CHR=gsub("(.+):(.+)_(.+)_(.+)","\\1",data$VARIANT)
	data$BP=gsub("(.+):(.+)_(.+)_(.+)","\\2",data$VARIANT)
	setcolorder(data, c("VARIANT", "ID", "CHR", "CHR_num", "BP", "MA_FREQ","MLOG10P_MAX", "INFO"))
        
	data$GWSIG=data$MLOG10P_MAX>9-log10(8.31)
}
write.table(data, file=paste(OUT_DIR, "/all_min.tsv", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
print("done")
