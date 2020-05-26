library(tidyverse)
library(data.table)

# Read all files and filter for PP01
all <- lapply(list.files("../../output/variant_predictions_svm/",full.names=T,pattern="-try1_FMall*."),function(name){
  df <- fread(cmd=paste0("zcat < ",name)) %>% setNames(.,c("ID","deltaSVM")) %>% unique()
  df$celltype <- gsub("../../output/variant_predictions_svm/","",name) %>% gsub("-try1_FMall.out.gz","",.) %>% gsub("-try1_FMall_v2.out.gz","",.)
  
  splits <- str_split_fixed(df$ID,"-",4)
  df$chr <- str_split_fixed(splits[,1],":",2)[,1]
  df$var <- paste0(gsub("chr","",splits[,1]),"_",splits[,2])
  df$trait <- splits[,3]
  df$PP <- splits[,4]
  df %>% dplyr::select(var,chr,celltype,deltaSVM) %>% unique() 
}) %>% bind_rows()