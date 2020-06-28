#!/software/R-3.3.0/bin/Rscript
library(ggplot2)
library(tools)
source(paste0(Sys.getenv('SCRIPT_DIR'), "/general_functions.R"))
OUT_DIR <- Sys.getenv('OUT_DIR')
traits=as.vector(read.table(Sys.getenv('PHENO_NAMES'), header=F, stringsAsFactors=FALSE)$V1)
disease_table=read.csv(Sys.getenv('DISEASE_NAMES'), header=T, stringsAsFactors=FALSE)


final_output <- read.csv(paste0(OUT_DIR, "/BCX_final_output_raw.csv"), header=T, stringsAsFactors=F)
final_output$novel[final_output$novel == "True"] = TRUE
final_output$novel[final_output$novel == "False"] = FALSE
final_output$novel = as.logical(final_output$novel)
scale_fill_Publication <- function(...){
  library(scales)
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  values = col_vector
  # just want to switch them because i prefer the order this way round (not important for actual function)
  values[2] = col_vector[3]
  values[3] = col_vector[2]
  discrete_scale("fill","Publication",manual_pal(values = values), ...)
  # c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
}

# first plot is per trait class and how many loci for each trait class colocalise
df_disease = c()
df_loci = c()
df_type = c()
for (loci in unique(final_output$clump_id)) {
  final_output.loci = final_output[final_output$clump_id == loci,]

  loci_type = NULL
  if (sum(!as.logical(final_output.loci$novel)) > 0) {
    loci_type = "not_novel"
  } else {
    loci_type = "novel"
  }
  
  loci_diseases = paste(final_output.loci$disease_coloc, collapse=",")
  for (disease in disease_table$full_name) {
    if (grepl(disease, loci_diseases) == TRUE) {
      df_loci = c(df_loci, loci)
      df_disease = c(df_disease, convert_disease_name(disease))
      df_type = c(df_type, loci_type)
    }
  }
}

traits.df <- data.frame(loci=df_loci, disease = df_disease, type=df_type, stringsAsFactors=F)
traits.df$type <- toTitleCase(gsub("_", " ", traits.df$type))
traits.df$disease <- toTitleCase(gsub("_", " ", traits.df$disease))
traits.df$disease <- factor(traits.df$disease, levels = names(table(traits.df$disease))[order(-as.vector(table(traits.df$disease)))])

p <- ggplot(data = traits.df, aes(x = disease, fill=type)) +
  geom_bar()+
  scale_fill_Publication() +
  theme_publication() +
  theme(legend.title=element_blank()) + 
  labs(y = "Associated Loci", x="Disease") +
  coord_flip()

ggsave(paste0(OUT_DIR, "/paper_figures/11_loci_coloc.svg"), plot=p, height=10, width=8)
ggsave(paste0(OUT_DIR, "/paper_figures/11_loci_coloc.png"), plot=p, height=10, width=8)
ggsave(paste0(OUT_DIR, "/paper_figures/11_loci_coloc.pdf"), plot=p, height=10, width=8)
