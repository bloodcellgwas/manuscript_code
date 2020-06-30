#!/software/R-3.2.2/bin/Rscript
library(ggplot2)
OUT_DIR <- Sys.getenv("OUT_DIR")
PHENO_NAMES <- Sys.getenv("PHENO_NAMES")
traits <- as.vector(read.table(PHENO_NAMES, header=F, stringsAsFactors=F)$V1)
# load user friendly file name function
# and user friendly trait name function
source("general_functions.R")

final_output = read.csv(paste0(OUT_DIR, "/final_output_raw.csv"), header=T, stringsAsFactors=F)
eqtl_coloc = read.csv(paste0(OUT_DIR, "/external_data/2019_07_10_ld_signif_eqtl_coloc_results.csv"), header=T, stringsAsFactors=F)

final_output = final_output[, c("VARIANT", "VEP_GENE_SYMBOL", "trait")]
#final_output$dup = paste(final_output$VARIANT, final_output$VEP_GENE_SYMBOL)
#final_output = final_output[!duplicated(final_output$dup),]

eqtl_coloc = eqtl_coloc[, c("condsig_var", "hgnc", "eqtl_celltype")]
eqtl_coloc$dup = paste(eqtl_coloc$condsig_var, eqtl_coloc$hgnc)
eqtl_coloc = eqtl_coloc[!duplicated(eqtl_coloc$dup),]

eqtl_tested = read.csv(paste0(OUT_DIR, "/external_data/2019_07_20_genes_tested.csv"))
final_output = final_output[get_traits_in_eqtl(final_output$trait),]

## does the eQTL annotation disagree with VEP annotation?
agree = 0
disagree = 0
final_output$eqtl_agree = NA
final_output$eqtl_hgnc = NA
for (rowi in 1:nrow(final_output)) {
    final_output_s = final_output[rowi,]
    ## check if gene was tested in eQTL data
    eqtl_coloc_s = eqtl_coloc[eqtl_coloc$condsig_var == final_output_s$VARIANT & eqtl_coloc$eqtl_celltype %in% get_trait_eqtl_type(final_output_s$trait),]
    #if (final_output_s$VEP_GENE_SYMBOL %in% eqtl_tested[eqtl_tested$eqtl_celltype %in% get_trait_eqtl_type(final_output_s$trait),]$hgnc  == FALSE) {next}
    if (nrow(eqtl_coloc_s) == 0) {next}
    final_output$eqtl_hgnc[rowi] = paste0(eqtl_coloc_s$hgnc, collapse=',')
    if (final_output_s$VEP_GENE_SYMBOL %in% eqtl_coloc_s$hgnc) {
        agree = 1 + agree
        final_output[rowi,]$eqtl_agree = TRUE
    } else {
        disagree = disagree + 1
        final_output[rowi,]$eqtl_agree = FALSE
    }
}
final_output = final_output[!is.na(final_output$eqtl_hgnc),]
final_output = final_output[!is.na(final_output$VEP_GENE_SYMBOL),]
final_output = final_output[order(final_output$VEP_GENE_SYMBOL),]
write.csv(final_output, paste0(OUT_DIR, "/paper_figures/02_eqtl_vs_vep.csv"), row.names=F)
stop('a')
scale_fill_Publication <- function(...){
  library(scales)
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  values = col_vector
  # just want to switch them because i prefer the order this way round (not important for actual function)
  values[1] = "#91bfdbff"
  values[2] = "#a50026ff"
  discrete_scale("fill","Publication",manual_pal(values = values), ...)
  # c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
  
}

theme_publication_specific <- function(base_size=14, base_family="helvetica") {
        library(grid)
              library(ggthemes)
              (theme_foundation(base_size=base_size, base_family=base_family)
                      + theme(plot.title = element_text(face = "bold", size = rel(0.6), hjust = 0.5),
                                             text = element_text(),
                                             panel.background = element_rect(colour = NA),
                                             plot.background = element_rect(colour = NA),
                                             panel.border = element_rect(colour = NA),
                                             axis.title = element_text(face = "bold",size = rel(1.2)),
                                             axis.title.y = element_text(angle=90,vjust =2),
                                             axis.title.x = element_text(vjust = -0.2),
                                             axis.text.x = element_text(size = rel(1.2)),
                                             axis.text.y = element_text(size = rel(1.2)),
                                             axis.line = element_line(colour="black"),
                                             axis.ticks = element_line(),
                                             panel.grid.major = element_line(colour="#f0f0f0"),
                                             panel.grid.minor = element_blank(),
                                             legend.key = element_rect(colour = NA),
                                             legend.position = "bottom",
                                             legend.direction = "horizontal",
                                             legend.key.size= unit(rel(0.5), "cm"),
                                             legend.margin = unit(0, "cm"),
                                             legend.title = element_text(face="italic"),
                                             legend.text=element_text(size=rel(1)),
                                             plot.margin=unit(c(10,5,5,5),"mm"),
                                             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                             strip.text = element_text(face="bold")
                                        ))

        }
final_output$trait = user_friendly_trait_name(final_output$trait)
##final_output$trait = as.factor(final_output$trait)
final_output = final_output[!is.na(final_output$eqtl_agree),]

## add percentages
for (trait in unique(final_output$trait)) {
    final_output_s = final_output[final_output$trait == trait,]
    pct = sum(final_output_s$eqtl_agree)/nrow(final_output_s)*100
    #final_output$trait[final_output$trait == trait] = rep(paste0(trait, " (", signif(pct, 2), "%)"), sum(final_output$trait == trait))
}

final_output$eqtl_agree[final_output$eqtl_agree == TRUE] = " eQTL Agree "
final_output$eqtl_agree[final_output$eqtl_agree == FALSE] = " eQTL Disagree "
p <- ggplot(data = final_output, aes(x = trait, fill = eqtl_agree)) +
  geom_bar(width=0.4) +
  scale_fill_Publication() +
  theme_publication_specific() +
  theme(legend.title=element_blank()) +
  labs(y = "Associated Signals", x="Trait") +
  coord_flip()

#svglite(file=paste0(OUT_DIR, "/paper_figures/01_bar_chart.svg"), width=6,height=8)
png(paste0(OUT_DIR, "/paper_figures/02_eqtl_vs_vep.png"), width=432, height=576)
print(p)
dev.off()

