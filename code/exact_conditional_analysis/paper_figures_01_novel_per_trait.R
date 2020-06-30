library(data.table)
library(ggplot2)
library(svglite)
source(paste0(Sys.getenv('SCRIPT_DIR'), "/general_functions.R"))
OUT_DIR = Sys.getenv('OUT_DIR')
traits = read.csv(Sys.getenv('PHENO_NAMES'), stringsAsFactors=F, header=F)$V1
final_output = read.csv(paste0(OUT_DIR, "/final_output_raw.csv"), header=T, stringsAsFactors=F)

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
final_output$novel[final_output$novel == "True"] = " Novel"
final_output$novel[final_output$novel == "False"] = " Not Novel  "
final_output$trait = factor(final_output$trait, levels=get_trait_cell_type_order(traits))
final_output$trait = user_friendly_trait_name(final_output$trait)
##final_output$trait = as.factor(final_output$trait)
p <- ggplot(data = final_output, aes(x = trait, fill = novel)) +
  geom_bar(width=0.4) +
  scale_fill_Publication() +
  theme_publication_specific() +
  theme(legend.title=element_blank()) + 
  labs(y = "Associated Signals", x="Trait")# +
#  coord_flip()

#svglite(file=paste0(OUT_DIR, "/paper_figures/01_bar_chart.svg"), width=6,height=8)
png(paste0(OUT_DIR, "/paper_figures/01_bar_chart.png"), width=432, height=576)
print(p)
dev.off()

