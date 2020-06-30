library(data.table)
library(ggplot2)
library(qqman)
library(svglite)
source(paste0(Sys.getenv('SCRIPT_DIR'), "/general_functions.R"))
OUT_DIR = Sys.getenv('OUT_DIR')
traits = read.csv(Sys.getenv('PHENO_NAMES'), stringsAsFactors=F, header=F)$V1

final_output = read.csv(paste0(OUT_DIR, "/final_output_raw.csv"), header=T, stringsAsFactors=F)
first_batch = c("ne_ssc_ch", "ne_sfl_ch", "ne_fsc_ch", "ig_10_9_l", "DF_EOSI_X", "DF_EOSI_Y", "DF_EOSI_Z", "WN_BASO_X", "WN_BASO_Y", "mo_x_ch", "mo_y_ch", "mo_z_ch", "RE_LYMP103uL", "ly_x_ch", "ly_y_ch", "ly_z_ch", "rpi", "hfr_pct", "lfr_pct")
second_batch = c("RET_Xch", "hyper_he_pct", "irf_y_ch", "macror_pct", "micror_pct", "rbc_he_pg", "ret_he_pg", "RET_RBC_Zch", "RET_RBC_Xch", "ret_rbc_y_ch", "ret_y_ch", "ipfx_10_9_l", "p_lcr_pct", "PLT_F_Zch", "PLT_F_Xch", "PLT_F_Ych", "ne_ssc_ch", "ne_sfl_ch", "ne_fsc_ch")
if (sum(!is.element(first_batch, final_output$trait)) > 0) {stop("asfs")}
final_output = final_output[final_output$trait %in% second_batch,]

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
                                             legend.key.size= unit(rel(1), "cm"),
                                             legend.margin = unit(0, "cm"),
                                             legend.title = element_text(face="italic"),
                                             legend.text=element_text(size=rel(1.2)),
                                             plot.margin=unit(c(10,5,5,5),"mm"),
                                             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                             strip.text = element_text(face="bold")
                                        ))

        }
final_output$novel_clump[final_output$novel_clump == 1] = "Novel"
final_output$novel_clump[final_output$novel_clump == 0] = "Not Novel"
final_output$novel_clump = as.factor(final_output$novel_clump)
final_output$trait = gsub("WDF-", "", gsub("WN-", "", user_friendly_trait_name(final_output$trait)))
final_output$trait = as.factor(final_output$trait)
p <- ggplot(data = final_output, aes(x = trait, fill = novel_clump)) +
  geom_bar(width=0.4) +
  scale_fill_Publication() +
  theme_publication_specific() +
  theme(legend.title=element_blank()) + 
  labs(y = "Associated Loci", x="Cell Type") +
  coord_flip()

svglite(file=paste0(OUT_DIR, "/paper_figures/24_bar_chart.svg"), width=6,height=8)
print(p)
dev.off()

