library(BuenColors)
library(dplyr)
library(data.table)
library(plotROC)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")

ff <- list.files("../../output/kernel_svm/", pattern = "cvPred.out")
lapply(ff[18:25],
       function(x) fread(paste0("../../output/kernel_svm/", x)) %>% data.frame() %>% mutate(cell = gsub("[.]", "-", gsub(".cvPred.out", "", x)))) %>% rbindlist() %>% 
  data.frame() -> dt


p1 <- ggplot(dt, aes(d = V3, m = V2, color = cell)) +
  geom_roc(n.cuts = 0) + pretty_plot() + L_border() +
  labs(x = "1 - Specificity", y = "Sensitivity", color = "Celltype") +
  scale_color_manual(values = eryth_color_maps) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

cowplot::ggsave2(p1, file = "../../output/plots_svm/erythroid-ROCmegaPlot.pdf", height = 5, width = 6)
