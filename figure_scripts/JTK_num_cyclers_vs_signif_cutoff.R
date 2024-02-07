library(tidyverse)
library(gridExtra)
setwd("~/Box Sync/Henry_stuff/AD_project/mouse_data/Dec13_2023_batch_corrected_data/JTK_output/")
ipbulk = read_csv("JTKresult_bulk_WT_data.csv")

get_num_genes_helper = function(df, threshold, amp_ratio_cutoff = 0.1){
  df %>% dplyr::filter(BH.Q < threshold & ampRatio > amp_ratio_cutoff) %>% nrow
}

get_num_genes = function(filename){
  file = read_csv(filename)
  file$mean = rowMeans(file[,-c(1:6)])
  file$ampRatio = file$AMP/file$mean
  trace = sapply(seq(0,1,by = 0.001), get_num_genes_helper, df = file)
}

WT_bulk = get_num_genes("JTKresult_bulk_WT_data.csv")
WT_astro = get_num_genes("JTKresult_astro_wt_data.csv")
WT_mglia = get_num_genes("JTKresult_micro_WT_data.csv")
APP_bulk = get_num_genes("JTKresult_bulk_APP_data.csv")
APP_astro = get_num_genes("JTKresult_astro_APP_data.csv")
APP_mglia = get_num_genes("JTKresult_micro_APP_data.csv")
plot_df = data.frame(WT_bulk = WT_bulk, WT_astro = WT_astro, WT_mglia = WT_mglia,
                     APP_bulk = APP_bulk, APP_astro = APP_astro, APP_mglia = APP_mglia,
                     axis = -log10(seq(0,1,by = 0.001)))

colors <- c("Bulk" = "darkgreen", "Astrocyte" = "blue", "Microglia" = "darkred")

p1 = ggplot(plot_df)+
  geom_line(size = 1, aes(axis, WT_bulk, color = "Bulk"))+
  geom_line(size = 1, aes(axis, WT_astro, color = "Astrocyte"))+
  geom_line(size = 1, aes(axis, WT_mglia, color = "Microglia"))+
  xlab("-Log10(significance cutoff)")+
  ylab("Number of Cycling Transcripts")+
  theme_minimal()+
  scale_color_manual(name = "Cell Type", values = colors)+
  ylim(c(0, 4e3))
type = "dotted"
p2 = ggplot(plot_df)+
  geom_line(size = 1, linetype = type, aes(axis, APP_bulk, color = "Bulk"))+
  geom_line(size = 1, linetype = type, aes(axis, APP_astro, color = "Astrocyte"))+
  geom_line(size = 1, linetype = type, aes(axis, APP_mglia, color = "Microglia"))+
  xlab("-Log10(significance cutoff)")+
  ylab("Number of Cycling Transcripts")+
  theme_minimal()+
  scale_color_manual(name = "Cell Type", values = colors)+
  ylim(c(0, 4e3))

grid.arrange(p1, p2, nrow = 1)
