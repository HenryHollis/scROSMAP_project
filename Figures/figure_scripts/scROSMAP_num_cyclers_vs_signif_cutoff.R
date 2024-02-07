library(tidyverse)
library(gridExtra)
library(pracma)
setwd("~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/TMMs_w_batch/Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_3EGdefault_noTransferFit/")

get_num_genes_helper = function(df, threshold){
  df %>% dplyr::filter(BHQ < threshold) %>% nrow
}

xaxis = logseq(1e-30, 1.001, n = 1000)
get_num_genes = function(filename){
  file = read_csv(filename)
  trace = sapply(xaxis, get_num_genes_helper, df = file)
}
CTL_bulk = get_num_genes("downstream_output_pseudobulkAll/cosinor_results_CTL.csv")
CTL_neuro = get_num_genes("downstream_output_Exc3_5/cosinor_results_CTL.csv")
CTL_astro = get_num_genes("downstream_output_Ast1_2_5_6_7_8/cosinor_results_CTL.csv")
CTL_mglia = get_num_genes("downstream_output_Micro4_5_6_7_8_10_11_13_15_16/cosinor_results_CTL.csv")
AD_bulk = get_num_genes("downstream_output_pseudobulkAll/cosinor_results_AD.csv")
AD_neuro = get_num_genes("downstream_output_Exc3_5/cosinor_results_AD.csv")
AD_astro = get_num_genes("downstream_output_Ast1_2_5_6_7_8/cosinor_results_AD.csv")
AD_mglia = get_num_genes("downstream_output_Micro4_5_6_7_8_10_11_13_15_16/cosinor_results_AD.csv")

plot_df = data.frame(CTL_bulk = CTL_bulk, CTL_neuro = CTL_neuro, CTL_astro = CTL_astro, CTL_mglia = CTL_mglia,
                     AD_bulk = AD_bulk, AD_neuro = AD_neuro, AD_astro = AD_astro, AD_mglia = AD_mglia,
                     axis = -log10(xaxis))

colors <- c("Bulk" = "darkgreen", "Astrocyte" = "blue", "Microglia" = "darkred", "Excitatory Neurons" = "darkorange")

p1 = ggplot(plot_df)+
  geom_line(size = 1, aes(axis, CTL_bulk, color = "Bulk"))+
  geom_line(size = 1, aes(axis, CTL_astro, color = "Astrocyte"))+
  geom_line(size = 1, aes(axis, CTL_mglia, color = "Microglia"))+
  geom_line(size = 1, aes(axis, CTL_neuro, color = "Excitatory Neurons"))+
  xlab("-log10(significance cutoff)")+
  ylab("Number of Cycling Transcripts")+
  theme_minimal()+
  scale_color_manual(values = colors)+
  ylim(c(0, 18e3))
type = "dotted"
p2 = ggplot(plot_df)+
  geom_line(size = 1, linetype = type, aes(axis, AD_bulk, color = "Bulk"))+
  geom_line(size = 1, linetype = type, aes(axis, AD_astro, color = "Astrocyte"))+
  geom_line(size = 1, linetype = type, aes(axis, AD_mglia, color = "Microglia"))+
  geom_line(size = 1, aes(axis, AD_neuro, color = "Excitatory Neurons"))+
  xlab("-log10(significance cutoff)")+
  ylab("Number of Cycling Transcripts")+
  theme_minimal()+
  scale_color_manual(values = colors)+
  ylim(c(0, 18e3))

grid.arrange(p1, p2, nrow = 1)
