library(tidyverse)
library(gridExtra)
ubiquity_cyclers = readxl::read_excel("~/Box Sync/Henry_stuff/AD_project/scROSMAP/Rscripts/figure_scripts/mouse_ubiquity_cyclers.xlsx", col_names = F) %>%
    unname %>% unlist
setwd("~/Box Sync/Henry_stuff/AD_project/mouse_data/Dec13_2023_batch_corrected_data/JTK_output/")

WT_bulk = read_csv("JTKresult_bulk_WT_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>%
  filter(Gene_Symbols %in% ubiquity_cyclers)
WT_astro = read_csv("JTKresult_astro_wt_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>%
  filter(Gene_Symbols %in% ubiquity_cyclers)
WT_mglia = read_csv("JTKresult_micro_WT_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>%
  filter(Gene_Symbols %in% ubiquity_cyclers)

#only keep the genes present in all cell types, otherwise, comment out 'drop_na'
df = data.frame(gene = ubiquity_cyclers, 
           Bulk = WT_bulk$AMP[match(ubiquity_cyclers, WT_bulk$Gene_Symbols)],
           Astrocyte = WT_astro$AMP[match(ubiquity_cyclers, WT_astro$Gene_Symbols)],
           Microglia = WT_mglia$AMP[match(ubiquity_cyclers, WT_mglia$Gene_Symbols)]) %>%
          drop_na()
df = df %>% pivot_longer(!gene, names_to = "tissue", values_to = "Amplitude")

p1 = ggplot(df, aes(x = Amplitude, fill = tissue))+
  geom_histogram(alpha = .7, position= "dodge")+
  theme(legend.position = "top")+
  scale_fill_manual(values = c( "#1F77B4","#2CA02C", "#D62728"), name = "Cell Population:" )+
  xlab("JTK Amplitude")+
  ylab("Count")+
  ylim(c(0, 20))+ 
  xlim(0,100)+
  ggtitle("WT")


APP_bulk = read_csv("JTKresult_bulk_APP_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>%
  filter(Gene_Symbols %in% ubiquity_cyclers)
APP_astro = read_csv("JTKresult_astro_APP_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>%
  filter(Gene_Symbols %in% ubiquity_cyclers)
APP_mglia = read_csv("JTKresult_micro_APP_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>%
  filter(Gene_Symbols %in% ubiquity_cyclers)

#only keep the genes present in all cell types, otherwise, comment out 'drop_na'
df = data.frame(gene = ubiquity_cyclers, 
                Bulk = APP_bulk$AMP[match(ubiquity_cyclers, APP_bulk$Gene_Symbols)],
                Astrocyte = APP_astro$AMP[match(ubiquity_cyclers, APP_astro$Gene_Symbols)],
                Microglia = APP_mglia$AMP[match(ubiquity_cyclers, APP_mglia$Gene_Symbols)]) %>%
                drop_na()
df = df %>% pivot_longer(!gene, names_to = "tissue", values_to = "Amplitude")

p2 = ggplot(df, aes(x = Amplitude, fill = tissue))+
  geom_histogram(alpha = .7, position= "dodge")+
  theme(legend.position = "top")+
  scale_fill_manual(values = c( "#1F77B4","#2CA02C", "#D62728"), name = "Cell Population:" )+
  xlab("JTK Amplitude")+
  ylab("Count")+
  ylim(c(0, 20))+
  xlim(0,100)+
  ggtitle("APP")

grid.arrange(p1, p2, nrow = 1)

#############################################
## What about the cycling genes in common? ##
#############################################

WT_bulk = read_csv("JTKresult_bulk_WT_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez")) %>% filter(BH.Q < 0.2)
WT_astro = read_csv("JTKresult_astro_wt_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez")) %>% filter(BH.Q < 0.2)
WT_mglia = read_csv("JTKresult_micro_WT_data.csv") %>%
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez")) %>% filter(BH.Q < 0.2)

#only keep the genes present in all cell types, otherwise, comment out 'drop_na'
merged = inner_join(WT_bulk, WT_astro, by = "Gene_Symbols", suffix =c("_Bulk", "_Astro")) %>%
  inner_join(WT_mglia, by = "Gene_Symbols") %>% select(Gene_Symbols, AMP_Bulk, AMP_Astro, AMP)
colnames(merged) = c("Gene_Symbols", "Bulk", "Astrocyte", "Microglia")
df = merged %>% pivot_longer(!Gene_Symbols, names_to = "tissue", values_to = "Amplitude")

p1 = ggplot(df, aes(x = Amplitude, fill = tissue))+
  geom_histogram(alpha = .7, position= "dodge")+
  theme(legend.position = "top")+
  scale_fill_manual(values = c( "#1F77B4","#2CA02C", "#D62728"), name = "Cell Population:" )+
  xlab("JTK Amplitude")+
  ylab("Count")+
  ylim(c(0, 20))+ 
  xlim(0,100)+
  ggtitle("WT")

