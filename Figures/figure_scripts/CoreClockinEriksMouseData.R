#Jan 2 2024, this script takes JTK results from Erik's mouse data and plots
# the acrophases of (core clock) genes in WT and APP mice. These are not the clock
# face, radial plots used in the figure. That comes from the python/jupyter notebook file
# called plot_mouse_core_clock_genes.ipynb
library(gridExtra)
library(readxl)
core_clock = c("Ciart", "Arntl", "Npas2", "Clock", "Cry1",  "Cry2",  "Nr1d1", "Nr1d2", "Per1",  "Per2",  "Per3"  ,"Dbp",  "Tef" )
#setwd("~/Box Sync/Henry_stuff/AD_project/mouse_data/custom_JTK_24period/")
setwd("~/Box Sync/Henry_stuff/AD_project/mouse_data/Dec13_2023_batch_corrected_data/JTK_output/")
ipbulk = read_csv("JTKresult_bulk_WT_data.csv") %>% 
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez")) %>% 
  filter(BH.Q < .2)
astros = read_csv("JTKresult_astro_wt_data.csv")%>% 
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>% 
  filter(BH.Q < .2)
mglia = read_csv("JTKresult_micro_WT_data.csv")%>% 
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez"))%>% 
  filter(BH.Q < .2)

ast_peaks = astros$LAG[match(core_clock,astros$Gene_Symbols)]
ip_peaks = ipbulk$LAG[match(core_clock, ipbulk$Gene_Symbols)]
mglia_peaks = mglia$LAG[match(core_clock, mglia$Gene_Symbols)]

df = data.frame(genenames = core_clock, Astrocyte = ast_peaks, Microglia = mglia_peaks,
           Bulk = ip_peaks)
df = df[-c(3, 4, 6, 8),]
WT_core_clock_df = df
plot_df = df %>% pivot_longer( c(2, 3, 4), names_to = "Cell Type", values_to = "Acrophase")
p1 = ggplot(plot_df, aes(x = genenames, y = Acrophase, fill =`Cell Type` ))+
  geom_point(colour="black",pch=21, size = 5,stroke = 1, position = position_jitter(width = 0.15))+
  theme_minimal()+
  xlab("Core Clock Genes")+
  ylab("Acrophase (hr)")
  

all_sig_cyclers = union(union(ipbulk$CycID, astros$CycID) , mglia$CycID)
ast_peaks = astros$LAG[match(all_sig_cyclers, astros$CycID)]
ip_peaks = ipbulk$LAG[match(all_sig_cyclers, ipbulk$CycID)]
mglia_peaks = mglia$LAG[match(all_sig_cyclers, mglia$CycID)]
df2 = data.frame(genenames = all_sig_cyclers, Astrocyte = ast_peaks, Microglia = mglia_peaks,
                Bulk = ip_peaks)
keep = which(rowSums(!is.na(df2[,-1])) > 1) #keep only genes which cycle in >1 cell type
df2 = df2[keep, ]
df2 = df2[-which(df2$genenames %in% df$genenames), ]
plot_df2 = df2 %>% pivot_longer( c(2, 3, 4), names_to = "Cell Type", values_to = "Acrophase")
p2 = ggplot(plot_df2, aes(x = genenames, y = Acrophase, fill =`Cell Type` ))+
  geom_point(colour="black",pch=21,size = 5, position = position_jitter(width = 0.15))+
  theme_minimal()+
  xlab("Rhythmic Genes")+
  ylab("Acrophase (hr)")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))

grid.arrange(p1, p2)

############# APP mice ########3
ipbulk = read_csv("JTKresult_bulk_APP_data.csv") %>% 
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez")) %>% 
  filter(BH.Q < .2)
astros =  read_csv("JTKresult_astro_APP_data.csv") %>% 
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez")) %>% 
  filter(BH.Q < .2)
mglia =  read_csv("JTKresult_micro_APP_data.csv") %>% 
  separate_wider_delim(1, "_", names = c("Ensembl", "Gene_Symbols", "Entrez")) %>% 
  filter(BH.Q < .2)

ast_peaks = astros$LAG[match(core_clock,astros$Gene_Symbols)]
ip_peaks = ipbulk$LAG[match(core_clock, ipbulk$Gene_Symbols)]
mglia_peaks = mglia$LAG[match(core_clock, mglia$Gene_Symbols)]

df = data.frame(genenames = core_clock, Astrocyte = ast_peaks, Microglia = mglia_peaks,
                Bulk = ip_peaks)
df = df[-c(3, 4,6,11),]
APP_core_clock_df = df
plot_df = df %>% pivot_longer( c(2, 3, 4), names_to = "Cell Type", values_to = "Acrophase")
p1 = ggplot(plot_df, aes(x = genenames, y = Acrophase, fill =`Cell Type` ))+
  geom_point(colour="black",pch=21, size = 5,stroke = 1, position = position_jitter(width = 0.15))+
  theme_minimal()+
  xlab("Core Clock Genes")+
  ylab("Acrophase (hr)")
WT_core_clock_df$treatment = 0
APP_core_clock_df$treatment = 1
write.table(rbind(WT_core_clock_df, APP_core_clock_df), "~/Desktop/WT_APP_core_clock.csv",sep =',',row.names = F)
        