#####
# Tues Sept 5 2023
# Script for taking a df containing TF_names (from pscan but also enrichR results),
# and checking if they are found in DE, DE, DM, cycling files. 
library(ggrepel)
library(tidyverse)
augment_tf_file = function(TF_filename, Diff_expr_filename, isCyclingBHQCutoff_str, abs_path_to_cyclops_ordering){
  current_dir = getwd()
  TF_file = read.csv(TF_filename)
  if (colnames(TF_file)[1]== "Rank"){ #if file comes from enrichR, make it look like pscan
    TF_file$TF_NAME = TF_file$Term.name
    TF_file$FDR = TF_file$Adjusted.p.value
    TF_file = TF_file[,-1]
    TF_file$TF_NAME = str_replace(TF_file$TF_NAME , " \\(human\\)", "")
    TF_file$TF_NAME = str_replace(TF_file$TF_NAME , " \\(mouse\\)", "")
    
  }
  edgeR_toptags = read.csv(paste(abs_path_to_cyclops_ordering, Diff_expr_filename, sep = '/'))
  DR_AR1_mthd2 = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio1.csv"))
  DR_AR1 = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio1.csv"))
  cycling_CTL = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/cosinor_results_CTL.csv"))
  cycling_AD = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/cosinor_results_AD.csv"))
  mesor_file = list.files(path = paste0(abs_path_to_cyclops_ordering, "/downstream_output") ,pattern = "\\.*mesor.*.csv$")
  Diff_mesor = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/", mesor_file))
    
  #Look for the TF name in the files above
  TF_file$cycling_in_CTL_BHQ = cycling_CTL$BHQ[match(toupper(TF_file$TF_NAME), cycling_CTL$Gene_Symbols)]
  TF_file$cycling_in_AD_BHQ = cycling_AD$BHQ[match(toupper(TF_file$TF_NAME), cycling_AD$Gene_Symbols)]
  TF_file$DR_AR1_BHQ = DR_AR1$BHQ[match(toupper(TF_file$TF_NAME), DR_AR1$Gene_Symbols)]
  TF_file$DR_AR1_mthd2_BHQ = DR_AR1_mthd2$BHQ[match(toupper(TF_file$TF_NAME), DR_AR1_mthd2$Gene_Symbols)]
  TF_file$DR_logAmpRatio = DR_AR1$Log_AD_CTL_ampRatio[match(toupper(TF_file$TF_NAME), DR_AR1$Gene_Symbols)] 
  TF_file$diff_mesor_BHQ = Diff_mesor$BHQ[match(toupper(TF_file$TF_NAME), Diff_mesor$Gene_Symbols)]
  TF_file$edgeR_DE_BHQ = edgeR_toptags$FDR[match(toupper(TF_file$TF_NAME), edgeR_toptags$X)] 
  
  write.table(TF_file, TF_filename, row.names = F, col.names = T, sep = ',')
  DR_tfs = dplyr::filter(TF_file, FDR < 0.1 & DR_AR1_BHQ < 0.3) 
  cycling_CTL_tfs = dplyr::filter(TF_file, FDR < 0.1 & cycling_in_CTL_BHQ < 0.1) 
  cycling_AD_tfs = dplyr::filter(TF_file, FDR < 0.1 & cycling_in_AD_BHQ < 0.1)
  DE_tfs = dplyr::filter(TF_file, FDR < 0.1 & edgeR_DE_BHQ < 0.1) 
  DM_tfs = dplyr::filter(TF_file, FDR < 0.1 & diff_mesor_BHQ < 0.1) 
  
  if(!(dir.exists(paste0(tools::file_path_sans_ext(TF_filename), "_plots")))){
    dir.create(paste0(tools::file_path_sans_ext(TF_filename), "_plots"))
  }
  setwd(paste0(tools::file_path_sans_ext(TF_filename), "_plots"))
  #Check diff rhythms for TFs:
  if (nrow(DR_tfs) > 0){
      p1 = ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(DR_AR1_BHQ)))+
        geom_label_repel(data = DR_tfs,
                         aes(x = -log10(FDR), y = -log10(DR_AR1_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.35, point.padding = 0.5)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(DR BH.q value) (cycling AR>0.1, BH.q<0.3)")+
        ggtitle("Differentially Rhythmic TF enriched in gene list")
      ggsave("DR_TF_enriched_in_list.png", p1, width = 6, height = 5, units = "in")
  }
  #Check genes cycling in CTL for TF:
  if (nrow(cycling_CTL_tfs)>0){
      p2 = ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(cycling_in_CTL_BHQ)))+
        geom_label_repel(data = cycling_CTL_tfs,
                         aes(x = -log10(FDR), y = -log10(cycling_in_CTL_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.35, point.padding = 0.5)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(Cycling in CTL BH.q value)")+
        ggtitle("Cycling in CTL TF enriched in gene list")
      ggsave("cycling_CTL_TF_enriched_in_list.png", p2, width = 6, height = 5, units = "in")
  }
  #Check genes cycling in AD for TF:
  if (nrow(cycling_AD_tfs)>0){
      p3 =  ggplot(TF_file)+
          geom_point(mapping = aes(x = -log10(FDR), y = -log10(cycling_in_AD_BHQ)))+
          geom_label_repel(data = cycling_AD_tfs,
                           aes(x = -log10(FDR), y = -log10(cycling_in_AD_BHQ),label = TF_NAME),
                           nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                           box.padding = 0.35, point.padding = 0.5)+
          xlab("-Log(Pscan/enrichR FDR)")+
          ylab("-Log(Cycling in AD BH.q value)")+
          ggtitle("Cycling in AD TF enriched in gene list")
      ggsave("cycling_AD_TF_enriched_in_list.png", p3, width = 6, height = 5, units = "in")
  }
  #Check Diff Expressed genes for TF:
  if(nrow(DE_tfs)>0){
    p4 =  ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(edgeR_DE_BHQ)))+
        geom_label_repel(data = DE_tfs,
                         aes(x = -log10(FDR), y = -log10(edgeR_DE_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue", 
                         box.padding = 0.35, point.padding = 0.5)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(edgeR Diff. Expr. BH.q value)")+
        ggtitle("EdgeR diff expr TF enriched in gene list")
    ggsave("DE_TF_enriched_in_list.png", p4, width = 6, height = 5, units = "in")
  }
  #Check genes with different mesor for TF:
  if(nrow(DM_tfs)>0){
    p5 =  ggplot(TF_file)+
      geom_point(mapping = aes(x = -log10(FDR), y = -log10(diff_mesor_BHQ)))+
      geom_label_repel(data = DM_tfs,
                       aes(x = -log10(FDR), y = -log10(diff_mesor_BHQ),label = TF_NAME),
                       nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                       box.padding = 0.35, point.padding = 0.5)+
      xlab("-Log(Pscan/enrichR FDR)")+
      ylab("-Log(Diff Mesor BH.q value)")+
      ggtitle("Diff Mesor TF enriched in gene list")
    ggsave("DM_TF_enriched_in_list.png", p5, width = 6, height = 5, units = "in")
  }
  setwd(current_dir)
}
