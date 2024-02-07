library(ggplot2)
fgsea_dotplot <- function(object, show = 15, 
                          font.size=12, my_title = "GSEA", enrichment = "all", BHQ_cutoff = 0.2) {
  object = as_tibble(object)
  if (enrichment == "pos"){
    print("Including only pos enriched pathways")
    object = dplyr::filter(object, NES > 0)
  }else if (enrichment == "neg"){
    print("Including only neg enriched pathways")
    object = dplyr::filter(object, NES < 0)
  }else{
    print("Using positive and negatively enriched pathways")
  }
  data = object %>% top_n(-1*show, wt = padj)%>% arrange(padj) %>%#-1 means "bottom" n pvalues
    dplyr::filter(padj < BHQ_cutoff)
  #fixed_labels <- str_replace_all(data$pathway, "_", " ")
  data$fixed_labels = str_replace_all(data$pathway, "_", " ")
  ggplot(data, aes(x=reorder(fixed_labels,padj, decreasing = T), y= abs(NES), color = padj, size = size)) +
    geom_point(alpha = 0.9, ) +
    coord_flip()+
    scale_color_continuous(low="red", high="blue", name = "BHQ",
                           guide=guide_colorbar(reverse=T)) +
    #scale_size_continuous( name = "Size", range  =c(3,8)) +
    scale_x_discrete(name = "Pathway", labels = function(x) str_wrap(x, width = 30)) +
    ylab("abs(NES)") +
    labs(title = str_wrap(my_title, 20)) #  # theme_dose(12)
  
}
setwd("~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/TMMs_w_batch/Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_04ContrVar4EGdefault_noTransferFit")
setwd("downstream_output_Exc3_5/")
# Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_3EGdefault_noTransferFit
cycling_ctl_pranked = read_csv("fGSEA/fGSEA_results/CTL_cyclers_minusLogPRanked.csv")
fgsea_dotplot(cycling_ctl_pranked, enrichment = "pos", my_title = "Pathways Enriched for Genes Cycling in CTL", show = 10)
cycling_ad_pranked = read_csv("fGSEA/fGSEA_results/AD_cyclers_minusLogPRanked.csv")
fgsea_dotplot(cycling_ad_pranked, enrichment = "pos", my_title = "Pathways Enriched for Genes Cycling in AD", BHQ_cutoff = 0.2, show = 10)
log_ad_ctl_25 = read_csv("fGSEA/fGSEA_results/DRgenesAmpRatio25_Log(AD-CTL)ranked.csv")
fgsea_dotplot(log_ad_ctl_25, enrichment = "neg", my_title = "Pathways Enriched for Genes Losing Rhyth in AD", BHQ_cutoff = 0.2)
fgsea_dotplot(log_ad_ctl_25, enrichment = "pos", my_title = "Pathways Enriched for Genes Gain Rhyth in AD", BHQ_cutoff = 0.3)
log_ad_ctl_25_mthd2 = read_csv("fGSEA/fGSEA_results/DRgenesAmpRatio25_Log(AD-CTL)ranked_method2.csv")
fgsea_dotplot(log_ad_ctl_25_mthd2, enrichment = "neg", my_title = "Pathways Enriched for Genes Losing Rhyth in AD", BHQ_cutoff = 0.2)
fgsea_dotplot(log_ad_ctl_25_mthd2, enrichment = "pos", my_title = "Pathways Enriched for Genes Gain Rhyth in AD", BHQ_cutoff = 0.3)

# Other cell types
setwd("downstream_output_Astro_1_2_5_6_7_8/")
cycling_ctl_pranked = read_csv("fGSEA/fGSEA_results/CTL_cyclers_minusLogPRanked.csv")
fgsea_dotplot(cycling_ctl_pranked, enrichment = "pos", my_title = "Pathways Enriched for Genes Cycling in CTL", show = 10)
cycling_ad_pranked = read_csv("fGSEA/fGSEA_results/AD_cyclers_minusLogPRanked.csv")
fgsea_dotplot(cycling_ad_pranked, enrichment = "pos", my_title = "Pathways Enriched for Genes Cycling in AD", BHQ_cutoff = 0.2, show = 10)
log_ad_ctl_1 = read_csv("fGSEA/fGSEA_results/DRgenesAmpRatio1_Log(AD-CTL)ranked.csv")
fgsea_dotplot(log_ad_ctl_1, enrichment = "neg", my_title = "Pathways Enriched for Genes Losing Rhyth in AD", BHQ_cutoff = 0.2)
fgsea_dotplot(log_ad_ctl_1, enrichment = "pos", my_title = "Pathways Enriched for Genes Gain Rhyth in AD", BHQ_cutoff = 0.2)
log_ad_ctl_25 = read_csv("fGSEA/fGSEA_results/DRgenesAmpRatio25_Log(AD-CTL)ranked.csv")
fgsea_dotplot(log_ad_ctl_25, enrichment = "neg", my_title = "Pathways Enriched for Genes Losing Rhyth in AD", BHQ_cutoff = 0.2)
fgsea_dotplot(log_ad_ctl_25, enrichment = "pos", my_title = "Pathways Enriched for Genes Gain Rhyth in AD", BHQ_cutoff = 0.2)
log_ad_ctl_25_mthd2 = read_csv("fGSEA/fGSEA_results/DRgenesAmpRatio25_Log(AD-CTL)ranked_method2.csv")
fgsea_dotplot(log_ad_ctl_25_mthd2, enrichment = "neg", my_title = "Pathways Enriched for Genes Losing Rhyth in AD", BHQ_cutoff = 0.2)
fgsea_dotplot(log_ad_ctl_25_mthd2, enrichment = "pos", my_title = "Pathways Enriched for Genes Gain Rhyth in AD", BHQ_cutoff = 0.3)
