
library(fgsea)
library(rstudioapi)
library(data.table)
library(stringr)
library(purrr)
run_fgsea = function(files, gene_dict, pathways,possible_human_gene_names, gsea_param = 1,min_size = 10,max_size = 500, dir = "fGSEA_results"){
  sapply(files, function(file){
  file_name = paste0("../", dir, "/", str_replace(file, ".rnk", "") , ".csv")
  print(file_name)
  score_type =  ifelse(grepl("minusLogPRanked", file_name, ignore.case = T) , "pos", "std")
  print(paste("Used", score_type, "fGSEA score type"))
  plot_name = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , ".png")
  
  ranks <- read.table(file, header=F, colClasses = c("character", "numeric"))
  #here I handle remapping my gene symbols to the "chip" file from MsigDB
  ranks_merged = merge(ranks, gene_dict, by.x = "V1", by.y = "Probe.Set.ID", all.y = F, all.x = T)
  #I don't want to remap if original gene name is in my data:
  ranks_merged = ranks_merged %>%
    mutate(unified_name = case_when(
      V1 == Gene.Symbol ~ V1, # If name1 and name2 agree, take either
      V1 != Gene.Symbol & V1 %in% possible_human_gene_names ~ V1, # If name1 is in ROSMAP, take name1
      V1 != Gene.Symbol & Gene.Symbol %in% possible_human_gene_names ~ Gene.Symbol, # If name2 is in ROSMAP, take name2
      TRUE ~ NA_character_ # If neither is in true_names, set unified to NA
    )) %>% arrange(V2)
  
  if(!is_empty(which(duplicated(ranks_merged$V1)))){
    ranks_merged = ranks_merged[-which(duplicated(ranks_merged$V1)),]   #remove duplicates if any
  }
  
  ranks_vec = ranks_merged[,2]
  names(ranks_vec) = ranks_merged[,1]
  set.seed(42)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks_vec,
                    minSize  = min_size,
                    #eps = 0.0,
                    gseaParam = gsea_param,
                    maxSize  = max_size,
                    scoreType = score_type,
                    nPermSimple = 10000)
  
  fgseaRes = fgseaRes %>% arrange(pval)
  fwrite(fgseaRes, file_name)
  
  topPathwaysUp <- fgseaRes[ES > 0 & !is.na(pval)][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0 & !is.na(pval)][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, topPathwaysDown)

  p = plotGseaTable(pathways[topPathways], ranks_vec, fgseaRes, 
                gseaParam=gsea_param)
  
  ggsave(plot_name, p, bg = 'white')
  
  
  # #if the rnk was ranked by -log(p), only the positive enriched pathways are important:
  # if (grepl("minusLogPRanked", file_name, ignore.case = T)){
  #   p1 = fgsea_dotplot(fgseaRes,enrichment = "pos" , my_title = paste(tools::file_path_sans_ext(file), "pos enriched pathways"))
  #   plot_name2 = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , "_dotplot",".png")
  #   ggsave(plot_name2, p1, bg = 'white')
  # 
  #   # if the rnk was ranked by log(ad_amp/ctl_amp), we want both pathways
  # }else{
  #   if (!purrr::is_empty(topPathwaysUp)){
  #     p2 = fgsea_dotplot(fgseaRes,enrichment = "pos" ,my_title = paste(tools::file_path_sans_ext(file), "pos enriched pathways"))
  #     plot_name2 = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , "_posEnrichment_dotplot",".png")
  #     ggsave(plot_name2, p2, bg = 'white')
  #   }
  #   if (!purrr::is_empty(topPathwaysDown)){
  #     p3 = fgsea_dotplot(fgseaRes,enrichment = "neg" ,my_title = paste(tools::file_path_sans_ext(file), "neg enriched pathways"))
  #     plot_name3 = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , "_negEnrichment_dotplot",".png")
  #     ggsave(plot_name3, p3, bg = 'white')
  #   }
  # }
  # 
})
}

