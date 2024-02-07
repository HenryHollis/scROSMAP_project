library(tidyverse)
library(VennDiagram)
setwd("~/Box Sync/Henry_stuff/AD_project/scROSMAP/figures/")
setwd("~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/TMMs_w_batch/Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_04ContrVar4EGdefault_noTransferFit_Feb5Redo/")
####Excitatory Neurons
Diff_mesor_EN = read_csv("downstream_output_Exc3_5/differential_mesor_all_genes.csv") %>%
  filter(BHQ < 0.1)

## Differential Expression
setwd("~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/Exc3_5/")
Exact_test_EN = read_csv("Exc3_5_filtByExpr_combatSeq_edgeR_DE_exactTest.csv") %>% 
  filter(FDR < 0.1)

#intersect(Diff_mesor_EN$Gene_Symbols, Exact_test_EN$...1)
DM_not_DE = setdiff(Diff_mesor_EN$Gene_Symbols,Exact_test_EN$...1 )

x <- list(
  `Diff.  Mesor` = Diff_mesor_EN$Gene_Symbols, 
  `Diff. Expr.` = Exact_test_EN$...1
)
library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),auto_scale = T,
  stroke_size = 0.02, set_name_size = 4
)+ggtitle("Excitatory Neurons 3 & 5")

###Inhibitory Neurons
setwd("~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/TMMs_w_batch/Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_04ContrVar4EGdefault_noTransferFit_Feb5Redo/")
Diff_mesor_IN = read_csv("downstream_output_Inh_2_8_11_12_13/differential_mesor_all_genes.csv") %>%
  filter(BHQ < 0.1)

## Differential Expression
setwd("~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/Inhib_2_8_11_12_13/")
Exact_test_IN = read_csv("Inhib_2_8_11_12_13_filtByExpr_combatSeq_edgeR_DE_exactTest.csv") %>% 
  filter(FDR < 0.1)

x <- list(
  `Diff.  Mesor` = Diff_mesor_IN$Gene_Symbols, 
  `Diff. Expr.` = Exact_test_IN$...1
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),auto_scale = T,
  stroke_size = 0.02, set_name_size = 4
)+ggtitle("Inhibitory Neurons 2,8,11,12,13")


###Astrocytes Neurons
setwd("~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/TMMs_w_batch/Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_04ContrVar4EGdefault_noTransferFit_Feb5Redo/")
Diff_mesor_AS = read_csv("downstream_output_Ast_1_2_5_6_7_8/differential_mesor_all_genes.csv") %>%
  filter(BHQ < 0.1)

## Differential Expression
setwd("~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/Ast_1_2_5_6_7_8/")
Exact_test_AS = read_csv("Ast_1_2_5_6_7_8_filtByExpr_combatSeq_edgeR_DE_exactTest.csv") %>% 
  filter(FDR < 0.1)

x <- list(
  `Diff.  Mesor` = Diff_mesor_AS$Gene_Symbols, 
  `Diff. Expr.` = Exact_test_AS$...1
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),auto_scale = T,
  stroke_size = 0.02, set_name_size = 4
)+ggtitle("Astrocytes 1, 2, 5, 6, 7, 8")


###Microglia Neurons
setwd("~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/TMMs_w_batch/Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_04ContrVar4EGdefault_noTransferFit_Feb5Redo/")
Diff_mesor_Mglia = read_csv("downstream_output_Mglia_4_5_6_7_8_10_11_13_15_16/differential_mesor_all_genes.csv") %>%
  filter(BHQ < 0.1)

## Differential Expression
setwd("~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/Mglia_4_5_6_7_8_10_11_13_15_16/")
Exact_test_Mglia = read_csv("Microglia_4_5_6_7_8_10_11_13_15_16_filtByExpr_combatSeq_edgeR_DE_exactTest.csv") %>% 
  filter(FDR < 0.1)

x <- list(
  `Diff.  Mesor` = Diff_mesor_Mglia$Gene_Symbols, 
  `Diff. Expr.` = Exact_test_Mglia$...1
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),auto_scale = T,
  stroke_size = 0.02, set_name_size = 4
)+ggtitle("Microglia 4,5,6,7,8,10,11,13,15,16")
