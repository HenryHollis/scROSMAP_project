meta = readRDS("~/Downloads/green_et_al_annotation_minimum.rds")
meta$cell.type = as.factor(meta$cell.type)
library(ggplot2)

ggplot(meta, aes(x = cell.type, y = nCount_RNA, fill = cell.type)) +
  geom_point(position = position_jitter(seed = 1, width =0.2), size = 0.3, alpha = 0.5) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  theme_minimal() +
  labs(x = "Cell Type", y = "Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
