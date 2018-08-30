# consensus DEGs from MAST and zinbwave
# Tiago 30/08/2018


library(data.table)
require(tibble)
library(dplyr)

degs_MAST_train <- fread('int/results/DE_MAST_denge_train.csv')
degs_deseq_zinbwave_train <-
  fread('int/results/DE_zinbwave_DESeq2_top1000var_denge_train.csv')

top_deseq_degs <-  filter_deseq_by_adjusted_P(degs_deseq_zinbwave_train)
top_mast_degs <- filter_mast_by_adjusted_P(degs_MAST_train)


top_intersect_degs <- top_deseq_degs[top_deseq_degs$V1 %in% top_mast_degs$primerid,]

write.csv(top_intersect_degs, './int/results/intersect_train_MAST_zinbwave.csv')


filter_deseq_by_adjusted_P <-  function(data,
                                        side = 'up',
                                        logFC_threshold = 0.3222 ,
                                        readjustedP_threshold = 0.05) {
  data <- rownames_to_column(data, var = 'gene')
  if (side == 'up') {
    data %>%
      filter(log2FoldChange > logFC_threshold)  %>%
      filter(padj < readjustedP_threshold)
  } else if (side == 'down') {
    data %>%
      filter(log2FoldChange < 0) %>%
      filter(abs(log2FoldChange) < logFC_threshold)  %>%
      filter(padj < readjustedP_threshold)
  }
}

filter_mast_by_adjusted_P <-  function(data,
                                       side = 'up',
                                       readjustedP_threshold = 0.05) {
  data <- rownames_to_column(data, var = 'gene')
  if (side == 'up') {
    data %>%
      filter(fdr < readjustedP_threshold)
  } else if (side == 'down') {
    data %>%
      filter(log2FoldChange < 0) %>%
      filter(fdr < readjustedP_threshold)
  }
}
