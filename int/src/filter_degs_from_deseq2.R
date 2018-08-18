# this script filters by adjusted-p and logfc a deseq 2 table
# written by tiago lubiana 17/08 

library(tibble)
library(dplyr)
filter_by_adjusted_P <-
  function(data,
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
        filter(log2FoldChange < 0 ) %>%
        filter(abs(log2FoldChange) < logFC_threshold)  %>%
        filter(padj < readjustedP_threshold)
    }
  }