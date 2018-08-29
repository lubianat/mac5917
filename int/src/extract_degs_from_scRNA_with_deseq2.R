# This script computes DEGs from the same cells in two conditions in scRNA-Seq
# Written by Tiago Lubiana, 17/08/2018

library("DESeq2")

# load files ----
load('./int/data/train_test_sets_dengue.Rdata')
expr_perturbed_train[1:10,1:10]

expr_control_test_count_matrix <- t(expr_control_test)
expr_perturbed_test_count_matrix <- t(expr_perturbed_test)

# ERROR: dependencies ‘png’, ‘hdf5r’ are not available for package ‘Seurat’
# Cleaning data for DESeq2 -----

expr_control_test_count_matrix[1:10, 1:10]

cts <-  merge(expr_control_test_count_matrix,
              expr_perturbed_test_count_matrix,
              by = "row.names",
              all.x = TRUE)

rownames(cts) <- cts[,1] 
cts[,1]  <- NULL

coldata <- as.data.frame(colnames(cts))  
names(coldata) <- c('cells')
coldata$condition <- ifelse(coldata$cells %in% colnames(expr_control_test_count_matrix), 'control', 'perturbed')

rownames(coldata) <- coldata$cells

# Running DESeq2 -----


# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = coldata,
#                               design = ~ condition)

# dds <- DESeq(dds)
# ERROR: every gene contains at least one zero, cannot compute log geometric means
# 1 will be added to every 0 arbitrarily. 
# WARNING: This may limit the validity of the results!

cts <- cts + 1


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)


# extracting DEGs from DESeq2 results

res_df <- as.data.frame(res)

source('./int/src/filter_degs_from_deseq2.R')

up_regulated <- filter_by_adjusted_P(res_df, side = 'up')
down_regulated <- filter_by_adjusted_P(res_df, side = 'down')

write.csv(x = up_regulated, file = './int/data/dengue_test_upregulated')
