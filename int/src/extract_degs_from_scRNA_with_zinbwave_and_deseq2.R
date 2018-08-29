# This script computes DEGs from the same cells in two conditions in scRNA-Seq
# using zinbwave and DESEQ2
# Written by Tiago Lubiana, 24/08/2018

library(data.table)
library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)

load(file = './int/data/train_test_sets_dengue.Rdata')
full_dengue_train <- rbind(expr_control_train, expr_perturbed_train)
annotation <- fread('./int/data/annot_dengue.csv')
annotation <-  annotation[annotation$title %in% rownames(full_dengue_train)]


full_dengue_train <- full_dengue_train[annotation$title,]

# assign counts to a single cell experiment object
counts <- as.matrix(full_dengue_train)
sce <- SingleCellExperiment(assays = list(counts = t(counts)))

# filter by genes with at least 5 counts in at least 5 cells
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce <- sce[filter,]


# get only top 1000 variable genes to build model
# assay(sce) %>% log1p %>% rowVars -> vars
# names(vars) <- rownames(sce)
# vars <- sort(vars, decreasing = TRUE)
# 
# sce <- sce[names(vars)[1:1000],]
assayNames(sce)[1] <- "counts"

zinb <- zinbFit(sce, K=2, epsilon=100000)

sce_zinb <- zinbwave(sce, fitted_model = zinb, K = 2, epsilon=1000)

# weights <- assay(sce, "weights")
# Error in assay(sce, "weights") : 
#   'assay(<SingleCellExperiment>, i="character", ...)' invalid subscript 'i'
# 'i' not in names(assays(<SingleCellExperiment>))

library(DESeq2)
assay(sce) <- as.matrix(assay(sce)) + 1
colData(sce_zinb)$Biological_Condition = as.factor(ifelse(annotation$`moi:ch1` == 0, 'ctrl', 'infected'))

dds <- DESeqDataSet(sce_zinb, design = ~ Biological_Condition)

dds <- DESeq(dds)

res <- lfcShrink(dds, contrast=c("Biological_Condition", "ctrl", "infected"),
                 type = "normal")
head(res)

results_zinbwave <- as.data.frame(res)

write.csv(results_zinbwave, file = './int/results/DE_zinbwave_DESeq2_top1000var_denge_train.csv')
