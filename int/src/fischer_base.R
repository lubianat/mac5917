#this script compares predicted genes with DEGs via Fischer's exact P test

rm(list=ls())


library(data.table)
degs_test <- fread('./int/data/dengue_test_upregulated')
degs_train <- fread('./int/data/dengue_train_upregulated')
predicted_degs <- read.csv('./int/data/gene_candidates_top_dengue')
load('./int/data/train_test_sets_dengue.Rdata')

total <- expr_control_test 

source('./int/src/ensg_to_gene_symbol.R')

# change from ENSG to gene symbol
degs_test <- ensg_to_gene_symbol(degs_test$gene)
degs_train <- ensg_to_gene_symbol(degs_train$gene)
total_with_gene_symbol <- ensg_to_gene_symbol(colnames(total))

#Find DEGS that are unique for this dataset
unique_degs <- degs_test[!(degs_test$hgnc_symbol %in% degs_train$hgnc_symbol),]

#Find which of the DEGs could be present in orginal data (because name space of DisGeNET is bigger than expr)
predicted_degs_in_original_data <- predicted_degs$gene[predicted_degs$gene %in% total_with_gene_symbol$hgnc_symbol]
predicted_degs_in_original_data <- droplevels(predicted_degs_in_original_data)

predicted_degs_in_original_data[predicted_degs_in_original_data %in% unique_degs$hgnc_symbol]

# create contingency table

degy_topy <- length(predicted_degs_in_original_data[predicted_degs_in_original_data %in% unique_degs$hgnc_symbol])
degn_topy <- length(predicted_degs_in_original_data) - degy_topy
degy_topn <- length(unique_degs$hgnc_symbol) - degy_topy
degn_topn <- 11999- (degy_topn+degn_topy+degy_topy)


# total is 11999 so the contingency table is 
ctg_table = rbind(c(degy_topy,degy_topn),
              c(degn_topy,degn_topn));

fisher.test(ctg_table)
