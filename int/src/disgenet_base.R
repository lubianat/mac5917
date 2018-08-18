# This script integrates a gene list with the genes retrieved from disgenet
# written by tiago lubiana, 18/08

df <- read.csv('./int/data/dengue_train_upregulated_list')

#library(devtools)
#install_bitbucket("ibi_group/disgenet2r")

library(disgenet2r)
library(data.table)
library(tidyr)

# search for dengue (code C0011311) related genes
# dq <- disgenetDisease(disease = "umls:C0011311",
#                      database = "ALL",
#                      score = c('>', 0.001))

# Warning: package is not working yet! Manually retrieving genes from database
dengue_disgenet <- fread('./data/tab8b.tsv')
dengue_disgenet_list <- dengue_disgenet$Gene


merged_list <- union(dengue_disgenet_list, df$hgnc_symbol)

write.table(x = merged_list , file = './int/data/dengue_train_upregulated_list_merged_with_disgenet.csv', sep =',', row.names = F, quote=FALSE)
