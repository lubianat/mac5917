# Generates a network from stringdb
#written by tiago lubiana 17/08/18
S

library(data.table)

source("https://bioconductor.org/biocLite.R")
biocLite("STRINGdb")
library("STRINGdb")
string_db <- STRINGdb::STRINGdb$new(version = "10", species = 9606, score_threshold = 0, input_directory = "" )
# bla <- fread('./int/data/dengue_train_upregulated_list_merged_with_disgenet')
names(bla) <- c('index', 'gene')
data(diff_exp_example1)

string_id <- string_db$map(data.frame(gene = moduleGenes), "gene",
                           takeFirst = FALSE) 

bla_mapped <- string_db$map(bla, "gene", removeUnmappedRows = TRUE )


# I was not able to use the package. THe internet browser version was used instead.
#ACtually, I have to convert gene symbols into protein symbols
#thats why I used network analyst for now

