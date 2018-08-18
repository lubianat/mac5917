# this script cleans the scrapes GSE110496 for exprs and annotation
# written by tiago lubiana, 16/08/2018

setwd("~/Documents/Disciplines/Fujita")

library(data.table)
library(GEOquery)
library(dplyr)
expr_dengue <- fread('./data/counts_dengue.tsv')
metadata_dengue <- fread('./data/cell_metadata_dengue.tsv')

get_gse <- function(gse_id) {
  gse <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)
  return(gse)
}

if (FALSE) {
  'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5826272/'
}

GSE110496 <- get_gse('GSE110496')
gse <- GSE110496

platform <- Biobase::annotation(gse[[1]])
expr <- as.data.frame(Biobase::exprs(gse[[1]]))

annot_dengue <- gse$GSE110496_series_matrix.txt.gz@phenoData@data %>%
  filter(`virus:ch1` == 'dengue') 

annot_dengue_infected <- gse$GSE110496_series_matrix.txt.gz@phenoData@data %>%
                filter(`virus:ch1` == 'dengue' & `moi:ch1` > 0) 

annot_dengue_control <- gse$GSE110496_series_matrix.txt.gz@phenoData@data %>%
  filter(`virus:ch1` == 'dengue' & `moi:ch1` == 0) 


write.csv(x = expr_dengue, './int/data/expr_dengue.csv')
write.csv(x = annot_dengue, './int/data/annot_dengue.csv')
write.csv(x = annot_dengue_infected, './int/data/annot_dengue_infected.csv')
write.csv(x = annot_dengue_control, './int/data/annot_dengue_control.csv')
