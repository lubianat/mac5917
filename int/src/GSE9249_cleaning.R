# this script cleans GEO dataset GSE92495
# written by tiago lubiana, 14/08/2018
library(GEOquery)
get_gse <- function(gse_id){
  gse <- GEOquery::getGEO(gse_id, GSEMatrix=TRUE, getGPL=TRUE)
  return(gse)
}

if (f){
'https://www.ncbi.nlm.nih.gov/pubmed/28192419'
}

GSE92495 <- get_gse('GSE92495')
gse <- GSE92495

platform <- Biobase::annotation(gse[[1]])
expr <- as.data.frame(Biobase::exprs(gse[[1]]))
