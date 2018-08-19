# this script runs biomart to chanse ENSG to gene symbols
# WRITTEN BY TIAGO lUBIANA 18/08/2018

#df <- read.csv('./int/data/dengue_train_upregulated')

library('biomaRt')
ensg_to_gene_symbol <- function(ensg_vector){
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-ensg_vector
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "entrezgene", "hgnc_symbol", "description"),values=genes,mart= mart)
glist <- as.data.frame(G_list)

return(glist)
}


# write.csv(x = glist, file = './int/data/dengue_train_upregulated_list')
