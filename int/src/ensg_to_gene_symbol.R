# this script runs biomart to chanse ENSG to gene symbols
# WRITTEN BY TIAGO lUBIANA 18/08/2018

#df <- read.csv('./int/data/dengue_train_upregulated')

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
chr1_genes <- getBM(attributes=c('ensembl_gene_id',
                                   'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters =
                        'chromosome_name', values ="1", mart = ensembl)

head(chr1_gene)

# write.csv(x = glist, file = './int/data/dengue_train_upregulated_list')
