# This script outlines an R pipeline for the discovery
#of genes related to DEGs in a dataset
#
# Tiago 08/09/18


disgenet_score <- c(0.5, 0.1, 0.05)
disgenet_articles <- c(1, 2, 3)
string_score <- c(0.4, 0.6)
encode_db <- c('direct', 'common_tf')
x2k <- c('same_regulation')
network_measure <- c('strength', 'eigen_centrality')
top_selected <- c(0.1, 0.5, 1.0)


x2k <- c(x2k, 'blank')
encode_db <- c(encode_db, 'blank')
string_score <- c(string_score, 'blank')
disgenet <-
  c('score_0.5',
    'score_0.1',
    'score_0.05',
    '1_article',
    '2_articles',
    '3_articles')

search_space <- expand.grid(
  disgenet = disgenet,
  string_score = string_score,
  encode_db = encode_db,
  x2k = x2k,
  network_measure = network_measure,
  top_selected = top_selected
)

library(dplyr)


search_space <- filter(search_space,
                       encode_db != 'blank' |
                         string_score != 'blank' |
                         x2k !=   'blank')


## pseudofunction for the exploration of  the pipeline -----------
function explore_options(
  vector_of_degs,
  disease_name,
  search_space,
  vector_of_target_genes,
  original_dataset_length
  number_of_runs = 5,
  start_n = 1,

) {
  set.seed(3)
  order <- sample(1:nrow(search_space))
  
  best_index = 0
  best_pvalue = 1
  
  
  for (run in start_n:(start_n + number_of_runs - 1)) {
    status_now <- search_space[order[run], ]
    disgenet_genes <-
      filter_disgenet_table(status_now$disgenet, disease_name)
    
    joint_genes <- c(disgenet_genes, vector_of_degs)
    
    edgelist_from_string <-
      get_STRING_edges(joint_genes, status_now$string_score)
    edgelist_from_encode <-
      get_encode_edges(joint_genes, status_now$encode_db)
    edgelist_from_kinases <-
      get_kinase_edges(joint_genes, status_now$x2k)
    
    total_edges <-
      rbind(edgelist_from_string,
            edgelist_from_encode,
            edgelist_from_kinases)
    
    total_edges$weights <-
      ifelse(todal_edges[, 1] %in% vector_of_degs &
               todal_edges[, 2] %in% vector_of_degs,
             4,
             1)
    
    
    total_edges$weights[xor(
      total_edges$weights[, 1] %in% vector_of_degs,
      total_edges$weights[, 2] %in% vector_of_degs
    )] <- 2
    
    library(igraph)
    
    g <- graph.data.frame(total_edges)
    
    gene_candidates <- get_candidates(g, status_now$network_measure, status_now$top_selected)
    
    
    unique_targets <- setdiff(vector_of_target_genes,vector_of_degs)
    
    n_candidates_and_targets_matched <- length(gene_candidates[gene_candidates %in% unique_targets])
    n_candidates_not_matched <- length(gene_candidates) - candidate_match
    n_targets_not_matched <- length(unique_degs$hgnc_symbol) - degy_topy
    n_of_irrelevant_genes <- original_dataset_length - (n_candidates_and_targets_matched+n_candidates_not_matched+n_targets_not_matched)
    
    # so the contingency table is 
    ctg_table = rbind(c(n_candidates_and_targets_matched,n_candidates_not_matched),
                      c(n_targets_not_matched,n_of_irrelevant_genes));
    
    fisher_result <- fisher.test(ctg_table)
    if (fisher_result$p.value < best_pvalue){
      best_pvalue <- fisher_result$p.value
      best_index <- order[run]
    }

    
    
    
    
  }
  
  
  
}






}

## disgenet filter pseudofunction -------

filter_disgenet_table <-
  function(status_now$disgenet, disease_name) {
    # read table from disgenet
    # subset table depending on criteria
    # return vector of genes
  }



##get_string_edges pseudofunction ----

get_string_edges <-
  function(vector_of_genes,
           status_now$string_score)
    # query string for edges
    # return dataframe of  edges
    
    
    ##get_encode_edges pseudofunction ----

get_encode_edges <-
  function(vector_of_genes, status_now$encode_db)
    # query encode for edges
    # return dataframe of  edges
    
##get_kinase_edges pseudofunction ----
    get_kinase_edges <- function(vector_of_genes, status_now$x2k)
      # query x2k for edges
      # return dataframe of  edges

##get_candidates pseudofunction ----
    
get_candidates <-  function(igraph_object,
                            status_now
  # uses status_now$top_selected and status_now$network_measure to select top hits
  # return vector of characters containing the candidates
          
          