# this script gets ctrl x perturbed scRNA seqs and divides into train and test (50% each)
# written by tiago lubiana, 16/08/2018

library(data.table)
library(dplyr)
library(tibble)

# Load data ------------
expr <- fread('./int/data/expr_dengue')
annot_perturbed <- fread('./int/data/annot_dengue_infected.csv')
annot_control <- fread ('./int/data/annot_dengue_control.csv')

# Define functions ---------
# extracts a training set from a expr file and a annotation
extract_train <- function(expr, annot) {
  expr_perturbed <- expr %>%
    select(c('EnsemblID', annot$title))
  
  expr_perturbed_matrix <- expr_perturbed[,-1]
  expr_perturbed_matrix <- t(expr_perturbed_matrix)
  
  expr_perturbed_df <- data.frame(expr_perturbed_matrix)
  colnames(expr_perturbed_df) <- expr_perturbed$EnsemblID
  
  expr_perturbed_train <- expr_perturbed_df %>%
    sample_frac(.5)
  
}

# Run script ---------------

set.seed('33')

# Select samples that are perturbed and split
expr_perturbed <- expr %>%
  select(c('EnsemblID', annot_perturbed$title))

expr_perturbed_matrix <- expr_perturbed[,-1]
expr_perturbed_matrix <- t(expr_perturbed_matrix)

expr_perturbed_df <- data.frame(expr_perturbed_matrix)
rownames(expr_perturbed_df) <- colnames(expr_perturbed[, -1])
colnames(expr_perturbed_df) <- expr_perturbed$EnsemblID

expr_perturbed_train <- extract_train(expr, annot_perturbed)

expr_perturbed_test <- expr_perturbed_df %>%
  rownames_to_column('sample') %>%
  filter(!(sample %in% rownames(expr_perturbed_train))) %>%
  column_to_rownames('sample')

expr_perturbed_test[1:10, 1:10]
# Select samples that are control and split

expr_control <- expr %>%
  select(c('EnsemblID', annot_control$title))

expr_control_matrix <- expr_control[,-1]
expr_control_matrix <- t(expr_control_matrix)

expr_control_df <- data.frame(expr_control_matrix)
rownames(expr_control_df) <- colnames(expr_control[, -1])
colnames(expr_control_df) <- expr_control$EnsemblID

expr_control_train <- extract_train(expr, annot_control)

expr_control_test <- expr_control_df %>%
  rownames_to_column('sample') %>%
  filter(!(sample %in% rownames(expr_control_train))) %>%
  column_to_rownames('sample')

expr_control_test[1:10, 1:10]

save(
  list =  c(
    "expr_control_test",
    "expr_control_train",
    "expr_perturbed_test",
    "expr_perturbed_train"
  ),
  file =  './int/data/train_test_sets_dengue.Rdata'
)
