# this script queries GEO and disgenet to find scRNA-Seq studies about diseases
# written by tiago lubiana, 14/08/2018


library(rentrez)

# exploratory searches ------
entrez_db_searchable('geoprofiles')


bla <- entrez_search(db = 'gds',
                     term = 'GPL[ETYP]',
                     use_history = F)

ble <- entrez_summary(db = 'gds', id = bla$ids)

gds_search <- rentrez::entrez_search(db = "gds",
                                     term = "GPL[ETYP]",
                                     use_history = TRUE)
study_count <- gds_search$count

id_vec <- c()

for (seq_start in seq(0, study_count, 50)) {
  message(seq_start)
  search_res <-
    rentrez::entrez_search(
      db = "gds",
      term = paste0("GPL[ETYP]"),
      use_history = TRUE,
      retmax = 50,
      retstart = seq_start
    )
  print(search_res$ids)
  id_vec <- c(id_vec, search_res$ids)

}
id_vec <- unlist(id_vec)
res <- entrez_summary(db="gds", id=id_vec[1:10])

