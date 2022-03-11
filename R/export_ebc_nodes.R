export_ebc_nodes <- function (network.data) {
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")
  
  all_nodes <- list()
  for (i in tax.levels) {
    all_nodes[[i]] <- purrr::reduce(network.data[[i]]$nodes, bind_rows)
    all_nodes[[i]] <- all_nodes[[i]] %>%
      dplyr::mutate(tax_level = i)
  }
  b <- purrr::map(all_nodes, 
                  ~distinct(., label, edge_btwn_cluster, group, .keep_all = T) %>%
                    dplyr::rename(degree=value,
                                  taxon_=label))
  b
}