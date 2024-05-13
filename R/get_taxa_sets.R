get_taxa_sets <- function (binned.taxonomy, metadata, set.group) {
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")
  
  binned.meta = list()
  for (i in tax.levels) {
    binned.meta[[i]] <- binned.taxonomy[[i]] %>%
      left_join(., metadata, "sample-id") %>%
      filter(rel_abun_binned > 0) %>%
      group_by(taxon_, .data[[set.group]]) %>%
      dplyr::mutate(pres_abs = 1) %>%
      distinct(taxon_, .data[[set.group]], .keep_all = TRUE) %>%
      ungroup() %>%
      select(phylum, taxon_, .data[[set.group]], pres_abs) %>%
      pivot_wider(values_from = "pres_abs", names_from = .data[[set.group]]) %>%
      replace(is.na(.), 0)
  }
  
  binned.meta
  
}