calc_ebc_comps <- function(ebc.nodes, binned.taxonomy, metadata.cols) {

  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")
  ebc.nodes.abun <- list()

  for (i in tax.levels) {
    ebc.nodes.abun[[i]] <- ebc.nodes[[i]] %>%
      left_join(select(binned.taxonomy[[i]],
                       `sample-id`, taxon_, rel_abun_binned),
                by = "taxon_") %>%
      left_join(., select(metadata, `sample-id`, all_of(subset_column), all_of(metadata.cols)),
                by = c("sample-id")) %>%
      filter(group == .data[[subset_column]]) %>%
      ungroup() %>% group_by(`sample-id`, edge_btwn_cluster, group) %>%
      distinct(`sample-id`, taxon_, edge_btwn_cluster, group, .keep_all = T) %>%
      ungroup() %>% group_by(edge_btwn_cluster, group) %>%
      dplyr::mutate(ebc_abun_sum = rel_abun_binned/sum(rel_abun_binned, na.rm=TRUE))
    # ungroup() %>%
    # distinct(`sample-id`, taxon_, edge_btwn_cluster, group, .keep_all = T)
  }

  ebc.nodes.abun

}
