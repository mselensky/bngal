import_network_rds <- function (rds.path) {
  
  ebc.groups <- readRDS(rds.path)
  
  tax_levels <- c("phylum", "class", "order", "family", "genus", "asv")
  ebc_df <- list()
  for (i in tax_levels) {
    ebc_df[[i]] <- ebc.groups[[i]] %>%
      select(., `sample-id`, taxon_, rel_abun_binned, edge_btwn_cluster, 
             degree, core, taxa_per_cluster,
             edge_btwn_cluster_color, tax_level, group)
  }
  
  ebc_df_long <- ebc_df %>%
    bind_rows() %>%
    group_by(`sample-id`, edge_btwn_cluster, group, tax_level) %>%
    dplyr::mutate(rel_abun_sum = sum(rel_abun_binned, na.rm = T))

}