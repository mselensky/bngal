#' Export EBC and taxonomic abundance summary data to CSV file
#'
#' @param binned.taxonomy Output from [`bngal::bin_taxonomy()`]
#' @param ebc.nodes.abun Output from [`bngal::ebc_compositions()`]
#' @param tax.level Taxonomic level at which to summarize data.
#' @param out.dr Output directory. Results will be saved to the `"network-summary-tables"` subfolder.
#'
#' @return
#' @export
#'
#' @examples
export_ebc_taxa_summary <- function(binned.taxonomy, ebc.nodes.abun, tax.level, out.dr) {

  for (x in names(ebc_comps[[tax.level]])) {

    communities = unique(ebc.nodes.abun[[tax.level]][[x]]$`sample-id`)
    n_samples = length(communities)

    output <- binned.taxonomy[[tax.level]] %>%
      filter(`sample-id` %in% communities) %>%
      left_join(select(ebc.nodes.abun[[tax.level]][[x]], `sample-id`, taxon_, edge_btwn_cluster, taxa_per_cluster, ebc_abun_sum), by = c("sample-id", "taxon_")) %>%
      group_by(edge_btwn_cluster) %>%
      dplyr::mutate(taxa_per_cluster = if_else(edge_btwn_cluster == 0,
                                               as.numeric(length(unique(taxon_))),
                                               as.numeric(taxa_per_cluster)))

    a <- output %>%
      dplyr::mutate(pres = if_else(rel_abun_binned == 0, 0, 1)) %>%
      distinct(`sample-id`, edge_btwn_cluster, pres) %>%
      group_by(edge_btwn_cluster) %>%
      dplyr::summarize(ebc_prev = sum(pres)/n_samples) %>%
      dplyr::arrange(desc(ebc_prev)) %>%
      dplyr::mutate(sub_comm = x)
    b <- output %>%
      distinct(`sample-id`, taxon_, edge_btwn_cluster, rel_abun_binned) %>%
      group_by(taxon_, edge_btwn_cluster) %>%
      dplyr::summarize(mean_rel_abun = mean(rel_abun_binned),
                       sd_rel_abun = sd(rel_abun_binned),
                       median_rel_abun = median(rel_abun_binned),
                       max_rel_abun = max(rel_abun_binned),
                       min_rel_abun = min(rel_abun_binned),
                       n_obs = n(),
                       .groups = "drop_last") %>%
      dplyr::arrange(desc(n_obs), edge_btwn_cluster) %>%
      dplyr::mutate(sub_comm = x)
    c <- output %>%
      distinct(`sample-id`, edge_btwn_cluster, ebc_abun_sum)
    d <- c %>%
      group_by(edge_btwn_cluster) %>%
      dplyr::summarize(mean_rel_abun = mean(ebc_abun_sum),
                       sd_rel_abun = sd(ebc_abun_sum),
                       median_rel_abun = median(ebc_abun_sum),
                       max_rel_abun = max(ebc_abun_sum),
                       min_rel_abun = min(ebc_abun_sum),
                       n_obs = n()) %>%
      left_join(a, by = "edge_btwn_cluster") %>%
      select(edge_btwn_cluster, ebc_prev, everything()) %>%
      dplyr::arrange(desc(n_obs), edge_btwn_cluster)

    out.dr.sum = file.path(out.dr, "network-summary-tables", tax.level)
    if (!dir.exists(out.dr.sum)) dir.create(out.dr.sum, recursive = TRUE)

    write_csv(d, file.path(out.dr.sum, paste0(x, "_ebc_distribution.csv")))
    write_csv(b, file.path(out.dr.sum, paste0(x, "_tax_spread.csv")))
    write_csv(c, file.path(out.dr.sum, paste0(x, "_ebc_abundance_per_sample.csv")))

  }

}
