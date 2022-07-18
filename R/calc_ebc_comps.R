#' Calculate compositions of edge between clusters based on defined metadata column(s)
#'
#' @param ebc.nodes
#' @param binned.taxonomy
#' @param tax.level
#' @param metadata.cols
#'
#' @return
#' @export
#'
#' @examples
calc_ebc_comps <- function(ebc.nodes, binned.taxonomy, tax.level, metadata.cols) {

  ebc.nodes.abun <- ebc.nodes %>%
    filter(tax_level %in% tax.level) %>%
    left_join(select(binned.taxonomy[[tax.level]],
                     `sample-id`, taxon_, rel_abun_binned, binned_count),
              by = "taxon_")

  # this ensures ebc relative abundance is calculated from full dataset
  # regardless of remove.singletons option from bngal::bin_taxonomy()
  full.data <- binned.taxonomy[[tax.level]] %>%
    select(`sample-id`, taxon_, binned_count) %>%
    left_join(select(ebc.nodes.abun, -binned_count), by = c("sample-id", "taxon_")) %>%
    ungroup() %>%
    select(`sample-id`, taxon_, binned_count) %>%
    pivot_wider(names_from = "taxon_", values_from = "binned_count") %>%
    filter(!is.na(`sample-id`)) %>%
    pivot_longer(cols = 2:ncol(.), names_to = "taxon_", values_to = "binned_count") %>%
    dplyr::mutate(binned_count = if_else(is.na(binned_count), 0, binned_count))

  full.data.ebc <- full.data %>%
    left_join(select(binned.taxonomy[[tax.level]], -binned_count), by = c("sample-id", "taxon_")) %>%
    left_join(select(ebc.nodes.abun, -binned_count, -rel_abun_binned)) %>%
    left_join(., select(metadata, `sample-id`, any_of(metadata.cols)),
              by = c("sample-id")) %>%
    distinct(`sample-id`, taxon_, edge_btwn_cluster, .keep_all = T) %>%
    ungroup() %>% group_by(`sample-id`, edge_btwn_cluster) %>%
    dplyr::mutate(ebc_count = sum(binned_count, na.rm=TRUE))

  dat.out <- full.data.ebc %>%
    distinct(`sample-id`, edge_btwn_cluster, ebc_count) %>%
    dplyr::mutate(edge_btwn_cluster = if_else(is.na(edge_btwn_cluster), "no_cluster", edge_btwn_cluster)) %>%
    group_by(`sample-id`) %>%
    dplyr::mutate(ebc_abun_sum = ebc_count/sum(ebc_count))

  phylum_colors <- dplyr::rename(bngal:::phylum_colors, color_order = order)

  full.data.ebc %>%
    select(-ebc_count) %>%
    left_join(dat.out, by = c("sample-id", "edge_btwn_cluster")) %>%
    left_join(phylum_colors2, by = "phylum") %>%
    dplyr::mutate(tax_level = if_else(is.na(tax_level), tax.level, tax_level))
  #dplyr::mutate(edge_btwn_cluster = if_else(is.na(edge_btwn_cluster), "no_cluster", as.character(edge_btwn_cluster)))

}
