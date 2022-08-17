#' Calculate compositions of edge between clusters based on defined metadata column(s)
#'
#' @param ebc.nodes Output from [`bngal::extract_node_data()`]
#' @param binned.taxonomy Output from [`bngal::bin_taxonomy()`]
#' @param tax.level Taxonomic level of classification at which to calculate EBC compositions
#' @param metadata.cols
#' @param sub.comms *(Optional)* Metadata column by which to split data into subcommunities
#'
#' @return
#' @export
#'
#' @examples
ebc_compositions <- function(ebc.nodes, binned.taxonomy, alpha.div, tax.level, metadata, metadata.cols, sub.comms) {

  if (missing(sub.comms) | is.null(sub.comms)) {
    sub.comms = "all_communities"
    metadata[[sub.comms]] = rep("all", nrow(metadata))
  }

  ebc.nodes <- split(ebc.nodes, ebc.nodes$sub_comm)
  out=list()
  for (x in unique(metadata[[sub.comms]])) {
    communities = metadata %>%
      filter(.data[[sub.comms]] %in% x) %>%
      pull(`sample-id`)

    full_abun_data <- binned.taxonomy[[tax.level]] %>%
      filter(`sample-id` %in% communities)

    ebc.nodes.abun <- ebc.nodes[[x]] %>%
      filter(tax_level %in% tax.level) %>%
      left_join(select(full_abun_data,
                       `sample-id`, taxon_, rel_abun_binned, binned_count),
                by = c("taxon_"))

    # this ensures ebc relative abundance is calculated from full dataset
    # regardless of remove.singletons option from bngal::bin_taxonomy()
    full.data <- full_abun_data %>%
      select(`sample-id`, taxon_, binned_count) %>%
      pivot_wider(names_from = "taxon_", values_from = "binned_count") %>%
      filter(!is.na(`sample-id`)) %>%
      pivot_longer(cols = 2:ncol(.), names_to = "taxon_", values_to = "binned_count") %>%
      dplyr::mutate(binned_count = if_else(is.na(binned_count), 0, binned_count))

    tax.levels = c("phylum", "class", "order", "family", "genus", "asv")
    select.by.tax = names(select(as.data.frame(t(data.frame(row.names = tax.levels))), phylum:.data[[tax.level]]))

    full.data.ebc <- full.data %>%
      left_join(select(full_abun_data, -binned_count), by = c("sample-id", "taxon_")) %>%
      left_join(select(ebc.nodes.abun, -binned_count, -rel_abun_binned),
                by = c("sample-id", "taxon_", "domain", all_of(select.by.tax))) %>%
                #by = c("sample-id", "taxon_", "domain", "phylum", "class", "order", "family", "genus", "asv")
      dplyr::mutate(rel_abun_binned = if_else(binned_count == 0, 0, rel_abun_binned),
                    edge_btwn_cluster = if_else(is.na(edge_btwn_cluster), 0, as.numeric(edge_btwn_cluster))) %>%
      left_join(., select(metadata, `sample-id`, any_of(metadata.cols)),
                by = c("sample-id")) %>%
      group_by(`sample-id`, edge_btwn_cluster) %>%
      dplyr::mutate(ebc_count = sum(binned_count, na.rm=TRUE))


    full.data.ebc.summ <- full.data.ebc %>%
      distinct(`sample-id`, edge_btwn_cluster, ebc_count) %>%
      group_by(`sample-id`) %>%
      dplyr::mutate(ebc_abun_sum = ebc_count/sum(ebc_count))

    phylum_colors2 <- dplyr::rename(bngal:::phylum_colors, color_order = order)

    dat.out <- full.data.ebc %>%
      select(-ebc_count) %>%
      left_join(full.data.ebc.summ, by = c("sample-id", "edge_btwn_cluster")) %>%
      left_join(phylum_colors2, by = "phylum") %>%
      dplyr::mutate(tax_level = if_else(is.na(tax_level), tax.level, tax_level))

    # merge alpha diversity data (shannon only for now)
    alpha.div <- alpha.div %>%
      filter(tax_level %in% tax.level & index == "shannon")

    # ensure no mismatches between subcommunities and samples (when '--subcommunities' option !is.null)
    out[[x]] <- dat.out %>%
      left_join(., select(metadata, `sample-id`, all_of(sub.comms)), by = "sample-id") %>%
      left_join(alpha.div, by = c("sample-id", "tax_level")) %>%
      dplyr::mutate(remove = if_else(.data[[sub.comms]] == .data$sub_comm, F, T),
                    remove = if_else(is.na(sub_comm), F, remove)) %>%
      filter(remove == FALSE) %>%
      select(-remove)
  }

  out

}
