#' Color network nodes by edge between cluster ID and phylum
#'
#' This function colors network nodes by phylum and by edge between cluster
#' using the defaul `bngal` color schemes. A custom list of phylum colors can
#' optionally be provided (see `phylum.colors` option).
#'
#' @param binned.tax Output from [`bngal::bin_taxonomy()`]
#' @param clusters.to.color Output from [`bngal::get_ebc_clusters()`]
#' @param phylum.colors *Optional* Dataframe with hex color codes assigned to
#' each phylum in the Silva 16S rRNA gene database (v. 138). Default color scheme
#' provided, but if a custom color scheme is desired, columns must be the
#' following:
#' * `phylum`: name of phylum in Silva database
#' * `hex.color`: color hexcode, including `#`
#' * `order`: desired order for downstream plots
#'
#' @return
#' @export
#'
#' @examples
color_nodes <- function(binned.tax, clusters.to.color, phylum.colors) {

  if (missing(phylum.colors)) {
    # R/sysdata.rda contains default color scheme for phyla
    phylum.colors = bngal:::phylum_colors %>%
      dplyr::rename(phylum_order = order)
  } else {
    phylum.colors <- phylum.colors %>%
      dplyr::rename(phylum_order = order)
  }


  full.data <- binned.tax %>%
    left_join(clusters.to.color$nodes, by = c("taxon_" = "label")) %>%
    ungroup() %>%
    select(`sample-id`, taxon_, binned_count) %>%
    pivot_wider(names_from = "taxon_", values_from = "binned_count") %>%
    filter(!is.na(`sample-id`)) %>%
    pivot_longer(cols = 2:ncol(.), names_to = "taxon_", values_to = "binned_count") %>%
    dplyr::mutate(binned_count = if_else(is.na(binned_count), 0, binned_count))

  full.data.ebc <- full.data %>%
    left_join(clusters.to.color$nodes, by = c("taxon_" = "label")) %>%
    dplyr::mutate(group = "all",
                  edge_btwn_cluster = if_else(is.na(edge_btwn_cluster), "no_cluster", as.character(edge_btwn_cluster))) %>%
    distinct(`sample-id`, taxon_, edge_btwn_cluster, group, .keep_all = T) %>%
    ungroup() %>% group_by(`sample-id`, edge_btwn_cluster, group) %>%
    dplyr::mutate(ebc_count = sum(binned_count, na.rm=TRUE))

  dat.out <- full.data.ebc %>%
    distinct(`sample-id`, group, edge_btwn_cluster, ebc_count) %>%
    group_by(`sample-id`) %>%
    dplyr::mutate(ebc_abun_sum = ebc_count/sum(ebc_count),
                  group = if_else(is.na(group), "all", group))

  full.data <- dat.out %>%
    left_join(select(full.data.ebc, -ebc_count), by = c("sample-id", "edge_btwn_cluster", "group"))

  rm(full.data.ebc)

  # this will arrange filled bars by summed EBC abundance
  ebc_arranged <- full.data %>%
    dplyr::arrange(ebc_abun_sum)

  # manually color the top 20 most abundant ebc nodes. greyscale the rest
  top_20_ebcs <- ebc_arranged %>%
    group_by(edge_btwn_cluster) %>%
    dplyr::summarize(sum_ebc_abun = sum(ebc_abun_sum)) %>%
    dplyr::arrange(desc(sum_ebc_abun)) %>%
    filter(!edge_btwn_cluster == "no_cluster") %>%
    slice_max(n = 20, order_by = sum_ebc_abun)

  # default EBC color scheme from bngal
  ebc.colors <- bngal:::ebc_colors

  no_clust <- ebc.colors %>%
    filter(color_name == "no_cluster")
  top_20_ebcs$color_name = as.character(seq(1:nrow(top_20_ebcs)))

  colorz <- top_20_ebcs %>%
    left_join(ebc.colors, by = "color_name") %>%
    full_join(no_clust, by = c("color_name", "hex.code")) %>%
    dplyr::mutate(edge_btwn_cluster = if_else(is.na(edge_btwn_cluster), "no_cluster", edge_btwn_cluster))

  nodes.colored <- clusters.to.color$nodes %>%
    dplyr::mutate(edge_btwn_cluster = as.character(edge_btwn_cluster),
                  shape = if_else(node_type == "env_var",
                                  "square", "dot")) %>%
    left_join(colorz, by = "edge_btwn_cluster") %>%
    dplyr::rename(edge_btwn_cluster_color = hex.code) %>%
    dplyr::mutate(edge_btwn_cluster_color = if_else(is.na(edge_btwn_cluster_color), # color all other taxa nodes black
                                                    "#000000", edge_btwn_cluster_color)) %>%
    left_join(., phylum.colors, by = "phylum") %>%
    dplyr::rename(phylum_color = hex.color)

  list(nodes = nodes.colored, edges = clusters.to.color$edges)
}
