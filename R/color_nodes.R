#' Color network nodes by edge between cluster ID and phylum
#'
#' @param clusters.to.color Output from [`bngal::get_ebc_clusters()`]
#' @param phylum.colors *Optional* Dataframe with hex color codes assigned to
#' each phylum in the Silva 16SrRNA gene database (v. 138). Default color scheme
#' provided, but if a custom color scheme is required, columns must be the
#' following:
#' * `phylum`: name of phylum in Silva database
#' * `hex.color`: color hexcode, including `#`
#' * `order`: desired order for downstream plots
#'
#' @return
#' @export
#'
#' @examples
color_nodes <- function(clusters.to.color, phylum.colors) {

  if (missing(phylum.colors)) {
    # R/sysdata.rda contains default color scheme
    phylum.colors = phylum_colors %>%
      dplyr::rename(phylum_order = order)
  } else {
  phylum.colors <- phylum.colors %>%
    dplyr::rename(phylum_order = order)
  }

  get_ebc_colors <- function(clusters.to.color) {
    clusters.to.color.filt <- clusters.to.color %>%
      filter(taxa_per_cluster > 1) %>%
      distinct(edge_btwn_cluster, taxa_per_cluster) %>%
      arrange(desc(taxa_per_cluster)) %>%
      head(n = 12)

    edge_btwn_color <- brewer.pal(nrow(clusters.to.color.filt), "Paired")

    names(edge_btwn_color) <- clusters.to.color.filt$edge_btwn_cluster

    edge_btwn_color <- edge_btwn_color %>%
      as.data.frame() %>%
      rownames_to_column("edge_btwn_cluster") %>%
      dplyr::rename(edge_btwn_cluster_color = ".") %>%
      dplyr::mutate(edge_btwn_cluster = as.numeric(edge_btwn_cluster))

    clusters.to.color %>%
      left_join(edge_btwn_color, by = "edge_btwn_cluster") %>%
      dplyr::mutate(edge_btwn_cluster_color = as.character(edge_btwn_cluster_color),
                    edge_btwn_cluster_color = if_else(is.na(edge_btwn_cluster_color) == TRUE,
                                                      true = "#f2f2f2",
                                                      false = edge_btwn_cluster_color),
                    shape = if_else(node_type == "env_var",
                                    "square", "dot")) %>%
      left_join(., phylum.colors, by = "phylum") %>%
      dplyr::rename(phylum_color = hex.color)

  }

  if (!is.null(nrow(clusters.to.color$nodes))) {

    edge_btwn_color <- get_ebc_colors(clusters.to.color$nodes)

    list(nodes = edge_btwn_color,
         edges = clusters.to.color$edges)



  } else if (class(clusters.to.color$nodes) == "list") {

    edge_btwn_color <- parallel::mclapply(X = clusters.to.color$nodes,
                                          FUN = get_ebc_colors,
                                          mc.cores = as.numeric(NCORES))

    list(nodes = edge_btwn_color,
         edges = clusters.to.color$edges)

  } else {
    stop("\n | [", Sys.time(), "] Unexpected input data for bngal::color_nodes()\n",
    " |   ->Requires output from bngal::get_ebc_clusters()")
  }

}
