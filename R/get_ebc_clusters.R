#' Return taxonomic compositions of each edge between cluster
#'
#' @param prepro.data Output from [`bngal::prepro_net_features()`]
#' @param ebc.member.ids Output from [`bngal::get_ebc_member_ids()`]
#' @param igraph.obj Output from [`bngal::get_igraph()`]
#' @param core.ness See `mode` argument from [`igraph::coreness()`]. Default = `"all"`
#'
#' @return
#' @export
#'
#' @examples get_ebc_clusters()
get_ebc_clusters <- function(prepro.data, ebc.member.ids, igraph.obj, core.ness="all") {

  extract_ebc <- function (prepro.data, ebc.member.ids, igraph.obj, core.ness="all") {

    node.clusters. <- prepro.data$nodes %>%
      left_join(., ebc.member.ids, by = "id") %>%
      dplyr::mutate(core = coreness(igraph.obj, core.ness)) %>%
      ungroup() %>%
      group_by(edge_btwn_cluster) %>%
      dplyr::mutate(taxa_per_cluster = n())

    edges.only <- prepro.data$edges

    list("nodes" = node.clusters.,
         "edges" = edges.only)

  }

  input.data.class = c("tbl_df", "tbl", "data.frame")
  if (any(class(prepro.data$edges) %in% input.data.class) &
      any(class(prepro.data$nodes) %in% input.data.class)) {

    igraph.obj = list(igraph.obj)

    dat.out <- extract_ebc(prepro.data, ebc.member.ids[[1]], igraph.obj[[1]], core.ness)

  } else if (class(prepro.data$edges) %in% "list" &
             class(prepro.data$nodes) %in% "list") {

    # might provide multicore support here if needed, but for now probably not worth the overhead:
    node.clusters. <- list()
    edges.only <- list()
    for (i in names(prepro.data$nodes)) {
      node.clusters.[[i]] <- prepro.data$nodes[[i]] %>%
        left_join(., ebc.member.ids[[i]], by = "id") %>%
        dplyr::mutate(core = coreness(igraph.obj[[i]], core.ness)) %>%
        ungroup() %>%
        group_by(edge_btwn_cluster) %>%
        dplyr::mutate(taxa_per_cluster = n())

      edges.only[[i]] <- prepro.data$edges[[i]]
    }

    dat.out <- list("nodes" = node.clusters.,
                    "edges" = edges.only)

  }

  dat.out
}
