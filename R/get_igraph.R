#' Create igraph object
#'
#' This is function is a wrapper around [`igraph::graph_from_data_frame()`] to
#' generate igraph objects for each subcommunity defined earlier in the
#' `bngal` pipeline.
#'
#' @param prepro.data Output from [`bngal::prepro_net_features`]
#' @param directed. See `directed` option from [`igraph::graph_from_data_frame`]
#'
#' @return
#' @export
#'
#' @examples get_igraph(prepped_nodes_edges)
get_igraph <- function(prepro.data, directed. = FALSE){

  input.data.class = c("data.frame")
  if (any(class(prepro.data$nodes) %in% input.data.class)) {
    igraph.list <- igraph::graph_from_data_frame(d = prepro.data$edges,
                                                 vertices = prepro.data$nodes,
                                                 directed = FALSE)
  } else if (class(prepro.data$nodes) == "list") {
    igraph.list = list()
    for (i in names(prepro.data$edges)) {
      igraph.list[[i]] <- igraph::graph_from_data_frame(d = prepro.data$edges[[i]],
                                                        vertices = prepro.data$nodes[[i]],
                                                        directed = FALSE)
    }
  }

  igraph.list
}
