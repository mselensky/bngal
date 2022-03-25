#' Get edge between cluster IDs
#'
#' Wrapper function around [`igraph::membership()`]
#'
#' @param betweenness. Output from [`bngal::get_edge_betweenness()`]
#'
#' @return
#' @export
#'
#' @examples
#' ebcs <- bngal::get_edge_betweenness(igraph_obj)
#' bngal::get_ebc_member_ids(ebcs)
get_ebc_member_ids <- function(betweenness.) {

  #input.data.class = c("tbl_df", "tbl", "data.frame")
  if (class(betweenness.) == "communities") {

    member.list <- igraph::membership(betweenness.)
    member.ids <- data.frame(id = as.character(names(member.list)),
                             edge_btwn_cluster = as.numeric(member.list))

  } else if (class(betweenness.) %in% "list") {

    # no multicore here, probably not worth the overhead
    member.list <- lapply(X = betweenness.,
                          FUN = igraph::membership)

    member.ids <- lapply(X = member.list,
                         FUN = function(i) {
                           data.frame(id = as.character(names(i)),
                                      edge_btwn_cluster = as.numeric(i)
                           )
                         }
    )

  } else {
    stop("
  | Unexpected input data for bngal::get_ebc_member_ids()
  |   ->Requires output from bngal::cluster_edge_betweenness()
         ")
  }

  member.ids

}
