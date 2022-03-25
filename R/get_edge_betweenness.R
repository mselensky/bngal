#' Return edge between cluster data
#'
#' This is a wrapper around [`igraph::cluster_edge_betweenness()`]. This can
#' be a memory- and time-intensive operation, especially for lower taxonomic
#' levels. Edge betweenness is processed in parallel for each subcommunity
#' by default to ease computational burden. For larger matrices of a single
#' community, consider using [`parallel::mcmapply()`] to increase computation
#' speed.
#'
#' @param igraph.list Output from [`bngal::get_igraph()`]
#'
#' @return
#' @export
#'
#' @examples
#' get_igraph_output <- bngal::get_igraph(prepro_data)
#' get_edge_betweenness(get_igraph_output)
get_edge_betweenness <- function(igraph.list) {
  # this is formatted for multicore processing on a SLURM-directed HPC system,
  # but any *nix-like machine can multithread here as well. otherwise
  # this will run on a single core.
  if (Sys.getenv("SLURM_NTASKS") > 1) {
    NCORES = Sys.getenv("SLURM_NTASKS")
  } else if (parallel::detectCores() > 2) {
    NCORES = parallel::detectCores()-1
  } else {
    NCORES = 1
  }

  cluster_ebc <- function (igraph.list) {
    igraph::cluster_edge_betweenness(igraph.list)
  }

  input.data.class = c("tbl_df", "tbl", "data.frame")
  if (class(igraph.list) == "list") {
    edge_betweenness <- mclapply(X = igraph.list,
                                 FUN = cluster_ebc,
                                 mc.cores = NCORES)
  } else if (class(igraph.list) %in% input.data.class) {
    edge_betweenness <- cluster_ebc(igraph.list)
    list(all=edge_betweenness)
  } else {
    stop("
    | Unexpected input data.
    |   ->Requires output from bngal::get_igraph()
           ")
  }

}