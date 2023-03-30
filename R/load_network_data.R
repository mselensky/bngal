#' Load exported network data
#'
#' @param network.dir Full path to network data output from either `[bngal::export_network_data]`
#' (if running interactively) or bngal-build-networks.R (if using bngal CLI)
#' @param tax.level Taxonomic level of network data output from either `[bngal::export_network_data]`
#' (if running interactively) or bngal-build-networks.R (if using bngal CLI)
#'
#' @return
#' @export
#'
#' @examples
load_network_data <- function (network.dir, tax.level) {

  network_data=list()
  out.data=list()
  network_data[[tax.level]] <- readRDS(paste0(
    network.dir, "/",
    list.files(network.dir, pattern = paste0("^.*", tax.level, ".*\\.rds$"))
  ))

  out.data$nodes[[tax.level]] = network_data[[tax.level]]$nodes
  out.data$edges[[tax.level]] = network_data[[tax.level]]$edges
  out.data

  # tax.levels <- c("asv", "genus", "family", "order", "class", "phylum")
  # network_data <- list()
  # a=list()
  # for (i in tax.levels) {
  #   network_data[[i]] <- readRDS(
  #     paste0(
  #       network.dir, "/",
  #       list.files(network.dir, pattern = paste0("^.*", i, ".*\\.rds$"))
  #     )
  #   )
  #
  #   a$nodes[[i]] = network_data[[i]]$nodes
  #   a$edges[[i]] = network_data[[i]]$edges
  # }
  # a
}
