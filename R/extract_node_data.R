#' Import network data
#'
#' @param network.data Output from Output from [`bngal::export_network_data()`]
#'
#' @return
#' @export
#'
#' @examples
extract_node_data <- function(network.data) {

  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")

  if (!is.null(nrow(network.data$nodes$asv))) {
    # list nodes and edges each under "all" for downstream syntax
    network.data.list=list()
    for (i in tax.levels) {
      network.data.list$nodes[[i]][["all"]] = network.data$nodes[[i]]
      network.data.list$edges[[i]][["all"]] = network.data$edges[[i]]
    }
  } else {
    network.data.list=network.data
  }
  rm(network.data)

  # a horrid double for loop to add "tax_level" to each node dataframe for each subnetwork

  network.nodes = list()
  for (i in tax.levels) {
    for (x in names(network.data.list$nodes[[i]])) {
      network.data.list$nodes[[i]][[x]] <- network.data.list$nodes[[i]][[x]] %>%
        dplyr::mutate(tax_level = i,
                      sub_comm = x)

    }

    network.nodes[[i]] <- Reduce(bind_rows, network.data.list$nodes[[i]])

    network.nodes[[i]] <- network.nodes[[i]] %>%
      distinct(label, edge_btwn_cluster, .keep_all = T) %>%
      dplyr::rename(degree=value,
                    taxon_=label)

  }

  Reduce(bind_rows, network.nodes)

}
