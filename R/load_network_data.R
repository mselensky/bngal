load_network_data <- function (network.dir) {
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")
  network_data <- list()
  for (i in tax.levels) {
    network_data[[i]] <- readRDS(
      paste0(
        network.dir, "/",
        list.files(network.dir, pattern = paste0("^", i, ".*\\.rds$"))
      )
    )
  }
  network_data
}