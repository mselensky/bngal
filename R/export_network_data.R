#' Export network node and ebc cluster data
#'
#' @param node.color.data Output from [`bngal::color_nodes()`]
#' @param tax.level See [`bngal::prepare_network_data`]
#' @param out.dr Data output directory (does not have to exists)
#'
#' @return
#' @export
#'
#' @examples
export_network_data <- function(node.color.data, tax.level, subset.column, out.dr){

  if (missing(subset.column)) {
    subset.column = "all"
  }

  write_data <- function (node.color.data, subset.column) {

      filename = paste0(subset.column, "_", tax.level, "_network-data.rds")

      saveRDS(node.color.data, file = filename)
  }

  work.dr <- getwd()
  setwd(out.dr)
  if (!is.null(nrow(node.color.data$nodes))) {
    subset.column = "all"
    write_data(node.color.data, subset.column)
  } else if (class(node.color.data$nodes) == "list") {
    subset.column = subset.column
    write_data(node.color.data, subset.column)
  } else {
    setwd(work.dr)
    stop("\n | [", Sys.time(), "] Unexpected input data for bngal::export_network_data() \n",
         " |   ->Requires output from bngal::color_nodes()
         ")
  }
  setwd(work.dr)

}
