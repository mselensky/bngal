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
export_network_data <- function(node.color.data, tax.level, out.dr){

  write_data <- function (node.color.data, subset.values) {
    for (i in subset.values) {

      if (subset.values == "all") {
        node.color.data$nodes$group = as.character(i)
        node.color.data$edges$group = as.character(i)
      } else {
        node.color.data$nodes[[i]]$group = as.character(i)
        node.color.data$edges[[i]]$group = as.character(i)
      }

      filename = paste0(i, "-", tax.level,
                        "-network-data.rds")

      saveRDS(node.color.data, file = filename)
    }
  }

  work.dr <- getwd()
  setwd(out.dr)
  if (!is.null(nrow(node.color.data$nodes))) {
    subset.values = "all"
    write_data(node.color.data, subset.values)
  } else if (class(node.color.data$nodes) == "list") {
    subset.values = names(node.color.data$nodes)
    write_data(node.color.data, subset.values)
  } else {
    setwd(work.dr)
    stop("
 | Unexpected input data for bngal::export_network_data()
 |   ->Requires output from bngal::color_nodes()
         ")
  }
  setwd(work.dr)

}
