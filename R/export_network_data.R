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
export_network_data <- function(node.color.data, tax.level, out.dr) {

  work.dr <- getwd()
  setwd(out.dr)
  filename = paste0(tax.level, "_network-data.rds")
  saveRDS(node.color.data, file = filename)
  setwd(work.dr)

}
