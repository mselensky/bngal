#' Make matrix from long-form microbiome data
#'
#' @param long_data Output from [`bngal::bin_taxonomy()`]
#' @param count_column_name name of count column
#' @param row_names column by which to name table rows (default = `"sample-id"`)
#'
#' @return
#' @export
#'
#' @examples make_matrix(long_data, "rel_abun", )
make_matrix <- function(long_data, count_column_name, row_names) {
  long_data_filtered <- long_data %>%
    select(all_of(row_names), taxon_, all_of(count_column_name))
  long_data_filtered %>%
    pivot_wider(names_from = "taxon_", values_from = all_of(count_column_name)) %>%
    column_to_rownames(var = row_names) %>%
    as.matrix() %>%
    replace(is.na(.), 0)
}
