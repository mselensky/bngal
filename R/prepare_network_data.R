#' Prepare ASV count data for correlation network analysis
#'
#' This function prepares data for correlation network analysis at a specified
#' level of taxonomic classification. Singleton ASVs are first removed from a
#' provided ASV table, then taxa are binned at a defined taxonomic level.
#' Binned data (from [`bngal::bin_taxonomy()`]) are filtered either above or below
#' a given abundance cutoff value based on user input. Corresponding metadata columns
#' (such as a geochemical parameter for each sample) may be included in the
#' preprocessed network data by specifying the `corr.cols` option.
#' Microbial abundance and environmental metadata (if applicable) are normalized
#' on a scale of 0-1 for downstream correlation matrix construction.
#'
#' @param binned.tax Output from [`bngal::bin_taxonomy()`]. Must be an absolute abundance ASV table.
#' @param meta.data See [`bngal::bin_taxonomy()`]
#' @param tax.level Level of taxonomic classification at which to bin absolute abundance.
#' @param direction Direction for abundance cutoff for data filtering (`"lessThan"`, `"greaterThan"` accepted)
#' @param cutoff.val Cutoff value for binned taxa that comprise more or less than this fraction of a sample's community (values `0` to `1` accepted)
#' @param corr.cols Metadata variables (such as sample environmental data) to include in network correlations. Variables must be supplied as a concatenated string of column names from `meta.data` (e.g., `c('var1','var2')`). Default = `NULL`
#'
#' @return
#' @export
#'
#' @examples
#'
prepare_network_data <- function(binned.tax, meta.data, corr.cols) {

  tax.level = names(binned.tax[ncol(binned.tax)-1])

  if (missing(corr.cols)) {

    long_norm_binned <- binned.tax %>%
      group_by(`sample-id`) %>%
      # normalize seq count data (0 to 1)
      dplyr::mutate(norm_vals = (binned_count-min(binned_count))/(max(binned_count)-min(binned_count)))

  } else {

    long_norm_binned <- binned.tax %>%
      pivot_wider(names_from = "taxon_", values_from = "binned_count") %>%
      left_join(., select(metadata, `sample-id`, all_of(corr.cols)), by = "sample-id") %>%
      # pivot binned count + correlated metadata columns only
      pivot_longer(cols = (which(names(.) == tax.level) + 1):ncol(.), names_to = "taxon_", values_to = "binned_count") %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>%
      filter(!is.na(binned_count)) %>%
      group_by(`sample-id`) %>%
      # normalize environmental and seq count data (0 to 1)
      dplyr::mutate(norm_vals = (binned_count-min(binned_count))/(max(binned_count)-min(binned_count)))

  }

  # rename empty "taxon_" columns to corr_col names
  long_norm_binned <- long_norm_binned %>%
    dplyr::mutate(
      node_type = if_else(is.na(taxon_),
                          true = "env_var", false = "taxon"),
      taxon_ = if_else(is.na(taxon_),
                       true = domain, false = taxon_),
      phylum = if_else(node_type == "env_var",
                       true = "env_var", false = phylum)) %>%
    dplyr::select(`sample-id`, node_type, everything())

  list(taxonomic_level = tax_level,
       data = long_norm_binned,
       metadata = select(meta.data, `sample-id`, everything()))

}
