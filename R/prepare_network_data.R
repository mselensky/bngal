#' Prepare ASV count data for correlation network analysis
#'
#' This function prepares data for correlation network analysis at a specified
#' level of taxonomic classification. Singleton ASVs are first removed from a
#' provided ASV table, then taxa are binned at a defined taxonomic level.
#' Binned data are then filtered either above or below a given abundance cutoff
#' value based on user input.
#'
#' @param asv.table See [`bngal::bin_taxonomy()`]. Must be an absolute abundance ASV table.
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
#' The following options results in a long-form dataframe of counts from `asv_table` binned by
#' prepare_network_data(asv_table, metadata, "family", "greaterThan", 0.01)
#'
prepare_network_data <- function(asv.table, meta.data, tax.level, direction, cutoff.val, corr.cols=NULL) {

  # long-form data and remove singletons for downstream analyses
  if (is.null(corr.cols) == TRUE) {
    asv.long <- asv.table %>%
      pivot_longer(cols = 2:ncol(.), names_to = "taxon", values_to = "count") %>%
      group_by(taxon) %>%
      dplyr::mutate(n_taxa_obs = sum(count)) %>% # id singletons
      filter(n_taxa_obs > 1) %>% # remove singletons
      ungroup() %>%
      group_by(`sample-id`) %>%
      dplyr::mutate(rel_abun = count/sum(count)) %>% # calculate relative abundance
      filter(rel_abun > 0 & !is.na(rel_abun))
  } else {
    # if optional corr.cols argument is supplied, add the following metadata columns to the asv table:
    asv.long <- asv.table %>%
      left_join(., select(metadata, `sample-id`, all_of(corr.cols)), by = "sample-id") %>%
      pivot_longer(cols = 2:ncol(.), names_to = "taxon", values_to = "count") %>%
      dplyr::filter(!is.na(count)) %>%
      group_by(taxon) %>%
      dplyr::mutate(n_taxa_obs = sum(count)) %>% # id singletons
      filter(n_taxa_obs > 1) %>% # remove singletons
      ungroup() %>%
      group_by(`sample-id`) %>%
      dplyr::mutate(rel_abun = (count-min(count))/(max(count)-min(count))) %>% # calculate relative abundance
      dplyr::filter(rel_abun > 0 & !is.na(rel_abun))
  }

  # split taxonomy into different levels
  long_rel_full_tax <- suppressWarnings(asv.long %>%
                                          dplyr::mutate(taxonomy = str_remove_all(taxon, "d__|p__|c__|o__|f__|g__|s__")) %>%
                                          separate(taxonomy,
                                                   c("domain", "phylum", "class", "order", "family", "genus", "asv"),
                                                   sep=';') %>%
                                          dplyr::mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum))
  )

  # bin by defined taxonomic level (tax.level)

  if (tax.level == "asv"){
    # long_rel_binned <- long_rel_full_tax %>%
    #   group_by(`sample-id`, taxon) %>%
    #   dplyr::mutate(rel_abun_binned = rel_abun) %>%
    #   select(`sample-id`, rel_abun_binned, taxon, domain:asv) %>%
    #   distinct() %>%
    #   ungroup() %>%
    #   dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, family, genus, asv, sep = ";"))
    long_rel_binned <- long_rel_full_tax %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, family, genus, asv, sep = ";")) %>%
      group_by(`sample-id`, taxon_) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon_, domain:asv) %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>% ungroup()

  } else if (tax.level == "genus"){
    long_rel_binned <- long_rel_full_tax %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, family, genus, sep = ";")) %>%
      group_by(`sample-id`, taxon_) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon_, domain:genus) %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>% ungroup()
  } else if (tax.level == "family"){
    long_rel_binned <- long_rel_full_tax %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, family, sep = ";")) %>%
      group_by(`sample-id`, taxon_) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon_, domain:family) %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>% ungroup()
  } else if (tax.level == "order"){
    long_rel_binned <- long_rel_full_tax %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, sep = ";")) %>%
      group_by(`sample-id`, taxon_) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon_, domain:order) %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>% ungroup()
  } else if (tax.level == "class"){
    long_rel_binned <- long_rel_full_tax %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, sep = ";")) %>%
      group_by(`sample-id`, taxon_) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon_, domain:class) %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>% ungroup()
  } else if (tax.level == "phylum"){
    long_rel_binned <- long_rel_full_tax %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, sep = ";")) %>%
      group_by(`sample-id`, taxon_) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon_, domain:phylum) %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>% ungroup()
  } else if (tax.level == "domain"){
    long_rel_binned <- long_rel_full_tax %>%
      dplyr::mutate(taxon_ = str_c(domain, sep = ";")) %>%
      group_by(`sample-id`, taxon_) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon_, domain) %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>% ungroup()
  } else {
    stop("\n | [", Sys.time(), "] prepare_data():\n",
    " | Must choose one of the following taxonomic levels:\n",
    " |     'asv', 'genus', 'family', 'order', 'class', 'phylum', 'domain'")

  }

  # rename empty "taxon_" columns to corr_col names
  long_rel_binned <- long_rel_binned %>%
    dplyr::mutate(
      node_type = if_else(is.na(taxon_),
                          true = "env_var", false = "taxon"),
      taxon_ = if_else(is.na(taxon_),
                       true = domain, false = taxon_),
      phylum = if_else(node_type == "env_var",
                       true = "env_var", false = phylum)
    ) %>%
    dplyr::select(`sample-id`, node_type, everything())

  # filter for binned taxa that comprise more or less than x% of the community
  # defined by direction and cutoff.val
  if (direction == "lessThan"){
    long_rel_binned_cutoff <- long_rel_binned %>%
      filter(rel_abun_binned <= cutoff.val & rel_abun_binned > 0)
  } else if (direction == "greaterThan") {
    long_rel_binned_cutoff <- long_rel_binned %>%
      filter(rel_abun_binned >= cutoff.val & rel_abun_binned > 0)
  } else {
    long_rel_binned_cutoff <- long_rel_binned %>%
      filter(rel_abun_binned > 0)
  }

  list(taxonomic_level = tax_level,
       data = long_rel_binned_cutoff,
       metadata = select(meta.data, `sample-id`, everything()))

}
