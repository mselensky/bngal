#' Bin ASV count table at a specified level of taxonomic classification
#'
#' This is the primary data formatting function of most `bngal` pipelines,
#' including the first step of network analysis.
#'
#' @param asv.table ASV count table named by Silva-138 L7 taxonomies. Ideally rarefied and quality filtered as necessary, though this script automatically removes singletons in the dataset.
#'   First column must be named "`sample-id`" and must contain unique identifiers.
#' @param meta.data Sample metadata corresponding to asv.table.
#' @param tax.level Level of taxonomic classification at which to calculate relative abundance.
#' @param remove.singletons *Optional* Drop singleton ASV observations from the dataset. Default = `TRUE`
#' @param cutoff.val *Optional* Cutoff value for ASVs that comprise more or less than this fraction of a sample's community (values `0` to `1` accepted). Required if 'direction' specified.
#' @param direction *Optional* Return binned data `'lessThan'` or `'greaterThan'` than cutoff.val. Required if 'cutoff.val' specified.
#' @param compositional *Optional* Determines whether output returns relative abundance (`TRUE`) or binned count data (`FALSE`). Default = `TRUE`.
#'
#' @return
#' @export
#'
#' @examples
#' # return long-form class-level count data
#' bin_taxonomy(asv_table, metadata, "class", compositional = FALSE)
#'
#' # return long-form family-level relative abundance data
#' bin_taxonomy(asv_table, metadata, "family")
#'
bin_taxonomy <- function(asv.table, meta.data, tax.level, remove.singletons=TRUE, cutoff.val, direction, compositional = TRUE) {

  taxa_levels = c("domain", "phylum", "class", "order", "family", "genus", "asv")

  if (!tax.level %in% taxa_levels) {
    stop("\n | [", Sys.time(), "] Must choose one of the following taxonomic levels:\n",
         " |                       'asv', 'genus', 'family', 'order', 'class', 'phylum', 'domain';\n",
         " |                       not '", tax.level, "'")
  }

  # convert table to long-form dataframe, removing singletons
  asv.long <- asv.table %>%
    pivot_longer(cols = 2:ncol(.), names_to = "taxon", values_to = "count") %>%
    filter(count > 0) %>%
    group_by(taxon) %>%
    dplyr::mutate(n_taxa_obs = sum(count)) # id ASV singletons

  if (remove.singletons == TRUE) {

    asv.long <- asv.long %>%
      filter(n_taxa_obs > 1 & !is.na(count)) # remove ASV singletons

  } else {
    asv.long <- asv.long %>%
      filter(!is.na(count))
  }

  if (!missing(direction)) {

    # if direction is specified but cutoff.val isn't, stop
    if (missing(cutoff.val)){
      stop("\n | [", Sys.time(), "] Argument mismatch! `cutoff.val` missing\n",
           " |                         * `cutoff.val` is required if `direction` is specified.\n",
           " |                         * see `?bngal::bin_taxonomy` for more details.")
    }


    if (!missing(cutoff.val) & direction == "lessThan") {
      asv.long <- asv.long %>%
        group_by(`sample-id`) %>%
        filter(count/sum(count) <= cutoff.val & count > 0)
    } else if (!missing(cutoff.val) & direction == "greaterThan") {
      asv.long <- asv.long %>%
        group_by(`sample-id`) %>%
        filter(count/sum(count) >= cutoff.val & count > 0)
    }
  }

  # split taxonomy into different levels
  taxa.levels = data.frame(tax_name = taxa_levels,
                           tax_rank = seq(1:length(taxa_levels)),
                           prefix = c("d__", "p__", "c__", "o__", "f__", "g__", "s__"))
  input_rank = taxa.levels[taxa.levels$tax_name == tax.level,]$tax_rank
  taxa.levels <- taxa.levels %>%
    dplyr::slice(., 1:input_rank)
  long_rel_full_tax <- asv.long %>%
    dplyr::mutate(taxonomy = str_remove_all(taxon, paste0(taxa.levels$prefix, collapse = "|"))) %>%
    tidyr::separate(taxonomy, sep = ";", into = taxa.levels$tax_name, extra = "drop", fill = "right")

  if (!tax.level %in% c('domain', 'phylum')) {
    long_rel_full_tax <- long_rel_full_tax %>%
      dplyr::mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum)) # we are opinionated
  }

  # default = calculate relative abundance; else bin counts (for e.g. prepare_network_data())
  long_binned <- long_rel_full_tax %>%
    group_by(`sample-id`, across(domain:.data[[tax.level]])) %>%
    dplyr::mutate(binned_count = sum(count)) %>%
    select(`sample-id`, binned_count, domain:.data[[tax.level]]) %>%
    distinct(`sample-id`, across(domain:.data[[tax.level]]), .keep_all = TRUE) %>%
    ungroup()

  if (compositional == TRUE) {
    long_binned <- long_binned %>%
      group_by(`sample-id`) %>%
      dplyr::mutate(rel_abun_binned = binned_count/sum(binned_count))
  } else if (compositional == FALSE) {
    long_binned <- long_binned
  } else {
    stop("\n | [", Sys.time(), "] Invalid selection for argument 'compositional': TRUE or FALSE accepted")
  }

  # bin taxonomic names for provided level of classification
  taxa_cols=list()
  for (i in taxa_levels) {
    taxa_cols[[i]] <- paste(long_binned[[i]], sep=";")
    if (i == tax.level) break
  }
  taxon_ = Reduce(paste, taxa_cols)
  taxon_ = gsub(" ", ";", taxon_)
  long_binned$taxon_ = taxon_

  long_binned

}
