#' Bin ASV count table at specified level of taxonomic classification.
#'
#' @param asv.table ASV count table named by Silva-138 L7 taxonomies. Ideally rarefied and filtered as necessary.
#'   First column must be named "`sample-id`" and must contain unique identifiers.
#' @param meta.data Sample metadata corresponding to asv.table.
#' @param tax.level Level of taxonomic classification at which to calculate relative abundance.
#'
#' @return
#' @export
#'
#' @examples bin_taxonomy(asv_table, metadata, "class")
bin_taxonomy <- function(asv.table, meta.data, tax.level) {
  # convert table to long-form dataframe, removing singletons
  asv.long <- asv.table %>%
    pivot_longer(cols = 2:ncol(.), names_to = "taxon", values_to = "count") %>%
    group_by(taxon) %>%
    filter(count > 0) %>%
    dplyr::mutate(n_taxa_obs = sum(count)) %>% # id singletons
    filter(n_taxa_obs > 1) %>% # remove singletons
    ungroup() %>%
    group_by(`sample-id`) %>%
    dplyr::mutate(rel_abun = count/sum(count)) %>% # calculate relative abundance
    filter(rel_abun > 0 & !is.na(rel_abun))

  # split taxonomy into different levels
  long_rel_full_tax <- asv.long %>%
    dplyr::mutate(taxonomy = str_remove_all(taxon, "d__|p__|c__|o__|f__|g__|s__")) %>%
    separate(taxonomy, sep=';',
             c("domain", "phylum", "class", "order", "family", "genus", "species")) %>%
    dplyr::mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum))

  # bin by defined taxonomic level (tax_level)
  if (tax.level == "asv"){
    long_rel_binned <- long_rel_full_tax %>%
      group_by(`sample-id`, taxon) %>%
      dplyr::mutate(rel_abun_binned = rel_abun) %>%
      select(`sample-id`, rel_abun_binned, taxon, domain:species) %>%
      distinct(`sample-id`, taxon, .keep_all = TRUE) %>%
      ungroup() %>%
      dplyr::rename(taxon_ = taxon)
  } else if (tax.level == "genus"){
    long_rel_binned <- long_rel_full_tax %>%
      group_by(`sample-id`, domain, phylum, class, order, family, genus) %>%
      dplyr::mutate(rel_abun_binned = sum(rel_abun)) %>%
      select(`sample-id`, rel_abun_binned, domain:genus) %>%
      distinct(`sample-id`, domain, phylum, class, order, family, genus, .keep_all = TRUE) %>%
      ungroup() %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, family, genus, sep = ";"))
  } else if (tax.level == "family"){
    long_rel_binned <- long_rel_full_tax %>%
      group_by(`sample-id`, domain, phylum, class, order, family) %>%
      dplyr::mutate(rel_abun_binned = sum(rel_abun)) %>%
      select(`sample-id`, rel_abun_binned, domain:family) %>%
      distinct(`sample-id`, domain, phylum, class, order, family, .keep_all = TRUE) %>%
      ungroup() %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, family, sep = ";"))
  } else if (tax.level == "order"){
    long_rel_binned <- long_rel_full_tax %>%
      group_by(`sample-id`, domain, phylum, class, order) %>%
      dplyr::mutate(rel_abun_binned = sum(rel_abun)) %>%
      select(`sample-id`, rel_abun_binned, domain:order) %>%
      distinct(`sample-id`, domain, phylum, class, order, .keep_all = TRUE) %>%
      ungroup() %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, order, sep = ";"))
  } else if (tax.level == "class"){
    long_rel_binned <- long_rel_full_tax %>%
      group_by(`sample-id`, domain, phylum, class) %>%
      dplyr::mutate(rel_abun_binned = sum(rel_abun)) %>%
      select(`sample-id`, rel_abun_binned, domain:class) %>%
      distinct(`sample-id`, domain, phylum, class, .keep_all = TRUE) %>%
      ungroup() %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, class, sep = ";"))
  } else if (tax.level == "phylum"){
    long_rel_binned <- long_rel_full_tax %>%
      group_by(`sample-id`, domain, phylum) %>%
      dplyr::mutate(rel_abun_binned = sum(rel_abun)) %>%
      select(`sample-id`, rel_abun_binned, domain:phylum) %>%
      distinct(`sample-id`, domain, phylum, .keep_all = TRUE) %>%
      ungroup() %>%
      dplyr::mutate(taxon_ = str_c(domain, phylum, sep = ";"))
  } else if (tax.level == "domain"){
    long_rel_binned <- long_rel_full_tax %>%
      group_by(`sample-id`, domain) %>%
      dplyr::mutate(rel_abun_binned = sum(rel_abun)) %>%
      select(`sample-id`, rel_abun_binned, domain) %>%
      distinct(`sample-id`, domain, .keep_all = TRUE) %>%
      ungroup() %>%
      dplyr::mutate(taxon_ = str_c(domain, sep = ";"),
                    phylum = taxon_)
  } else {
    message("
----------
ERROR, bin_taxonomy(): Must choose one of the following taxonomic levels:
                       'asv', 'genus', 'family', 'order', 'class', 'phylum', 'domain'
----------
          ")
  }
}
