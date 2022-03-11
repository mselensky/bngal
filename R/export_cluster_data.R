#' Export filtered edge between cluster abundance data
#'
#' This function returns a list of nodes (i.e., taxa) at each level of
#' taxonomic classification, split by their statistical distributions into
#' four categories listed below (defined by the option ```quantile.cutoff```) :
#'   * `high.core_high.abun` : Nodes with high coreness and high abundance
#'   * `high.core_low.abun`  : Nodes with high coreness and low abundance
#'   * `low.core_high.abun`  : Nodes with low coreness and high abundance
#'   * `low.core_low.abun`   : Nodes with low coreness and low abundance
#'
#' @param ebc.comp.data Output from `bngal::calc_ebc_comps()`
#' @param quantile.cutoff (numeric) Quantile distribution cutoff for coreness/abundance data:
#' \itemize{
#'  \item 1: 0-25%
#'  \item 2: 25-50%
#'  \item 3: 50-75%
#'  \item 4: 75-100%
#'}
#'
#' @return
#' @export
#'
#' @examples export_cluster_data(ebc_composition_data, 4)
export_cluster_data <- function(ebc.comp.data, quantile.cutoff) {
  high.core_high.abun <- list()
  high.core_low.abun <- list()
  low.core_high.abun <- list()
  low.core_low.abun <- list()
  split_df <- list()
  sub.comms <- list()
  core.distr = list()
  abun.distr = list()
  tax_abun.distr = list()

  taxa.of.interest <- function(ebc.comp.data) {
    for (i in tax.levels) {

      split_df[[i]] <- split(ebc.comp.data[[i]], ebc.comp.data[[i]][["group"]])

      for (sub.comm in names(split_df[[i]])) {
        # calculate coreness and abundance distributions of each taxon in each
        # network
        core.distr[[i]][[sub.comm]] <- quantile(split_df[[i]][[sub.comm]][["core"]], na.rm = T)
        abun.distr[[i]][[sub.comm]] <- quantile(split_df[[i]][[sub.comm]][["ebc_abun_sum"]], na.rm = T)
        tax_abun.distr[[i]][[sub.comm]] <- quantile(split_df[[i]][[sub.comm]][["rel_abun_binned"]], na.rm = T)
      }
    }
    # list data, core quantiles, and abundance quantiles
    unfiltered.data <- list(data = split_df, core.quantile = core.distr, abun.quantile = abun.distr)
  }
  # filter data based on defined quantile cutoffs by user
  filter.ebc.data <- function (unfiltered.data, sub.comm, quantile.cutoff) {
    for (i in tax.levels) {
      high.core_high.abun[[i]] <- unfiltered.data[["data"]][[i]][[sub.comm]] %>%
        filter(core >= unfiltered.data[["core.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]
               & ebc_abun_sum >= unfiltered.data[["abun.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]) %>%
        distinct(taxon_, tax_level, core, degree, edge_btwn_cluster) %>%
        dplyr::mutate(taxa_group = "high.core_high.abun")

      high.core_low.abun[[i]] <- unfiltered.data[["data"]][[i]][[sub.comm]] %>%
        filter(core >= unfiltered.data[["core.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]
               & ebc_abun_sum <= unfiltered.data[["abun.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]) %>%
        distinct(taxon_, tax_level, core, degree, edge_btwn_cluster) %>%
        dplyr::mutate(taxa_group = "high.core_low.abun")

      low.core_high.abun[[i]] <- unfiltered.data[["data"]][[i]][[sub.comm]] %>%
        filter(core <= unfiltered.data[["core.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]
               & ebc_abun_sum >= unfiltered.data[["abun.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]) %>%
        distinct(taxon_, tax_level, core, degree, edge_btwn_cluster) %>%
        dplyr::mutate(taxa_group = "low.core_high.abun")

      low.core_low.abun[[i]] <- unfiltered.data[["data"]][[i]][[sub.comm]] %>%
        filter(core <= unfiltered.data[["core.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]
               & ebc_abun_sum <= unfiltered.data[["abun.quantile"]][[i]][[sub.comm]][[quantile.cutoff]]) %>%
        distinct(taxon_, tax_level, core, degree, edge_btwn_cluster) %>%
        dplyr::mutate(taxa_group = "low.core_low.abun")

    }
    list(high.core_high.abun = high.core_high.abun,
         high.core_low.abun = high.core_low.abun,
         low.core_high.abun = low.core_high.abun,
         low.core_low.abun = low.core_low.abun)
  }

  split.data <- taxa.of.interest(ebc.comp.data)

  filtered.data = list()
  for (sub.comm in names(split.data$data[[1]])) {
    filtered.data[[sub.comm]] <- split.data %>%
      filter.ebc.data(., sub.comm, quantile.cutoff)
  }
  filtered.data

  join.dfs <- function (filtered.data, sub.comm) {
    taxa.groups = names(filtered.data[[1]])

    # collapse taxonomic levels into one list for each taxa group in a given
    # subcommunity
    joined.data = list()
    for (taxa.group in taxa.groups) {

      # full.join <- function (x) {
      #   full_join(x, y, by = )
      # }
      #
      joined.data[[taxa.group]] <- purrr::reduce(filtered.data[[sub.comm]][[taxa.group]], full_join)
    }

    joined.data

  }

  data.out = list()
  for (sub.comm in names(filtered.data)) {
    data.out[[sub.comm]] <- join.dfs(filtered.data, sub.comm)
  }

  data.out

}
