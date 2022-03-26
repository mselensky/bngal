#' Split prepared data for separate networks
#'
#' This function produces a list of dataframes split by a categorical metadata variable
#' defined by `subset.column`. Comparisons within each `subset.column` group are dropped
#' if the number of observations is less than 5. Each split dataframe represents a
#' subcommunity from which separate networks will be constructed in downstream functions.
#'
#' @param prepared.data Output from [`bngal::prepare_network_data()`].
#' @param subset.column Metadata variable by which `prepared.data` is split into separate dataframes for separate network analyses.
#'
#' @return
#' @export
#'
#' @examples split_network_data(prepared_data, "zone")
split_network_data <- function(prepared.data, subset.column) {

  if (missing(subset.column)) {

    return(prepared.data)
    warning(" | [", Sys.time(), "] No subcommunities defined; returning original data.\n",
            " |   ->you can define subcommunities via subset.column!\n",
            " |   ->if you want to analyze the entire community, consider using the\n",
            " |     output of prepare_network_data() downstream instead.\n")
    stop("\n | Execution halted; no subcommunities defined. You can still use\n",
         " | this output for downstream analysis on the full community.")

  } else {

    prepared.data.df = Reduce(left_join, prepared.data[-1])

    # filter samples where n_obs for a given subset.column group > 4 (for spearman calcs)
    network.groups <- prepared.data.df %>%
      dplyr::distinct(`sample-id`, .data[[subset.column]]) %>%
      group_by(.data[[subset.column]]) %>%
      dplyr::summarise(n_obs = n()) %>%
      filter(n_obs > 4) %>% ungroup()

    # filter dataframe for samples in groups with >4 observations and group split
    unsplit.data <- prepared.data.df %>%
      filter(.data[[subset.column]] %in% pull(network.groups, .data[[subset.column]])) %>%
      group_by(.data[[subset.column]]) %>%
      select(`sample-id`, node_type, rel_abun_binned, taxon_, domain:.data[[prepared.data$taxonomic_level]],
             .data[[subset.column]])

      out.data = split(unsplit.data, f = as.factor(unsplit.data[[subset.column]]))

      out.data = parallel::mclapply(X = out.data,
                                    FUN = function(i){
                                      i = i[1:(length(i)-1)]
                                      }
                                    )

        message(" | Subcommunities defined for:")
      for (i in names(out.data)) {
        message(" |   --", i)
      }

      list(data = out.data,
           metadata = prepared.data[["metadata"]]
          )
  }

}
