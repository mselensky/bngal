#' Prepare ASV count data for correlation network analysis
#'
#' This function prepares data for correlation network analysis at a specified
#' level of taxonomic classification. Singleton ASVs are first removed from a
#' provided binned data (from [`bngal::bin_taxonomy()`]) that are filtered either above or below
#' a given abundance cutoff value based on user input. Corresponding metadata columns
#' (such as a geochemical parameter for each sample) may be included in the
#' preprocessed network data by specifying the `corr.cols` option.
#' Microbial abundance and environmental metadata (if applicable) are normalized
#' on a scale of 0-1 for downstream correlation matrix construction.
#'
#' @param binned.tax Output from [`bngal::bin_taxonomy()`]. Must be an absolute abundance ASV table.
#' @param meta.data See [`bngal::bin_taxonomy()`]
#' @param corr.cols Metadata variables (such as sample environmental data) to include in network correlations. Variables must be supplied as a concatenated string of column names from `meta.data` (e.g., `c('var1','var2')`). Default = `NULL`
#'
#' @return
#' @export
#'
#' @examples
#'
prepare_network_data <- function(binned.tax, meta.data, corr.cols, sub.comms) {

  # define subfunction for optional parallelization

  prepare.data <- function(binned.tax, meta.data, corr.cols) {

  if (missing(corr.cols) | is.null(corr.cols)) {

    long_norm_binned <- binned.tax %>%
      group_by(`sample-id`) %>%
      # normalize seq count data (0 to 1)
      dplyr::mutate(norm_vals = (binned_count-min(binned_count))/(max(binned_count)-min(binned_count)))

  } else {

    env.corrs <- select(meta.data, `sample-id`, all_of(corr.cols))
    # env.corrs <- na.omit(env.corrs)
    suppressWarnings(
      env.corrs <- env.corrs %>%
        filter(`sample-id` %in% unique(binned.tax$`sample-id`)) %>%
        dplyr::mutate(across(.cols = all_of(corr.cols), as.numeric)) %>%
        pivot_longer(cols = 2:ncol(.), names_to = "env_var", values_to = "val") %>%
        filter(!is.na(val)) %>%
        group_by(env_var) %>%
        dplyr::mutate(norm_vals = (val-min(val))/(max(val)-min(val))) %>%
        select(-val) %>%
        pivot_wider(names_from = "env_var", values_from = "norm_vals")
    )

    long_norm_binned <- binned.tax %>%
      select(`sample-id`, taxon_, binned_count) %>%
      pivot_wider(names_from = "taxon_", values_from = "binned_count") %>%
      inner_join(., env.corrs, by = "sample-id") %>%
      pivot_longer(cols = 2:ncol(.), names_to = "taxon_", values_to = "binned_count") %>%
      distinct(`sample-id`, taxon_, .keep_all = TRUE) %>%
      filter(!is.na(binned_count)) %>%
      group_by(`sample-id`) %>%
      # normalize environmental and seq count data (0 to 1)
      dplyr::mutate(norm_vals = (binned_count-min(binned_count))/(max(binned_count)-min(binned_count)))

    # rejoin separated taxonomies
    taxonomy_list <- unique(binned.tax[3:ncol(binned.tax)])

    long_norm_binned <- long_norm_binned %>%
      left_join(taxonomy_list, by = "taxon_")

    # warn user that the following samples dropped due to lack of metadata
    # (this is for possible future support, currently the code keeps all samples regardless of metadata density)
    all.samples = tibble(`sample-id` = unique(binned.tax$`sample-id`))
    removed.samples = anti_join(all.samples, env.corrs, by = "sample-id")
    if (nrow(removed.samples) > 0) {
      message(" | [", Sys.time(), "] The following samples were dropped from the network analysis due to a lack of metadata:")
      for (i in unique(removed.samples$`sample-id`)) {
        message(" | -- ", i)
      }
    }

  }

  # rename empty "taxon_" columns to corr_col names
  long_norm_binned <- long_norm_binned %>%
    dplyr::mutate(
      node_type = if_else(is.na(domain),
                          true = "env_var", false = "taxon"),
      taxon_ = if_else(is.na(taxon_),
                       true = domain, false = taxon_),
      phylum = if_else(node_type == "env_var",
                       true = "env_var", false = phylum)) %>%
    dplyr::select(`sample-id`, node_type, everything())

  list(taxonomic_level = tax_level,
       data = long_norm_binned,
       metadata = meta.data %>%
         filter(`sample-id` %in% unique(long_norm_binned$`sample-id`)))
  }


  # optionally split dataframe into list of subcommunities based on provided metadata column
  if (!is.null(sub.comms)) {

    # this is formatted for multithreading on a SLURM-directed HPC system,
    # but any *nix-like machine can multithread here as well. otherwise
    # this runs on a single core.
    if (Sys.getenv("SLURM_NTASKS") >= 1) {
      NCORES = Sys.getenv("SLURM_NTASKS")
    } else if (parallel::detectCores() > 2) {
      NCORES = parallel::detectCores()-1
    } else {
      NCORES = 1
    }

    df <- binned.tax %>%
      left_join(., select(meta.data, `sample-id`, all_of(sub.comms)), by = "sample-id")
    split.df <- split(df, df[[sub.comms]])
    for (i in names(split.df)) {
      split.df[[i]] <- split.df[[i]] %>%
        dplyr::select(-all_of(sub.comms))
    }

    parallel::mclapply(X = split.df,
                       FUN = function(i){prepare.data(i, meta.data, corr.cols)},
                       mc.cores = NCORES)

  } else { # otherwise run on single core
    prepare.data(binned.tax, meta.data, corr.cols)
  }



}
