#' Filter preprocessed data for correlation matrix
#'
#' This function filters input data and returns a normalized abundance matrix.
#'
#' The parameter `obs.cutoff` refers to the minimum number of observations
#' per pairwise relationship required for inclusion in the output matrix.
#'
#' @param prepared.data Output from [`bngal::prepare_network_data`] or [`bngal::split_network_data`]
#' @param transformation *Optional* Numeric transformation to apply to data (default = none). `"log10"` accepted.
#' @param obs.cutoff *Optional* Minimum number of observations required per pairwise relationship to be included in correlation matrix (default = `10`).
#'
#' @return
#' @export
#'
#' @examples
prepare_corr_data <- function(prepared.data, obs.cutoff, transformation) {

  # this is formatted for multithreading on a SLURM-directed HPC system,
  # but any *nix-like machine can multithread here as well. otherwise
  # this runs on a single core.
  if (Sys.getenv("SLURM_NTASKS") > 1) {
    NCORES = Sys.getenv("SLURM_NTASKS")
  } else if (parallel::detectCores() > 2) {
    NCORES = parallel::detectCores()-1
  } else {
    NCORES = 1
  }

  comp_corr <- function(prepared_data, transformation, obs.cutoff) {
    if (missing(transformation)) {
      matrix. <- prepared_data %>%
        dplyr::select(`sample-id`, taxon_, norm_vals) %>%
        pivot_wider(names_from = "taxon_",
                    values_from = "norm_vals") %>%
        column_to_rownames("sample-id") %>%
        as.matrix() %>%
        replace_na(., 0)
    } else if (transformation == "log10") {
      matrix. <- prepared_data %>%
        dplyr::select(`sample-id`, taxon_, norm_vals) %>%
        dplyr::mutate(transf = log10(norm_vals+1)) %>%
        dplyr::select(-norm_vals) %>% ungroup() %>%
        pivot_wider(names_from = "taxon_", values_from = "transf") %>%
        column_to_rownames("sample-id") %>%
        as.matrix() %>%
        replace_na(., 0)
    }

    matrix.l <- as.data.frame(matrix.)
    matrix.l <- matrix.l %>%
      rownames_to_column("sample-id") %>%
      pivot_longer(2:ncol(.), names_to = "taxa1", values_to = "norm_vals") %>%
      filter(norm_vals > 0)

    samp_tabs <- split(matrix.l, matrix.l$`sample-id`)
    pw = list()
    message(" | [", Sys.time(), "] Filtering data for correlation matrix...\n",
            " | Possible pairwise '", tax_level, "'-level relationships per sample:\n",
            " | --------------------------------------------------------------------")
    for (sample.name in names(samp_tabs)) {
      pw[[sample.name]] = data.frame(t(combn(samp_tabs[[sample.name]][["taxa1"]], m = 2)))
      names(pw[[sample.name]]) = c("taxa1", "taxa2")
      pw[[sample.name]]$taxa_pair = paste0(pw[[sample.name]]$taxa1, '~<>~', pw[[sample.name]]$taxa2)
      pw[[sample.name]]$`sample-id` = sample.name

      # test3[[sample.name]] = semi_join(pw, samp_tabs[[sample.name]], by =c("taxa1"))
      # test3[[sample.name]]$`sample-id` = sample.name
      # test3[[sample.name]] <- test3[[sample.name]] %>% select(`sample-id`, taxa_pair)
      message(" |   * '", sample.name, "': ", nrow(pw[[sample.name]]))
    }
    message(" | --------------------------------------------------------------------")

    pw_full = suppressMessages(Reduce(full_join, pw))
    message(" | ", length(unique(pw_full$taxa_pair)), " total possible pairwise '", tax_level, "'-level relationships observed across the dataset.")

    pw_joined <- left_join(pw_full, matrix.l, by = c("sample-id", "taxa1")) %>%
      rename(taxa1_norm_val = norm_vals) %>%
      left_join(matrix.l, by = c("sample-id", "taxa2"="taxa1")) %>%
      rename(taxa2_norm_val = norm_vals)

    pw_counts <- pw_joined %>%
      group_by(taxa_pair) %>%
      dplyr::summarize(n_pairs = n()) %>%
      filter(n_pairs >= obs.cutoff)

    pw_filtered <- pw_joined %>%
      semi_join(pw_counts, by ="taxa_pair") %>%
      distinct()

    matrix.out <- pw_filtered %>%
      select(`sample-id`, taxa1, taxa1_norm_val) %>%
      distinct() %>%
      dplyr::mutate(taxa1_norm_val = as.numeric(taxa1_norm_val)) %>%
      pivot_wider(names_from = "taxa1", values_from = "taxa1_norm_val") %>%
      replace(is.na(.), 0) %>%
      column_to_rownames("sample-id")

    message(" | [", Sys.time(), "] Filtered data for correlation matrix:\n",
            " | * Minimum observation threshold: ", obs.cutoff, "\n",
            " | * # of unique pairwise '", tax_level, "'-level relationships observed: ", nrow(pw_counts), "\n",
            " | * # of unique '", tax_level, "'-level taxa involved: ", ncol(matrix.out), "\n",
            " | * # of total observations included: ", sum(pw_counts$n_pairs))

    as.matrix(matrix.out)
  }

  if (!is.null(nrow(prepared.data$data))) {
    comp_corr(prepared.data$data, transformation, obs.cutoff)

  } else {
    dat.in = list()
    for (i in names(prepared.data$data)) {
      message("\n\n | [", Sys.time(), "] Preparing network data for subcommunity '", i, "' ...")
      dat.in[[i]] <- comp_corr(prepared.data$data[[i]], transformation, obs.cutoff)
    }
  }
}
