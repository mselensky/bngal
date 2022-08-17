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
#' @param out.dr Required. Output directory for pairwise summary data.
#'
#' @return
#' @export
#'
#' @examples prepare_corr_data(prepared.data, obs.cutoff, transformation, out.dr)
prepare_corr_data <- function(prepared.data, obs.cutoff, transformation, out.dr) {

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

  comp_corr <- function(prepared_data, transformation, obs.cutoff, out.dr) {
    if (missing(transformation) | is.null(transformation)) {
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

      message(" |   * '", sample.name, "': ", nrow(pw[[sample.name]]))
    }
    message(" | --------------------------------------------------------------------")

    pw_full = suppressMessages(Reduce(full_join, pw))
    message(" | [", Sys.time(), "] ", length(unique(pw_full$taxa_pair)), " unique '", tax_level, "'-level pairwise relationships observed across the dataset.")

    pw_joined <- left_join(pw_full, matrix.l, by = c("sample-id", "taxa1")) %>%
      dplyr::rename(taxa1_norm_val = norm_vals) %>%
      left_join(matrix.l, by = c("sample-id", "taxa2"="taxa1")) %>%
      dplyr::rename(taxa2_norm_val = norm_vals)

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
            " |   * Minimum observation threshold: ", obs.cutoff, "\n",
            " |   * # of unique pairwise '", tax_level, "'-level relationships passing threshold: ", nrow(pw_counts), "\n",
            " |   * # of unique '", tax_level, "'-level taxa involved: ", ncol(matrix.out), "\n",
            " |   * # of total pairwise observations included: ", sum(pw_counts$n_pairs))

    # export per-sample pairwise summary data to csv file
    all_pw <- pw_full %>%
      group_by(`sample-id`) %>%
      dplyr::summarize(possible_pairwise = n())
    passed_qc <- pw_filtered %>%
      group_by(`sample-id`) %>%
      dplyr::summarize(prefiltered_pairwise = n())

    summ.out <- left_join(all_pw, passed_qc, by = "sample-id")
    summ.out$tax_level = tax_level

    if (!is.null(nrow(prepared.data$data))) {
      write_csv(summ.out, file.path(out.dr, paste0("pairwise_summary_", tax_level, "all.csv")))
    } else {
      write_csv(summ.out, file.path(out.dr, paste0("pairwise_summary_", tax_level, "-", i, ".csv")))
    }


    # final output matrix
    as.matrix(matrix.out)

  }

  if (!is.null(nrow(prepared.data$data))) {
    dat.in <- comp_corr(prepared.data$data, transformation, obs.cutoff, out.dr)

  } else {
    # will add multicore support sometime in the future:
    # parallel::mclapply(X = prepared.data,
    #                    FUN = function(i){comp_corr(i$data, transformation, obs.cutoff, out.dr)},
    #                    mc.cores = NCORES)
    dat.in = list()
    for (i in names(prepared.data)) {
      message("\n | [", Sys.time(), "] Preparing network data for subcommunity '", i, "' ...")
      dat.in[[i]] <- comp_corr(prepared.data[[i]]$data, transformation, obs.cutoff, out.dr)
      message(" |   * Subcommunity analyzed: ", i)
      message(" | --------------------------------------------------------------------")
    }

  }

  for (i in names(dat.in)) {
    if (length(dat.in[[i]]) == 0) {
      dat.in[[i]] <- NULL
      message(" |   * WARNING: subcommunity '", i, "' removed from '", tax_level, "'-level network analysis due to lack of data after quality filtering.")
    }
  }
  dat.in

}
