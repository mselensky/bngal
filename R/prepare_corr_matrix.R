#' Convert preprocessed abundance data into pairwise correlation matrix
#'
#' This function calculates pairwise correlation coefficients and associated p-values from pre-processed abundance data. This is essentially a convenience wrapper around [`Hmisc::rcorr()`].
#'
#' @param prepared.data Output from [`bngal::prepare_network_data`] or [`bngal::split_network_data`]
#' @param transformation *Optional* Numeric transformation to apply to relative abundance data (default = none). `"log10"` accepted.
#' @param correlation Metric for pairwise comparisons (see [`Hmisc::rcorr()`])
#'
#' @return
#' @export
#'
#' @examples
prepare_corr_matrix <- function(prepared.data, transformation, correlation) {

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

  comp_corr <- function(prepared_data, transformation) {
    if (missing(transformation)) {
      matrix. <- prepared_data %>%
        dplyr::select(`sample-id`, taxon_, rel_abun_binned) %>%
        pivot_wider(names_from = "taxon_",
                    values_from = "rel_abun_binned") %>%
        column_to_rownames("sample-id") %>%
        as.matrix() %>%
        replace_na(., 0)
    } else if (transformation == "log10") {
      matrix. <- prepared_data %>%
        dplyr::select(`sample-id`, taxon_, rel_abun_binned) %>%
        group_by(`sample-id`) %>%
        dplyr::mutate(transf = log10(rel_abun_binned+1)) %>%
        dplyr::select(-rel_abun_binned) %>%
        pivot_wider(names_from = "taxon_", values_from = "transf") %>%
        column_to_rownames("sample-id") %>%
        as.matrix() %>%
        replace_na(., 0)
    }
  }

  # pairwise dissim. and p-values
  # did prepared.data get split into subcommunities?
  # if so, break subcommunity data preparation steps across NCORES to multithread
  if (!is.null(nrow(prepared.data$data))) {
    dat.in <- comp_corr(prepared.data$data, transformation)
    Hmisc::rcorr(dat.in, type = correlation)

  } else {
    dat.in <- parallel::mclapply(X = prepared.data$data,
                                 FUN = function(i){comp_corr(i, transformation)},
                                 mc.cores = NCORES)
    parallel::mclapply(X = dat.in,
                       FUN = function(i){Hmisc::rcorr(i, type = correlation)},
                       mc.cores = NCORES)

          #x
    # pairwise dissim. and p-values
    # parallel::mclapply(X = dat.in,
    #                    FUN = function(i){Hmisc::rcorr(i, correlation)},
    #                    mc.cores = NCORES)
    # map(dat.in, ~Hmisc::rcorr(., type = correlation)) # pairwise dissim. and p-values
  }

  # corr_matrix <- map(rel_matrices, ~Hmisc::rcorr(., type = correlation)) # pairwise dissim. and p-values

}
