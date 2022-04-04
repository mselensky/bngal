#' Compute correlation matrix
#'
#' Calculate correlation matrix from filtered data. This is essentially a
#' convenience wrapper around `[Hmisc::rcorr()]`.
#'
#' @param filtered.matrix Output from `[bngal::prepare_corr_data]`
#' @param correlation See `[Hmisc::rcorr()]`
#'
#' @return
#' @export
#'
#' @examples
corr_matrix <- function(filtered.matrix, correlation) {

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

  # pairwise dissim. and p-values
  # did prepared.data get split into subcommunities?
  # if so, break subcommunity data preparation steps across NCORES to multithread
  if (!is.null(nrow(filtered.matrix))) {
    Hmisc::rcorr(filtered.matrix, type = correlation)

  } else {
    parallel::mclapply(X = filtered.matrix,
                       FUN = function(i){Hmisc::rcorr(i, type = correlation)},
                       mc.cores = NCORES)

    message(" | [", Sys.time(), "] Correlation matrices computed for the following subcommunities:")
    for (i in names(filtered.matrix)) {
      message(" |   --", i)
    }
  }

}
