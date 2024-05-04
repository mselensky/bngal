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
corr_matrix <- function(filtered.matrix, correlation, cores=NULL) {

  NCORES = bngal::check_cores(cores)

  # pairwise dissim. and p-values
  # did prepared.data get split into subcommunities?
  # if so, break subcommunity data preparation steps across NCORES
  if (!is.null(nrow(filtered.matrix$output_matrix))) {
    out <- Hmisc::rcorr(filtered.matrix, type = correlation)
  } else {
    out <- parallel::mclapply(X = filtered.matrix,
                              FUN = function(i){Hmisc::rcorr(i$output_matrix, type = correlation)},
                              mc.cores = NCORES)
    message(" | [", Sys.time(), "] Correlation matrices computed for the following subcommunities:")
    for (i in names(filtered.matrix)) {
      message(" |   --", i)
    }
  }

  out

}
