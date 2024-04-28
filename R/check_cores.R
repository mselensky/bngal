#' Check number of available CPUs
#'
#' This function checks for the number of available CPUs (cores) available
#' on the user's machine. This is formatted for parallel processing on a
#' SLURM-directed HPC system, but any *nix-like machine can multithread here
#' as well. If a Windows OS is detected, this function automatically returns `1`
#' to run on a single core.
#'
#' It is recommended to save this value to be called from multithreaded bngal
#' functions such as `bngal::color_nodes`
#'
#' @return
#' @export
#' @param num.cores User-provided number of cores to use. If the variable is left as NULL (or missing entirely), bngal will automatically detect all available cores on the (non-Windows) machine. If a Windows platform is detected, num.cores is set to 1.
#'
#' @examples n.cores <- bngal::check_cores()
check_cores <- function(num.cores = NULL) {

  if (missing(num.cores) | is.null(num.cores)) {
    if (grepl(Sys.info()["sysname"], "windows", ignore.case = TRUE)) {
      NCORES = 1
    } else {

      if (Sys.getenv("SLURM_NTASKS") > 1) {
        NCORES = Sys.getenv("SLURM_NTASKS")
      } else if (parallel::detectCores() > 1) {
        NCORES = parallel::detectCores()
      } else {
        NCORES = 1
      }
    }
  } else {
    NCORES = num.cores
  }

  NCORES

}
