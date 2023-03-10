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
#'
#' @examples n.cores <- bngal::check_cores()
check_cores <- function() {

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

  NCORES

}
