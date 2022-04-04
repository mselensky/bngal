#' Get network node IDs
#'
#' This function yields a dataframe containing each unique, quality-filtered
#' node to be used in constructing a correlation network.
#'
#' @param prepared.data Output from [`bngal::prepare_network_data`] or [`bngal::split_network_data`]
#' @param corr.matrix Output from [`bngal::prepare_corr_data()`]
#' @return
#' @export
#'
#' @examples
get_node_ids <- function(prepared.data, corr.matrix){

  prepared_data_df = prepared.data$data

  # this is formatted for multicore processing on a SLURM-directed HPC system,
  # but any *nix-like machine can multithread here as well. otherwise
  # this will run on a single core.
  if (Sys.getenv("SLURM_NTASKS") > 1) {
    NCORES = Sys.getenv("SLURM_NTASKS")
  } else if (parallel::detectCores() > 2) {
    NCORES = parallel::detectCores()-1
  } else {
    NCORES = 1
  }

  get_ids <- function(prepared.data.df, corr.cols) {
    nodes.1 <- prepared.data.df %>%
      as.data.frame() %>%
      distinct(taxon_, .keep_all = TRUE)

    # tax_level corresponds to last column name
    # from bngal::prepare_network_data
    tax_level = names(prepared.data.df[,ncol(prepared.data.df)])
    message(" | [", Sys.time(), "] Taxonomic level detected: '", tax_level, "'")

    taxa.levels = c("domain", "phylum", "class", "order", "family", "genus", "asv")

    if (tax_level %in% taxa.levels) {
      nodes. <- nodes.1 %>%
        select(taxon_, domain:.data[[tax_level]], node_type)
    } else {
      stop("\n | [", Sys.time(), "] get_node_ids(): Must choose one of the following taxonomic levels:\n",
           " | 'asv', 'genus', 'family', 'order', 'class', 'phylum', 'domain'\n")
    }

    nodes. %>%
      dplyr::rename(label = taxon_) %>%
      tibble::rowid_to_column("id") %>%
      dplyr::mutate(id = as.character(id)) %>%
      filter(label %in% rownames(corr.matrix$P))
  }

  if (any(class(prepared_data_df) %in% c("tbl_df", "tbl", "data.frame"))) {
    dat.out = get_ids(prepared_data_df)
    message(
      " |   --Extracting node IDs from object '", deparse(substitute(prepared.data)), "'\n"
    )
  } else if (class(prepared_data_df) == "list"){
    dat.out = parallel::mclapply(X = prepared_data_df,
                                 FUN = function(i){get_ids(i)},
                                 mc.cores = NCORES)
    message(
      " |   --Extracting node IDs from object '", deparse(substitute(prepared.data)), "'\n"
    )

  } else {

    stop("\n | [", Sys.time(), "] Unexpected input data.\n",
         " |   ->Is prepared.data the output from bngal::prepare_network_data() or bngal::split_network_data() ?\n",
         " |   ->Is corr.matrix the output from bngal::prepare_corr_data() ?")
  }

  dat.out

}
