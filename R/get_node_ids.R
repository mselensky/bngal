#' Get network node IDs
#'
#' This function yields a dataframe containing each unique node to be used in
#' constructing a correlation network.
#'
#' @param prepared.data Output from [`bngal::prepare_network_data`] or [`bngal::split_network_data`]
#' @param corr.cols See [`bngal::prepare_network_data()`]
#' @return
#' @export
#'
#' @examples
get_node_ids <- function(prepared.data, corr.cols){

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
    message(" | Taxonomic level detected: '", tax_level, "'")

    if(tax_level == "domain") {
      nodes. <- nodes.1 %>%
        dplyr::select(taxon_, domain)
    } else if(tax_level == "phylum"){
      nodes. <- nodes.1 %>%
        dplyr::select(taxon_, domain:phylum)
    } else if(tax_level == "class"){
      nodes. <- nodes.1 %>%
        dplyr::select(taxon_, domain:class)
    } else if(tax_level == "order"){
      nodes. <- nodes.1 %>%
        dplyr::select(taxon_, domain:order)
    } else if(tax_level == "family"){
      nodes. <- nodes.1 %>%
        dplyr::select(taxon_, domain:family)
    } else if(tax_level == "genus"){
      nodes. <- nodes.1 %>%
        dplyr::select(taxon_, domain:genus)
    } else if(tax_level == "asv"){
      nodes. <- nodes.1 %>%
        dplyr::select(taxon_, domain:asv)
    } else {
      stop()
      message(" | ERROR, get_node_ids(): Must choose one of the following taxonomic levels:\n",
              " | 'asv', 'genus', 'family', 'order', 'class', 'phylum', 'domain'\n"
      )
    }

    # join "node_type" column if optional corr.cols argument is supplied
    if (!is.null(corr.cols)) {
      nodes. <- nodes. %>%
        left_join(., dplyr::select(nodes.1, taxon_, node_type), by = "taxon_")
    }

    nodes. %>%
      dplyr::rename(label = taxon_) %>%
      tibble::rowid_to_column("id") %>%
      dplyr::mutate(id = as.character(id))
  }

  if (any(class(prepared_data_df) %in% c("tbl_df", "tbl", "data.frame"))) {
    dat.out = get_ids(prepared_data_df, corr.cols)
    message(
      " |   --Extracting node IDs from object '", deparse(substitute(prepared.data)), "'\n"
      )
    } else if (class(prepared_data_df) == "list"){
      dat.out = parallel::mclapply(X = prepared_data_df,
                                   FUN = function(i){get_ids(i, corr.cols)},
                                   mc.cores = NCORES)
      message(
        " |   --Extracting node IDs from object '", deparse(substitute(prepared.data)), "'\n"
        )

  } else {

    stop("
  | Unexpected input data.
  |   ->Is prepared.data the output from bngal::prepare_network_data() or bngal::split_network_data() ?
         ")
  }

  dat.out

}
