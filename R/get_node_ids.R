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

  if (any(names(prepared.data) %in% c("taxonomic_level", "data", "metadata"))) {
    prepared.data.df = prepared.data$data
    tax_level = prepared.data[["taxonomic_level"]]
  } else {
    prepared.data.df = list()
    for (i in names(prepared.data)) {
      prepared.data.df[[i]] <- prepared.data[[i]]$data
      tax_level = prepared.data[[i]][["taxonomic_level"]]
    }
  }

  get_ids <- function(prepared.data.df, tax_level) {

    nodes.1 <- prepared.data.df %>%
      as.data.frame() %>%
      distinct(taxon_, .keep_all = TRUE)

    taxa.levels = c("domain", "phylum", "class", "order", "family", "genus", "asv")

    if (tax_level %in% taxa.levels) {
      nodes. <- nodes.1 %>%
        select(taxon_, domain:.data[[tax_level]], node_type)
    } else {
      stop("\n | [", Sys.time(), "] get_node_ids(): Must choose one of the following taxonomic levels:\n",
           " | 'asv', 'genus', 'family', 'order', 'class', 'phylum', 'domain'\n")
    }

    nodes. <- nodes. %>%
      dplyr::rename(label = taxon_) %>%
      tibble::rowid_to_column("id") %>%
      dplyr::mutate(id = as.character(id))

    if (any(names(prepared.data) %in% c("taxonomic_level", "data", "metadata"))) {
      nodes. <- nodes. %>%
        filter(label %in% rownames(corr.matrix$P))
    } else {
      nodes. <- nodes. %>%
        filter(label %in% rownames(corr.matrix[[i]]$P))
    }

    nodes.
  }

  #message(" |   --Extracting node IDs from object '", deparse(substitute(prepared.data)), "'\n")

  if (any(names(prepared.data) %in% c("taxonomic_level", "data", "metadata"))) {
    dat.out = get_ids(prepared.data, tax_level)
  } else if (any(names(prepared.data[[1]]) %in% c("taxonomic_level", "data", "metadata"))) {
    dat.out = parallel::mclapply(X = prepared.data.df,
                                 FUN = function(i){get_ids(i, tax_level)},
                                 mc.cores = NCORES)

  } else {

    stop("\n | [", Sys.time(), "] Unexpected input data.\n",
         " |   ->Is prepared.data the output from bngal::prepare_network_data() or bngal::split_network_data() ?\n",
         " |   ->Is corr.matrix the output from bngal::prepare_corr_data() ?")
  }

  dat.out

}
