#' Get network node IDs
#'
#' This function yields a dataframe containing each unique, quality-filtered
#' node to be used in constructing a correlation network.
#'
#' @param prepared.data Output from [`bngal::prepare_network_data`]
#' @param corr.matrix Output from [`bngal::prepare_corr_data()`]
#' @return
#' @export
#'
#' @examples
get_node_ids <- function(prepared.data, corr.matrix){

  NCORES <- bngal::check_cores()

  if (any(names(prepared.data) %in% c("taxonomic_level", "data", "metadata"))) {
    prepared.data.df = prepared.data$data
    sub.comms = NULL
    tax_level = prepared.data[["taxonomic_level"]]
  } else {
    prepared.data.df = list()
    sub.comms = names(prepared.data)
    for (i in sub.comms) {
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
        select(taxon_, domain:all_of(tax_level), node_type)
    } else {
      stop("\n | [", Sys.time(), "] get_node_ids(): Must choose one of the following taxonomic levels:\n",
           " | 'asv', 'genus', 'family', 'order', 'class', 'phylum', 'domain'\n")
    }

    nodes. <- nodes. %>%
      dplyr::rename(label = taxon_) %>%
      tibble::rowid_to_column("id") %>%
      dplyr::mutate(id = as.character(id))

    # if (any(names(prepared.data) %in% c("taxonomic_level", "data", "metadata"))) {
    #   nodes. <- nodes. %>%
    #     filter(label %in% rownames(corr.matrix$P))
    # } else {
    #   nodes. <- nodes. %>%
    #     filter(label %in% rownames(corr.matrix[[i]]$P))
    # }

    # nodes. <- nodes. %>%
    #   filter(label %in% rownames(corr.matrix$P))

    nodes.
  }

  #message(" |   --Extracting node IDs from object '", deparse(substitute(prepared.data)), "'\n")

  if (any(names(prepared.data) %in% c("taxonomic_level", "data", "metadata"))) {
    dat.out = get_ids(prepared.data.df, tax_level)
  } else if (any(names(prepared.data[[1]]) %in% c("taxonomic_level", "data", "metadata"))) {
    dat.out = parallel::mclapply(X = prepared.data.df,
                                 FUN = function(i){get_ids(i, tax_level)},
                                 mc.cores = NCORES)

    for (i in sub.comms) {
      dat.out[[i]] <- dat.out[[i]] %>%
        filter(label %in% rownames(corr.matrix[[i]]$P))
    }


  } else {

    stop("\n | [", Sys.time(), "] Unexpected input data.\n",
         " |   ->Is prepared.data the output from bngal::prepare_network_data() or bngal::split_network_data() ?\n",
         " |   ->Is corr.matrix the output from bngal::prepare_corr_data() ?")
  }

  dat.out

}
