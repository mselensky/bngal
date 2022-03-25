#' Generate network edges
#'
#' @param corr.matrix Output from [`bngal::prepare_corr_matrix()`]
#' @param correlation Metric for pairwise comparisions (see [`Hmisc::rcorr()`])
#' @param node.ids Output from [`bngal::get_node_ids()`]
#'
#' @return
#' @export
#'
#' @examples generate_edges(corr_matrix, "spearman", nodes)
generate_edges <- function(corr.matrix, correlation, node.ids) {

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

  # define function to extract rho and pval for given correlation matrix
  get_rho_pval <- function(corr.matrix){
    rho.pval <- list(
      "r" = corr.matrix[["r"]] %>% # class data.frame w/ pairwise correlation values
        as.data.frame() %>%
        tibble::rownames_to_column(var = "taxon1") %>%
        tidyr::pivot_longer(cols = 2:ncol(.), names_to = "taxon2", values_to = correlation) %>%
        dplyr::mutate(test = (taxon1 == taxon2)) %>%
        dplyr::filter(!test == "TRUE") %>% # remove if ID of taxon1 = taxon2
        dplyr::select(-test),

      "P" = corr.matrix[["P"]] %>% # class data.frame w/ pairwise p-values
        as.data.frame() %>%
        tibble::rownames_to_column(var = "taxon1") %>%
        tidyr::pivot_longer(cols = 2:ncol(.), names_to = "taxon2", values_to = "p_value") %>%
        dplyr::mutate(is_na = is.na(.$p_value)) %>%
        dplyr::filter(!is_na == "TRUE") %>% # remove if pairwise comparison of same taxon
        dplyr::select(-is_na)
    )

    purrr::reduce(rho.pval, left_join, by = c("taxon1", "taxon2"))

  }
  join_ids <- function(rho.pval, node.ids) {
    rho.pval %>%
      dplyr::left_join(., select(node.ids, id, label), by = c("taxon1" = "label")) %>%
      dplyr::left_join(., select(node.ids, id, label), by = c("taxon2" = "label")) %>%
      dplyr::rename(to = id.x,
                    from = id.y)
  }

  # if (class(corr.matrix) != "list") {
  #   corr_matrix <- list(corr.matrix)
  # }

  # extract rho and pvals and join with node ids to identify edges in network
  # rho_pval <- map(corr_matrix, get_rho_pval)
  # edges. <- map(rho_pval, ~join_ids(., node_ids))
  corr.matrix=corr_matrix
  correlation="spearman"
  if (class(corr.matrix) != "list") {
    rho_pval <- get_rho_pval(corr.matrix)
    edges. <- join_ids(rho_pval, node_ids)
  } else if (class(corr.matrix) == "list") {
    rho_pval <- parallel::mclapply(X = corr.matrix,
                                   FUN = get_rho_pval,
                                   mc.cores = NCORES)
    edges. = list()
    for (i in names(rho_pval)) {
      edges.[[i]] <- join_ids(rho_pval[[i]], node_ids[[i]])
    }
  }
  edges.

}
