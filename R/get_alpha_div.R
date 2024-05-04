#' Calculate alpha diversity for a defined level of taxonomic classification
#'
#' @param binned.taxonomy Output from [`bngal::bin_taxonomy()`]
#' @param tax.level Taxonomic level of network data output from either `[bngal::export_network_data]`
#' (if running interactively) or bngal-build-networks.R (if using bngal CLI)
#'
#' @return This function can be applied directly to the output from
#' [bngal::bin_taxonomy()].
#' @export
#'
get_alpha.div <- function (binned.taxonomy, tax.level) {

  rel_mats <- bngal::make_matrix(binned.taxonomy, "rel_abun_binned", "sample-id")

  # define alpha diversity indices to be calculated
  indices = c("shannon", "simpson", "invsimpson")

  # wrap tidying commands
  alpha_diversity <- function(abundance_matrix, index = ".") {
      alpha_diversity_result <- abundance_matrix %>%
        vegan::diversity(., index = index) %>%
        as.data.frame() %>%
        rownames_to_column(var = "sample-id")
      names(alpha_diversity_result) <- c("sample-id", paste0(tax.level, "_", index))
    alpha_diversity_result
  }

  alpha_res_list <- list()
  for (index in indices) {
    alpha_res_list[[index]] <- alpha_diversity(rel_mats, index)
  }

  suppressMessages(
    alpha_res_list %>%
      purrr::reduce(., left_join) %>%
      pivot_longer(2:ncol(.), names_to = "alpha", values_to = "value") %>%
      separate(alpha, into = c("tax_level", "index"), sep = "_")
  )
  # tax_levels = c("phylum", "class", "order", "family", "genus", "asv")
  # rel_mats=list()
  # for (i in tax_levels) {
  #   rel_mats[[i]] <- bngal::make_matrix(binned.taxonomy[[i]], "rel_abun_binned", "sample-id")
  # }
  #
  #
  # # define alpha diversity indices to be calculated
  # indices = c("shannon", "simpson", "invsimpson")
  #
  # # wrap tidying commands
  # alpha_diversity <- function(abundance_matrix, index = ".") {
  #   alpha_diversity_result <- list()
  #   for (i in tax_levels) {
  #     alpha_diversity_result[[i]] <- abundance_matrix[[i]] %>%
  #       vegan::diversity(., index = index) %>%
  #       as.data.frame() %>%
  #       rownames_to_column(var = "sample-id")
  #     names(alpha_diversity_result[[i]]) <- c("sample-id", paste0(i, "_", index))
  #   }
  #   alpha_diversity_result
  # }
  #
  # alpha_res_list <- list()
  # for (index in indices) {
  #   alpha_res_list[[index]] <- alpha_diversity(rel_mats, index)
  # }
  #
  # alpha_res <- list()
  # suppressMessages(
  #   for (i in indices) {
  #     alpha_res[[i]] <- alpha_res_list[[i]] %>%
  #       purrr::reduce(left_join)
  #   }
  # )
  #
  # suppressMessages(
  #   alpha_res_joined <- alpha_res %>%
  #     purrr::reduce(., left_join) %>%
  #     pivot_longer(2:ncol(.), names_to = "alpha", values_to = "value") %>%
  #     separate(alpha, into = c("tax_level", "index"), sep = "_")
  # )

}
