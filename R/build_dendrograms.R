#' Build dendrograms
#'
#' This function performs hierarchical clustering on microbial relative
#' abundance data binned at a specified level of taxonomic classification.
#'
#' @param binned.taxonomy Output from [bngal::bin_taxonomy()]
#' @param metadata Sample metadata
#' @param color.by Metadata column by which to color
#' @param trans Transformation to apply to relative abundance data (default = none). Can be one of `"log10"`, `"log"`, `"sqrt"`, or `"none"`
#'
#' @return This function can be applied directly to the output from
#' [bngal::bin_taxonomy()].
#' A list is returned:
#'   * tip_data: underlying dendrogram data; can be used in [bngal::plot_core_abun()]
#'   * hclust_plots: hierarchical clustering plot\cr
#'   * legends: legend corresponding to hclust_plot\cr
#' @export
#'
#' @examples build_dendrograms(binned_tax, metadata, "sample_type", "sqrt")
build_dendrograms <- function(binned.taxonomy, metadata, color.by, trans="log10") {
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")
  # create relative abundance matrices from binned_tax
  rel_abun_mats <- list()
  for (i in tax.levels) {
    rel_abun_mats[[i]] <- binned.taxonomy[[i]] %>%
      make_matrix("rel_abun_binned", "sample-id") %>%
      replace_na(replace = 0)
  }
  rel_abun_dist <- list()
  hclust_res <- list()
  for (i in tax.levels) {
    # calculate manhattan due to abundance of 0s
    if (trans=="log10") {
      rel_abun_dist[[i]] <- dist(log10(rel_abun_mats[[i]]+1), method = "manhattan")
      } else if (trans=="sqrt") {
        rel_abun_dist[[i]] <- dist(sqrt(rel_abun_mats[[i]]+1), method = "manhattan")
      } else if (trans=="log") {
        rel_abun_dist[[i]] <- dist(log(rel_abun_mats[[i]]+1), method = "manhattan")
      } else if (trans=="none") {
        rel_abun_dist[[i]] <- dist(rel_abun_mats[[i]], method = "manhattan")
      } else {
        message("
        -----ERROR-| build_dendrogram():
                     'trans' variable must be missing or one of
                     'log10', 'log', 'sqrt', 'none'.
                ")}

    # generate dendrogram from distance matrix
    hclust_res[[i]] <- hclust(rel_abun_dist[[i]], method = "complete")
    }

  ggt <- list()
  ggt_df <- list()
  for (i in tax.levels) {
    # use ggtree to manipulate dendrogram as needed for aesthetics
    ggt[[i]] <- ggtree(hclust_res[[i]])

    ggt_df[[i]] <- get.tree(ggt[[i]])$tip.label %>%
      as.data.frame() %>%
      rename(`sample-id` = ".") %>%
      left_join(metadata)
  }

  hclust_plots <- list()
  for (i in tax.levels) {
    hclust_plots[[i]] <- ggt[[i]] %<+% ggt_df[[i]] +
      #geom_tiplab(aes(text = label, angle = 90)) +
      geom_tippoint(aes(color = .data[[color.by]])) +
      coord_flip() +
      scale_y_reverse()
  }

  # export hclust_plot legend as separate object
  legends <- list()
  for (i in tax.levels) {
    legends[[i]] <- hclust_plots[[i]] %>%
      ggpubr::get_legend()
  }

  ordered_names <- list()
  for (i in tax.levels) {
    # ordered names for dendrogram
    ordered_names[[i]] <- tibble(`sample-id` = get_taxa_name(hclust_plots[[i]]))
  }

  plot_label_data <- list()
  for (i in tax.levels) {
    # merge `ordered_names` with ggt[[i]]$data to join `y` values,
    # which is what we need for ordering figure
    plot_label_data[[i]] <- ggt[[i]]$data %>%
      rename(`sample-id` = label) %>%
      filter(!is.na(`sample-id`))
  }
  merged_labs <- list()
  for (i in tax.levels) {
    merged_labs[[i]] <- ordered_names[[i]] %>%
      left_join(., plot_label_data[[i]]) %>%
      rownames_to_column("plot_order") %>%
      select(`sample-id`, plot_order, node) %>%
      dplyr::mutate(plot_order = as.numeric(plot_order))
  }

  list(tip_data = merged_labs,
       hclust_plots = hclust_plots,
       legends = legends)

}
