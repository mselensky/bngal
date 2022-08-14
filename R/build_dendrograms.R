#' Build dendrograms
#'
#' This function performs hierarchical clustering on microbial relative
#' abundance data binned at a specified level of taxonomic classification.
#'
#' @param binned.taxonomy Output from [`bngal::bin_taxonomy()`]
#' @param metadata Sample metadata
#' @param color.by Metadata column by which to color
#' @param trans Transformation to apply to relative abundance data (default = none). Can be one of `"log10"`, `"log"`, `"sqrt"`, or `"none"`
#' @param sub.comms *(Optional)* Metadata column by which to split dendrograms into subcommunities
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
build_dendrograms <- function(binned.taxonomy, metadata, color.by, trans="log10", sub.comms) {
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")

  # split binned.taxonomy by optional sub.comms
  test = list()
  threads = list()
  for (i in tax.levels) {

    # if no subcommunity column is defined, refer to full community as "all"
    if (missing(sub.comms) | is.null(sub.comms)) {
      test[[i]] <- list("all" = binned.taxonomy[[i]])
      threads = "all"

    } else {

      binned.taxonomy[[i]] <- binned.taxonomy[[i]] %>%
        left_join(select(metadata, `sample-id`, .data[[sub.comms]]), by = "sample-id")
      test[[i]] <- split(binned.taxonomy[[i]], binned.taxonomy[[i]][[sub.comms]])
      threads[[i]] = names(test[[i]])

    }

  }
  message(" | [", Sys.time(), "] Formatting complete.")

  # create relative abundance matrices from split binned.taxonomy
  rel_abun_mats=list()
  rel_abun_dist=list()
  hclust_res=list()
  for (i in tax.levels) {
    rel_abun_mats[[i]] <- parallel::mclapply(X = threads[[i]],
                                             FUN = function(x){
                                               bngal::make_matrix(test[[i]][[x]],
                                                                  "rel_abun_binned",
                                                                  "sample-id") %>%
                                                 replace_na(replace = 0)
                                             },
                                             mc.cores = NCORES)
    names(rel_abun_mats[[i]]) = threads[[i]]

    # calculate manhattan due to abundance of 0s
    if (trans=="log10") {
      rel_abun_dist[[i]] <- parallel::mclapply(X = threads[[i]],
                                               FUN = function(x) {
                                                 dist(log10(rel_abun_mats[[i]][[x]]+1), method = "manhattan")
                                               },
                                               mc.cores = NCORES)
    } else if (trans=="sqrt") {
      rel_abun_dist[[i]] <- parallel::mclapply(X = threads[[i]],
                                               FUN = function(x) {
                                                 dist(sqrt(rel_abun_mats[[i]][[x]]+1), method = "manhattan")
                                               },
                                               mc.cores = NCORES)
    } else if (trans=="log") {
      rel_abun_dist[[i]] <- parallel::mclapply(X = threads[[i]],
                                               FUN = function(x) {
                                                 dist(log(rel_abun_mats[[i]][[x]]+1), method = "manhattan")
                                               },
                                               mc.cores = NCORES)
    } else if (trans=="none") {
      rel_abun_dist[[i]] <- parallel::mclapply(X = threads[[i]],
                                               FUN = function(x) {
                                                 dist((rel_abun_mats[[i]][[x]]+1), method = "manhattan")
                                               },
                                               mc.cores = NCORES)
    } else {
      message("
        -----ERROR-| build_dendrogram():
                     'trans' variable must be missing or one of
                     'log10', 'log', 'sqrt', 'none'.
                ")}
    names(rel_abun_dist[[i]]) = threads[[i]]

    # generate dendrogram from distance matrix
    hclust_res[[i]] <- parallel::mclapply(X = threads[[i]],
                                          FUN = function(x) {
                                            hclust(rel_abun_dist[[i]][[x]], method = "complete")
                                          },
                                          mc.cores = NCORES)
    names(hclust_res[[i]]) = threads[[i]]

  }
  message(" | [", Sys.time(), "] Relative abundance matrices created.")

  ggt <- list()
  ggt_df <- list()
  hclust_plots <- list()
  for (i in tax.levels) {
    # use ggtree to manipulate dendrogram as needed for aesthetics
    # ggt[[i]] <- ggtree(hclust_res[[i]])

    ggt[[i]] <- parallel::mclapply(X = threads[[i]],
                                   FUN = function(x) {
                                     ggtree(hclust_res[[i]][[x]])
                                   },
                                   mc.cores = NCORES)
    names(ggt[[i]]) = threads[[i]]

    ggt_df[[i]] <- parallel::mclapply(X = threads[[i]],
                                      FUN = function(x) {
                                        get.tree(ggt[[i]][[x]])$tip.label %>%
                                          as.data.frame() %>%
                                          rename(`sample-id` = ".") %>%
                                          left_join(metadata, by = "sample-id")
                                      },
                                      mc.cores = NCORES)
    names(ggt_df[[i]]) = threads[[i]]

    hclust_plots[[i]] <- parallel::mclapply(X = threads[[i]],
                                            FUN = function(x) {
                                              suppressMessages(
                                                ggt[[i]][[x]] %<+% ggt_df[[i]][[x]] +
                                                  #geom_tiplab(aes(text = label, angle = 90)) +
                                                  coord_flip() +
                                                  scale_y_reverse()
                                              )
                                            },
                                            mc.cores = NCORES)
    names(hclust_plots[[i]]) = threads[[i]]

    if (missing(color.by) | is.null(color.by)) {

      hclust_plots[[i]] <- parallel::mclapply(X = threads[[i]],
                                              FUN = function(x) {
                                                hclust_plots[[i]][[x]] +
                                                  geom_tippoint()
                                              },
                                              mc.cores = NCORES)
    } else {
      hclust_plots[[i]] <- parallel::mclapply(X = threads[[i]],
                                              FUN = function(x) {
                                                hclust_plots[[i]][[x]] +
                                                  geom_tippoint(aes(color = .data[[color.by]]))
                                              },
                                              mc.cores = NCORES)
    }
    names(hclust_plots[[i]]) = threads[[i]]

  }
  message(" | [", Sys.time(), "] Raw dendrogram plots created.")

  legends <- list()
  ordered_names <- list()
  for (i in tax.levels) {
    for (x in names(hclust_plots[[i]])) {

      # export hclust_plot legend as separate object
      legends[[i]][[x]] <- hclust_plots[[i]][[x]] %>%
        ggpubr::get_legend()

      # ordered names for dendrogram
      ordered_names[[i]][[x]] <- tibble(`sample-id` = get_taxa_name(hclust_plots[[i]][[x]]))

    }
  }
  message(" | [", Sys.time(), "] Dendrogram legends created.")


  plot_label_data <- list()
  merged_labs <- list()
  for (i in tax.levels) {
    # merge `ordered_names` with ggt[[i]][[x]]$data to join `y` values,
    # which is what we need for ordering figure

    plot_label_data[[i]] <- parallel::mclapply(X = threads[[i]],
                                               FUN = function(x) {
                                                 ggt[[i]][[x]]$data %>%
                                                   rename(`sample-id` = label) %>%
                                                   filter(!is.na(`sample-id`))
                                               },
                                               mc.cores = NCORES)
    names(plot_label_data[[i]]) = threads[[i]]

    merged_labs[[i]] <- parallel::mclapply(X = threads[[i]],
                                           FUN = function(x) {
                                             ordered_names[[i]][[x]] %>%
                                               left_join(., plot_label_data[[i]][[x]], by = "sample-id") %>%
                                               rownames_to_column("plot_order") %>%
                                               select(`sample-id`, plot_order, node) %>%
                                               dplyr::mutate(plot_order = as.numeric(plot_order))
                                           },
                                           mc.cores = NCORES)
    names(merged_labs[[i]]) = threads[[i]]

    message(" | [", Sys.time(), "] Final dendrograms constructed at the '", i, "' level.")
  }

  list(tip_data = merged_labs,
       hclust_plots = hclust_plots,
       legends = legends)

}
