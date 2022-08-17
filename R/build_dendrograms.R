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
#' @param cores *(Optional)* Number of CPUs (default = 1)
#'
#' @return This function can be applied directly to the output from
#' [bngal::bin_taxonomy()].
#' A list is returned:
#'   * tip_data: underlying dendrogram data; can be used in [bngal::build_taxa.barplot()]
#'   * hclust_plots: hierarchical clustering plot\cr
#'   * legends: legend corresponding to hclust_plot\cr
#' @export
#'
#' @examples build_dendrograms(binned_tax, metadata, "sample_type", "sqrt")
build_dendrograms <- function(binned.taxonomy, metadata, color.by, trans="log10", sub.comms, cores = 1) {
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")

  NCORES=as.numeric(cores)

  # split binned.taxonomy by optional sub.comms
  test = list()
  threads = list()
  for (i in tax.levels) {

    # if no subcommunity column is defined, refer to full community as "all"
    if (missing(sub.comms) | is.null(sub.comms)) {
      test[[i]] <- list("all" = binned.taxonomy[[i]])
      threads[[i]] = "all"

    } else {

      binned.taxonomy[[i]] <- binned.taxonomy[[i]] %>%
        left_join(select(metadata, `sample-id`, .data[[sub.comms]]), by = "sample-id")
      test[[i]] <- split(binned.taxonomy[[i]], binned.taxonomy[[i]][[sub.comms]])
      threads[[i]] = names(test[[i]])

    }

  }
  # message(" | [", Sys.time(), "] Formatting complete.")
  gc()

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
  # message(" | [", Sys.time(), "] Relative abundance matrices created.")
  gc()

  ggt <- list()
  ggt_df <- list()
  hclust_plots <- list()
  for (i in tax.levels) {

    # use ggtree to manipulate dendrogram as needed for aesthetics
    for (x in names(hclust_res[[i]])) {

      ggt[[i]][[x]] <- ggtree(hclust_res[[i]][[x]])
      message(" | [", Sys.time(), "] ggtree for '", x, "' completed at ", i,"' level.")
      ggt_df[[i]][[x]] <- get.tree(ggt[[i]][[x]])$tip.label %>%
        as.data.frame() %>%
        rename(`sample-id` = ".") %>%
        left_join(metadata, by = "sample-id")

      hclust_plots[[i]][[x]] <- suppressMessages(
        ggt[[i]][[x]] %<+% ggt_df[[i]][[x]] +
          #geom_tiplab(aes(text = label, angle = 90)) +
          coord_flip() +
          scale_y_reverse()
      )

      if (missing(color.by) | is.null(color.by)) {
        hclust_plots[[i]][[x]] <- hclust_plots[[i]][[x]] +
          geom_tippoint()
      } else {
        hclust_plots[[i]][[x]] <- hclust_plots[[i]][[x]] +
          geom_tippoint(aes(color = .data[[color.by]]))
      }
    }

    # message(" | [", Sys.time(), "] Base dendrograms constructed at the '", i,"' level.")

  }

    # ggt[[i]] <- parallel::mclapply(X = threads[[i]],
    #                                FUN = function(x) {
    #                                  ggtree(hclust_res[[i]][[x]])
    #                                },
    #                                mc.cores = NCORES)
    # names(ggt[[i]]) = threads[[i]]
    #
    # gc()

    # ggt_df[[i]] <- parallel::mclapply(X = threads[[i]],
    #                                   FUN = function(x) {
    #                                     get.tree(ggt[[i]][[x]])$tip.label %>%
    #                                       as.data.frame() %>%
    #                                       rename(`sample-id` = ".") %>%
    #                                       left_join(metadata, by = "sample-id")
    #                                   },
    #                                   mc.cores = NCORES)
    # names(ggt_df[[i]]) = threads[[i]]

    # gc()
    #
    # hclust_plots[[i]] <- parallel::mclapply(X = threads[[i]],
    #                                         FUN = function(x) {
    #                                           suppressMessages(
    #                                             ggt[[i]][[x]] %<+% ggt_df[[i]][[x]] +
    #                                               #geom_tiplab(aes(text = label, angle = 90)) +
    #                                               coord_flip() +
    #                                               scale_y_reverse()
    #                                           )
    #                                         },
    #                                         mc.cores = NCORES)
    # names(hclust_plots[[i]]) = threads[[i]]
    #
    # gc()
#
#     if (missing(color.by) | is.null(color.by)) {
#
#       hclust_plots[[i]] <- parallel::mclapply(X = threads[[i]],
#                                               FUN = function(x) {
#                                                 hclust_plots[[i]][[x]] +
#                                                   geom_tippoint()
#                                               },
#                                               mc.cores = NCORES)
#       gc()
#     } else {
#       hclust_plots[[i]] <- parallel::mclapply(X = threads[[i]],
#                                               FUN = function(x) {
#                                                 hclust_plots[[i]][[x]] +
#                                                   geom_tippoint(aes(color = .data[[color.by]]))
#                                               },
#                                               mc.cores = NCORES)
#       gc()
#     }
#     names(hclust_plots[[i]]) = threads[[i]]

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
  # message(" | [", Sys.time(), "] Dendrogram legends created.")
  gc()


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

    gc()

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
    gc()

    message(" | [", Sys.time(), "] Final dendrograms constructed at the '", i, "' level.")
  }
  gc()

  list(tip_data = merged_labs,
       hclust_plots = hclust_plots,
       legends = legends)

}
