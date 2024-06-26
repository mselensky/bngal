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
build_dendrograms <- function(binned.taxonomy, metadata, color.by, shape.by, trans="log10", sub.comms, cores = 1) {
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")

  NCORES=as.numeric(cores)

  # split binned.taxonomy by optional sub.comms
  community.counts = list()
  # if no subcommunity column is defined, refer to full community as "all"
  if (missing(sub.comms) | is.null(sub.comms)) {
    community.counts <- list("all" = binned.taxonomy)
    threads = "all"

  } else {
    binned.taxonomy <- binned.taxonomy %>%
      left_join(select(metadata, `sample-id`, .data[[sub.comms]]), by = "sample-id")
    community.counts <- split(binned.taxonomy, binned.taxonomy[[sub.comms]])

  }

  # create relative abundance matrices from split binned.taxonomy
  rel_abun_mats=list()
  rel_abun_dist=list()
  hclust_res=list()
  rel_abun_mats <- parallel::mclapply(X = community.counts,
                                      FUN = function(x){
                                        bngal::make_matrix(x,
                                                           "rel_abun_binned",
                                                           "sample-id") %>%
                                          replace_na(replace = 0)
                                      },
                                      mc.cores = NCORES)



  # calculate manhattan due to abundance of 0s
  if (trans=="log10") {
    rel_abun_dist <- parallel::mclapply(X = rel_abun_mats,
                                        FUN = function(x) {
                                          dist(log10(x+1), method = "manhattan")
                                        },
                                        mc.cores = NCORES)
  } else if (trans=="sqrt") {
    rel_abun_dist <- parallel::mclapply(X = rel_abun_mats,
                                        FUN = function(x) {
                                          dist(sqrt(x+1), method = "manhattan")
                                        },
                                        mc.cores = NCORES)
  } else if (trans=="log") {
    rel_abun_dist <- parallel::mclapply(X = rel_abun_mats,
                                        FUN = function(x) {
                                          dist(sqrt(x+1), method = "manhattan")
                                        },
                                        mc.cores = NCORES)
  } else if (trans=="none") {
    rel_abun_dist <- parallel::mclapply(X = rel_abun_mats,
                                        FUN = function(x) {
                                          dist(x, method = "manhattan")
                                        },
                                        mc.cores = NCORES)
  } else {
    message("
-----ERROR-| build_dendrogram():
             'trans' variable must be missing, or one of
             'log10', 'log', 'sqrt', or 'none', not '", trans, "'")
  }

  # generate dendrogram from distance matrix
  hclust_res <- parallel::mclapply(X = rel_abun_dist,
                                        FUN = function(x) {
                                          hclust(x, method = "complete")
                                        },
                                        mc.cores = NCORES)

  hclust_plots <- list()
  ordered_names <- list()

    # use ggtree to manipulate dendrogram as needed for aesthetics
  for (x in names(hclust_res)) {

    # extract dendrogram segment data
    dendrogram. = as.dendrogram(hclust_res[[x]])
    dendrogram_data <- ggdendro::dendro_data(dendrogram.)
    dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data

    # get terminal dendrogram segments
    dendrogram_ends <- dendrogram_segments %>%
      filter(yend == 0) %>% # filter for terminal dendrogram ends
      left_join(dendrogram_data$labels, by = "x") %>%
      rename(`sample-id` = label) %>%
      left_join(metadata, by = "sample-id")

    # avoid overplotting by removing ends from dendrogram segments
    dendrogram_segments <- dendrogram_segments %>%
      filter(!yend == 0)

    # dendrogram_ends also contains hc-ordered sample info, export for later plotting
    ordered_names[[x]] <- dendrogram_ends %>%
      select(`sample-id`, x) %>%
      dplyr::rename(hc.order=x)

    # plot dendrogram
    hclust_plots[[x]] <- ggplot() +
      geom_segment(data = dendrogram_segments,
                   aes(x=x, y=y, xend=xend, yend=yend)) +
      # scale_y_reverse() +
      # scale_x_reverse() +
      theme_minimal() +
      theme(legend.position = "top",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid=element_blank()) +
      ylab("Manhattan distance")

    # user-defined plot variables
    
    # color not defined, without dendrogram tip shape
    if (missing(color.by) && missing(shape.by) | is.null(color.by) && is.null(shape.by)) {
      hclust_plots[[x]] <- hclust_plots[[x]] +
        geom_segment(data = dendrogram_ends,
                     aes(x=x, y=y.x, xend=xend, yend=yend)) +
        geom_point(data = dendrogram_ends,
                   aes(x=x, y=0)) +
        theme(legend.box = "vertical")
      #guides(guide_legend(direction = "vertical"))
    # color not defined, with dendrogram tip shape
    } else if (missing(color.by) && !missing(shape.by) | is.null(color.by) && !is.null(shape.by)) {
      hclust_plots[[x]] <- hclust_plots[[x]] +
        geom_segment(data = dendrogram_ends,
                     aes(x=x, y=y.x, xend=xend, yend=yend)) +
        geom_point(data = dendrogram_ends,
                   shape = .data[[shape.by]],
                   aes(x=x, y=0)) +
        theme(legend.box = "vertical")
      #guides(guide_legend(direction = "vertical"))
    # color defined, without dendrogram tip shape
    } else if (!missing(color.by) && missing(shape.by) | !is.null(color.by) && is.null(shape.by)) {
      hclust_plots[[x]] <- hclust_plots[[x]] +
        geom_segment(data = dendrogram_ends,
                     aes(x=x, y=y.x, xend=xend, yend=yend,
                         color = .data[[color.by]])) +
        geom_point(data = dendrogram_ends,
                   aes(x=x, y=0, 
                       color = .data[[color.by]])) +
        theme(legend.box = "vertical")
      #guides(guide_legend(direction = "vertical"))
    # color defined, with dendrogram tip shape
    } else if (!missing(color.by) && !missing(shape.by) | !is.null(color.by) && !is.null(shape.by)) {
      hclust_plots[[x]] <- hclust_plots[[x]] +
        geom_segment(data = dendrogram_ends,
                     aes(x=x, y=y.x, xend=xend, yend=yend,
                         color = .data[[color.by]])) +
        geom_point(data = dendrogram_ends,
                   aes(x=x, y=0, 
                       color = .data[[color.by]],
                       shape = .data[[shape.by]])) +
        theme(legend.box = "vertical")
        #guides(guide_legend(direction = "vertical"))
    } else {
      stop("[", Sys.time(), "] build_dendrograms(): I'm confused.")
    }
  }

  list(ordered_names = ordered_names,
       hclust_plots = hclust_plots)

}

