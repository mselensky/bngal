#' Plot correlation networks
#'
#' @param node.color.data Output from [`bngal::color_nodes()`]
#' @param selected.By Drop-down option for interactive network plot. Values of `'phylum'`, `'edge_btwn_cluster'`, or `'other'` accepted.
#' If `'other'` is chosen, the desired custom color scheme must be saved as a column labeled `hex.color` (containing hex codes) in `node.color.data[[nodes]]`
#' @param graph.layout igraph-style layout option for interactive network plot. Refer to the [igraph documentation](https://igraph.org/r/html/latest/layout_.html) for the full list of options. See also [visNetwork::visIgraphLayout()]
#' @param out.dr Figure output directory (does not have to exist in parent directory)
#' @param sign *Optional* See [`bngal::prepro_net_features()`]. For plot/output file name, if desired.
#' @param direction *Optional* See [`bngal::prepare_network_data()`]. For plot/output file name, if desired.
#' @param cutoff.val *Optional* See [`bngal::prepare_network_data()`]. For plot/output file name, if desired.
#' @param pval.cutoff *Optional* See [`bngal::prepro_net_features()`]. For plot/output file name, if desired.
#' @param other.variable *Optional* Only required if `selected.By == 'other'`. Name of nodes column by which to select node groupings.
#'
#' @return
#' @export
#'
#' @examples
plot_networks <- function (node.color.data, selected.By, graph.layout, out.dr,
                           sign, direction, cutoff.val, pval.cutoff, other.variable) {

  nodes. = node.color.data$nodes
  edges. = node.color.data$edges

  plot_vis_nets <- function(nodes., edges., subset.values, selected.By, graph.layout,
                            sign, tax_level, correlation, pval.cutoff, direction, other.variable) {

    if (selected.By == "phylum") {
      nodes.$color = nodes.$phylum_color
    } else if (selected.By == "edge_btwn_cluster") {
      nodes.$color = nodes.$edge_btwn_cluster_color
    } else if (selected.By == "other") {
      if (missing(other.variable)) {
        stop("\n If selected.By == 'other', then you must define 'other.variable' as the column by which to select nodes!")
      } else nodes.$color = nodes.$hex.color
    } else {
      stop("\n | Incorrect 'selected_By value.")
    }

    # save absolute value of correlation column as "value" (for edges width)
    edges.[["value"]] = abs(edges.[[3]])

    network.plot <- visNetwork(nodes = nodes.,
               edges = edges.,
               main = paste0(sign,
                             " co-occurrences of ",
                             tax_level,
                             "-level taxa (", correlation, " rho; p < ", pval.cutoff, ")"),
               submain = paste0(subset.values,
                                " communities (", direction, " ", cutoff.val*100,
                                "% relative abundance)")) %>%
      visIgraphLayout(layout = graph.layout, randomSeed = 99)

    if (selected.By == "other") {
      network.plot <- network.plot %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           hover = TRUE),
                   selectedBy = list(
                     multiple = TRUE,
                     variable = c(other.variable))) %>%
        visEvents(type = "on", doubleClick = "networkOpenCluster")
    } else {
      network.plot <- network.plot %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           hover = TRUE),
                   selectedBy = list(
                     multiple = TRUE,
                     variable = c(selected.By))) %>%
        visEvents(type = "on", doubleClick = "networkOpenCluster")
    }
    network.plot
  }
  export_net_plot <- function(plot.out, selected.By){

    for (i in subset.values) {
      filename = paste0(i, "-", tax_level,
                        "-network-coloredBy-", selected.By, ".html")
      htmlwidgets::saveWidget(plot.out[[i]],
                              file=filename,
                              selfcontained = T)

    }
  }

  if (!is.null(nrow(nodes.))) {
    subset.values = "all"
    plot.out <- plot_vis_nets(nodes., edges., subset.values, selected.By,
                              graph.layout,
                              sign, tax_level, correlation, pval.cutoff, direction, other.variable)
    plot.out = list(all=plot.out)

  } else if (class(nodes.) == "list") {
    subset.values = names(node.color.data$nodes)

    plot.out = list()
    for (i in subset.values) {
      plot.out[[i]] <- plot_vis_nets(nodes.[[i]],
                                     edges.[[i]],
                                     i,
                                     selected_By,
                                     graph.layout,
                                     sign, tax_level, correlation, pval.cutoff, direction)
    }

  } else {
    stop("\n | [", Sys.time(), "] Unexpected input data for bngal::plot_networks()\n",
    "  |   ->Requires output from bngal::color_nodes()")
  }

  if (!dir.exists(out.dr)) {
    dir.create(out.dr)
  }

  work.dir = getwd()
  setwd(out.dr)
  if (selected.By == "other") {
    export_net_plot(plot.out, other.variable)
  } else {
    export_net_plot(plot.out, selected.By)
  }
  # test <- mclapply(X = plot.out,
  #                  FUN = function(i){export_net_plot(i, "phylum")},
  #                  mc.cores = NCORES)
  setwd(work.dir)

}
