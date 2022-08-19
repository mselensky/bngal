#' Plot correlation networks
#'
#' @param node.color.data Output from [`bngal::color_nodes()`]
#' @param filled.by Value by which to color network nodes. Also defines the drop-down option for interactive network plot. Values of `'phylum'`, `'edge_btwn_cluster'`, or `'other'` accepted.
#' If `'other'` is chosen, the desired custom color scheme must be saved as a column labeled `hex.color` (containing hex codes) in `node.color.data[[nodes]]`
#' @param graph.layout igraph-style layout option for interactive network plot. Refer to the [igraph documentation](https://igraph.org/r/html/latest/layout_.html) for the full list of options. See also [visNetwork::visIgraphLayout()]
#' @param out.dr Figure output directory (does not have to exist in parent directory)
#' @param sign *Optional* See [`bngal::prepro_net_features()`]. For plot/output file name, if desired.
#' @param direction *Optional* See [`bngal::prepare_network_data()`]. For plot/output file name, if desired.
#' @param cutoff.val *Optional* See [`bngal::prepare_network_data()`]. For plot/output file name, if desired.
#' @param pval.cutoff *Optional* See [`bngal::prepro_net_features()`]. For plot/output file name, if desired.
#' @param other.variable *Optional* Only required if `filled.by == 'other'`. Name of nodes column by which to select node groupings.
#'
#' @return
#' @export
#'
#' @examples
plot_networks <- function (node.color.data, filled.by, graph.layout, out.dr,
                           sign, direction, cutoff.val, pval.cutoff, other.variable) {

  if (filled.by == "other") {
    plot.type = other.variable
  } else {
    plot.type = filled.by
  }

  pdf.path = file.path(out.dr, "network-plots", "pdfs", plot.type, tax_level)
  html.path = file.path(out.dr, "network-plots", "html", plot.type, tax_level)

  if (!dir.exists(pdf.path)) dir.create(pdf.path, recursive = TRUE)
  if (!dir.exists(html.path)) dir.create(html.path, recursive = TRUE)

  plot_nets <- function(nodes., edges., subset.values, filled.by, graph.layout,
                        sign, tax_level, correlation, pval.cutoff, direction, other.variable, pdf_path = pdf.path) {

    if (filled.by == "phylum") {
      nodes.[[subset.values]]$color = nodes.[[subset.values]]$phylum_color
    } else if (filled.by == "edge_btwn_cluster") {
      nodes.[[subset.values]]$color = nodes.[[subset.values]]$edge_btwn_cluster_color
    } else if (filled.by == "other") {
      if (missing(other.variable)) {
        stop("\n | [", Sys.time(), "] If filled.by='other', then you must define 'other.variable' as the column by which to select your nodes!")
      } else nodes.[[subset.values]]$color = nodes.[[subset.values]]$hex.color
    } else {
      stop("\n | [", Sys.time(), "] Incorrect 'filled.by' value.")
    }

    nodes.[[subset.values]] <- nodes.[[subset.values]] %>%
      rename(shape_visN = shape) %>%
      dplyr::mutate(shape = if_else(shape_visN == "dot", "circle", shape_visN),
                    edge_btwn_cluster = as.numeric(edge_btwn_cluster))

    # save absolute value of correlation column as "value" (for edges width)
    edges.[[subset.values]][["value"]] = edges.[[subset.values]][[3]]

    # define plot titles
    if (direction == "greaterThan") {
      direction. = ">"
    } else if (direction == "lessThan") {
      direction. = "<"
    }
    plot_title = paste0(sign, " co-occurrences of ", tax_level, "-level taxa (", correlation, direction., correlation_cutoff, "; p<", pval.cutoff, ")")
    plot_subtitle = paste0(subset.values, " communities (", direction., cutoff.val*100, "% relative abundance)")

    # create igraph object from nodes. and edges.
    igraph.obj <- igraph::graph_from_data_frame(d = edges.[[subset.values]],
                                                vertices = nodes.[[subset.values]],
                                                directed = T)

    # scale edge widths
    E(igraph.obj)$arrow.mode <- 0
    edge_widths = abs(E(igraph.obj)$value)
    edge_widths_orig = E(igraph.obj)$value # save this for visNetwork edge labels
    norm_widths = 5*(0.2+(edge_widths-min(edge_widths))/(max(edge_widths)-min(edge_widths)))
    E(igraph.obj)$value <- norm_widths
    n = 3
    size_vec = seq_len(n)
    edge_widths_cut = cut(edge_widths, n, dig.lab = min(edge_widths))


    # scale node diameters
    core_scale = abs(V(igraph.obj)$core)
    norm_core = 5*(0.2+(core_scale-min(core_scale))/(max(core_scale)-min(core_scale)))
    V(igraph.obj)$norm_core <- norm_core

    node_size_cut = cut(core_scale, n, dig.lab = min(core_scale))
    vertex_shapes = V(igraph.obj)$shape
    edge_color  = E(igraph.obj)$color

    set.seed(1)
    layout = layout_nicely(igraph.obj)

    if (filled.by == "other") {
      pdf(file = file.path(pdf_path, paste0(subset.values, "-", tax_level, "-by-grouping.pdf")))
    } else {
      pdf(file = file.path(pdf_path, paste0(subset.values, "-", tax_level, "-by-", filled.by, ".pdf")))
    }
    plot(igraph.obj,
         vertex.label = NA,
         vertex.size = norm_core,
         vertex.shapes = vertex_shapes,
         edge.color = edge_color,
         edge.width	= norm_widths,
         main = plot_title,
         sub = plot_subtitle,
         layout = layout,
         margin = 0
    )
    legend("topleft",
           title = expression(paste("| ", rho, " |")),
           legend = levels(edge_widths_cut),
           lty = 1,
           bty = "n",
           lwd = c(min(norm_widths),
                   max(norm_widths) - 2*((max(norm_widths)-min(norm_widths))/n),
                   max(norm_widths)),
           inset = -0.05,
           y.intersp = 0.7
    )
    legend("bottomleft",
           legend = c(min(core_scale),
                      (max(core_scale)-min(core_scale))/2,
                      max(core_scale)),
           title = "Coreness",
           pch = 21,
           bty = "n",
           pt.cex = c(min(norm_core),
                      (max(norm_core)-min(norm_core))/2,
                      max(norm_core))/2.5,
           inset = -0.05,
           y.intersp = 0.7)
    if (filled.by == "edge_btwn_cluster") {

      ebc.leg.cols <- nodes.[[subset.values]] %>%
        distinct(edge_btwn_cluster, color)

      legend("topright",
             #inset=c(-0.2,0),
             bty = "n",
             title = "EBC",
             legend = ebc.leg.cols$edge_btwn_cluster,
             col = ebc.leg.cols$color,
             pch = 19,
             pt.cex = (max(norm_core)-min(norm_core))/2,
             y.intersp = 0.9,
             x.intersp = 0.6,
             #ncol = 2
      )
    }
    suppressMessages(
      dev.off()
    )


    ## plot interactive version with visNetwork
    # visNetwork inverses y axis layout for some reason...
    layout2 <- as.matrix(data.frame(layout[,1], -1*layout[,2]))
    V(igraph.obj)$shape=V(igraph.obj)$shape_visN
    E(igraph.obj)$title = paste0(correlation, ": ", round(edge_widths_orig, 2), " ||| ",
                                 "p-value: ", round(E(igraph.obj)$p_value, 3))
    V(igraph.obj)$title = paste0("Phylum: ", V(igraph.obj)$phylum, " ||| ",
                                 "EBC: ", V(igraph.obj)$edge_btwn_cluster, " ||| ",
                                 "Coreness: ", V(igraph.obj)$core)

    visnet.nodes <- igraph::as_data_frame(igraph.obj, what = "vertices")
    visnet.nodes$id = visnet.nodes$name
    visnet.edges <- igraph::as_data_frame(igraph.obj, what = "edges")

    network.plot <- visNetwork(nodes = visnet.nodes,
                               edges = visnet.edges,
                               main = plot_title,
                               submain = plot_subtitle) %>%
      visIgraphLayout(layout = "layout.norm",
                      layoutMatrix = layout2) %>%
      visEvents(type = "on", doubleClick = "networkOpenCluster")

    if (filled.by == "other") {
      network.plot <- network.plot %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           hover = TRUE),
                   selectedBy = list(
                     multiple = TRUE,
                     variable = other.variable))
    } else {
      network.plot <- network.plot %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           hover = TRUE),
                   selectedBy = list(
                     multiple = TRUE,
                     variable = filled.by))
    }
    network.plot
  }

  export_visnet_plot <- function(plot.out, filled.by){

    for (i in subset.values) {
      filename = paste0(i, "-", tax_level,
                        "-network-coloredBy-", filled.by)

      plot.out[[i]] <- plot.out[[i]] %>%
        visExport(type = "pdf",
                  name = filename)

      htmlwidgets::saveWidget(plot.out[[i]],
                              file=paste0(filename, ".html"),
                              selfcontained = T)

    }

  }

  # import functional grouping data
  func.group.data <- system.file("data", "16S_families.csv", package = "bngal")
  func.groups <- read_csv(func.group.data, col_types = cols())

  # rename for clarity
  nodes. = node.color.data$nodes
  edges. = node.color.data$edges


  if (!is.null(nrow(nodes.))) {

    if (tax_level %in% c("family", "genus", "asv") & filled.by == "other") {
      nodes. <- nodes. %>%
        left_join(., func.groups, by = c("phylum", "class", "order", "family")) %>%
        dplyr::rename(hex.color=hex.code)
    }

    subset.values = "all"
    plot.out <- plot_nets(nodes., edges., subset.values, filled.by,
                          graph.layout,
                          sign, tax_level, correlation, pval.cutoff, direction, other.variable)
    plot.out = list(all=plot.out)

  } else if (class(nodes.) == "list") {
    subset.values = names(nodes.)

    plot.out = list()
    for (i in subset.values) {

      if (tax_level %in% c("family", "genus", "asv") & filled.by == "other") {
        nodes.[[i]] <- nodes.[[i]] %>%
          left_join(., func.groups, by = c("phylum", "class", "order", "family")) %>%
          dplyr::rename(hex.color=hex.code)
      }

      plot.out[[i]] <- plot_nets(nodes.,
                                 edges.,
                                 i,
                                 filled.by,
                                 graph.layout,
                                 sign, tax_level, correlation, pval.cutoff, direction, other.variable)
    }

  } else {
    stop("\n | [", Sys.time(), "] Unexpected input data for bngal::plot_networks()\n",
         "  |   ->Requires output from bngal::color_nodes()")
  }

  work.dir = getwd()
  setwd(html.path)
  if (filled.by == "other") {
    export_visnet_plot(plot.out, other.variable)
  } else {
    export_visnet_plot(plot.out, filled.by)
  }
  # test <- mclapply(X = plot.out,
  #                  FUN = function(i){export_net_plot(i, "phylum")},
  #                  mc.cores = NCORES)
  setwd(work.dir)

}
