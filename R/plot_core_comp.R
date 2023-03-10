#' Plot edge between cluster compositions against node coreness
#'
#' @param ebc.nodes.abun Output from [`bngal::ebc_compositions()`]
#' @param tax.level Level of taxonomic classification at which to create dendrogram
#' @param metadata See [`bngal::bin_taxonomy()`]
#' @param fill.by *Optional* Metadata column by which to fill EBC composition plots
#'
#' @return
#' @export
#'
#' @examples
plot_core_comp <- function(ebc.nodes.abun, tax.level, metadata, fill.by) {

  # build network coreness and cluster abundance plots for each defined
  # community at each level of taxonomic classification:
  # ebc.nodes.abun = output of merged EBC node & relative abundance data
  # sub.comm = sub-community network to be summarized (defined by '--subcommunity' option --> "all" if NULL)

  split.comms = ebc.nodes.abun[[tax.level]]

  core.abun.grobs <- list()
  for (i in names(split.comms)) {
    ebc.comp.plot <- split.comms[[i]] %>%
      dplyr::mutate(edge_btwn_cluster = if_else(edge_btwn_cluster == 0, "no_cluster", as.character(edge_btwn_cluster))) %>%
      ungroup() %>%
      filter(phylum != "env_var") %>%
      filter(tax_level %in% tax.level) %>%
      #filter(.data[["group"]] %in% sub.comm) %>%
      distinct(`sample-id`, taxon_, edge_btwn_cluster, .keep_all = TRUE) %>%
      ggplot() +
      geom_bar(aes(y = reorder(edge_btwn_cluster, core),
                   x = ebc_abun_sum,
                   fill = .data[[fill.by]]), position="fill", stat="identity") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 6),
            legend.direction = "horizontal") +
      labs(fill = fill.by) +
      xlab("EBC composition") + ylab("Network cluster")

    ebc.comp.plot.legend <- ggpubr::get_legend(ebc.comp.plot)

    core_data <- split.comms[[i]] %>%
      ungroup() %>%
      distinct(taxon_, core, sub_comm) %>%
      filter(!is.na(core))

    core.plot <- split.comms[[i]] %>%
      ungroup() %>%
      filter(tax_level %in% tax.level) %>%
      distinct(taxon_, edge_btwn_cluster, .keep_all = TRUE) %>%
      ggplot(aes(fill = phylum)) +
      geom_vline(xintercept = quantile(core_data$core)) +
      geom_point(aes(x = core,
                     y = reorder(edge_btwn_cluster, core),
                     shape = node_type
      ), size = 3, position = "jitter") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.y = element_blank(), axis.title.y = element_blank()) +
      scale_fill_manual(values = pull(bngal:::phylum_colors_tol, phylum_color, Silva_phylum)) +
      scale_shape_manual(values = c("env_var" = 24, "taxon" = 21)) +
      xlab("Coreness centrality") + ylab("Network cluster")

    suppressWarnings(
    core.abun.grobs[[i]] <-
      ggpubr::ggarrange(ebc.comp.plot +
                          theme(legend.position = "none"),
                        core.plot,
                        ebc.comp.plot.legend,
                        ncol=2, nrow=2,
                        heights = c(3,.5)) %>%
      ggpubr::annotate_figure(text_grob(paste0(tax.level, "-level node summaries for '", i, "' communities"),
                                        face = "italic"))
    )
  }

    gridExtra::marrangeGrob(core.abun.grobs, nrow=1, ncol=1)
}
