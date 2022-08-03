#' Plot edge between cluster compositions against node coreness
#'
#' @param ebc.nodes.abun
#' @param tax.level
#' @param metadata
#' @param fill.by
#'
#' @return
#' @export
#'
#' @examples
plot_core_comp <- function(ebc.nodes.abun, tax.level, metadata, fill.by) {

  # build network coreness and cluster abundance plots for each defined
  # community at each level of taxonomic classification:
  # ebc.nodes.abun = output of merged EBC node & relative abundance data
  # sub.comm = sub-community network to be summarized (defined by subset_column global var)

  ebc.nodes.abun <- ebc.nodes.abun %>%
    filter(!is.na(edge_btwn_cluster)) %>%
    filter(!is.na(core))

  ebc.comp.plot <- ebc.nodes.abun %>%
    ungroup() %>%
    filter(phylum != "env_var") %>%
    filter(tax_level %in% tax.level) %>%
    #filter(.data[["group"]] %in% sub.comm) %>%
    distinct(`sample-id`, taxon_, edge_btwn_cluster, .keep_all = TRUE) %>%
    ggplot() +
    geom_bar(aes(y = reorder(edge_btwn_cluster, core),
                 #y = reorder(as.factor(edge_btwn_cluster), edge_btwn_cluster),
                 x = ebc_abun_sum,
                 fill = .data[[fill.by]]), position="fill", stat="identity") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6),
          legend.direction = "horizontal") +
    labs(fill = fill.by) +
    xlab("EBC composition") + ylab("Network cluster")

  ebc.comp.plot.legend <- ggpubr::get_legend(ebc.comp.plot)

  core_data <- ebc.nodes.abun %>%
    ungroup() %>%
    distinct(taxon_, core) %>%
    filter(!is.na(core))

  core.plot <- ebc.nodes.abun %>%
    ungroup() %>%
    filter(tax_level %in% tax.level) %>%
    #filter(.data[["group"]] %in% sub.comm) %>%
    distinct(taxon_, edge_btwn_cluster, .keep_all = TRUE) %>%
    ggplot(aes(fill = phylum)) +
    geom_vline(xintercept = quantile(core_data$core)) +
    geom_point(aes(x = core,
                   y = reorder(edge_btwn_cluster, core),
                   #y = reorder(as.factor(edge_btwn_cluster), edge_btwn_cluster),
                   shape = node_type
    ), size = 3, position = "jitter") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_blank(), axis.title.y = element_blank()) +
    scale_fill_manual(values = pull(phylum_colors, hex.color, phylum)) +
    scale_shape_manual(values = c("env_var" = 24, "taxon" = 21)) +
    xlab("Coreness centrality") + ylab("Network cluster")

  core.abun.grobs <-
    ggpubr::ggarrange(ebc.comp.plot +
                        theme(legend.position = "none"),
                      core.plot,
                      ebc.comp.plot.legend,
                      ncol=2, nrow=2,
                      heights = c(3,.5)) %>%
    ggpubr::annotate_figure(text_grob(paste0(tax.level, "-level node summaries"),
                                      face = "italic"))

  core.abun.grobs

}
