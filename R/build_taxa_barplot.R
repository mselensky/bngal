#' Build summary taxa barplots from network abundance data
#'
#' @param plotdata Output from [`bngal::ebc_compositions`]
#' @param tax.level Taxonomic level
#' @param dendrogram Output from [`bngal::build_dendrograms`]
#' @param fill.by Fill option for summary barplot. Can be one of `"phylum"` (to fill by phylum abundance), `"ebc"` (to fill by edge between cluster abundance), or `"grouping"` (to fill by family-level functional groupings). `"grouping"` is only available for family-level and below.
#' @param interactive If `TRUE`, plot exported as interactive HTML file. If `FALSE`, plot exported as PDF. Default = `TRUE`.
#' @param out.dr Output directory
#' @param metadata.cols
#' @param other.variable
#'
#' @return
#' @export
#'
#' @examples
build_taxa.barplot <- function(plotdata, tax.level, dendrogram, fill.by="phylum", interactive=TRUE, out.dr, metadata.cols, other.variable) {

  dendro_names <- names(dendrogram[["ordered_names"]][[tax.level]])
  out.dr.taxa.bp = file.path(out.dr, "taxa-barplots", fill.by)
  if (!dir.exists(out.dr.taxa.bp)) dir.create(out.dr.taxa.bp, recursive = TRUE)

  legend.out = list()
  for (x in dendro_names) {

    dat.in = plotdata[[tax.level]][[x]] %>%
      left_join(dendrogram[["ordered_names"]][[tax.level]][[x]], by = "sample-id")
    hc.order = pull(dendrogram[["ordered_names"]][[tax.level]][[x]], `sample-id`, hc.order)
    communities = pull(dendrogram[["ordered_names"]][[tax.level]][[x]], `sample-id`)

    # this will arrange filled bars by summed EBC abundance
    ebc_arranged <- dat.in %>%
      filter(`sample-id` %in% communities) %>%
      distinct(`sample-id`, edge_btwn_cluster, hc.order, ebc_abun_sum) %>%
      dplyr::arrange(as.numeric(ebc_abun_sum))
    ebc_arranged["bar_order"] = seq(1:nrow(ebc_arranged))

    # order ebc legend by overall ebc abundance
    ebc_legend_order <- ebc_arranged %>%
      group_by(edge_btwn_cluster) %>%
      dplyr::summarise(sum_ebc_abun = sum(ebc_abun_sum)/as.numeric(length(unique(ebc_arranged$`sample-id`)))) %>%
      distinct(edge_btwn_cluster, sum_ebc_abun) %>%
      arrange(desc(sum_ebc_abun))

    # extract alpha diversity data for downstream plots
    alpha.div.data <- dat.in %>%
      filter(`sample-id` %in% communities) %>%
      distinct(`sample-id`, edge_btwn_cluster, index, value) %>%
      group_by(edge_btwn_cluster) %>%
      dplyr::mutate(ebc_median_alpha = median(value)) %>%
      arrange(ebc_median_alpha)
    # alpha.div.data %>%
    #   ggplot(aes(value, edge_btwn_cluster)) +
    #   geom_boxplot() +
    #   xlab("Shannon diversity")


    if (fill.by == "phylum") {
      taxa_barplot.d <- dat.in %>%
        filter(`sample-id` %in% communities) %>%
        left_join(., select(ebc_arranged, `sample-id`, edge_btwn_cluster, bar_order), by = c("sample-id", "edge_btwn_cluster")) %>%
        group_by('sample-id') %>%
        dplyr::arrange(phylum, .by_group = TRUE)

    } else if (fill.by == "ebc") {

      plotdata2 <- dat.in %>%
        filter(`sample-id` %in% communities) %>%
        dplyr::select(`sample-id`, edge_btwn_cluster, edge_btwn_cluster_color, any_of(metadata.cols), index, value) %>%
        distinct()

      taxa_barplot.d <- ebc_arranged %>%
        left_join(., plotdata2, by = c("sample-id", "edge_btwn_cluster")) %>%
        dplyr::mutate(edge_btwn_cluster_color = if_else(edge_btwn_cluster == 0, "#f0f0f0", edge_btwn_cluster_color))

      ebc.color.order <- taxa_barplot.d %>%
        ungroup() %>% distinct(edge_btwn_cluster, edge_btwn_cluster_color) %>%
        filter(edge_btwn_cluster %in% ebc_legend_order$edge_btwn_cluster) %>%
        dplyr::mutate(edge_btwn_cluster = if_else(edge_btwn_cluster == 0,
                                                  9999,
                                                  as.numeric(edge_btwn_cluster))) %>%
        dplyr::arrange(edge_btwn_cluster) %>%
        dplyr::mutate(edge_btwn_cluster = if_else(edge_btwn_cluster == 9999,
                                                  "none",
                                                  as.character(edge_btwn_cluster))) %>%
        dplyr::mutate(edge_btwn_cluster.plot = if_else(edge_btwn_cluster_color == "#000000",
                                                       "other",
                                                       edge_btwn_cluster))
      ebc.color.order.plot <- ebc.color.order %>%
        distinct(edge_btwn_cluster.plot, edge_btwn_cluster_color) %>%
        filter(edge_btwn_cluster.plot != "other") %>%
        filter(edge_btwn_cluster.plot != "none") %>%
        add_row(edge_btwn_cluster.plot = c("other", "none"), edge_btwn_cluster_color = c("#000000", "#f0f0f0"))

      ebc.colors <- ebc.color.order.plot %>%
        pull(edge_btwn_cluster_color, edge_btwn_cluster.plot)

      taxa_barplot.d <- taxa_barplot.d %>%
        dplyr::mutate(edge_btwn_cluster = if_else(edge_btwn_cluster == 0,
                                                  "none",
                                                  as.character(edge_btwn_cluster))) %>%
        left_join(select(ebc.color.order, -edge_btwn_cluster_color), by = "edge_btwn_cluster")


    } else if (fill.by == "other") {

      if (missing(other.variable)) stop(" | [", Sys.time(), "] If fill.by == 'other', then other.variable == name of column by which to fill bars!" )

      taxa_barplot.d <- dat.in %>%
        filter(`sample-id` %in% communities) %>%
        filter(!is.na(.data[[other.variable]])) %>%
        dplyr::mutate(hex.code = if_else(is.na(hex.code), "#000000", hex.code)) %>%
        dplyr::arrange(desc(as.numeric(edge_btwn_cluster)))

      custom.colorz <- taxa_barplot.d %>%
        ungroup() %>%
        distinct(.data[[other.variable]], hex.code) %>%
        pull(hex.code, .data[[other.variable]])
    } else if (fill.by == "grouping") {

      if (!tax.level %in% c("family", "genus", "asv")) {
        stop("\n | [", Sys.time(), "] build_taxa.barplot: \n | Functional groupings can only be visualized at family level or lower.")
      }

      func.groups <- read_csv(system.file("data", "16S_families.csv", package = "bngal"), col_types = cols())
      func.groups.key <- read_csv(system.file("data", "groupings.csv", package = "bngal"), col_types = cols())

      taxa_barplot.d <- dat.in %>%
        filter(`sample-id` %in% communities) %>%
        left_join(., func.groups, by = c("phylum", "class", "order", "family")) %>%
        left_join(., func.groups.key, by = "grouping") %>%
        dplyr::arrange(desc(as.numeric(edge_btwn_cluster)))

      taxa.group.colors <- taxa_barplot.d %>%
        ungroup() %>%
        distinct(grouping_ID, hex.code) %>%
        filter(!is.na(grouping_ID) & !is.na(hex.code))

      taxa.group.colors$grouping_ID = factor(taxa.group.colors$grouping_ID,
                                             levels = func.groups.key$grouping_ID)

      group.colorz = pull(taxa.group.colors, hex.code, grouping_ID)

    } else {
      stop("\n | [", Sys.time(), "] 'fill.by' must be one of 'phylum', 'ebc', 'grouping', or 'other'")
    }

    # initialize plot
    taxa_barplot <- taxa_barplot.d %>%
      ggplot() +
      theme_minimal()

    if (fill.by == "phylum") {

      phylum.colors = bngal:::phylum_colors %>%
        filter(phylum != "env_var") %>%
        pull(hex.color, phylum)

      taxa_barplot <- taxa_barplot +
        geom_bar(aes(hc.order, rel_abun_binned * 100,
                     fill = phylum,
                     text = paste("<br>sample: ", `sample-id`,
                                  "<br>shannon diversity: ", round(value,3),
                                  "<br>-------------------",
                                  "<br>p: ", phylum,
                                  "<br>full taxon ID: ", taxon_,
                                  "<br>relative abundance: ", round(rel_abun_binned,3)*100,"%",
                                  "<br>-------------------",
                                  "<br>edge between cluster: ", edge_btwn_cluster,
                                  "<br>cluster abundance: ", round(ebc_abun_sum,3)*100,"%"
                     )),
                 stat = "identity") +
        scale_fill_manual(values = phylum.colors) +
        labs(fill = "Phylum")
    } else if (fill.by == "ebc") {
      taxa_barplot <- taxa_barplot +
        geom_bar(aes(hc.order, ebc_abun_sum*100,
                     fill = as.factor(edge_btwn_cluster.plot),
                     text = paste0("<br>sample: ", `sample-id`,
                                   "<br>shannon diversity: ", round(value,3),
                                   "<br>-------------------",
                                   "<br>edge between cluster: ", edge_btwn_cluster,
                                   "<br>cluster abundance: ", ebc_abun_sum)
        ),
        stat = "identity") +
        scale_fill_manual(values = ebc.colors, na.value = "#000000") +
        labs(fill = "EBC")
    } else if (fill.by == "grouping") {
      taxa_barplot <- taxa_barplot +
        geom_bar(aes(as.numeric(hc.order), rel_abun_binned * 100,
                     fill = as.factor(grouping_ID),
                     text = paste("<br>sample: ", `sample-id`,
                                  "<br>shannon diversity: ", round(value,3),
                                  "<br>-------------------",
                                  "<br>p: ", phylum,
                                  "<br>full taxon ID: ", taxon_,
                                  "<br>relative abundance: ", round(rel_abun_binned,3)*100,"%",
                                  "<br>-------------------",
                                  "<br>node group: ", .data[[fill.by]],
                                  "<br>edge between cluster: ", edge_btwn_cluster
                     )),
                 stat = "identity") +
        scale_fill_manual(values = group.colorz, na.value = "#000000") +
        theme(legend.title = element_blank())
        #labs(fill = "Grouping")
    } else if (fill.by == "other") {
      taxa_barplot <- taxa_barplot +
        geom_bar(aes(hc.order, rel_abun_binned * 100,
                     fill = as.factor(.data[[other.variable]]),
                     text = paste("<br>sample: ", `sample-id`,
                                  "<br>shannon diversity: ", round(value,3),
                                  "<br>-------------------",
                                  "<br>p: ", phylum,
                                  "<br>full taxon ID: ", taxon_,
                                  "<br>relative abundance: ", round(rel_abun_binned,3)*100,"%",
                                  "<br>-------------------",
                                  "<br>node group: ", .data[[other.variable]],
                                  "<br>edge between cluster: ", edge_btwn_cluster
                     )),
                 stat = "identity") +
        scale_fill_manual(values = custom.colorz)
    } else {
      stop("\n | [", Sys.time(), "] 'fill.by' values not inherited from taxa_barplot.d")
    }


    if (interactive==TRUE){
      out.plot <- plotly::ggplotly(
        taxa_barplot +
          theme(legend.position = "none",
                axis.text.x = element_blank(),
                axis.title.x = element_blank()),
        tooltip = "text"
      )


      work.dir = getwd()
      setwd(out.dr.taxa.bp)
      filename = paste0(tax.level, "-", x, "-clustered-barplot-", fill.by, ".html")
      htmlwidgets::saveWidget(out.plot,
                              file=filename,
                              selfcontained = T)
      setwd(work.dir)

    } else if (interactive==FALSE) {

        if (fill.by == "phylum") {
          out.plot <- taxa_barplot +
            theme(legend.position = "none",
                  axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  panel.grid=element_blank())

          taxa_barplot_legend <- taxa_barplot +
            guides(fill=guide_legend(ncol=6))

          out.legend <- ggpubr::get_legend(taxa_barplot_legend) %>%
            ggpubr::as_ggplot() #+
            #guides(fill=guide_legend(ncol=6))

          filename = paste0(x, "-clustered-barplot", fill.by)
          out.dr.taxa.bp.p = file.path(out.dr, "taxa-barplots")
          legends.path = file.path(out.dr.taxa.bp.p, "legends")

          legend.out[[x]][['legend']] = out.legend
          legend.out[[x]][['path']] = legends.path

          # if (!dir.exists(legends.path)) dir.create(legends.path)
          # suppressMessages(
          #   ggplot2::ggsave(filename = file.path(legends.path, paste0(filename, "-legend.pdf")),
          #                   plot = out.legend,
          #                   width = 11, height = 8.5, units = "in")
          # )

        } else {
          out.plot <- taxa_barplot +
            guides(fill=guide_legend(nrow = 3, byrow = FALSE)) +
            theme(legend.position = "bottom",
                  axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  panel.grid=element_blank())
        }

      dendro_to_plot = dendrogram[["hclust_plots"]][[tax.level]][[x]]
      n_samples = nrow(dendrogram[["ordered_names"]][[tax.level]][[x]])

      out.plot.joined <- ggpubr::ggarrange(dendro_to_plot +
                                             coord_cartesian(xlim = c(1,n_samples)),
                                           out.plot +
                                             ylab("Relative abundance (%)") +
                                             coord_cartesian(xlim = c(1,n_samples)),
                                           heights = c(1,2),
                                           align = "v",
                                           ncol = 1, hjust = -0.1,
                                           labels = paste0("'", x, "' communities clustered at the '", tax.level, "' level"))

      work.dir = getwd()
      setwd(out.dr.taxa.bp)
      filename = paste0(tax.level, "-", x, "-clustered-barplot-", fill.by)
      suppressMessages(
        ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                        plot = out.plot.joined)
      )

      setwd(work.dir)

    }

  }

  # export phylum legend color key
  if (fill.by == "phylum"){
    out.legend = legend.out[[1]][['legend']]
    legends.path = legend.out[[1]][['path']]

    if (!dir.exists(legends.path)) dir.create(legends.path)
    suppressMessages(
      ggplot2::ggsave(filename = file.path(legends.path, paste0("phylum-legend.pdf")),
                      plot = out.legend,
                      width = 11, height = 8.5, units = "in")
    )
  }


}
