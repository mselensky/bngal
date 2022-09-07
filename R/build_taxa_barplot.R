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

  dendro_names <- names(dendrogram[["tip_data"]][[tax.level]])
  out.dr.taxa.bp = file.path(out.dr, "taxa-barplots", fill.by)
  if (!dir.exists(out.dr.taxa.bp)) dir.create(out.dr.taxa.bp, recursive = TRUE)

  for (x in dendro_names) {

    dat.in = plotdata[[tax.level]][[x]]
    hc.order = pull(dendrogram[["tip_data"]][[tax.level]][[x]], `sample-id`, node)
    communities = pull(dendrogram[["tip_data"]][[tax.level]][[x]], `sample-id`)

    # this will arrange filled bars by summed EBC abundance
    ebc_arranged <- dat.in %>%
      filter(`sample-id` %in% communities) %>%
      distinct(`sample-id`, edge_btwn_cluster, ebc_abun_sum) %>%
      dplyr::arrange(as.numeric(ebc_abun_sum))
    ebc_arranged["bar_order"] = seq(1:nrow(ebc_arranged))

    # rename sum_ebc_abun to rel_abun_ebc_tot
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
        dplyr::arrange(edge_btwn_cluster) %>%
        dplyr::mutate(edge_btwn_cluster = if_else(edge_btwn_cluster == 0,
                                                  "no_cluster",
                                                  as.character(edge_btwn_cluster))) %>%
        dplyr::add_row(edge_btwn_cluster = "other_cluster", edge_btwn_cluster_color = "#000000")
      ebc.colors <- ebc.color.order %>%
        pull(edge_btwn_cluster_color, edge_btwn_cluster)

      taxa_barplot.d <- taxa_barplot.d %>%
        dplyr::mutate(edge_btwn_cluster = if_else(edge_btwn_cluster == 0,
                                                  "no_cluster",
                                                  as.character(edge_btwn_cluster)))


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

      taxa_barplot.d <- dat.in %>%
        filter(`sample-id` %in% communities) %>%
        left_join(., func.groups, by = c("phylum", "class", "order", "family")) %>%
        dplyr::arrange(desc(as.numeric(edge_btwn_cluster)))

      taxa.group.colors <- taxa_barplot.d %>%
        ungroup() %>%
        distinct(.data[[fill.by]], hex.code) %>%
        filter(!is.na(.data[[fill.by]]))

      group.colorz = pull(taxa.group.colors, hex.code, .data[[fill.by]])

    } else {
      stop("\n | [", Sys.time(), "] 'fill.by' must be one of 'phylum', 'ebc', 'grouping', or 'other'")
    }

    # initialize plot
    taxa_barplot <- taxa_barplot.d %>%
      ggplot() +
      theme_bw()

    if (fill.by == "phylum") {

      phylum.colors = bngal:::phylum_colors %>%
        filter(phylum != "env_var") %>%
        pull(hex.color, phylum)

      taxa_barplot <- taxa_barplot +
        geom_bar(aes(factor(`sample-id`, levels = hc.order), rel_abun_binned,
                     fill = phylum,
                     text = paste("<br>sample: ", `sample-id`,
                                  "<br>shannon diversity: ", round(value,3),
                                  # "<br>cave: ", cave_name,
                                  # "<br>region: ", region,
                                  # "<br>water column zone: ", zone,
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
        labs(fill = "Phylum") +
        ggtitle(paste0("'", x, "' communities filled by phylum"))
    } else if (fill.by == "ebc") {
      taxa_barplot <- taxa_barplot +
        geom_bar(aes(factor(`sample-id`, levels = hc.order), ebc_abun_sum,
                     fill = as.factor(edge_btwn_cluster),
                     text = paste0("<br>sample: ", `sample-id`,
                                   "<br>shannon diversity: ", round(value,3),
                                   # "<br>zone: ", zone,
                                   # "<br>cave: ", cave_name,
                                   # "<br>region: ", region,
                                   "<br>-------------------",
                                   "<br>edge between cluster: ", edge_btwn_cluster,
                                   "<br>cluster abundance: ", ebc_abun_sum)
        ),
        stat = "identity") +
        scale_fill_manual(values = ebc.colors, na.value = "#000000") +
        labs(fill = "Edge between cluster (EBC)") +
        ggtitle(paste0("'", x, "' communities filled by edge between cluster ID"))
    } else if (fill.by == "grouping") {
      taxa_barplot <- taxa_barplot +
        geom_bar(aes(factor(`sample-id`, levels = hc.order), rel_abun_binned,
                     fill = as.factor(.data[[fill.by]]),
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
        guides(fill = guide_legend(title="Grouping")) +
        scale_fill_manual(values = group.colorz) +
        ggtitle(paste0("'", x, "' communities filled by family functional groupings"))
    } else if (fill.by == "other") {
      taxa_barplot <- taxa_barplot +
        geom_bar(aes(factor(`sample-id`, levels = hc.order), rel_abun_binned,
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
        scale_fill_manual(values = custom.colorz) +
        ggtitle(label = paste0("'", x, "' communities filled by ", other.variable))
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
      out.plot <- taxa_barplot +
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.title.x = element_blank())

      out.legend <- ggpubr::get_legend(taxa_barplot) %>%
        ggpubr::as_ggplot()

      work.dir = getwd()
      setwd(out.dr.taxa.bp)
      filename = paste0(tax.level, "-", x, "-clustered-barplot-", fill.by)
      suppressMessages(
        ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                        plot = out.plot)
      )
      if (!dir.exists("legends")) dir.create("legends")
      suppressMessages(
        ggplot2::ggsave(filename = file.path("legends", paste0(filename, "-legend.pdf")),
                        plot = out.legend)
      )
      setwd(work.dir)

    }

    #message(" | [", Sys.time(), "] Exported summary barplots to:\n |   ", file.path(out.dr.taxa.bp, tax.level))

  }

}
