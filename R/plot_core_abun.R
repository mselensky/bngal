plot_core_abun <- function(ebc.nodes.abun, ebc.nodes, metadata, dendrograms) {
  
  tax.levels <- c("phylum", "class", "order", "family", "genus", "asv")
  # build network coreness and cluster abundance plots for each defined 
  # community at each level of taxonomic classification:
  # ebc.nodes.abun = output of merged EBC node & relative abundance data
  # sub.comm = sub-community network to be summarized (defined by subset_column global var)
  core.abun.plots <- function(ebc.nodes.abun, sub.comm) {
    core.abun.grobs <- list()
    #core.plots <- list()
    # for loop to create individual plots at a given taxonomic level
    for (i in tax.levels) {
      ebc.nodes.abun[[i]] <- ebc.nodes.abun[[i]] %>%
        left_join(., select(dendrograms[["tip_data"]][[i]], `sample-id`, plot_order), by = "sample-id") %>%
        ungroup() %>%
        distinct(`sample-id`, taxon_, edge_btwn_cluster, group, .keep_all = T)
      
      mean.abun.plot <- ebc.nodes.abun[[i]] %>%
        filter(phylum != "env_var") %>%
        filter(.data[["group"]] %in% sub.comm) %>%
        ggplot() +
        geom_tile(aes(y = reorder(as.factor(edge_btwn_cluster), edge_btwn_cluster),
                      #x = factor(`sample-id`, levels = hc.order[[i]]),
                      x = plot_order,
                      fill = log10(ebc_abun_sum*100))) +
        theme_bw() +
        theme(axis.text.y = element_text(size = 6)) +
        scale_fill_viridis(scales::pseudo_log_trans(sigma = 0.001)) +
        labs(fill = "Relative\nabundance (%)") +
        xlab("Sample") + ylab("Network cluster") 
      
      core.plot <- ebc.nodes[[i]] %>%
        #filter(phylum != "env_var") %>%
        filter(.data[["group"]] %in% sub.comm) %>%
        ggplot(aes(fill = phylum)) +
        geom_point(aes(x = core, 
                       y = reorder(as.factor(edge_btwn_cluster), edge_btwn_cluster),
                       shape = node_type
        ), size = 3, position = "jitter") +
        theme_bw() +
        theme(legend.position = "none", 
              axis.text.y = element_blank(), axis.title.y = element_blank()) +
        scale_fill_manual(values = pull(phylum_colors, hex.color, phylum)) +
        scale_shape_manual(values = c("env_var" = 24, "taxon" = 21)) +
        xlab("Coreness centrality") + ylab("Network cluster") 
      
      
      core.abun.grobs[[i]] <- 
        ggpubr::ggarrange(mean.abun.plot + theme(legend.position = "none"),
                          core.plot, 
                          dendrograms[["hclust_plots"]][[i]] + theme(legend.position = "none"),
                          dendrograms[["legends"]][[i]],
                          ncol=2, nrow=2,
                          heights = c(3,1)) %>%
        ggpubr::annotate_figure(text_grob(paste0(i),
                                          face = "italic"))
      
    }
    core.abun.grobs
  }
  
  # determine unique sub-communities represented in network data:
  sub.comms <- unique(ebc.nodes.abun[["asv"]][["group"]])
  plotz <- list()
  for (sub.comm in sub.comms) {
    plotz[[sub.comm]] <- core.abun.plots(ebc.nodes.abun, sub.comm)
  }
  plotz
  
}