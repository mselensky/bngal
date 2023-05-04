#' Summarize co-occurrences across subnetworks
#'
#' @param nodes. Node output from `bngal::extract_node_data`
#' @param edges. Network edge output data from `bngal::load_network_data`
#' @param tax.level Level of taxonomic classification
#'
#' @return
#' @export
#'
#' @examples
summarize_cooccurrence <- function(nodes., edges., tax.level) {

  # primary connections
  tmp.joined <- nodes. %>%
    select(id, taxon_, phylum:all_of(tax.level), edge_btwn_cluster, edge_btwn_cluster_color) %>%
    dplyr::rename_with(~paste0("from_", .x), .cols = everything()) %>%
    right_join(edges., by = c("from_id" = "from"))

  tmp.full.data <- nodes. %>%
    select(id, taxon_, phylum:all_of(tax.level), edge_btwn_cluster, edge_btwn_cluster_color) %>%
    dplyr::rename_with(~paste0("to_", .x), .cols = everything()) %>%
    right_join(tmp.joined, by = c("to_id" = "to")) %>%
    dplyr::mutate(same_ebc = if_else(from_edge_btwn_cluster == to_edge_btwn_cluster,
                                     TRUE, FALSE),
                  same_phy = if_else(from_phylum == to_phylum,
                                     TRUE, FALSE))

  if (tax.level %in% c("family", "genus", "asv")) {
    tmp.full.data <- tmp.full.data %>%
      dplyr::mutate(same_fam = if_else(from_family == to_family,
                                       TRUE, FALSE))
  }

  # number of nodes classified to each phylum
  n.nodes.per.phy <- nodes. %>%
    group_by(phylum) %>%
    dplyr::summarize(n_taxa_per_phy = n())

  # total number of connections per node
  k1_total_xions <- tmp.full.data %>%
    group_by(from_taxon_) %>%
    distinct(from_taxon_, to_taxon_) %>%
    dplyr::summarize(total_xions = n())

  # number of inter-ebc connections per node
  k1_inter_ebc <- tmp.full.data %>%
    group_by(from_taxon_) %>%
    filter(same_ebc == FALSE) %>%
    dplyr::summarize(inter_ebc_xions = n())

  # number of inter-phylum
  k1_inter_phy <- tmp.full.data %>%
    group_by(from_taxon_) %>%
    filter(same_phy == FALSE) %>%
    dplyr::summarize(inter_phy_xions = n())

  # number of unique phyla, ebcs, and (if applicable) families involved in each xion
  unique.phyla.k1 <- tmp.full.data %>%
    ungroup() %>%
    distinct(from_taxon_, to_phylum) %>%
    group_by(from_taxon_) %>%
    dplyr::summarize(unique_inter_phy_xions = n())
  unique.ebcs.k1 <- tmp.full.data %>%
    ungroup() %>%
    distinct(from_taxon_, to_edge_btwn_cluster) %>%
    group_by(from_taxon_) %>%
    dplyr::summarize(unique_inter_ebc_xions = n())

  node_summaries <- left_join(k1_total_xions, k1_inter_ebc, by = "from_taxon_") %>%
    left_join(k1_inter_phy, by = "from_taxon_") %>%
    left_join(unique.phyla.k1, by = "from_taxon_") %>%
    left_join(unique.ebcs.k1, by = "from_taxon_") %>%
    left_join(select(nodes., edge_btwn_cluster, taxon_, phylum), by = c("from_taxon_" = "taxon_")) %>%
    left_join(n.nodes.per.phy, by = "phylum") %>%
    replace(is.na(.), 0) %>%
    dplyr::rename(taxon_ = from_taxon_)



  # # number of inter-family (if applicable)
  # if (tax.level %in% c("family", "genus", "asv")) {
  #   k1_inter_fam <- tmp.full.data %>%
  #     group_by(from_taxon_) %>%
  #     filter(same_fam == FALSE) %>%
  #     dplyr::summarize(inter_fam_xions = n())
  #
  #   node_summaries <- node_summaries %>%
  #     left_join(k1_inter_fam, by = c("taxon_" = "from_taxon_"))
  # }
  # ebc legend colors
  ebc.leg.cols <- nodes. %>%
    distinct(edge_btwn_cluster, edge_btwn_cluster_color) %>%
    filter(edge_btwn_cluster_color != "#000000") %>%
    dplyr::mutate(edge_btwn_cluster = as.character(edge_btwn_cluster))
  ebc.leg.cols <- rbind(ebc.leg.cols, c("edge_btwn_cluster" = "other", "edge_btwn_cluster_color" = "#000000"))
  ebc.colors <- ebc.leg.cols$edge_btwn_cluster_color
  names(ebc.colors) <- as.character(ebc.leg.cols$edge_btwn_cluster)


  # secondary connections
  a <- split(edges.$to, edges.$from)
  sec_neighbors = list()
  k2.out = list()
  for (i in names(a)) {

    from_names = i
    to_names <- a[[from_names]]

    # name of taxon for which we are analyzing secondary connections
    taxon.name <- nodes.[nodes.$id == from_names,]$taxon_

    # message(from_names, " <-> ", paste0(to_names, collapse = "-"))
    # from_names: list of node IDs
    # to_names: vectors of node IDs of primary xions for each from_names element

    prim_and_sec_neighbors = a[names(a) %in% to_names]
    frmt_data <- function(i) {
      data.frame(prim_sec = i)
    }
    prim_and_sec_neighbors_l <- Reduce(rbind,
                                       lapply(prim_and_sec_neighbors, frmt_data))

    # filter for secondary connections only for a given from_names
    sec_neighbors[[from_names]] = unique(prim_and_sec_neighbors_l[prim_and_sec_neighbors_l$prim_sec != from_names,])
    sec_neighbors[[from_names]] = sec_neighbors[[from_names]][sec_neighbors[[from_names]] %in% to_names == FALSE]

    k2s <- tmp.full.data[tmp.full.data$from_id %in% sec_neighbors[[from_names]],]

    see <- k2s %>%
      ungroup() %>%
      distinct(from_id, .keep_all = TRUE)

    k2_total_xions <- length(sec_neighbors[[from_names]])

    # save ebc of node being compared (node.ebc) and sum non-matches for
    # number of k2-inter-ebc connections
    node.ebc <- nodes.[nodes.$id == from_names,]$edge_btwn_cluster
    k2_inter_ebc <- nrow(see[!see$from_edge_btwn_cluster %in% node.ebc,])

    # do the same for inter-phylum and (if applicable) inter-family k2 xions
    node.phy <- nodes.[nodes.$id == from_names,]$phylum
    k2_inter_phy <- nrow(see[!see$from_phylum %in% node.phy,])

    k2.df <- data.frame(id = from_names,
                        taxon_ = taxon.name,
                        total_xions_k2 = k2_total_xions,
                        inter_ebc_k2 = k2_inter_ebc,
                        inter_phy_k2 = k2_inter_phy)

    # if (tax.level %in% c("family", "genus", "asv")) {
    #   node.fam <- nodes.[nodes.$id == from_names,]$family
    #   k2_inter_fam <- nrow(see[!see$from_family %in% node.fam,])
    #   k2.df$inter_fam_k2 = k2_inter_fam
    # }

    k2.out[[from_names]] = k2.df

  }

  node_summaries_k2 <- Reduce(rbind, k2.out)

  # number of positive connections per node
  pos.xions <- edges. %>%
    dplyr::mutate(pos = if_else(spearman > 0, 1, 0)) %>%
    group_by(from) %>%
    dplyr::summarize(pos_xions = sum(pos))

  # join k1 and k2 summaries for output
  node_summaries_out <- merge(node_summaries, node_summaries_k2) %>%
    left_join(pos.xions, by = c("id" = "from")) %>%
    select(taxon_, phylum, id, edge_btwn_cluster, total_xions, pos_xions, everything())

  list("node_summaries" = node_summaries_out,
       "edge_details" = tmp.full.data)

}
