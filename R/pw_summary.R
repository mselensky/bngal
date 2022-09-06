#' Summarize and export pairwise data for each sample after quality control
#'
#' @param corr.data Output from [`bngal::prepare_corr_data()`]
#' @param preprocessed.features Required. Output from [`bngal::prepare_net_features()`]
#' @param tax.level Taxonomic level at which to summarize pairwise data.
#' @param out.dr Required. Output directory for pairwise summary data. Should be the same path as defined in [`bngal::prepare_corr_data()`]
#' @param cores *(Optional)* Number of CPUs. Default = 1
#'
#' @return
#' @export
#'
#' @examples
pw_summary <- function(corr.data, preprocessed.features, tax.level, out.dr, cores=1){

  # look in "pairwise-summaries" folder for input data
  pw.out.dr = file.path(out.dr, "pairwise-summaries")

  # define number of CPUs
  NCORES = cores

  # define taxonomic level of classification
  tax_level = tax.level

  input.data.class = c("tbl_df", "tbl", "data.frame")
  if (class(preprocessed.features$nodes) %in%  input.data.class) {
    preprocessed.features$edges = list("all" = preprocessed.features$edges)
    preprocessed.features$nodes = list("all" = preprocessed.features$nodes)
    corr.data = list("all" = corr.data)
  }

  message(" | [", Sys.time(), "] Number of unique nodes passing final quality filtering: ")
  for (x in names(preprocessed.features$nodes)) {
    # return number of nodes that passed bngal::prepare_net_features()
    qc.nodes <- preprocessed.features$nodes[[x]] %>%
      distinct(label) %>%
      pull()
    message(" |   ** '", x, "': ", length(qc.nodes))


    # pairwise input data
    pw.in <- read_csv(file.path(pw.out.dr, paste0("pairwise_summary_", tax_level, "-", x, ".csv")),
                      col_types = cols())

    keep.samples <- pw.in %>%
      filter(post_obs_filt.pairwise > 0) %>%
      pull(`sample-id`)

    corr.data.sub <- corr.data[[x]] %>%
      as.data.frame() %>%
      rownames_to_column("sample-id") %>%
      filter(`sample-id` %in% keep.samples) %>%
      pivot_longer(2:ncol(.), names_to = "taxon_", values_to = "norm_vals") %>%
      left_join(select(preprocessed.features$nodes[[x]], label, node_type),
                by = c("taxon_" = "label"))
    # if is.na(corr.data.sub$node_type) == TRUE, then node was removed during
    # bngal::prepare_net_features()
    corr.data.sub <- corr.data.sub %>%
      filter(!is.na(node_type)) %>%
      filter(case_when(node_type == "taxon" ~ norm_vals > 0,
                       node_type == "env_var" ~ norm_vals >= 0))

    corr.data.sub.split <- split(corr.data.sub, corr.data.sub$`sample-id`)

    b <- parallel::mclapply(corr.data.sub.split,
                            function(i){choose(nrow(i), 2)},
                            mc.cores = NCORES)

    # determine final number of quality-controlled nodes + pairwise relationships
    # per sample and save to input csv as new columns "QC.nodes" and "QC.pairwise"
    qc_nodes <- corr.data.sub %>%
      group_by(`sample-id`) %>%
      dplyr::summarize(QC.nodes = n())
    qc_pw <- Reduce(rbind, b)
    qc_pw <- data.frame(qc_pw, row.names = names(b)) %>%
      rownames_to_column("sample-id")
    names(qc_pw)[names(qc_pw) == "qc_pw"] = "QC.pairwise"

    pw.out <- left_join(qc_nodes, qc_pw, by = "sample-id")

    dat.out <- pw.in %>%
      left_join(pw.out, by = "sample-id") %>%
      select(`sample-id`, tax_level, everything())

    write_csv(dat.out, file.path(pw.out.dr, paste0("pairwise_summary_", tax_level, "-", x, ".csv")))

  }


}
