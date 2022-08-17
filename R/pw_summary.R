#' Summarize and export pairwise data for each sample after quality control
#'
#' @param binned.taxonomy Output from [`bngal::bin_taxonomy()`]
#' @param preprocessed.features Required. Output from [`bngal::prepro_net_features()`]
#' @param out.dr Required. Output directory for pairwise summary data. Should be the same path as defined in [`bngal::prepare_corr_data()`]
#' @param cores *(Optional)* Number of CPUs. Default = 1
#'
#' @return
#' @export
#'
#' @examples
pw_summary <- function(binned.taxonomy, preprocessed.features, out.dr, cores=1){

  pw.out.dr = file.path(out.dr, "pairwise-summaries")
  NCORES = cores

  input.data.class = c("tbl_df", "tbl", "data.frame")
  if (class(preprocessed.features$nodes) %in%  input.data.class) {
    preprocessed.features$edges = list("all" = preprocessed.features$edges)
    preprocessed.features$nodes = list("all" = preprocessed.features$nodes)
  }

  tax_level = names(binned.taxonomy[ncol(binned.taxonomy)-1])

  for (x in names(preprocessed.features$nodes)) {

    pw.in <- read_csv(file.path(pw.out.dr, paste0("pairwise_summary_", tax_level, "-", x, ".csv")),
                      col_types = cols())

    # subset binned taxonomy based on subcommunity
    binned.taxonomy.sub <- binned.taxonomy %>%
      filter(`sample-id` %in% pw.in[["sample-id"]])

    # filter for nodes that passed all QC steps
    filtered.node.names <- binned.taxonomy.sub %>%
      dplyr::semi_join(preprocessed.features$nodes[[x]], by = c("taxon_" = "label")) %>%
      left_join(select(preprocessed.features$nodes[[x]],
                       label, id), by = c("taxon_" = "label"))

    # summarize number of unique taxa in each sample
    unique.taxa <- filtered.node.names %>%
      distinct(`sample-id`, taxon_) %>%
      group_by(`sample-id`) %>%
      dplyr::summarize(n_unique_taxa = n())

    # get list of each unique taxon for each sample
    qc.pairwise <- filtered.node.names %>%
      left_join(., preprocessed.features$edges[[x]], by = c("id" = "from")) %>%
      distinct(`sample-id`, taxon_)

    qc.pairwise.split <- split(qc.pairwise, qc.pairwise$`sample-id`)

    # calculate number of final QC'd pairwise relationships included in network
    calc.qc.pw <- function(y) {
      n = nrow(y)
      r = 2
      final.qc.pw = factorial(n)/(factorial(r)*factorial(n-r))
    }
    final.pw <- parallel::mclapply(X = qc.pairwise.split,
                                   FUN = calc.qc.pw,
                                   mc.cores = NCORES)

    final.pw = as_tibble(final.pw) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("sample-id") %>%
      rename("finalQC_pairwise"="V1")

    pw.out <- pw.in %>%
      left_join(final.pw, by = "sample-id") %>%
      left_join(unique.taxa, by = "sample-id") %>%
      select(`sample-id`, tax_level, everything())

    write_csv(pw.out, file.path(pw.out.dr, paste0("pairwise_summary_", tax_level, "-", x, ".csv")))

  }


}
