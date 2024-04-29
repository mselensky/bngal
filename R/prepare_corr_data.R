#' Filter preprocessed data for correlation matrix
#'
#' This function filters input data and returns a normalized abundance matrix.
#'
#' The parameter `obs.cutoff` refers to the minimum number of observations
#' per pairwise relationship required for inclusion in the output matrix.
#'
#' @param prepared.data Output from [`bngal::prepare_network_data`] or [`bngal::split_network_data`]
#' @param transformation *Optional* Numeric transformation to apply to data (default = none). `"log10"` accepted.
#' @param obs.cutoff *Optional* Minimum number of observations required per pairwise relationship to be included in correlation matrix (default = `10`).
#' @param out.dr Required. Output directory for pairwise summary data.
#' @param num.cores See [`bngal::check_cores`]
#'
#' @return
#' @export
#'
#' @examples prepare_corr_data(prepared.data, obs.cutoff, transformation, out.dr)
prepare_corr_data <- function(prepared.data, obs.cutoff, transformation, out.dr, num.cores = NULL) {

  NCORES <- bngal::check_cores(num.cores)

  comp_corr <- function(prepared_data, transformation, obs.cutoff, out.dr) {
    if (missing(transformation) | is.null(transformation)) {
      matrix. <- prepared_data %>%
        dplyr::select(`sample-id`, taxon_, norm_vals) %>%
        pivot_wider(names_from = "taxon_",
                    values_from = "norm_vals") %>%
        column_to_rownames("sample-id") %>%
        as.matrix() %>%
        replace_na(., 0)
    } else if (transformation == "log10") {
      matrix. <- prepared_data %>%
        dplyr::select(`sample-id`, taxon_, norm_vals) %>%
        dplyr::mutate(transf = log10(norm_vals+1)) %>%
        dplyr::select(-norm_vals) %>% ungroup() %>%
        pivot_wider(names_from = "taxon_", values_from = "transf") %>%
        column_to_rownames("sample-id") %>%
        as.matrix() %>%
        replace_na(., 0)
    }

    # keep env.var '0' values (should be treated differently than taxa count '0' values)
    matrix.l <- as.data.frame(matrix.)
    matrix.l <- matrix.l %>%
      rownames_to_column("sample-id") %>%
      pivot_longer(2:ncol(.), names_to = "taxon_", values_to = "norm_vals") %>%
      left_join(., select(prepared_data, `sample-id`, taxon_, node_type), by = c("sample-id", "taxon_")) %>%
      filter(case_when(node_type == "taxon" ~ norm_vals > 0,
                       node_type == "env_var" ~ norm_vals >= 0))

    # apply pw.taxa() to determine all possible combinations
    z <- matrix.l %>%
      distinct(taxon_, node_type)

    # make taxon key to save memory on long taxa names
    taxon.key <- z %>%
      dplyr::mutate(taxon_id = seq(1:nrow(.))) %>%
      select(-node_type)

    pw.taxa <- function(z) {
      n = as.numeric(nrow(z))
      r = 2
      z = as.data.frame(z)
      pw.z = t(combn(z[["taxon_id"]], m = r))
      pw.z.df = data.frame("taxa1" = pw.z[,1], "taxa2" = pw.z[,2])
      taxa_pair = paste0(pw.z.df$taxa1, "~<>~", pw.z.df$taxa2)
      #`sample-id` = unique(z$`sample-id`)
      pw.z.df <- cbind(pw.z.df, taxa_pair)
      #pw.z.df <- cbind(pw.z.df, `sample-id`)
    }
    message(" | --------------------------------------------------------------------")
    message(" | [", Sys.time(), "] Filtering data for '", tax_level, "'-level correlation matrix...")
    pw <- pw.taxa(taxon.key)
    message(" | [", Sys.time(), "] ", length(unique(pw$taxa_pair)), " unique '", tax_level, "'-level pairwise relationships observed before observational threshold filtering.")

    # return number of pre-filtered observations per taxa pair and filter out
    # comparisons under observational threshold
    # this section requires a lot of memory for large matrices so we create a
    # unique numeric "taxon_id" for each "taxon_"

    pw_all <- matrix.l %>%
      left_join(taxon.key, by = "taxon_") %>%
      dplyr::select(-taxon_) %>%
      left_join(pw, by = c("taxon_id" = "taxa1"))

    pw_counts <- pw_all %>%
      group_by(taxa_pair) %>%
      dplyr::summarize(n_pairs = n()) %>%
      filter(n_pairs >= obs.cutoff)

    pw_filtered <- pw_all %>%
      semi_join(pw_counts, by =c("taxa_pair")) %>%
      select(-norm_vals, -node_type)

    # output matrix
    pw_filtered.mat <- pw_filtered %>%
      distinct(`sample-id`, taxon_id) %>%
      left_join(., taxon.key, by = "taxon_id") %>%
      left_join(., select(matrix.l, -node_type), by = c("sample-id", "taxon_")) %>%
      select(-taxon_id) %>%
      pivot_wider(names_from = "taxon_", values_from = "norm_vals") %>%
      replace(is.na(.), 0)

    # construct binary matrix
    pw.bin <- pw_filtered.mat %>%
      pivot_longer(2:ncol(.)) %>%
      dplyr::mutate(pres = if_else(value > 0, 1, 0)) %>%
      select(-value) %>%
      pivot_wider(names_from = "name", values_from = "pres") %>%
      column_to_rownames("sample-id")

    # return number of "post observational threshold-filtered"
    # (post_obs_filt) nodes per sample present
    post_obs_filt_nodes <- as.data.frame(rowSums(pw.bin))
    names(post_obs_filt_nodes) <- "post_obs_filt.nodes"

    # return number of "post observational threshold-filtered"
    # pairwise associations (edges) per sample present
    a <- split(post_obs_filt_nodes, rownames(post_obs_filt_nodes))
    b <- parallel::mclapply(a,
                            function(i){choose(i[,1], 2)},
                            mc.cores = NCORES)
    post_obs_filt_pw <- Reduce(rbind, b)
    post_obs_filt_pw <- data.frame(post_obs_filt_pw, row.names = names(b)) %>%
      rownames_to_column("sample-id")
    names(post_obs_filt_pw)[names(post_obs_filt_pw) == "post_obs_filt_pw"] = "post_obs_filt.pairwise"

    # return number of nodes/pairwise associations before filtering
    matrix.bin <- matrix. %>%
      as.data.frame() %>%
      rownames_to_column("sample-id") %>%
      pivot_longer(2:ncol(.)) %>%
      dplyr::mutate(pres = if_else(value > 0, 1, 0)) %>%
      select(-value) %>%
      pivot_wider(names_from = "name", values_from = "pres") %>%
      column_to_rownames("sample-id")

    pre_obs_filt_nodes <- as.data.frame(rowSums(matrix.bin, na.rm = TRUE))
    names(pre_obs_filt_nodes) <- "input.nodes"

    a <- split(pre_obs_filt_nodes, rownames(pre_obs_filt_nodes))
    b <- parallel::mclapply(a,
                            function(i){choose(i[,1], 2)},
                            mc.cores = NCORES)
    pre_obs_filt_pw <- Reduce(rbind, b)
    pre_obs_filt_pw <- data.frame(pre_obs_filt_pw, row.names = names(b)) %>%
      rownames_to_column("sample-id")
    names(pre_obs_filt_pw)[names(pre_obs_filt_pw) == "pre_obs_filt_pw"] = "input.pairwise"

    # join summary data for export
    pre.filter.data <- pre_obs_filt_nodes %>%
      rownames_to_column("sample-id") %>%
      left_join(pre_obs_filt_pw, by = "sample-id")
    post.filter.data <- post_obs_filt_nodes %>%
      rownames_to_column("sample-id") %>%
      left_join(post_obs_filt_pw, by = "sample-id")

    summ.out <- left_join(pre.filter.data, post.filter.data, by = "sample-id")
    summ.out$tax_level = tax_level


    # export per-sample pairwise summary data to csv file

    if (!is.null(nrow(prepared.data$data))) {
      write_csv(summ.out, file.path(out.dr, paste0("pairwise_summary_", tax_level, "-all.csv")))
    } else {
      write_csv(summ.out, file.path(out.dr, paste0("pairwise_summary_", tax_level, "-", i, ".csv")))
    }

    message(" | [", Sys.time(), "] Filtered data for correlation matrix:\n",
            " |   * Minimum observation threshold: ", obs.cutoff, "\n",
            " |   * # of unique pairwise '", tax_level, "'-level relationships passing threshold: ", nrow(pw_counts), "\n",
            # remove sample-id, then ncol(pw_filtered.mat)=number of unique node IDs:
            " |   * # of unique '", tax_level, "'-level nodes involved: ", ncol(pw_filtered.mat)-1, "\n",
            " |   * # of individual pairwise observations included: ", sum(pw_counts$n_pairs))


    if (!is.null(nrow(prepared.data$data))) {
      write_csv(summ.out, file.path(out.dr, paste0("pairwise_summary_", tax_level, "-all.csv")))
    } else {
      write_csv(summ.out, file.path(out.dr, paste0("pairwise_summary_", tax_level, "-", i, ".csv")))
    }


    # final output matrix
    pw_filtered.mat %>%
      column_to_rownames("sample-id") %>%
      as.matrix()

  }

  if (!is.null(nrow(prepared.data$data))) {
    dat.in <- comp_corr(prepared.data$data, transformation, obs.cutoff, out.dr)

  } else {
    dat.in = list()

    for (i in names(prepared.data)){
      prepared.data[[i]] = prepared.data[[i]]$data
    }
    dat.in = parallel::mclapply(X = prepared.data,
                                FUN = function(i){comp_corr(i, transformation, obs.cutoff, out.dr)},
                                mc.cores = NCORES)
  }

  # drop subcommunity if there is a lack of data
  for (i in names(dat.in)) {
    if (length(dat.in[[i]]) == 0) {
      dat.in[[i]] <- NULL
      message(" |   * WARNING: subcommunity '", i, "' removed from '", tax_level, "'-level network analysis due to lack of data after observational threshold filtering.")
    }
  }
  dat.in

}
