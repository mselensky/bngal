#' Preprocess network edge data
#'
#' This function preprocesses edge and node data, combining both outputs into a list object to be used in
#' downstream functions, such as [`bngal::get_igraph()`].
#'
#' @param edges. Output from [`bngal::generate_edges()`]
#' @param node.ids Output from [`bngal::get_node_ids()`]
#' @param p.val.cutoff Any pairwise comparison with a p-value greater than this value will be dropped
#' @param correlation Metric for pairwise comparisions (see [`Hmisc::rcorr()`])
#' @param correlation.cutoff Values from `0` to `1` accepted. Only relationships with an absolute value of this or greater are kept.
#' @param sign Direction(s) of pairwise relationship to include in network. Values from `"positive"`, `"negative"`, or `"all"` accepted.
#' @param num.cores See [`bngal::check_cores()`]
#'
#' @return
#' @export
#'
#' @examples
prepare_net_features <- function(edges., node.ids, p.val.cutoff, correlation, correlation.cutoff, sign, num.cores = NULL) {

 NCORES <- bngal::check_cores(num.cores)

  format_edges <- function(edges.) {

    if (sign == "positive") {
      edges_filt <- edges. %>%
        filter(.data[[correlation]] > abs(correlation.cutoff))
    } else if (sign == "negative") {
      edges_filt <- edges. %>%
        filter(.data[[correlation]] < -1*abs(correlation.cutoff))
    } else if (sign == "all") {
      edges_filt <- edges. %>%
        filter(.data[[correlation]] > abs(correlation.cutoff) | .data[[correlation]] < -1*abs(correlation.cutoff))
    } else {
      stop("\n | [", Sys.time(), "] prepare_net_features():\n",
      " |   ->sign for ",correlation," must be 'positive', 'negative', or 'all', not ", "'",sign,"'")
    }

    edges_filt <- edges_filt %>%
      dplyr::group_by(from) %>%
      dplyr::filter(p_value <= p.val.cutoff) %>%
      dplyr::mutate(degree = sum(n()),
                    color = if_else(.data[[correlation]] > 0,
                                    true = "#2f0f70", false = "#ed362f")) %>% # color edges by direction of relationship
      dplyr::select(from, to, all_of(correlation), p_value, degree, color) %>%
      filter(!is.na(from) & !is.na(to))
  }

  format_nodes <- function (node.ids, edges_filt) {
    nodes_filt <- node.ids %>%
      dplyr::left_join(., select(edges_filt, from, degree), by = c("id" = "from")) %>%
      dplyr::distinct(id, label, degree, .keep_all = TRUE) %>%
      dplyr::rename(value = degree) %>% # degree column must be renamed "value" to plot this column by size in visNetwork
      dplyr::filter(value != "NA")
  }


  input.data.class = c("tbl_df", "tbl", "data.frame")
  if (any(class(edges.) %in% input.data.class) &
      any(class(node.ids) %in% input.data.class)) {

    edges_f=format_edges(edges.)
    nodes_f=format_nodes(node.ids, edges_f)

    dat.out = list(edges = edges_f,
                   nodes = nodes_f)
  } else if (class(edges.) %in% "list" &
             class(node.ids) %in% "list") {
    # parallel approach is at least twice as fast but I can't figure out how to transform
    # the output data lol
    # t_init=Sys.time()
    # dat.out <- parallel::mcmapply(edges., node.ids,
    #                               FUN = function(i, j){preprocess_data(i, j)},
    #                               mc.preschedule = TRUE,
    #                               mc.cores = NCORES
    #                               )
    # t_par=Sys.time()
    edges_f <- parallel::mclapply(X = edges.,
                                  FUN = function(i){format_edges(i)},
                                  mc.cores = NCORES)
    nodes_f = list()
    for (i in names(edges_f)) {
      nodes_f[[i]] = format_nodes(node.ids[[i]], edges_f[[i]])
    }

    # standardize output: data$edges and .$nodes
    dat.out.seq = list(edges = edges_f,
                       nodes = nodes_f)

    for (i in names(dat.out.seq$nodes)) {
      if (nrow(dat.out.seq$nodes[[i]]) == 0) {
        message(" |   * WARNING: subcommunity '", i, "' removed from '", tax_level, "'-level network analysis due to lack of significant pairwise associations.")
        dat.out.seq$nodes[[i]] <- NULL
        dat.out.seq$edges[[i]] <- NULL
      }
    }
    dat.out.seq

  } else {
    stop("\n | [", Sys.time(), "] Unexpected input data.\n",
    " |   ->Do your edges match your nodes?")
  }

}
