% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alpha_div.R
\name{get_alpha.div}
\alias{get_alpha.div}
\title{Calculate alpha diversity at each level of taxonomic classification}
\usage{
get_alpha.div(binned.taxonomy, tax.level)
}
\arguments{
\item{binned.taxonomy}{Output from \code{\link[=bin_taxonomy]{bin_taxonomy()}}}

\item{tax.level}{Taxonomic level of network data output from either \verb{[bngal::export_network_data]}
(if running interactively) or bngal-build-networks.R (if using bngal CLI)}
}
\value{
This function can be applied directly to the output from
\code{\link[=bin_taxonomy]{bin_taxonomy()}}.
}
\description{
Calculate alpha diversity at each level of taxonomic classification
}
