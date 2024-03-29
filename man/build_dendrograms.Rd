% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_dendrograms.R
\name{build_dendrograms}
\alias{build_dendrograms}
\title{Build dendrograms}
\usage{
build_dendrograms(
  binned.taxonomy,
  metadata,
  color.by,
  trans = "log10",
  sub.comms,
  cores = 1
)
}
\arguments{
\item{binned.taxonomy}{Output from \code{\link[=bin_taxonomy]{bin_taxonomy()}}}

\item{metadata}{Sample metadata}

\item{color.by}{Metadata column by which to color}

\item{trans}{Transformation to apply to relative abundance data (default = none). Can be one of \code{"log10"}, \code{"log"}, \code{"sqrt"}, or \code{"none"}}

\item{sub.comms}{\emph{(Optional)} Metadata column by which to split dendrograms into subcommunities}

\item{cores}{\emph{(Optional)} Number of CPUs (default = 1)}
}
\value{
This function can be applied directly to the output from
\code{\link[=bin_taxonomy]{bin_taxonomy()}}.
A list is returned:
\itemize{
\item tip_data: underlying dendrogram data; can be used in \code{\link[=build_taxa.barplot]{build_taxa.barplot()}}
\item hclust_plots: hierarchical clustering plot\cr
\item legends: legend corresponding to hclust_plot\cr
}
}
\description{
This function performs hierarchical clustering on microbial relative
abundance data binned at a specified level of taxonomic classification.
}
\examples{
build_dendrograms(binned_tax, metadata, "sample_type", "sqrt")
}
