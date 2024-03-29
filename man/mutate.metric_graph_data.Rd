% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{mutate.metric_graph_data}
\alias{mutate.metric_graph_data}
\alias{mutate}
\title{A version of \code{dplyr::mutate()} function for datasets on metric graphs}
\usage{
\method{mutate}{metric_graph_data}(.data, ...)
}
\arguments{
\item{.data}{The data list or \code{tidyr::tibble} obtained from a metric graph object.}

\item{...}{Additional parameters to be passed to \code{dplyr::mutate()}.}
}
\value{
A \code{tidyr::tibble} with the resulting selected columns.
}
\description{
Applies \code{dplyr::mutate()} function for datasets obtained from a metric graph object.
}
