#'  Gaussian processes on metric graphs
#'
#' `GPGraph` is used for specifying Gaussian random fields on metric graphs.
#' The main models are the Whittle-Matérn fields, which are specified through the
#' fractional elliptic SPDE
#' \deqn{(\kappa^2 - \Delta)^{\alpha/2} (\tau u(s)) = W,}
#' \eqn{\kappa,\tau>0} and \eqn{\alpha>1/2} are parameters and \eqn{W} is Gaussian white noise.
#' The package also implements models based on graph Laplacians and isotropic covariance functions,
#' and contains various utility functions for specifying graphs, computing likelihoods, performing
#' prediction, simulating processes, and visualizing results and graphs.
#'
#' At the heart of the package is the `R6` class `[metric_graph()]`. This is used for specifying
#' metric graphs, and contains various utility functions which are needed for specifying Gaussian
#' processes on such spaces.
#'
#' Basic statistical operations such as likelihood evaluations (see
#' `[likelihood.graph.covariance], [likelihood.graph]`) and
#' predictions (see ...) are also implemented...
#'
#'
#' For a more detailed introduction to the package, see the GPGraph Vignettes.
#'
#' @docType package
#' @name GPGraph
#' @import Matrix
#' @export Cholesky t
#' @importFrom igraph graph distances
#' @importFrom rSPDE matern.covariance gg_df
#' @importFrom methods is slot
#' @export gg_df
#' @importFrom stats approx dist dnorm pnorm rnorm
#' @importFrom sp Line Lines SpatialPoints SpatialPointsDataFrame CRS Polygon Polygons SpatialPolygons coordinates
#' @importFrom ggplot2 ggplot geom_path aes geom_point coord_fixed labs scale_colour_gradientn guide_legend
#' @importFrom CB c_basis2
#' @importFrom rgeos gInterpolate gProject
#' @importFrom igraph E E<-
#' @importFrom plotly plot_ly add_trace
#' @importFrom magrittr %>%
#' @importFrom viridis scale_color_viridis viridis
NULL
