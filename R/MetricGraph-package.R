#'  Gaussian processes on metric graphs
#'
#' 'MetricGraph' is used for creation and manipulation of metric graphs, such as street or river networks.
#' It also has several functions thatfacilitates operations and visualizations of data on metric graphs,
#' and the creation of a large class of random fields and stochastic partial differential equations on
#' such spaces.
#' The main models are the Whittle-Matérn fields, which are specified through the
#' fractional elliptic SPDE
#' \deqn{(\kappa^2 - \Delta)^{\alpha/2} (\tau u(s)) = W,}
#' \eqn{\kappa,\tau>0} and \eqn{\alpha>1/2} are parameters and \eqn{W} is Gaussian white noise.
#' It contains exact implementations of the above model for \eqn{\alpha=1} and \eqn{\alpha=2},
#' and contains approximate implementations, via the finite element method, for any \eqn{\alpha > 0.5}.
#' It also implements models based on graph Laplacians and isotropic covariance functions.
#' Several utility functions for specifying graphs, computing likelihoods, performing
#' prediction, simulating processes, and visualizing results on metric graphs are provided.
#' In particular, linear mixed effects models including random field components can be fitted to
#' data based on computationally efficient sparse matrix representations. Interfaces to the R
#' packages 'INLA' and 'inlabru' are also provided, which facilitate working with Bayesian statistical
#' models on metric graphs.
#'
#' At the heart of the package is the `R6` class `[metric_graph()]`. This is used for specifying
#' metric graphs, and contains various utility functions which are needed for specifying Gaussian
#' processes on such spaces.
#'
#' Linear mixed effects models are provided (see
#' `[graph_lme]`) and perform
#' predictions (see `[predict.graph_lme]`). The package also has interfaces for
#' 'INLA' (see `[graph_spde]`), and it this interface also works with 'inlabru'.
#'
#' For a more detailed introduction to the package, see the 'MetricGraph' Vignettes.
#'
"_PACKAGE"
#' @import Matrix
#' @export Cholesky t
#' @importFrom igraph graph distances is_tree
#' @importFrom rSPDE matern.covariance gg_df
#' @importFrom methods is slot new
#' @export gg_df
#' @importFrom stats predict approx nobs deviance lm logLik na.omit dist sd dnorm pnorm rnorm var as.formula delete.response model.matrix optim terms simulate rpois runif profile qchisq 
#' @importFrom ggplot2 ggplot geom_path aes geom_point coord_fixed labs scale_colour_gradientn guide_legend scale_colour_discrete
#' @importFrom igraph E E<-
#' @importFrom Rcpp evalCpp
#' @importFrom R6 R6Class
#' @importFrom graphics layout
#' @importFrom lifecycle deprecated
#' @importFrom dplyr mutate filter select summarise
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom zoo na.approx
#' @importFrom ggnewscale new_scale_color
#' @export %>%
#' @export summarise
#' @export mutate
#' @export filter
#' @export select
#' @export drop_na
#' @importFrom broom glance augment
# @importFrom broom tidy
#' @export glance
#' @export augment
# @export tidy
#' @useDynLib MetricGraph, .registration = TRUE
NULL

