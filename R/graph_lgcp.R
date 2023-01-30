#' Simulation of log-Gaussian Cox process on metric graph
#' @param intercept mean value of the Gaussian process
#' @param kappa range parameter kappa
#' @param tau variance parameter
#' @param alpha smoothness parameter (1 or 2)
#' @param graph metric_graph object
#' @return List with Gaussian process sample and simulated points
#' @export
graph_lgcp <- function(intercept = 0, kappa, tau, alpha, graph) {

  if(is.null(graph$mesh)) {
    stop("No mesh provided")
  }

  if(is.null(graph$mesh$C)) {
    graph$compute_fem()
  }

  C <- Diagonal(dim(graph$mesh$C)[1], rowSums(graph$mesh$C))
  L <- kappa^2*C + graph$mesh$G

  if(alpha == 1) {
    Q <- L
  } else if(alpha == 2) {
    Q <- L%*%solve(C, L)
  } else {
    stop("not implemented yet")
  }

  R <- chol(Q)
  u <- intercept + solve(R, rnorm(dim(Q)[1]))

  lambda_max <- max(exp(u))
  domain_size <- sum(graph$edge_lengths)

  #simulate Poisson number of points
  N <- rpois(1, lambda_max*domain_size)

  #simulate locations of points from uniform distribution on edges
  p_edge <- graph$edge_lengths/domain_size
  edge_numbers <- sample(1:graph$nE,size = N, replace = TRUE, prob = p_edge)
  edge_loc <- runif(N)
  points <- cbind(edge_numbers, edge_loc)

  #Thin the sample
  lambda_loc <- exp(graph$mesh_A(points)%*%u)
  p_loc <- as.double(lambda_loc/lambda_max)
  ind_keep <- runif(N) < p_loc
  points <- points[ind_keep,]

  return(list(u = u, edge_numbers = points[,1],
              edge_loc = points[,2]))
}
