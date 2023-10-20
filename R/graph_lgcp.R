#' Simulation of log-Gaussian Cox processes driven by Whittle-Matérn
#' fields on metric graphs
#' @param n Number of samples.
#' @param intercept Mean value of the Gaussian process.
#' @param sigma Parameter for marginal standard deviations.
#' @param range Parameter for practical correlation range.
#' @param alpha Smoothness parameter (1 or 2).
#' @param graph A `metric_graph` object.
#' @return List with Gaussian process sample and simulated points.
#' @export
graph_lgcp <- function(n = 1, intercept = 0, sigma, range, alpha, graph) {

  if(is.null(graph$mesh)) {
    stop("No mesh provided")
  }

  if(n < 1 || n%%1 != 0){
    stop("n must be an integer")
  }

  if(is.null(graph$mesh$C)) {
    graph$compute_fem()
  }
  nu <- alpha - 1/2
  kappa <- sqrt(8*nu)/range
  tau <- sqrt(gamma(nu)/(gamma(alpha)*sqrt(4*pi)*kappa^(2*nu)*sigma^2))
  C <- Diagonal(dim(graph$mesh$C)[1], rowSums(graph$mesh$C))
  L <- kappa^2*C + graph$mesh$G

  if(alpha == 1) {
    Q <- tau^2*L
  } else if(alpha == 2) {
    Q <- tau^2*L%*%solve(C, L)
  } else {
    stop("not implemented yet")
  }
  R <- chol(Q)
  result <- list()
  for(i in 1:n){
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
    lambda_loc <- exp(graph$fem_basis(points)%*%u)
    p_loc <- as.double(lambda_loc/lambda_max)
    ind_keep <- runif(N) < p_loc
    edge_numbers <- edge_loc <- NULL
    if (length(ind_keep) > 0) {
      edge_numbers <- points[ind_keep,1]
      edge_loc <- points[ind_keep,2]
    }
    result[[i]] <- list(u = u, edge_numbers = edge_numbers, edge_loc = edge_loc)
  }

  if(n == 1){
    return(result[[1]])
  }

  return(result)
}
