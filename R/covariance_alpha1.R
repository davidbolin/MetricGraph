
#' The exponential covariance
#' @param D vector or matrix with distances
#' @param kappa parameter kappa
#' @param sigma parameter sigma
r_1 <- function(D, kappa, sigma) {
  return((sigma^2 / (2 * kappa)) * exp(-kappa * abs(D)))
}


#' the Whittle--Matern covariance with alpha=1 on a circle
#' @param t locations
#' @param l_e circle perimeter
#' @param kappa parameter kappa
#' @param sigma parameter sigma
r_1_circle <- function(t, l_e, kappa, sigma) {

  c <- sigma^2 / (2 * kappa * sinh(kappa *l_e/2))
  r <- matrix(0, length(t), length(t))
  r_0 <- cosh(-kappa * l_e / 2)
  if (length(t) == 1) {
    r[1, 1] <- c * r_0
    return(r)
  }

  for (i in 1:(length(t) - 1)) {

    for (ii in (i + 1):length(t)){
      r[i, ii] <- cosh(kappa * (abs(t[i] - t[ii]) - l_e / 2))
    }
  }
  r <- r + t(r)
  diag(r) <- r_0
  return(c * r)
}


#' Precision matrix for exponential covariance
#' @description Computes the precision matrix of observations on an interval
#' for a Gaussian process with an exponential covariance
#' \deqn{r(h) = \sigma^2\exp(-\kappa h)/(2\kappa)}{r(h) = sigma^2*(exp(-kappa*h)/(2*kappa)}
#' @param kappa parameter kappa
#' @param sigma parameter sigma
#' @param t (n x 1) position on the line
#' @param t_sorted (bool)
#' @export
precision_exp_line <- function(kappa, sigma, t,  t_sorted = FALSE) {

  l_t <- length(t)
  i_ <- j_ <- x_ <- rep(0, 4*(l_t-1) + 2)
  count <- 0
  if (t_sorted == FALSE) {
    order_t  <- order(t)
    t        <- t[order_t]
    P <- sparseMatrix(seq_along(order_t), order_t, x = 1)
  }

  for (i in 2:l_t) {

    c1 <- exp(-kappa * (t[i] - t[i - 1]))
    c2 <- c1^2
    one_m_c2 <- 1 - c2
    c_1 <- 0.5 + c2 / one_m_c2
    c_2 <- -c1 / one_m_c2

    i_[count + 1] <- i
    j_[count + 1] <- i - 1
    x_[count + 1] <- c_2

    i_[count + 2] <- i - 1
    j_[count + 2] <- i
    x_[count + 2] <- c_2


    i_[count + 3] <- i
    j_[count + 3] <- i
    x_[count + 3] <- c_1

    i_[count + 4] <- i - 1
    j_[count + 4] <- i - 1
    x_[count + 4] <- c_1
    count <- count + 4
  }
  #boundary term correction
  i_[count + 1] <- 1
  j_[count + 1] <- 1
  x_[count + 1] <- 0.5
  i_[count + 2] <- l_t
  j_[count + 2] <- l_t
  x_[count + 2] <- 0.5
  count <- count + 2
  Q <- Matrix::sparseMatrix(i = i_[1:count],
                            j = j_[1:count],
                            x = (2 * kappa / sigma^2) * x_[1:count],
                            dims = c(l_t, l_t))
  if (t_sorted == FALSE)
    Q <- Matrix::t(P) %*% Q %*% P

  return(Q)
}



#' Compute covariance of a point to the entire graph (discretized) for
#' alpha=1 model
#' @param P (2 x 1) point with edge number and normalized location on the edge
#' @param kappa parameter kappa
#' @param sigma parameter sigma
#' @param  graph metric_graph object
#' @param  n.p number of points to compute the covariance on each edge
#' @return C (n.p*number of edges x 3) `[,1]` edge number `[,2]` distance from
#' lower edge `[,3]` covariance
#' @param scale scale the covariance by 2*kappa so that sigma corresponds to
#' the marginal standard deviation (default FALSE)
#' @export
covariance_alpha1 <- function(P, kappa, sigma, graph, n.p = 50,
                              scale = FALSE){

  gpgraph_check_graph(graph)

  #compute covarains of the two edges of EP[1]
  Q <- spde_precision(kappa = kappa, sigma = sigma,
                      alpha = 1, graph = graph)
  R <- Cholesky(Q, LDL = FALSE, perm = TRUE)
  Vs <- graph$E[P[1],]
  Z <- matrix(0, nrow = dim(graph$V)[1], ncol = 2)
  Z[Vs[1], 1] = 1
  Z[Vs[2], 2] = 1
  V  <- solve(R, solve(R,Z,system = 'P'), system='L')
  CV <- solve(R, solve(R,V,system = 'Lt'), system='Pt')

  # compute covariance between two edges and the point
  t_norm <- P[2]
  l <- graph$edge_lengths[P[1]]
  Q_line <- as.matrix(precision_exp_line(kappa = kappa, sigma = sigma,
                                         t = c(0, l * t_norm,l)))
  Q_AB <- Q_line[2, c(1,3), drop = FALSE]
  Q_AA <- Q_line[2, 2]
  B <- -Q_AB / Q_AA

  # covariance of a point to an edge
  #COV[X,Y] = cov[Xtilde+BZ,Y] = B Cov[Z,Y]
  CV_P <- CV %*% t(B)
  C <- matrix(0, nrow = n.p * graph$nE, ncol = 3)
  for (i in 1:graph$nE) {
    l <- graph$edge_lengths[i]
    t_s <- seq(0, 1, length.out = n.p)
    if (i == P[1]) {
      D_matrix <- as.matrix(dist(c(0, l, l * t_norm, l * t_s)))
      S <- r_1(D_matrix, kappa = kappa, sigma = sigma)

      #covariance update
      E.ind <- c(1:2)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind, drop = FALSE],
                  S[E.ind, Obs.ind,drop=FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      C_P <- CV_P[graph$E[i, ]] %*% Bt[, -1] + Sigma_i[1, -1]
    } else {
      D_matrix <- as.matrix(dist(c(0, l, l * t_s)))
      S <- r_1(D_matrix, kappa = kappa, sigma = sigma)

      #covariance update
      E.ind <- c(1:2)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind, drop = FALSE],
                  S[E.ind, Obs.ind, drop = FALSE])
      C_P <- CV_P[graph$E[i, ]] %*% Bt
    }
    C[ (i-1) * n.p + (1:n.p), ] <- cbind(i, l * t_s, c(C_P))
  }
  if(scale){
    return(2*kappa*C)
  } else {
    return(C)
  }
}

#' Compute covariance of a point to the mesh points of the graph for
#' alpha=1 model
#' @param P (2 x 1) point with edge number and normalized location on the edge
#' @param kappa parameter kappa
#' @param sigma parameter sigma
#' @param  graph metric_graph object
#' @param scale scale the covariance by 2*kappa so that sigma corresponds to
#' the marginal standard deviation (default FALSE)
#' @return vector with covariance values (order Vertex of Graph then mesh$PtE)
#' @export
covariance_alpha1_mesh <- function(P, kappa, sigma, graph, scale = FALSE) {

  check <- gpgraph_check_graph(graph)

  if(!check$has.mesh) {
    stop("No mesh provided.")
  }

  #compute covarains of the two edges of EP[1]
  Q <- spde_precision(kappa = kappa, sigma = sigma,
                      alpha = 1, graph = graph)
  R <- Cholesky(Q, LDL = FALSE, perm = TRUE)
  Vs <- graph$E[P[1],]
  Z <- matrix(0, nrow = dim(graph$V)[1], ncol = 2)
  Z[Vs[1], 1] = 1
  Z[Vs[2], 2] = 1
  V  <- solve(R, solve(R,Z,system = 'P'), system='L')
  CV <- solve(R, solve(R,V,system = 'Lt'), system='Pt')

  # compute covariance between two edges and the point
  t_norm <- P[2]
  l <- graph$edge_lengths[P[1]]
  Q_line <- as.matrix(precision_exp_line(kappa = kappa, sigma = sigma,
                                         t = c(0, l * t_norm,l)))
  Q_AB <- Q_line[2, c(1,3), drop = FALSE]
  Q_AA <- Q_line[2, 2]
  B <- -Q_AB / Q_AA

  # covariance of a point to an edge
  #COV[X,Y] = cov[Xtilde+BZ,Y] = B Cov[Z,Y]
  CV_P <- CV %*% t(B)
  C <- c(as.vector(CV_P))
  inds_PtE <- sort(unique(graph$mesh$PtE[,1])) #inds
  for (i in inds_PtE) {
    l <- graph$edge_lengths[i]
    t_s <- graph$mesh$PtE[graph$mesh$PtE[,1] == i,2]
    if (i == P[1]) {
      D_matrix <- as.matrix(dist(c(0, l, l * t_norm, l * t_s)))
      S <- r_1(D_matrix, kappa = kappa, sigma = sigma)

      #covariance update
      E.ind <- c(1:2)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind, drop = FALSE],
                  S[E.ind, Obs.ind,drop=FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      C_P <- CV_P[graph$E[i, ]] %*% Bt[, -1] + Sigma_i[1, -1]
    } else {
      D_matrix <- as.matrix(dist(c(0, l, l * t_s)))
      S <- r_1(D_matrix, kappa = kappa, sigma = sigma)

      #covariance update
      E.ind <- c(1:2)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind, drop = FALSE],
                  S[E.ind, Obs.ind, drop = FALSE])
      C_P <- CV_P[graph$E[i, ]] %*% Bt
    }
    C <- c(C, C_P)
  }
  if(scale){
    return(2*kappa*C)
  } else {
    return(C)
  }
}

