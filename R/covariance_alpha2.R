#' The Matern covariance with alpha = 2
#' @param D vector or matrix with distances
#' @param sigma parameter sigma
#' @param kappa parameter kappa
#' @param deriv (0,1,2) no derivative, first, or second order
r_2 <- function(D, kappa, sigma, deriv = 0){
  aD <- abs(D)
  c <- ( sigma^2/(4 * kappa^3))

  R0 <-  exp( -kappa * aD)
  if (deriv == 0)
    return( c * (1 + kappa * aD) * R0)

  d1 <- -kappa^2 * c * D * R0

  if (deriv == 1)
    return( d1)
  if (deriv == 2)
    return(kappa^2 * c * ( kappa* aD - 1) * R0)
  stop("deriv must be either (0,1,2)")
}

#' plot the Matern alpha = 2
#' @param kappa parameter kappa
#' @param sigma parameter sigma
#' @param sigma_e parameter sigma_e
plot_r_2 <-function(kappa, sigma, sigma_e, t = NULL) {
  if (is.null(t)) {
    r_p <- 4.743859 / kappa
    t <- seq(0, r_p, length.out = 100)
  }
  r_  <- r_2(t,kappa, sigma)
  if (t[1] == 0)
    r_[1] = r_[1] + sigma_e
  plot(t, r_, type = "l")

}

#' Sample Gaussian process with alpha = 2 on a line given end points
#' @details Samples a Gaussian process \eqn{u(t)} with alpha = 2 on an
#' interval \eqn{(0,l_e)} conditionally on \eqn{u(0), u(l_e)}.
#' If `y` and `py` are supplied, the sampling is done conditionally on observations
#' \deqn{y_i = u(t_i) + sigma_e e_i}{y_i = u(t_i) + sigma_e e_i}
#' where \eqn{e_i} are independent standard Gaussian variables.
#' @param kappa parameter kappa
#' @param theta parameter theta
#' @param sigma_e parameter sigma_e
#' @param u_e  (2 x 1) the two end points
#' @param l_e (1 x 1) line length
#' @param  t (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py  (n x 1) observation locations
#' @param  y (n x 1) observations
#' @param  sample (bool) if true sample else return posterior mean
#' @export
sample_alpha2_line <-function(kappa, sigma, sigma_e,
                                 u_e, l_e, t=NULL, Line=NULL,
                                 nt=100,  py=NULL, y=NULL, sample=TRUE) {

  if(is.null(t)){
    t  = seq(0, 1, length.out = nt)
    t <- t * l_e
  }
  t_end <- c(0, l_e)
  t <- unique(t)
  t0 <- t
  if (is.null(y) == FALSE) {
    ind_remove = which(py %in% t_end)
    if (length(ind_remove) > 0) {
      y  <- y[-ind_remove]
      py <- py[-ind_remove]
    }

    ind_remove_t <- which(t %in% c(t_end,py))
    if(length(ind_remove_t) > 0)
      t <- t[-ind_remove_t]

    if (length(py) == 0)
      py <- NULL
  } else {
    ind_remove_t <- which(t %in% t_end)
    if (length(ind_remove_t)>0)
      t <- t[-ind_remove_t]
  }

  t <- c(t_end, t)
  if (is.null(py) == FALSE) {
    t <- c(py, t)
  }
  Sigma <- matrix(0, length(t) + 2, length(t) + 2)
  d.index <- c(1, 2)
  index_E <- 2 + length(py) + 1:2
  D <- outer (t, t, `-`)
  Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                   sigma = sigma, deriv = 0)
  Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l_e))),
                                  kappa = kappa, sigma = sigma, deriv = 2)
  Sigma[d.index, -d.index] <- -r_2(D[index_E-2,],kappa = kappa,
                                   sigma = sigma, deriv = 1)
  Sigma[-d.index,  d.index] <- t(Sigma[d.index,  -d.index])

  index_boundary <- c(d.index,index_E)
  u_e <- u_e[c(2, 4, 1, 3)]
  SinvS <- solve(Sigma[index_boundary, index_boundary],
                 Sigma[index_boundary, -index_boundary])
  Sigma_X <- Sigma[-index_boundary, -index_boundary] -
    Sigma[-index_boundary, index_boundary] %*% SinvS
  mu_X <- - t(SinvS) %*% (0-u_e)

  if(is.null(py) == FALSE){
    index_y <- 1:length(py)
    Sigma_Y <- Sigma_X[index_y, index_y, drop = FALSE]
    Matrix::diag(Sigma_Y) <-  Matrix::diag(Sigma_Y) + sigma_e^2

    SinvS <- solve(Sigma_Y, Sigma_X[index_y,,drop = FALSE])
    Sigma_X <- Sigma_X - Sigma_X[,index_y, drop = FALSE] %*% SinvS
    mu_X <- mu_X + t(SinvS) %*% (y- mu_X[index_y])
  }

  x <- rep(0, length(t))

  if(sample){
    R_X <- chol(Sigma_X)
    z <- rnorm(dim(Sigma_X)[1])
    x[-c(1:2)] <- mu_X + t(R_X) %*% z
    x[c(1:2)] <- u_e[c(3, 4)]
  }else{
    x[-c(1:2)] <- mu_X
    x[1:2] <- u_e[c(3, 4)]

  }
  x_out <- matrix(0, nrow = length(t0), 2)
  x_out[, 1] <- t0
  for(i in 1:length(t0))
  {
    ind <- which(t == t0[i])
    x_out[i, 2] <- x[ind]
  }
  return(x_out)
}

#' Generates samples of the entire graph
#' @param graph   - graphical object
#' @param theta  - (sigma_e, sigma, kappa)
#' @export
graph_posterior_mean_matern2 <- function(graph,  theta, sample=F){

  X <- c()
  V.post.mean <- posterior.mean.matern2(theta, graph)
  for(i in 1:dim(graph$EtV)[1]){
    V.i <-   V.post.mean[4*(i-1) + 1:4]

    ind <- which(graph$PtE[,1] == i)
    if(length(ind)==0){
      X.i   <- sample.line.matern2(theta,
                           V.i,
                           graph$edge_lengths[i],
                           nt = 100,
                           sample=sample)
    }else{
      X.i   <- sample.line.matern2(theta,
                           V.i,
                           graph$edge_lengths[i],
                           nt = 100,
                           y = graph$y[ind],
                           py =graph$PtE[ind,2],
                           sample=sample)

    }
    X.i[,1] <- X.i[,1]/graph$edge_lengths[i]
    X <- rbind(X, cbind(X.i,i))

  }
  return(X)

}

#' Compute covariance of a point to the entire graph (discretized) for
#' alpha=2 model
#' @param P (2 x 1) point with edge number and normalized location on the edge
#' @param kappa parameter kappa
#' @param sigma parameter sigma
#' @param  graph metric_graph object
#' @param  n.p number of points to compute the covariance on each edge
#' @return C (n.p*numer of edges x 3) [,1] edge number [,2] distance from
#' lower edge [,3] covariance
#' @export
covariance_alpha2 <- function(P, kappa, sigma, graph, n.p = 50){

  #compute covaraince of the two edges of P[1]
  Q <- spde_precision(kappa = kappa, sigma = sigma,
                      alpha = 2, graph = graph)
  if (is.null(graph$CBobj))
    graph$buildC(2, FALSE)

  n_const <- length(graph$CBobj$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CBobj$T[-ind.const, ]
  Q_mod <- Tc %*% Q %*% t(Tc)
  R <- Cholesky(Q_mod, LDL = FALSE, perm = TRUE)
  Vs <- graph$E[P[1], ]


  #COV[X,Y] = cov[Xtilde+BZ,Y] = B Cov[Z,Y]
  #the four Z are u(0),u'(0), u(T), u'(T) for each edge
  Z <- matrix(0, nrow = 4 * graph$nE, ncol = 4)
  Z[cbind(4 * (P[1] - 1) + 1:4, 1:4)] = 1
  TZ = Tc %*% Z
  V  <- Matrix::solve(R, Matrix::solve(R,TZ,system = 'P'),
                      system='L')
  TCV <- Matrix::solve(R,Matrix::solve(R,V,system = 'Lt'),
                       system='Pt')
  CZ <- t(Tc) %*% TCV

  #modfing  u(0),u'(0),u(T),u'(T)
  CZ <- CZ[,c(2, 4, 1, 3)]

  # compute covariance between
  # u(0),u'(0),u(T),u'(T) and the point u(p)
  t_norm <- P[2]
  l <- graph$edge_lengths[P[1]]
  Sigma <- matrix(0, 5, 5)
  t <- l * c(0, 1, t_norm)
  D <- outer (t, t, `-`)
  d.index <- c(1, 2)
  Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa, sigma = sigma)
  Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                  kappa = kappa, sigma = sigma,
                                  deriv = 2)
  Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa,
                                   sigma = sigma, deriv = 1)
  Sigma[-d.index, d.index] <- t(Sigma[d.index, -d.index])

  B <- Sigma[1:4, 5] %*% solve(Sigma[1:4, 1:4])
  CV_P <- CZ %*% t(B)
  C <- matrix(0, nrow = n.p * graph$nE, ncol = 3)
  for (i in 1:graph$nE) {
    l <- graph$edge_lengths[i]

    t_s <- seq(0, 1, length.out = n.p)
    if(i == P[1]){
      Sigma <- matrix(0, length(t_s) + 5, length(t_s) + 5)
      d.index <- c(1, 2)
      index_boundary <- c(d.index, 3:4)
      t <- l*c(0, 1, t_norm, t_s)
      D <- outer (t, t, `-`)
      Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa, sigma = sigma)
      Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                      kappa = kappa, sigma = sigma,
                                      deriv = 2)
      Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa, sigma = sigma,
                                       deriv = 1)
      Sigma[-d.index, d.index] <- t(Sigma[d.index,  -d.index])

      u_e <- CV_P[4 * (i-1) + (1:4)]
      u_e <- u_e[c(2, 4, 1, 3)]
      SinvS <- solve(Sigma[index_boundary, index_boundary],
                     Sigma[index_boundary, -index_boundary] )
      Sigma_X <- Sigma[-index_boundary,-index_boundary] -
        Sigma[-index_boundary, index_boundary]%*%SinvS
      C_P <- t(SinvS)[-1, ]%*%u_e + Sigma_X[1, -1]
    } else {
      Sigma <- matrix(0, length(t_s) + 4, length(t_s) + 4)
      d.index <- c(1, 2)
      index_boundary <- c(d.index, 3:4)
      t <- l*c(0, 1, t_s)
      D <- outer (t, t, `-`)
      Sigma[-d.index, -d.index] <- r_2(D,kappa = kappa, sigma = sigma)
      Sigma[ d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                       kappa = kappa, sigma = sigma,
                                       deriv = 2)
      Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa, sigma = sigma,
                                       deriv = 1)
      Sigma[-d.index, d.index] <- t(Sigma[d.index, -d.index])

      u_e <- CV_P[4 * (i-1) + (1:4)]
      u_e <- u_e[c(2, 4, 1, 3)]
      SinvS <- solve(Sigma[index_boundary, index_boundary],
                     Sigma[index_boundary, -index_boundary])
      C_P <-  t(SinvS) %*% u_e
    }
    C[ (i-1) * n.p + (1:n.p), ] <- cbind(rep(i,n.p), l*t_s, c(C_P))
  }
  return(C)
}

#' Compute covariance of a point to the mesh locations in a graph for
#' alpha=2 model
#' @param P (2 x 1) point with edge number and normalized location on the edge
#' @param kappa parameter kappa
#' @param sigma parameter sigma
#' @param  graph metric_graph object
#' @return a vector with covariance values
#' @export
covariance_alpha2_mesh <- function(P, kappa, sigma, graph){

  #compute covaraince of the two edges of P[1]
  Q <- spde_precision(kappa = kappa, sigma = sigma,
                      alpha = 2, graph = graph)
  if (is.null(graph$CBobj))
    graph$buildC(2, FALSE)

  n_const <- length(graph$CBobj$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CBobj$T[-ind.const, ]
  Q_mod <- Tc %*% Q %*% t(Tc)
  R <- Cholesky(Q_mod, LDL = FALSE, perm = TRUE)
  Vs <- graph$E[P[1], ]


  #COV[X,Y] = cov[Xtilde+BZ,Y] = B Cov[Z,Y]
  #the four Z are u(0),u'(0), u(T), u'(T) for each edge
  Z <- matrix(0, nrow = 4 * graph$nE, ncol = 4)
  Z[cbind(4 * (P[1] - 1) + 1:4, 1:4)] = 1
  TZ = Tc %*% Z
  V  <- Matrix::solve(R, Matrix::solve(R,TZ,system = 'P'),
                      system='L')
  TCV <- Matrix::solve(R,Matrix::solve(R,V,system = 'Lt'),
                       system='Pt')
  CZ <- t(Tc) %*% TCV

  #modfing  u(0),u'(0),u(T),u'(T)
  CZ <- CZ[,c(2, 4, 1, 3)]

  # compute covariance between
  # u(0),u'(0),u(T),u'(T) and the point u(p)
  t_norm <- P[2]
  l <- graph$edge_lengths[P[1]]
  Sigma <- matrix(0, 5, 5)
  t <- l * c(0, 1, t_norm)
  D <- outer (t, t, `-`)
  d.index <- c(1, 2)
  Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa, sigma = sigma)
  Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                  kappa = kappa, sigma = sigma,
                                  deriv = 2)
  Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa,
                                   sigma = sigma, deriv = 1)
  Sigma[-d.index, d.index] <- t(Sigma[d.index, -d.index])

  B <- Sigma[1:4, 5] %*% solve(Sigma[1:4, 1:4])
  CV_P <- CZ %*% t(B)
  C <- NULL
  for (i in 1:graph$nE) {
    l <- graph$edge_lengths[i]

    t_s <- graph$mesh$PtE[graph$mesh$PtE[,1] == i,2]
    if(i == P[1]){
      Sigma <- matrix(0, length(t_s) + 5, length(t_s) + 5)
      d.index <- c(1, 2)
      index_boundary <- c(d.index, 3:4)
      t <- l*c(0, 1, t_norm, t_s)
      D <- outer (t, t, `-`)
      Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa, sigma = sigma)
      Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                      kappa = kappa, sigma = sigma,
                                      deriv = 2)
      Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa, sigma = sigma,
                                       deriv = 1)
      Sigma[-d.index, d.index] <- t(Sigma[d.index,  -d.index])

      u_e <- CV_P[4 * (i-1) + (1:4)]
      u_e <- u_e[c(2, 4, 1, 3)]
      SinvS <- solve(Sigma[index_boundary, index_boundary],
                     Sigma[index_boundary, -index_boundary] )
      Sigma_X <- Sigma[-index_boundary,-index_boundary] -
        Sigma[-index_boundary, index_boundary]%*%SinvS
      C_P <- t(SinvS)[-1, ]%*%u_e + Sigma_X[1, -1]
    } else {
      Sigma <- matrix(0, length(t_s) + 4, length(t_s) + 4)
      d.index <- c(1, 2)
      index_boundary <- c(d.index, 3:4)
      t <- l*c(0, 1, t_s)
      D <- outer (t, t, `-`)
      Sigma[-d.index, -d.index] <- r_2(D,kappa = kappa, sigma = sigma)
      Sigma[ d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                       kappa = kappa, sigma = sigma,
                                       deriv = 2)
      Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa, sigma = sigma,
                                       deriv = 1)
      Sigma[-d.index, d.index] <- t(Sigma[d.index, -d.index])

      u_e <- CV_P[4 * (i-1) + (1:4)]
      u_e <- u_e[c(2, 4, 1, 3)]
      SinvS <- solve(Sigma[index_boundary, index_boundary],
                     Sigma[index_boundary, -index_boundary])
      C_P <-  t(SinvS) %*% u_e
    }
    C <- c(C,c(C_P))
  }
  return(C)
}
