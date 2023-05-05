#' The Matern covariance with alpha = 2
#' @param D vector or matrix with distances
#' @param sigma parameter sigma
#' @param kappa parameter kappa
#' @param deriv (0,1,2) no derivative, first, or second order
#' @noRd
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


#' Compute covariance of a point to the entire graph (discretized) for
#' alpha=2 model
#' @param P (2 x 1) point with edge number and normalized location on the edge
#' @param kappa parameter kappa
#' @param sigma parameter sigma
#' @param  graph metric_graph object
#' @param  n.p number of points to compute the covariance on each edge
#' @return C (n.p*numer of edges x 3) `[,1]` edge number `[,2]` distance from
#' lower edge `[,3]` covariance
#' @noRd
covariance_alpha2 <- function(P, kappa, sigma, graph, n.p = 50){

  check <- check_graph(graph)

  #compute covaraince of the two edges of P[1]
  Q <- spde_precision(kappa = kappa, sigma = sigma,
                      alpha = 2, graph = graph)
  if (is.null(graph$CoB))
    graph$buildC(2, FALSE)

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
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
#' @noRd
covariance_alpha2_mesh <- function(P, kappa, sigma, graph){

  check <- check_graph(graph)

  if(!check$has.mesh) {
    stop("No mesh provided.")
  }

  #compute covaraince of the two edges of P[1]
  Q <- spde_precision(kappa = kappa, sigma = sigma,
                      alpha = 2, graph = graph)
  if (is.null(graph$CoB))
    graph$buildC(2, FALSE)

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
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
  VtE <- graph$VtEfirst()
  C <- CV_P[ 4*(VtE[,1]-1)  + 1 + 2 * VtE[,2]]
  inds_PtE <- sort(unique(graph$mesh$PtE[,1])) #inds
  for (i in inds_PtE) {
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
