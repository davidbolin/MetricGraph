#' Covariance function for Whittle-Matérn fields
#' 
#' Computes the covariance function for a Whittle-Matérn field.
#' 
#' @param P Location (edge number and normalized location on the edge) for the
#' location to evaluate the covariance function at.
#' @param kappa Parameter kappa from the SPDE.
#' @param tau Parameter tau from the SPDE.
#' @param range Range parameter.
#' @param sigma Standard deviation parameter.
#' @param alpha Smoothness parameter (1 or 2).
#' @param graph A `metric_graph` object.
#' @details Compute the covariance function \eqn{\rho(P,s_i)}{\rho(P,s_i)} where
#' P is the provided location and \eqn{s_i}{s_i} are all locations in the mesh
#' of the graph.
#' @return Vector with the covariance function evaluate at the mesh locations.
#' @export
#'
spde_covariance <- function(P, kappa, tau, range, sigma, alpha, graph) {

  check <- check_graph(graph)

  if(!check$has.mesh) {
    stop("No mesh provided.")
  }

  if(alpha == 1) {
    if((missing(kappa) || missing(tau)) && (missing(sigma) || missing(range))){
      stop("You should either provide either kappa and tau, or sigma and range.")
    } else if(!missing(sigma) && !missing(range)){
      nu <- 1 - 0.5
      kappa <- sqrt(8 * nu) / range
      tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
                                 (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
    }

    #compute covarains of the two edges of EP[1]
    Q <- spde_precision(kappa = kappa, tau = tau,
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
    Q_line <- as.matrix(precision_exp_line(kappa = kappa, tau = tau,
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
        S <- r_1(D_matrix, kappa = kappa, tau = tau)

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
        S <- r_1(D_matrix, kappa = kappa, tau = tau)

        #covariance update
        E.ind <- c(1:2)
        Obs.ind <- -E.ind
        Bt <- solve(S[E.ind, E.ind, drop = FALSE],
                    S[E.ind, Obs.ind, drop = FALSE])
        C_P <- CV_P[graph$E[i, ]] %*% Bt
      }
      C <- c(C, C_P)
    }
  } else if (alpha == 2) {
    if((missing(kappa) || missing(tau)) && (missing(sigma) || missing(range))){
    stop("You should either provide either kappa and tau, or sigma and range.")
  } else if(!missing(sigma) && !missing(range)){
    nu <- 2 - 0.5
    kappa <- sqrt(8 * nu) / range
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
    (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
  }

  #compute covaraince of the two edges of P[1]
  Q <- spde_precision(kappa = kappa, tau = tau,
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
  Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa, tau = tau)
  Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                  kappa = kappa, tau = tau,
                                  deriv = 2)
  Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa,
                                   tau = tau, deriv = 1)
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
      Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa, tau = tau)
      Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                      kappa = kappa, tau = tau,
                                      deriv = 2)
      Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa, tau = tau,
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
      Sigma[-d.index, -d.index] <- r_2(D,kappa = kappa, tau = tau)
      Sigma[ d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                       kappa = kappa, tau = tau,
                                       deriv = 2)
      Sigma[d.index, -d.index] <- -r_2(D[3:4-2,], kappa = kappa, tau = tau,
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
  } else {
    stop("alpha should be 1 or 2.")
  }
  return(C)
}



#' The exponential covariance
#' @param D vector or matrix with distances
#' @param kappa parameter kappa from the SPDE
#' @param tau parameter tau from the SPDE
#' @noRd
r_1 <- function(D, kappa, tau) {
  return((1 / (2 * kappa * tau^2)) * exp(-kappa * abs(D)))
}


#' the Whittle--Matern covariance with alpha=1 on a circle
#' @param t locations
#' @param l_e circle perimeter
#' @param kappa parameter kappa
#' @param tau parameter tau
#' @noRd
r_1_circle <- function(t, l_e, kappa, tau) {

  c <- 1 / (2 * kappa * sinh(kappa *l_e/2) * tau^2)
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
#' @param kappa parameter kappa
#' @param tau parameter tau
#' @param t (n x 1) position on the line
#' @param t_sorted (bool)
#' @noRd
precision_exp_line <- function(kappa, tau, t,  t_sorted = FALSE) {

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
                            x = (2 * kappa * tau^2) * x_[1:count],
                            dims = c(l_t, l_t))
  if (t_sorted == FALSE)
    Q <- Matrix::t(P) %*% Q %*% P

  return(Q)
}

#' The Matern covariance with alpha = 2
#' @param D vector or matrix with distances
#' @param tau parameter tau
#' @param kappa parameter kappa
#' @param deriv (0,1,2) no derivative, first, or second order
#' @noRd
r_2 <- function(D, kappa, tau, deriv = 0){
  aD <- abs(D)
  c <- ( 1/(4 * (kappa^3) * (tau^2)))

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
