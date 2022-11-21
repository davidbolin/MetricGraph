#' Log-likelihood calculation for alpha=1 model
#'
#' @param theta (sigma_e, sigma, kappa)
#' @param graph metric_graph object
#' @param version if 1, the likelihood is computed by integrating out
#' the vertex locations, if 2, no integration is done
#' @return The log-likelihood
#' @export
likelihood_graph_spde <- function(theta, graph, alpha = 1, version = 1) {

  check <- gpgraph_check_graph(graph)

  if (alpha == 1 && version == 1) {
    return(likelihood_alpha1(theta, graph))
  } else if (alpha == 1 && version == 2) {
    return(likelihood_alpha1_v2(theta, graph))
  } else if (alpha == 2) {
    return(likelihood_alpha2(theta, graph))
  } else {
    stop("Only alpha = 1 and alpha = 2 implemented")
  }
}

#' Computes the log likelihood function fo theta for the graph object
#' @param theta parameters (sigma_e, sigma, kappa)
#' @param graph  metric_graph object
likelihood_alpha2 <- function(theta, graph) {
  if(is.null(graph$C)){
    graph$buildC(2)
  }
  sigma_e <- theta[1]
  sigma <- theta[2]
  kappa <- theta[3]
  #build Q

  n_const <- length(graph$CBobj$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CBobj$T[-ind.const, ]
  Q <- spde_precision(kappa = kappa, sigma = sigma,
                      alpha = 2, graph = graph)
  R <- Matrix::Cholesky(Tc%*%Q%*%t(Tc),
                        LDL = FALSE, perm = TRUE)
  loglik <- Matrix::determinant(R)$modulus[1]

  #build BSIGMAB
  Qpmu      <- rep(0, 4*nrow(graph$E))
  obs.edges <- unique(graph$PtE[,1])

  i_ <- j_ <- x_ <- rep(0, 16 * length(obs.edges))
  count <- 0
  for (e in obs.edges) {
    obs.id <- graph$PtE[,1] == e
    y_i <- graph$y[obs.id]
    l <- graph$edge_lengths[e]
    t <- c(0, l, l*graph$PtE[obs.id,2])

    D <- outer (t, t, `-`)
    S <- matrix(0, length(t) + 2, length(t) + 2)

    d.index <- c(1,2)
    S[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                 sigma = sigma, deriv = 0)
    S[ d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                 kappa = kappa, sigma = sigma,
                                 deriv = 2)
    S[d.index, -d.index] <- -r_2(D[1:2,], kappa = kappa,
                                 sigma = sigma, deriv = 1)
    S[-d.index, d.index] <- t(S[d.index, -d.index])

    #covariance update see Art p.17
    E.ind <- c(1:4)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind, drop = FALSE])
    Sigma_i <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

    R <- chol(Sigma_i, pivot = TRUE)
    if(attr(R, "rank") < dim(R)[1])
      return(-Inf)

    Sigma_iB <- t(Bt)
    Sigma_iB[attr(R,"pivot"),] <- forwardsolve(R,
                                        backsolve(R, t(Bt[, attr(R,"pivot")]),
                                        transpose = TRUE), upper.tri = TRUE)
    BtSinvB <- Bt %*% Sigma_iB

    E <- graph$E[e, ]
    if (E[1] == E[2]) {
      cat("Warning: circle not implemented\n")
    }
      BtSinvB <- BtSinvB[c(3,1,4,2), c(3,1,4,2)]
      Qpmu[4 * (e - 1) + 1:4] <- Qpmu[4 * (e - 1) + 1:4] +
        (t(Sigma_iB) %*% y_i)[c(3, 1, 4, 2)]

      #lower edge precision u
      i_[count + 1] <- 4 * (e - 1) + 1
      j_[count + 1] <- 4 * (e - 1) + 1
      x_[count + 1] <- BtSinvB[1, 1]

      #lower edge  u'
      i_[count + 2] <- 4 * (e - 1) + 2
      j_[count + 2] <- 4 * (e - 1) + 2
      x_[count + 2] <- BtSinvB[2, 2]

      #upper edge  u
      i_[count + 3] <- 4 * (e - 1) + 3
      j_[count + 3] <- 4 * (e - 1) + 3
      x_[count + 3] <- BtSinvB[3, 3]

      #upper edge  u'
      i_[count + 4] <- 4 * (e - 1) + 4
      j_[count + 4] <- 4 * (e - 1) + 4
      x_[count + 4] <- BtSinvB[4, 4]

      #lower edge  u, u'
      i_[count + 5] <- 4 * (e - 1) + 1
      j_[count + 5] <- 4 * (e - 1) + 2
      x_[count + 5] <- BtSinvB[1, 2]
      i_[count + 6] <- 4 * (e - 1) + 2
      j_[count + 6] <- 4 * (e - 1) + 1
      x_[count + 6] <- BtSinvB[1, 2]

      #upper edge  u, u'
      i_[count + 7] <- 4 * (e - 1) + 3
      j_[count + 7] <- 4 * (e - 1) + 4
      x_[count + 7] <- BtSinvB[3, 4]
      i_[count + 8] <- 4 * (e - 1) + 4
      j_[count + 8] <- 4 * (e - 1) + 3
      x_[count + 8] <- BtSinvB[3, 4]

      #lower edge  u, upper edge  u,
      i_[count + 9]  <- 4 * (e - 1) + 1
      j_[count + 9]  <- 4 * (e - 1) + 3
      x_[count + 9]  <- BtSinvB[1, 3]
      i_[count + 10] <- 4 * (e - 1) + 3
      j_[count + 10] <- 4 * (e - 1) + 1
      x_[count + 10] <- BtSinvB[1, 3]

      #lower edge  u, upper edge  u',
      i_[count + 11] <- 4 * (e - 1) + 1
      j_[count + 11] <- 4 * (e - 1) + 4
      x_[count + 11] <- BtSinvB[1, 4]
      i_[count + 12] <- 4 * (e - 1) + 4
      j_[count + 12] <- 4 * (e - 1) + 1
      x_[count + 12] <- BtSinvB[1, 4]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4 * (e - 1) + 2
      j_[count + 13] <- 4 * (e - 1) + 3
      x_[count + 13] <- BtSinvB[2, 3]
      i_[count + 14] <- 4 * (e - 1) + 3
      j_[count + 14] <- 4 * (e - 1) + 2
      x_[count + 14] <- BtSinvB[2, 3]

      #lower edge  u', upper edge  u',
      i_[count + 15] <- 4 * (e - 1) + 2
      j_[count + 15] <- 4 * (e - 1) + 4
      x_[count + 15] <- BtSinvB[2, 4]
      i_[count + 16] <- 4 * (e - 1) + 4
      j_[count + 16] <- 4 * (e - 1) + 2
      x_[count + 16] <- BtSinvB[2, 4]

      count <- count + 16

    loglik <- loglik - 0.5  * t(y_i)%*%solve(Sigma_i, y_i)
    loglik <- loglik - sum(log(diag(R)))
  }
  i_ <- i_[1:count]
  j_ <- j_[1:count]
  x_ <- x_[1:count]
  BtSB <- Matrix::sparseMatrix(i = i_,
                               j = j_,
                               x = x_,
                               dims = dim(Q))
  Qp <- Q + BtSB
  Qp <- Tc %*% Qp %*% t(Tc)
  R <- Matrix::Cholesky(forceSymmetric(Qp), LDL = FALSE, perm = TRUE)
  loglik <- loglik - Matrix::determinant(R)$modulus[1]

  v <- c(as.matrix(Matrix::solve(R, Matrix::solve(R, Tc%*%Qpmu, system = 'P'),
                                 system='L')))

  n.o <- length(graph$y)
  loglik <- loglik + 0.5  * t(v) %*% v  - 0.5 * n.o * log(2 * pi)

  return(loglik[1])
}



#' Log-likelihood calculation for alpha=1 model
#'
#' @param theta (sigma_e, sigma, kappa)
#' @param graph metric_graph object
#' @details This function computes the likelihood without integrating out
#' the vertices.
#' @return The log-likelihood
likelihood_alpha1_v2 <- function(theta, graph) {
  sigma_e <- theta[1]
  #build Q
  Q <- spde_precision(kappa = theta[3], sigma = theta[2],
                      alpha = 1, graph = graph)
  A <- Diagonal(graph$nV, rep(1, graph$nV))[graph$PtV, ]
  Q.p <- Q  + t(A) %*% A / sigma_e^2
  mu.p <- solve(Q.p, as.vector(t(A) %*% graph$y/sigma_e^2))

  R <- chol(Q)
  R.p <- chol(Q.p)

  n.o <- length(graph$y)
  l <- sum(log(diag(R))) - sum(log(diag(R.p))) - n.o * log(sigma_e)
  v <- graph$y  - A %*% mu.p
  l <- l - 0.5*(t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
    0.5 * n.o * log(2*pi)
  return(as.double(l))
}

#' Log-likelihood calculation for alpha=1 model
#' @param theta (sigma_e, sigma, kappa)
#' @param graph metric_graph object
likelihood_alpha1 <- function(theta, graph) {
  sigma_e <- theta[1]
  #build Q

  Q.list <- spde_precision(kappa = theta[3], sigma = theta[2], alpha = 1,
                           graph = graph, build = FALSE)

  Qp <- Matrix::sparseMatrix(i = Q.list$i,
                             j = Q.list$j,
                             x = Q.list$x,
                             dims = Q.list$dims)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  loglik <- Matrix::determinant(R)$modulus[1]

  #build BSIGMAB
  Qpmu <- rep(0, nrow(graph$V))
  obs.edges <- unique(graph$PtE[, 1])

  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))

  count <- 0
  for (e in obs.edges) {
    obs.id <- graph$PtE[,1] == e
    y_i <- graph$y[obs.id]
    l <- graph$edge_lengths[e]
    D_matrix <- as.matrix(dist(c(0, l, l*graph$PtE[obs.id, 2])))
    S <- r_1(D_matrix, kappa = theta[3], sigma = theta[2])

    #covariance update see Art p.17
    E.ind <- c(1:2)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE])
    Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
      S[Obs.ind, E.ind, drop = FALSE] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    R <- chol(Sigma_i)
    Sigma_iB <- solve(Sigma_i, t(Bt))
    BtSinvB <- Bt %*% Sigma_iB

    E <- graph$E[e, ]
    if (E[1] == E[2]) {
      Qpmu[E[1]] <- Qpmu[E[1]] + sum(t(Sigma_iB) %*% y_i)
      i_[count + 1] <- E[1]
      j_[count + 1] <- E[1]
      x_[count + 1] <- sum(Bt %*% Sigma_iB)
    } else {
      Qpmu[E] <- Qpmu[E] + t(Sigma_iB) %*% y_i
      i_[count + (1:4)] <- c(E[1], E[1], E[2], E[2])
      j_[count + (1:4)] <- c(E[1], E[2], E[1], E[2])
      x_[count + (1:4)] <- c(BtSinvB[1, 1], BtSinvB[1, 2],
                             BtSinvB[1, 2], BtSinvB[2, 2])
      count <- count + 4
    }

    loglik <- loglik - 0.5 * t(y_i) %*% solve(Sigma_i, y_i)
    loglik <- loglik - sum(log(diag(R)))
  }

  i_ <- c(Q.list$i, i_[1:count])
  j_ <- c(Q.list$j, j_[1:count])
  x_ <- c(Q.list$x, x_[1:count])
  Qp <- Matrix::sparseMatrix(i = i_,
                             j = j_,
                             x = x_,
                             dims = Q.list$dims)

  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  loglik <- loglik - Matrix::determinant(R)$modulus[1]
  n.o <- length(graph$y)
  v <- c(as.matrix(Matrix::solve(R, Matrix::solve(R, Qpmu,system = "P"),
                                 system = "L")))

  loglik <- loglik + 0.5  * t(v) %*% v - 0.5 * n.o * log(2*pi)
  return(loglik[1])
}


#' Likelihood evaluation not using sparsity
#'
#' @param theta parameters (sigma_e, sigma, kappa)
#' @param graph metric_graph object
#' @param model Model for Gaussian process, supported options are alpha1 (SPDE with alpha=1),
#' alpha2 (SPDE with alpha=2), GL (graph Laplacian model with alpha=1), GL2 (graph Laplacian
#' model with alpha=2) and isoExp (model with isotropic exponential covariance)
#' @return The log-likelihood
#' @export
likelihood_graph_covariance <- function(theta, graph, model = "alpha1") {

  check <- gpgraph_check_graph(graph)

  sigma_e <- theta[1]
  sigma <- theta[2]
  kappa <- theta[3]

  #build covariance matrix
  if (model == "alpha1") {

    Q <- spde_precision(kappa = kappa, sigma = sigma,
                        alpha = 1, graph = graph)
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]

  } else if (model == "alpha2") {

    n.c <- 1:length(graph$CBobj$S)
    Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 2,
                        graph = graph, BC = 1)
    Qtilde <- (graph$CBobj$T) %*% Q %*% t(graph$CBobj$T)
    Qtilde <- Qtilde[-n.c,-n.c]
    Sigma.overdetermined  = t(graph$CBobj$T[-n.c,]) %*% solve(Qtilde) %*%
      (graph$CBobj$T[-n.c,])
    index.obs <- 4 * (graph$PtE[,1] - 1) + 1.0 * (abs(graph$PtE[, 2]) < 1e-14) +
      3.0 * (abs(graph$PtE[, 2]) > 1e-14)
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])

  } else if (model == "GL"){

    Q <- (kappa^2 * Diagonal(graph$nV, 1) + graph$Laplacian) / sigma^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]

  } else if (model == "GL2"){

    Q <- (kappa^2 * Diagonal(graph$nV, 1) + graph$Laplacian) / sigma^2
    Q <- Q %*% Q
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]

  } else if (model == "isoExp"){

    if(is.null(graph$res.dist)){
      stop("You must first compute the resistance metric for the observations")
    }
    Sigma <- as.matrix(sigma^2 * exp(-kappa * graph$res.dist))
  } else {
    stop("wrong model choice.")
  }

  diag(Sigma) <- diag(Sigma) + sigma_e^2

  R <- chol(Sigma)
  return(as.double(-sum(log(diag(R))) - 0.5*t(graph$y)%*%solve(Sigma,graph$y) -
                     length(graph$y)*log(2*pi)/2))
}


#' Log-likelihood calculation for graph Laplacian model with alpha=1
#'
#' @param theta (sigma_e, sigma, kappa)
#' @param graph metric graph object
#' @param alpha integer values supported
#' @return The log-likelihood
#' @export
likelihood_graph_laplacian <- function(theta, graph, alpha) {

  check <- gpgraph_check_graph(graph)

  if(alpha%%1 != 0){
    stop("only integer values of alpha supported")
  }

  sigma_e <- theta[1]
  sigma <- theta[2]
  kappa <- theta[3]

  K <- kappa^2*Diagonal(graph$nV, 1) + graph$Laplacian
  Q <- K / sigma^2
  if (alpha>1) {
    for (k in 2:alpha) {
      Q <- Q %*% K
    }
  }

  Q.p <- Q  + t(graph$A) %*% graph$A/sigma_e^2
  mu.p <- solve(Q.p,as.vector(t(graph$A) %*% graph$y / sigma_e^2))

  R <- chol(Q)
  R.p <- chol(Q.p)

  n.o <- length(graph$y)
  l <- sum(log(diag(R))) - sum(log(diag(R.p))) - n.o*log(sigma_e)
  v <- graph$y - graph$A%*%mu.p
  l <- l - 0.5*(t(mu.p)%*%Q%*%mu.p + t(v)%*%v/sigma_e^2) - 0.5 * n.o*log(2*pi)
  return(as.double(l))
}
