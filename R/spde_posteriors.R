
#' Computes the posterior expectation for SPDE models
#' @param theta parameters (sigma_e, sigma, kappa)
#' @param graph  metric_graph object
#' @param alpha alpha parameter (1 or 2)
#' @param obs if TRUE, then the posterior mean is calculated at the observation locations, otherwise
#' it is computed at the vertices in the graph
#' @param leave_edge_out if obs=TRUE, then leave_edge_out = TRUE means that the posterior mean
#' is computed for each observation based on all observations which are not on that edge. If
#' obs=FALSE, then leave_edge_out = e means that the posterior mean is computed based on all 
#' observations except those on that particular edge. 
#' @export
spde_posterior_mean <- function(theta, graph, alpha = 1, obs = TRUE, leave_edge_out = FALSE) {
    if (alpha == 1 && obs == TRUE) {
        return(posterior.mean.obs.exp(theta = theta, graph = graph,
                                      leave.edge.out = leave_edge_out))
    } else if (alpha == 2 && obs == TRUE) {
        return(posterior.mean.obs.matern2(theta = theta, graph = graph,
                                          leave.edge.out = leave_edge_out))
    } else if(alpha == 1 && obs == FALSE) {
        return(posterior.mean.exp(theta = theta, graph = graph,
                                  rem.edge = leave_edge_out))
    } else if (alpha == 2 && obs == FALSE) {
        return(posterior.mean.matern2(theta = theta, graph = graph,
                                      rem.edge = leave_edge_out))
    } else {
        stop("alpha should be 1 or 2")
    }
}

#' Computes the posterior expectation for each observation in the graph
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph      - metric_graph object
#' @param leave.edge.out - compute the expectation of the graph if the observatrions are not on the edge
posterior.mean.obs.exp <- function(theta, graph, leave.edge.out = FALSE) {

  sigma_e = theta[1]

  Qp <- spde_precision(sigma= theta[2], kappa = theta[3], alpha = 1, graph = graph)
  if(leave.edge.out == FALSE)
    V.post <- posterior.mean.exp(theta, graph)

  y_hat <- rep(0, length(graph$y))
  Qpmu <- rep(0, nrow(graph$V))
  obs.edges <- unique(graph$PtE[,1])

  for (e in obs.edges) {

    if(leave.edge.out == TRUE)
      V.post <- posterior.mean.exp(theta, graph, rem.edge = e)

    obs.id <- graph$PtE[,1] == e
    y_i <- graph$y[obs.id]
    l <- graph$edge_lengths[e]
    D_matrix <- as.matrix(dist(c(0,l,graph$PtE[obs.id,2])))
    S <- r_1(D_matrix,theta[2:3])

    #covariance update see Art p.17
    E.ind <- c(1:2)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])

    E <- graph$E[e,]
    y_hat[obs.id] <- t(Bt)%*%V.post[E]
    if(leave.edge.out == FALSE){
      Sigma_i <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
      Sigma_noise <- Sigma_i
      diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2

      y_hat[obs.id] <- y_hat[obs.id] +  Sigma_i%*%solve(Sigma_noise,y_i-y_hat[obs.id] )
    }
  }
  return(y_hat)
}

#' Computes the posterior mean for each observation in the graph based
#' on the alpha=2 model
#' @param theta parameters (sigma_e, sigma, kappa)
#' @param graph metric_graph object
#' @param leave.edge.out compute the expectation of the graph if the
#' observatrions are not on the edge
posterior.mean.obs.matern2 <- function(theta, graph, leave.edge.out = FALSE) {

  sigma_e = theta[1]

  if(leave.edge.out == FALSE)
    E.post <- posterior.mean.matern2(theta, graph)

  y_hat <- rep(0, length(graph$y))
  obs.edges <- unique(graph$PtE[,1])

  for(e in obs.edges){

    if(leave.edge.out == TRUE)
      E.post <- posterior.mean.matern2(theta, graph, rem.edge = e)

    obs.id <- graph$PtE[,1] == e
    y_i <- graph$y[obs.id]
    l <- graph$edge_lengths[e]
    t <- c(0,l,graph$PtE[obs.id,2])
    D <- outer (t, t, `-`)
    S <- matrix(0, length(t) + 2, length(t) + 2)

    d.index <- c(1,2)
    S[-d.index, -d.index] <- r_2(D, theta[2:3])
    S[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))), theta[2:3], 2)
    S[d.index, -d.index] <- -r_2(D[1:2,], theta[2:3], 1)
    S[-d.index, d.index] <- t(S[d.index, -d.index])

    #covariance update see Art p.17
    E.ind <- c(1:4)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])

    u_e <- E.post[4 * (e - 1) + c(2, 4, 1, 3)]
    y_hat[obs.id] <- t(Bt) %*% u_e
    if (leave.edge.out == FALSE) {
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] - S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      Sigma_noise  <- Sigma_i
      diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2

      y_hat[obs.id] <- y_hat[obs.id] + Sigma_i%*%solve(Sigma_noise, y_i - y_hat[obs.id])
    }
  }
  return(y_hat)
}


#' Computes the posterior expectation for each node in the graph
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph - metric_graph object
#' @param rem.edge  - remove edge
posterior.mean.exp <- function(theta, graph, rem.edge = FALSE) {
  sigma_e <- theta[1]
  #build Q
  Qp.list <- spde_precision(kappa = theta[3], sigma = theta[2], alpha = 1,
                            graph = graph, build = FALSE)
  #build BSIGMAB
  Qpmu <- rep(0, graph$nV)

  obs.edges <- graph$nE
  if(is.logical(rem.edge) == FALSE)
    obs.edges <- setdiff(obs.edges, rem.edge)
  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))
  count <- 0
  for (e in obs.edges) {
    y_i <- graph$y[e]
    l <- graph$edge_lengths[e]
    D_matrix <- as.matrix(dist(c(0, l, graph$PtE[e, 2])))
    S <- r_1(D_matrix, theta[2:3])

    #covariance update see Art p.17
    E.ind <- c(1:2)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])
    Sigma_i <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    R <- chol(Sigma_i)
    Sigma_iB <- solve(Sigma_i, t(Bt))
    BtSinvB <- Bt %*% Sigma_iB

    E <- graph$E[e, ]
    if(E[1] == E[2]){
      Qpmu[E[1]] <- Qpmu[E[1]] + sum(t(Sigma_iB)%*%y_i)
      Qp[E[1],E[1]] <- Qp[E[1],E[1]] + sum(Bt %*% Sigma_iB)
      i_[count+1] <- E[1]
      j_[count+1] <- E[1]
      x_[count+1] <- sum(Bt %*% Sigma_iB)
      count <- count + 1
    }else{
      i_[count+(1:4)] <- c(E[1], E[1], E[2], E[2])
      j_[count+(1:4)] <- c(E[1], E[2], E[1], E[2])
      x_[count+(1:4)] <- c(BtSinvB[1,1], BtSinvB[1,2],
                           BtSinvB[1,2], BtSinvB[2,2])
      count <- count + 4
      Qpmu[E] <- Qpmu[E] + t(Sigma_iB) %*% y_i
    }
  }
  i_ <- c(Qp.list$i, i_[1:count])
  j_ <- c(Qp.list$j, j_[1:count])
  x_ <- c(Qp.list$x, x_[1:count])
  Qp <- Matrix::sparseMatrix(i = i_,
                             j = j_,
                             x = x_,
                             dims = Qp.list$dims)

  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Qpmu,system = 'P'), system='L')))
  Qpmu <- as.vector(Matrix::solve(R,Matrix::solve(R, v,system = 'Lt'), system='Pt'))

  return(Qpmu)

}

#' Computes the posterior mean for alpha = 2
#' @param theta parameters (sigma_e, sigma, kappa)
#' @param graph metric_graph object
posterior.mean.matern2 <- function(theta, graph, rem.edge = NULL) {
  sigma_e <- theta[1]
  
  n_const <- length(graph$CBobj$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CBobj$T[-ind.const,]
  Q <- spde_precision(kappa = theta[3], sigma = theta[2], alpha = 2, graph = graph)

  #build BSIGMAB
  Qpmu <- rep(0, 4*graph$nE)
  obs.edges <- unique(graph$PtE[,1])
  if(is.logical(rem.edge) == FALSE)
    obs.edges <- setdiff(obs.edges, rem.edge)

  i_ <- j_ <- x_ <- rep(0,16*length(obs.edges))
  count <- 0
  for(e in obs.edges){
    obs.id <- graph$PtE[,1] == e
    y_i <- graph$y[obs.id]
    l <- graph$edge_lengths[e]
    t <- c(0,l,graph$PtE[obs.id,2])

    D <- outer (t, t, `-`)
    S <- matrix(0, length(t) + 2, length(t) + 2)

    d.index <- c(1,2)
    S[-d.index, -d.index] <- r_2(D, theta[2:3])
    S[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))), theta[2:3], 2)
    S[d.index, -d.index] <- -r_2(D[1:2,], theta[2:3], 1)
    S[-d.index, d.index] <- t(S[d.index, -d.index])

    #covariance update see Art p.17
    E.ind <- c(1:4)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind,drop=F])
    Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] - S[Obs.ind, E.ind, drop = FALSE] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

    R <- chol(Sigma_i, pivot=T)
    Sigma_iB <- t(Bt)
    Sigma_iB[attr(R,"pivot"),] <- forwardsolve(R, backsolve(R, t(Bt[,attr(R,"pivot")]), transpose = TRUE), upper.tri = TRUE)
    BtSinvB <- Bt %*% Sigma_iB

    E <- graph$E[e,]
    if (E[1] == E[2]) {
      error("circle not implemented")
    } else {
      BtSinvB <- BtSinvB[c(3, 1, 4, 2), c(3, 1, 4, 2)]
      Qpmu[4 * (e - 1) + 1:4] <- Qpmu[4*(e-1)+1:4] + (t(Sigma_iB)%*%y_i)[c(3, 1, 4, 2)]

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
    }
  }
  i_ <- i_[1:count]
  j_ <- j_[1:count]
  x_ <- x_[1:count]
  BtSB <- Matrix::sparseMatrix(i = i_,
                               j = j_,
                               x = x_,
                               dims = dim(Q))
  Qp <- Q + BtSB
  Qp <- Tc%*%Qp%*%t(Tc)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Tc%*%Qpmu,system = 'P'), system='L')))
  Qpmu <- as.vector(Matrix::solve(R,Matrix::solve(R, v,system = 'Lt'), system='Pt'))

  return(t(Tc)%*%Qpmu)
}
