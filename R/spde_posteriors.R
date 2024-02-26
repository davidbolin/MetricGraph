
#' Computes the posterior expectation for SPDE models.
#' @param theta parameters (sigma_e, tau, kappa).
#' @param graph A `metric_graph` object.
#' @param alpha Smoothness parameter (1 or 2).
#' @param type Decides where to predict. Set to 'obs' for prediction at the
#' observation locations or to 'mesh' for prediction at the mesh locations.
#' @param leave_edge_out If `TRUE`, the posterior mean is computed for each
#' observation based on all observations which are not on that edge.
#' @return A vector with the posterior expectations.
#' @noRd
spde_posterior_mean <- function(theta,
                                graph,
                                alpha = 1,
                                type = "mesh",
                                leave_edge_out = FALSE) {

  check <- check_graph(graph)

  if (!(type %in% c("mesh", "obs"))) {
    stop("Type must be 'mesh' or 'obs'.")
  }
  if( type == "mesh" && !check$has.mesh) {
    stop("mesh must be provided")
  }
  if(!check$has.data){
    stop("The graph contains no data.")
  }
  if (alpha == 1) {
    return(posterior_mean_obs_alpha1(theta = theta, graph = graph,
                                     type = type,
                                     leave.edge.out = leave_edge_out))
  } else if (alpha == 2) {
    if(is.null(graph$CoB)){
      graph$buildC(2)
    }
    return(posterior_mean_obs_alpha2(theta = theta, graph = graph,
                                     type = type,
                                     leave.edge.out = leave_edge_out))
  } else {
    stop("alpha should be 1 or 2")
  }
}

#' Computes the posterior mean for the alpha=1 model
#' @param theta parameters (sigma_e, tau, kappa)
#' @param graph metric_graph object
#' @param type decides where to predict, 'obs' or 'mesh'.
#' @param leave.edge.out compute the mean of the graph if the observations
#' are not on the edge
#' @noRd
posterior_mean_obs_alpha1 <- function(theta,
                                      graph,
                                      resp, #resp must be in the graph's internal order
                                      PtE_resp,
                                      PtE_pred,
                                      type = "PtE",
                                      leave.edge.out = FALSE, no_nugget = FALSE) {

  sigma_e <- theta[1]
  tau <- theta[2]
  kappa <- theta[3]

  Qp <- spde_precision(tau= theta[2], kappa = theta[3],
                       alpha = 1, graph = graph)
  if(leave.edge.out == FALSE)
    V.post <- posterior_mean_alpha1(theta = theta, graph = graph,
                                    resp = resp, PtE_resp = PtE_resp, no_nugget = no_nugget)

  Qpmu <- rep(0, nrow(graph$V))
  if(type == "obs") {
    # y_hat <- rep(0, length(graph$y))
    y_hat <- rep(0, length(resp))
    obs.edges <- unique(PtE_resp[,1])
  }  else {
    y_hat <- rep(0, dim(PtE_pred)[1])
    obs.edges <- unique(PtE_pred[,1])
  }


  for (e in obs.edges) {
    if(leave.edge.out == TRUE)
      V.post <- posterior_mean_alpha1(theta = theta, graph = graph,
                                      rem.edge = e, resp = resp,
                                      PtE_resp = PtE_resp, no_nugget = no_nugget)

    obs.id <- which(PtE_resp[,1] == e)
    obs.loc <- PtE_resp[obs.id,2]
    y_i <- resp[obs.id]
    l <- graph$edge_lengths[e]
    if (type == "obs") {
      D <- as.matrix(dist(c(0,l, l*obs.loc)))
      S <- r_1(D,kappa = kappa, tau = tau)

      E.ind <- c(1:2)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind])

      y_hat[obs.id] <- t(Bt) %*% V.post[graph$E[e, ]]
      if(leave.edge.out == FALSE){
        Sigma_i <- S[Obs.ind, Obs.ind] - S[Obs.ind, E.ind] %*% Bt
        Sigma_noise <- Sigma_i
        if(!no_nugget){
          diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2  
        }

        y_hat[obs.id] <- y_hat[obs.id] + Sigma_i %*% solve(Sigma_noise,
                                                           y_i-y_hat[obs.id])
      }
    } else {
      pred.id <- PtE_pred[, 1] == e
      pred.loc <- PtE_pred[pred.id,2]
      D <- as.matrix(dist(c(0,l, l*obs.loc, l*pred.loc)))
      S <- r_1(D,kappa = kappa, tau = tau)
      E.ind <- c(1:2)
      Obs.ind <- 2 + seq_len(length(obs.loc))
      Pred.ind <- 2 + length(obs.loc) + seq_len(length(pred.loc))
      Bt_p <- solve(S[E.ind, E.ind], S[E.ind, Pred.ind])

      y_hat[pred.id] <- t(Bt_p) %*% V.post[graph$E[e, ]]
      if(leave.edge.out == FALSE && length(obs.loc)>0){
        Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind])
        Sigma_noise <- S[Obs.ind, Obs.ind] - S[Obs.ind, E.ind] %*% Bt
        Sigma_op <- S[Obs.ind, Pred.ind] - S[Obs.ind, E.ind] %*% Bt_p
        if(!no_nugget){
          diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2
        }
        y_hat_obs <- t(Bt) %*% V.post[graph$E[e, ]]

        y_hat[pred.id] <- y_hat[pred.id] + t(Sigma_op) %*% solve(Sigma_noise,
                                                           y_i-y_hat_obs)
      }
    }

  }
  if(type == "obs"){
    return(y_hat)
  } else {
    # return(c(V.post,y_hat))
    return(y_hat)
  }
}

#' Computes the posterior mean for the alpha=2 model
#' @param theta parameters (sigma_e, tau, kappa)
#' @param graph metric_graph object
#' @param leave.edge.out compute the expectation of the graph if the
#' @param type Set to 'obs' for computation at observation locations, or to
#' 'PtE' for computation at PtE locations.
#' @noRd
posterior_mean_obs_alpha2 <- function(theta,
                                      graph,
                                      resp, #resp must be in the graph's internal order
                                      PtE_resp,
                                      PtE_pred,
                                      type = "PtE",
                                      leave.edge.out = FALSE,
                                      no_nugget = FALSE) {


  if(type == "obs" && leave.edge.out) {
    stop("leave.edge.out only possible for type = 'obs'.")
  }
  sigma_e <- theta[1]
  tau <- theta[2]
  kappa <- theta[3]

  if(is.null(PtE_resp)){
    PtE <- graph$get_PtE()
  }
  if(is.null(resp)){
    stop("Please, provide 'resp'")
  }

  if(leave.edge.out == FALSE)
    E.post <- posterior_mean_alpha2(theta = theta, graph = graph,
                                    resp = resp, PtE_resp = PtE_resp,
                                    no_nugget = no_nugget)

  y_hat <- rep(0, length(resp))

  if (type == "obs") {
    y_hat <- rep(0, length(resp))
    obs.edges <- unique(PtE_resp[,1])
  }  else {
    y_hat <- rep(0, dim(PtE_pred)[1])
    obs.edges <- unique(PtE_pred[,1])
  }

  for(e in obs.edges){

    if(leave.edge.out == TRUE)
      E.post <- posterior_mean_alpha2(theta, graph, rem.edge = e, no_nugget = no_nugget)

    obs.id <- which(PtE_resp[, 1] == e)
    obs.loc <- PtE_resp[obs.id, 2]
    y_i <- resp[obs.id]
    l <- graph$edge_lengths[e]

    if(type == "obs") {
      t <- c(0, l, l * obs.loc)
      D <- outer (t, t, `-`)
      S <- matrix(0, length(t) + 2, length(t) + 2)

      d.index <- c(1, 2)
      S[-d.index, -d.index] <- r_2(D, kappa = kappa, tau = tau, deriv = 0)
      S[d.index, d.index] <- -r_2(as.matrix(dist(c(0, l))),
                                  kappa = kappa, tau = tau,
                                  deriv = 2)
      S[d.index, -d.index] <- -r_2(D[1:2, ], kappa = kappa,
                                   tau = tau, deriv = 1)
      S[-d.index, d.index] <- t(S[d.index, -d.index])

      #covariance update see Art p.17
      E.ind <- c(1:4)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind])

      u_e <- E.post[4 * (e - 1) + c(2, 4, 1, 3)]
      y_hat[obs.id] <- t(Bt) %*% u_e
      if (leave.edge.out == FALSE) {
        Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
          S[Obs.ind, E.ind, drop = FALSE] %*% Bt
        Sigma_noise  <- Sigma_i
        if(!no_nugget){
          diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2
        }

        y_hat[obs.id] <- y_hat[obs.id] + Sigma_i%*%solve(Sigma_noise,
                                                         y_i - y_hat[obs.id])
      }
    } else {
      pred.id <- PtE_pred[, 1] == e
      pred.loc <- PtE_pred[pred.id, 2]

      t <- c(0,l,l*obs.loc, l*pred.loc)
      D <- outer (t, t, `-`)
      S <- matrix(0, length(t) + 2, length(t) + 2)

      d.index <- c(1,2)
      S[-d.index, -d.index] <- r_2(D, kappa = kappa, tau = tau, deriv = 0)
      S[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                  kappa = kappa, tau = tau,
                                  deriv = 2)
      S[d.index, -d.index] <- -r_2(D[1:2,], kappa = kappa,
                                   tau = tau, deriv = 1)
      S[-d.index, d.index] <- t(S[d.index, -d.index])

      #covariance update see Art p.17
      E.ind <- c(1:4)
      Obs.ind <- 4 + seq_len(length(obs.loc))
      Pred.ind <- 4 + length(obs.loc) + seq_len(length(pred.loc))
      Bt_p <- solve(S[E.ind, E.ind],S[E.ind, Pred.ind])

      u_e <- E.post[4 * (e - 1) + c(2, 4, 1, 3)]
      u_e_tmp <- t(Bt_p) %*% u_e
      u_e_tmp <- u_e_tmp[,1]
      y_hat[pred.id] <- u_e_tmp

      if (leave.edge.out == FALSE && length(obs.loc)>0) {
        Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind])
        Sigma_noise <- S[Obs.ind, Obs.ind, drop = FALSE] -
          S[Obs.ind, E.ind, drop = FALSE] %*% Bt
        Sigma_op <- S[Obs.ind, Pred.ind] - S[Obs.ind, E.ind] %*% Bt_p
        if(!no_nugget){
          diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2
        }
        Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind])
        u_e <- E.post[4 * (e - 1) + c(2, 4, 1, 3)]
        y_hat_obs <- t(Bt) %*% u_e
        y_hat[pred.id] <- y_hat[pred.id] + t(Sigma_op) %*% solve(Sigma_noise,
                                                                 y_i - y_hat_obs)
      }
    }
  }
  return(y_hat)
}


#' Computes the posterior expectation for each node in the graph
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph - metric_graph object
#' @param rem.edge  - remove edge
#' @noRd
posterior_mean_alpha1 <- function(theta, graph, resp,
                                  PtE_resp, rem.edge = FALSE,
                                  no_nugget = FALSE) {

  sigma_e <- theta[1]
  tau <- theta[2]
  kappa <- theta[3]

  Qp.list <- spde_precision(kappa = theta[3], tau = theta[2], alpha = 1,
                            graph = graph, build = FALSE)
  #build BSIGMAB
  Qpmu <- rep(0, graph$nV)

  # obs.edges <- unique(graph$PtE[,1])
  obs.edges <- unique(PtE_resp[,1])
  if(is.logical(rem.edge) == FALSE)
    obs.edges <- setdiff(obs.edges, rem.edge)
  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))
  count <- 0
  for (e in obs.edges) {
    # obs.id <- graph$PtE[,1] == e
    obs.id <- PtE_resp[,1] == e
    # y_i <- graph$y[obs.id]
    y_i <- resp[obs.id]
    l <- graph$edge_lengths[e]
    # D_matrix <- as.matrix(dist(c(0, l, l*graph$PtE[obs.id, 2])))
    D_matrix <- as.matrix(dist(c(0, l, l*PtE_resp[obs.id, 2])))
    S <- r_1(D_matrix, kappa = kappa, tau = tau)

    #covariance update see Art p.17
    E.ind <- c(1:2)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])
    Sigma_i <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
    if(!no_nugget){
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    }
    R <- base::chol(Sigma_i)
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

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Qpmu,system = 'P'),
                                 system='L')))
  Qpmu <- as.vector(Matrix::solve(R,Matrix::solve(R, v,system = 'Lt'),
                                  system='Pt'))

  return(Qpmu)

}

#' Computes the posterior mean for alpha = 2
#' @param theta parameters (sigma_e, tau, kappa)
#' @param graph metric_graph object
#' @noRd
posterior_mean_alpha2 <- function(theta, graph, resp,
                                  PtE_resp, rem.edge = NULL,
                                  no_nugget = FALSE) {

  sigma_e <- theta[1]
  tau <- theta[2]
  kappa <- theta[3]
  if(is.null(PtE_resp)){
    PtE_resp <- graph$get_PtE()
  }

  PtE <- PtE_resp

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const,]

  Q <- spde_precision(kappa = theta[3], tau = theta[2],
                      alpha = 2, graph = graph)


  #build BSIGMAB
  Qpmu <- rep(0, 4 * graph$nE)
  obs.edges <- unique(PtE[, 1])
  if(is.logical(rem.edge) == FALSE)
    obs.edges <- setdiff(obs.edges, rem.edge)

  i_ <- j_ <- x_ <- rep(0, 16 * length(obs.edges))
  count <- 0

  for (e in obs.edges) {
    obs.id <- PtE[, 1] == e
    y_i <- resp[obs.id]
    l <- graph$edge_lengths[e]
    t <- c(0, l, l * PtE[obs.id, 2])

    D <- outer (t, t, `-`)
    S <- matrix(0, length(t) + 2, length(t) + 2)

    d.index <- c(1, 2)
    S[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                 tau = tau, deriv = 0)
    S[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                kappa = kappa, tau = tau,
                                deriv = 2)
    S[d.index, -d.index] <- -r_2(D[1:2,], kappa = kappa,
                                 tau = tau, deriv = 1)
    S[-d.index, d.index] <- t(S[d.index, -d.index])

    #covariance update see Art p.17
    E.ind <- c(1:4)
    Obs.ind <- -E.ind
    Bt <- solve(S[E.ind, E.ind],
                S[E.ind, Obs.ind, drop = FALSE])
    Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
      S[Obs.ind, E.ind, drop = FALSE] %*% Bt
    if(!no_nugget){
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    }

    R <- base::chol(Sigma_i, pivot=T)
    Sigma_iB <- t(Bt)
    Sigma_iB[attr(R,"pivot"),] <- base::forwardsolve(R,
                                               base::backsolve(R,
                                                         t(Bt[,attr(R,"pivot")]),
                                                         transpose = TRUE),
                                               upper.tri = TRUE)
    BtSinvB <- Bt %*% Sigma_iB

    E <- graph$E[e,]
    if (E[1] == E[2]) {
      stop("circle not implemented")
    } else {
      BtSinvB <- BtSinvB[c(3, 1, 4, 2), c(3, 1, 4, 2)]
      Qpmu[4 * (e - 1) + 1:4] <- Qpmu[4*(e-1)+1:4] +
        (t(Sigma_iB)%*%y_i)[c(3, 1, 4, 2)]

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

  v <- c(as.matrix(Matrix::solve(R,
                                 Matrix::solve(R,
                                               Tc%*%Qpmu,
                                               system = 'P'),
                                 system='L')))
  Qpmu <- as.vector(Matrix::solve(R,
                                  Matrix::solve(R,
                                                v,
                                                system = 'Lt'),
                                  system='Pt'))

  return(t(Tc)%*%Qpmu)
}

