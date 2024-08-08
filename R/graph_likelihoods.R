# #' Function factory for likelihood evaluation for the metric graph SPDE model
# #'
# #' @param graph metric_graph object
# #' @param alpha Order of the SPDE, should be either 1 or 2.
# #' @param data_name Name of the response variable
# #' @param covariates OBSOLETE
# #' @param log_scale Should the parameters `theta` of the returning function be
# #' given in log scale?
# #' @param version if 1, the likelihood is computed by integrating out
# #' @param maximize If `FALSE` the function will return minus the likelihood, so
# #' one can directly apply it to the `optim` function.
# #' @param BC which boundary condition to use (0,1) 0 is no adjustment on boundary point
# #'        1 is making the boundary condition stationary
# #' @return The log-likelihood function, which is returned as a function with
# #' parameter 'theta'.
# #' The parameter `theta` must be supplied as
# #' the vector `c(sigma_e, sigma, kappa)`.
# #'
# #' If `covariates` is `TRUE`, then the parameter `theta` must be supplied as the
# #' vector `c(sigma_e, sigma, kappa, beta[1], ..., beta[p])`,
# #' where `beta[1],...,beta[p]` are the coefficients and `p` is the number of
# #' covariates.
# #' @noRd

# likelihood_graph_spde <- function(graph,
#                                   alpha = 1,
#                                   covariates = FALSE,
#                                   data_name,
#                                   log_scale = TRUE,
#                                   maximize = FALSE,
#                                   version = 1,
#                                   repl=NULL,
#                                   BC = 1) {

#   check <- check_graph(graph)

#   if(!(alpha%in%c(1,2))){
#     stop("alpha must be either 1 or 2!")
#   }

#   loglik <- function(theta){
#         if(log_scale){
#           theta_spde <- exp(theta[1:3])
#           theta_spde <- c(theta_spde, theta[-c(1:3)])
#         } else{
#           theta_spde <- theta
#         }

#       switch(alpha,
#       "1" = {
#         if(version == 1){
#           loglik_val <- likelihood_alpha1(theta_spde, graph, data_name, covariates,BC=BC)
#         } else if(version == 2){
#           loglik_val <- likelihood_alpha1_v2(theta_spde, graph, X_cov, y, repl,BC=BC)
#         } else{
#           stop("Version should be either 1 or 2!")
#         }
#       },
#       "2" = {
#         loglik_val <- likelihood_alpha2(theta_spde, graph, data_name, covariates, BC=BC)
#       }
#       )
#       if(maximize){
#         return(loglik_val)
#       } else{
#         return(-loglik_val)
#       }
#   }
# }



#' Log-likelihood calculation for alpha=1 model directional
#' @param theta (sigma_e, reciprocal_tau, kappa)
#' @param graph metric_graph object
#' @param data_name name of the response variable
#' @param repl replicates
#' @param X_cov matrix of covariates
#' @noRd
likelihood_alpha1_directional <- function(theta,
                                          graph,
                                          data_name = NULL,
                                          manual_y = NULL,
                                          X_cov = NULL,
                                          repl = NULL,
                                          parameterization="matern") {
  #build Q
  if(is.null(graph$C)){
    graph$buildDirectionalConstraints(alpha = 1)
  } else if(graph$CoB$alpha == 2){
    graph$buildDirectionalConstraints(alpha = 1)
  }



  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  sigma_e <- exp(theta[1])
  if(parameterization == "matern"){
    kappa = sqrt(8 * 0.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  reciprocal_tau <- exp(theta[2])
  Q.list <- Qalpha1_edges(c( 1/reciprocal_tau,kappa),
                           graph,
                           w = 0,
                           BC=1, build=FALSE)
  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
  Q <- Matrix::sparseMatrix(i = Q.list$i,
                            j = Q.list$j,
                            x = Q.list$x,
                            dims = Q.list$dims)
  R <- Matrix::Cholesky(forceSymmetric(Tc%*%Q%*%t(Tc)),
                        LDL = FALSE, perm = TRUE)
  det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]

  #build BSIGMAB
  PtE <- graph$get_PtE()
  obs.edges <- unique(PtE[, 1])

  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))

  if(is.null(repl)){
    u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  ind_repl <- (graph$.__enclos_env__$private$data[[".group"]] %in% u_repl)

  loglik <- 0

  det_R_count <- NULL
  if(is.null(manual_y)){
    y_resp <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y_resp <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }
  n.o <- sum(ind_repl)

  for(repl_y in 1:length(u_repl)){
    loglik <- loglik + det_R
    count <- 0
    Qpmu <- rep(0, 2*nrow(graph$E))
    for (e in obs.edges) {
      ind_repl <- graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y]
      obs.id <- PtE[,1] == e
      y_i <- y_resp[ind_repl]
      y_i <- y_i[obs.id]
      idx_na <- is.na(y_i)
      y_i <- y_i[!idx_na]

      if(sum(!idx_na) == 0){
        next
      }

      if(!is.null(X_cov)){
        n_cov <- ncol(X_cov)
        if(n_cov == 0){
          X_cov_repl <- 0
        } else{
          X_cov_repl <- X_cov[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y], , drop=FALSE]
          X_cov_repl <- X_cov_repl[PtE[,1] == e, ,drop = FALSE]
          X_cov_repl <- X_cov_repl[!idx_na, , drop = FALSE]
          y_i <- y_i - X_cov_repl %*% theta[4:(3+n_cov)]
        }
      }

      l <- graph$edge_lengths[e]

      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]

      D_matrix <- as.matrix(dist(c(0, l, l*PtE_temp)))

      S <- r_1(D_matrix, kappa = kappa, tau = 1/reciprocal_tau)

      #covariance update see Art p.17
      E.ind <- c(1:2)
      Obs.ind <- -E.ind

      Bt <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R <- base::chol(Sigma_i)
      Sigma_iB <- solve(Sigma_i, t(Bt))
      BtSinvB <- Bt %*% Sigma_iB

      E <- graph$E[e, ]
      if (E[1] == E[2]) {
        Qpmu[2*(e-1)+1] <- Qpmu[2*(e-1)+1] + sum(t(Sigma_iB) %*% y_i)
        i_[count + 1] <- 2*(e-1)+1
        j_[count + 1] <- 2*(e-1)+1
        x_[count + 1] <- sum(Bt %*% Sigma_iB)
      } else {
        Qpmu[2*(e-1) + c(1, 2)] <- Qpmu[2*(e-1) + c(1, 2)] + t(Sigma_iB) %*% y_i
        i_[count + (1:4)] <- c(2*(e-1)+1, 2*(e-1)+1, 2*(e-1)+2, 2*(e-1)+2)
        j_[count + (1:4)] <- c(2*(e-1)+1, 2*(e-1)+2, 2*(e-1)+1, 2*(e-1)+2)
        x_[count + (1:4)] <- c(BtSinvB[1, 1], BtSinvB[1, 2],
                               BtSinvB[1, 2], BtSinvB[2, 2])
        count <- count + 4
      }

      loglik <- loglik - 0.5 * t(y_i) %*% solve(Sigma_i, y_i)
      loglik <- loglik - sum(log(diag(R)))

    }

    if(is.null(det_R_count)){
      i_ <- c(Q.list$i, i_[1:count])
      j_ <- c(Q.list$j, j_[1:count])
      x_ <- c(Q.list$x, x_[1:count])


      Qp <- Matrix::sparseMatrix(i = i_,
                                 j = j_,
                                 x = x_,
                                 dims = Q.list$dims)
      Qp <- Tc %*% Qp %*% t(Tc)
      R_count <- Matrix::Cholesky(forceSymmetric(Qp), LDL = FALSE, perm = TRUE)
      det_R_count <- Matrix::determinant(R_count, sqrt=TRUE)$modulus[1]
    }

    loglik <- loglik - det_R_count

    v <- c(as.matrix(Matrix::solve(R_count, Matrix::solve(R_count, Tc%*%Qpmu,
                                                          system = "P"),
                                   system = "L")))

    loglik <- loglik + 0.5  * t(v) %*% v - 0.5 * n.o * log(2*pi)
  }

  return(loglik[1])
}


#' Computes the log likelihood function fo theta for the graph object
#' @param theta parameters (sigma_e, reciprocal_tau, kappa)
#' @param graph  metric_graph object
#' @param data_name name of the response variable
#' @param BC which boundary condition to use (0,1)
#' @param covariates OBSOLETE
#' @noRd
likelihood_alpha2 <- function(theta, graph, data_name = NULL, manual_y = NULL,
                             X_cov = NULL, repl, BC, parameterization) {
  if(is.null(graph$C)){
    graph$buildC(2)
  } else if(graph$CoB$alpha == 1){
     graph$buildC(2)
  }

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  if(parameterization == "matern"){
    kappa = sqrt(8 * 1.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  if(is.null(repl)){
    u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  ind_repl <- (graph$.__enclos_env__$private$data[[".group"]] %in% u_repl)

  if(is.null(manual_y)){
    y <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }

  n.o <- sum(ind_repl)

  #build Q

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
  Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                      alpha = 2, graph = graph, BC=BC)
  R <- Matrix::Cholesky(forceSymmetric(Tc%*%Q%*%t(Tc)),
                        LDL = FALSE, perm = TRUE)

  PtE <- graph$get_PtE()
  obs.edges <- unique(PtE[, 1])

  loglik <- 0
  det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]
  # n.o <- length(graph$.__enclos_env__$private$data[[data_name]])
  # u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  det_R_count <- NULL

  # y <- graph$.__enclos_env__$private$data[[data_name]]



  for(repl_y in 1:length(u_repl)){
      loglik <- loglik + det_R
      Qpmu <- rep(0, 4 * nrow(graph$E))
      # y_rep <- graph$.__enclos_env__$private$data[[data_name]][graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y]]
      y_rep <- y[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y]]
      #build BSIGMAB

      i_ <- j_ <- x_ <- rep(0, 16 * length(obs.edges))
      count <- 0
      for (e in obs.edges) {
        obs.id <- PtE[,1] == e

        y_i <- y_rep[obs.id]

      if(!is.null(X_cov)){
          n_cov <- ncol(X_cov)
          if(n_cov == 0){
            X_cov_repl <- 0
          } else{
            X_cov_repl <- X_cov[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y], , drop=FALSE]
            X_cov_repl <- X_cov_repl[obs.id, ,drop = FALSE]
            y_i <- y_i - X_cov_repl %*% theta[4:(3+n_cov)]
          }
      }

        l <- graph$edge_lengths[e]
        t <- c(0, l, l*PtE[obs.id, 2])

        D <- outer (t, t, `-`)
        S <- matrix(0, length(t) + 2, length(t) + 2)

        d.index <- c(1,2)
        S[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                     tau = 1/reciprocal_tau, deriv = 0)
        S[ d.index, d.index] <- -r_2(as.matrix(dist(c(0,l))),
                                     kappa = kappa, tau = 1/reciprocal_tau,
                                     deriv = 2)
        S[d.index, -d.index] <- -r_2(D[1:2,], kappa = kappa,
                                     tau = 1/reciprocal_tau, deriv = 1)
        S[-d.index, d.index] <- t(S[d.index, -d.index])

        #covariance update see Art p.17
        E.ind <- c(1:4)
        Obs.ind <- -E.ind
        Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind, drop = FALSE])
        Sigma_i <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
        diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

        R <- base::chol(Sigma_i, pivot = TRUE)
        if(attr(R, "rank") < dim(R)[1])
          return(-Inf)

        Sigma_iB <- t(Bt)
        Sigma_iB[attr(R,"pivot"),] <- base::forwardsolve(R,
                                            base::backsolve(R, t(Bt[, attr(R,"pivot")]),
                                            transpose = TRUE), upper.tri = TRUE)


        # R <- Matrix::Cholesky(Sigma_i)
        # Sigma_iB <- solve(R, t(Bt), system = "A")

        BtSinvB <- Bt %*% Sigma_iB

        E <- graph$E[e, ]
        if (E[1] == E[2]) {
          warning("Circle not implemented")
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
        # loglik <- loglik - 0.5 * c(determinant(R, logarithm = TRUE)$modulus)

      }
      if(is.null(det_R_count)){
        i_ <- i_[1:count]
        j_ <- j_[1:count]
        x_ <- x_[1:count]
        BtSB <- Matrix::sparseMatrix(i = i_,
                                     j = j_,
                                     x = x_,
                                     dims = dim(Q))
        Qp <- Q + BtSB
        Qp <- Tc %*% Qp %*% t(Tc)
        R_count <- Matrix::Cholesky(forceSymmetric(Qp), LDL = FALSE, perm = TRUE)
        det_R_count <- Matrix::determinant(R_count, sqrt=TRUE)$modulus[1]
      }
      loglik <- loglik - det_R_count

      v <- c(as.matrix(Matrix::solve(R_count, Matrix::solve(R_count, Tc%*%Qpmu, system = 'P'),
                                     system='L')))

      loglik <- loglik + 0.5  * t(v) %*% v  - 0.5 * n.o * log(2 * pi)
  }

  return(loglik[1])
}



#' Log-likelihood calculation for alpha=1 model
#'
#' @param theta (sigma_e, reciprocal_tau, kappa)
#' @param graph metric_graph object
#' @param X_cov matrix of covariates
#' @param y response vector
#' @param repl replicates to be considered
#' @param BC boundary conditions
#' @param parameterization parameterization to be used.
#' @details This function computes the likelihood without integrating out
#' the vertices.
#' @return The log-likelihood
#' @noRd
likelihood_alpha1_v2 <- function(theta, graph, X_cov, y, repl, BC, parameterization) {

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  if(parameterization == "matern"){
    kappa = sqrt(8 * 0.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  #build Q
  Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                      alpha = 1, graph = graph, BC=BC)
  if(is.null(graph$PtV)){
    stop("No observation at the vertices! Run observation_to_vertex().")
  }

  # R <- chol(Q)
  # R <- Matrix::chol(Q)

  R <- Matrix::Cholesky(Q)

  l <- 0

  for(i in repl){
      A <- Matrix::Diagonal(graph$nV)[graph$PtV, ]
      ind_tmp <- (repl_vec %in% i)
      y_tmp <- y[ind_tmp]
      if(ncol(X_cov) == 0){
        X_cov_tmp <- 0
      } else {
        X_cov_tmp <- X_cov[ind_tmp,,drop=FALSE]
      }
      na_obs <- is.na(y_tmp)

      y_ <- y_tmp[!na_obs]
      n.o <- length(y_)
      Q.p <- Q  + t(A[!na_obs,]) %*% A[!na_obs,]/sigma_e^2
      # R.p <- Matrix::chol(Q.p)
      R.p <- Matrix::Cholesky(Q.p)

      # l <- l + sum(log(diag(R))) - sum(log(diag(R.p))) - n.o*log(sigma_e)

      l <- l + determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus - determinant(R.p, logarithm = TRUE, sqrt = TRUE)$modulus - n.o * log(sigma_e)

      v <- y_

      if(ncol(X_cov) != 0){
        X_cov_tmp <- X_cov_tmp[!na_obs, ]
        v <- v - X_cov_tmp %*% theta[4:(3+ncol(X_cov))]
      }

      # mu.p <- solve(Q.p,as.vector(t(A[!na_obs,]) %*% v / sigma_e^2))

      mu.p <- solve(R.p, as.vector(t(A[!na_obs,]) %*% v / sigma_e^2), system = "A")

      v <- v - A[!na_obs,]%*%mu.p

      l <- l - 0.5*(t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
        0.5 * n.o * log(2*pi)

  }

  return(as.double(l))
}

#' Log-likelihood calculation for alpha=1 model
#' @param theta (sigma_e, reciprocal_tau, kappa)
#' @param graph metric_graph object
#' @param data_name name of the response variable
#' @param repl replicates
#' @param X_cov matrix of covariates
#' @param BC. - which boundary condition to use (0,1)
#' @noRd
likelihood_alpha1 <- function(theta, graph, data_name = NULL, manual_y = NULL,
                             X_cov = NULL, repl, BC, parameterization) {
  sigma_e <- exp(theta[1])
  #build Q

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  if(parameterization == "matern"){
    kappa = sqrt(8 * 0.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  reciprocal_tau <- exp(theta[2])

  Q.list <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau, alpha = 1,
                           graph = graph, build = FALSE,BC=BC)

  Qp <- Matrix::sparseMatrix(i = Q.list$i,
                             j = Q.list$j,
                             x = Q.list$x,
                             dims = Q.list$dims)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]

  #build BSIGMAB
  PtE <- graph$get_PtE()
  obs.edges <- unique(PtE[, 1])

  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))

  if(is.null(repl)){
    u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  ind_repl <- (graph$.__enclos_env__$private$data[[".group"]] %in% u_repl)

  loglik <- 0

  det_R_count <- NULL
  if(is.null(manual_y)){
    y_resp <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y_resp <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }
  n.o <- sum(ind_repl)

  for(repl_y in 1:length(u_repl)){
    loglik <- loglik + det_R
    count <- 0
    Qpmu <- rep(0, nrow(graph$V))
    for (e in obs.edges) {
      ind_repl <- graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y]
      obs.id <- PtE[,1] == e
      y_i <- y_resp[ind_repl]
      y_i <- y_i[obs.id]
      idx_na <- is.na(y_i)
      y_i <- y_i[!idx_na]

      if(sum(!idx_na) == 0){
        next
      }

      #   if(covariates){ #obsolete
      #     n_cov <- ncol(graph$covariates[[1]])
      #     if(length(graph$covariates)==1){
      #     X_cov <- graph$covariates[[1]]
      #     } else if(length(graph$covariates) == ncol(graph$.__enclos_env__$private$data[[data_name]])){
      #       X_cov <- graph$covariates[[repl_y]]
      #     } else{
      #       stop("You should either have a common covariate for all the replicates, or one set of covariates for each replicate!")
      #     }
      #   X_cov <- X_cov[obs.id,]
      #   y_i <- y_i - X_cov %*% theta[4:(3+n_cov)]
      # }

      if(!is.null(X_cov)){
          n_cov <- ncol(X_cov)
          if(n_cov == 0){
            X_cov_repl <- 0
          } else{
            X_cov_repl <- X_cov[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y], , drop=FALSE]
            X_cov_repl <- X_cov_repl[PtE[,1] == e, ,drop = FALSE]
            X_cov_repl <- X_cov_repl[!idx_na, , drop = FALSE]
            y_i <- y_i - X_cov_repl %*% theta[4:(3+n_cov)]
          }
      }

      l <- graph$edge_lengths[e]

      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]

      D_matrix <- as.matrix(dist(c(0, l, l*PtE_temp)))

      S <- r_1(D_matrix, kappa = kappa, tau = 1/reciprocal_tau)

      #covariance update see Art p.17
      E.ind <- c(1:2)
      Obs.ind <- -E.ind

      Bt <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R <- base::chol(Sigma_i)
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

    if(is.null(det_R_count)){
      i_ <- c(Q.list$i, i_[1:count])
      j_ <- c(Q.list$j, j_[1:count])
      x_ <- c(Q.list$x, x_[1:count])


        Qp <- Matrix::sparseMatrix(i = i_,
                             j = j_,
                             x = x_,
                             dims = Q.list$dims)

        R_count <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

        det_R_count <- Matrix::determinant(R_count, sqrt=TRUE)$modulus[1]
    }

    loglik <- loglik - det_R_count

    v <- c(as.matrix(Matrix::solve(R_count, Matrix::solve(R_count, Qpmu,
                                                          system = "P"),
                                   system = "L")))

    loglik <- loglik + 0.5  * t(v) %*% v - 0.5 * n.o * log(2*pi)
  }

  return(loglik[1])
}



#' Function factory for likelihood evaluation not using sparsity
#' @param graph A `metric_graph` object.
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL1" gives
#' the model based on the graph Laplacian with smoothness 1, "GL2" gives the
#' model based on the graph Laplacian with smoothness 2, and "isoCov" gives a
#' model with isotropic covariance.
#' @param cov_function The covariance function to be used in case 'model' is
#' chosen as 'isoCov'. `cov_function` must be a function of `(h, theta_cov)`,
#' where `h` is a vector, or matrix, containing the distances to evaluate the
#' covariance function at, and `theta_cov` is the vector of parameters of the
#' covariance function `cov_function`.
#' @param y_graph Response vector given in the same order as the internal
#' locations from the graph.
#' @param X_cov Matrix with covariates. The order must be the same as the
#' internal order from the graph.
#' @param repl Vector with the replicates to be considered. If `NULL` all
#' replicates will be considered.
#' @param log_scale Should the parameters `theta` of the returning function be
#' given in log-scale?
#' @param maximize If `FALSE` the function will return minus the likelihood, so
#' one can directly apply it to the `optim` function.
#' @return The log-likelihood function.
#' @details The log-likelihood function that is returned is a function of a
#' parameter `theta`. For models 'alpha1', 'alpha2', 'GL1' and 'GL2', the
#' parameter `theta` must be supplied as the vector `c(sigma_e, sigma, kappa)`.
#'
#' For 'isoCov' model, theta must be a vector such that `theta[1]` is `sigma.e`
#' and the vector `theta[2:(q+1)]` is the input of `cov_function`, where `q` is
#' the number of parameters of the covariance function.
#'
#' If `covariates` is `TRUE`, then the parameter `theta` must be supplied as the
#' vector `c(sigma_e, theta[2], ..., theta[q+1], beta[1], ..., beta[p])`,
#' where `beta[1],...,beta[p]` are the coefficients and `p` is the number of
#' covariates.
#'
#' For the remaining models, if `covariates` is `TRUE`, then `theta` must be
#' supplied as the vector `c(sigma_e, sigma, kappa, beta[1], ..., beta[p])`,
#' where `beta[1],...,beta[p]` are the coefficients and `p` is the number of
#' covariates.
#' @noRd
likelihood_graph_covariance <- function(graph,
                                        model = "alpha1",
                                        y_graph,
                                        cov_function = NULL,
                                        X_cov = NULL,
                                        repl,
                                        log_scale = TRUE,
                                        maximize = FALSE,
                                        fix_vec = NULL,
                                        fix_v_val = NULL,
                                        check_euclidean = TRUE) {

  # check <- check_graph(graph)

  if(!(model%in%c("WM1", "WM2", "GL1", "GL2", "isoCov"))){
    stop("The available models are: 'WM1', 'WM2', 'GL1', 'GL2' and 'isoCov'!")
  }

  loglik <- function(theta){

      if(!is.null(X_cov)){
            n_cov <- ncol(X_cov)
      } else{
            n_cov <- 0
      }

    if(!is.null(fix_v_val)){
      # new_theta <- fix_v_val
      fix_v_val_full <- c(fix_v_val, rep(NA, n_cov))
      fix_vec_full <- c(fix_vec, rep(FALSE, n_cov))
      new_theta <- fix_v_val_full
      new_theta[!fix_vec_full] <- theta
    } else{
      new_theta <- theta
    }


      if(model == "isoCov"){
        if(log_scale){
          sigma_e <- exp(new_theta[1])
          theta_cov <- exp(new_theta[2:(length(new_theta)-n_cov)])
        } else{
          sigma_e <- new_theta[1]
          theta_cov <- new_theta[2:(length(new_theta)-n_cov)]
        }
      } else{
        if(log_scale){
          sigma_e <- exp(new_theta[1])
          reciprocal_tau <- exp(new_theta[2])
          kappa <- exp(new_theta[3])
        } else{
          sigma_e <- new_theta[1]
          reciprocal_tau <- new_theta[2]
          kappa <- new_theta[3]
        }
      }

      if(n_cov >0){
        theta_covariates <- new_theta[(length(new_theta)-n_cov+1):length(new_theta)]
      }

      if(is.null(graph$Laplacian) && (model %in% c("GL1", "GL2"))) {
        graph$compute_laplacian()
      }

      #build covariance matrix
      switch(model,
      WM1 = {

        Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                            alpha = 1, graph = graph)
        Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]

      },
      WM2 = {
        PtE <- graph$get_PtE()
        n.c <- 1:length(graph$CoB$S)
        Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau, alpha = 2,
                            graph = graph, BC = 1)
        Qtilde <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
        Qtilde <- Qtilde[-n.c,-n.c]
        Sigma.overdetermined  = t(graph$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
          (graph$CoB$T[-n.c,])
        index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-14) +
          3.0 * (abs(PtE[, 2]) > 1e-14)
        Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])

      }, GL1 = {
        Q <- (kappa^2 * Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]) / reciprocal_tau^2
        Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
      }, GL2 = {

        Q <- kappa^2 * Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]
        Q <- Q %*% Q / reciprocal_tau^2
        Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]

      }, isoCov = {
        if(is.null(cov_function)){
          stop("If model is 'isoCov' the covariance function must be supplied!")
        }
        if(!is.function(cov_function)){
          stop("'cov_function' must be a function!")
        }

        if(is.null(graph$res_dist)){
          graph$compute_resdist(full = TRUE, check_euclidean = check_euclidean)
        }

        Sigma <- as.matrix(cov_function(as.matrix(graph$res_dist[[".complete"]]), theta_cov))
        # nV <- nrow(graph$res_dist[[".complete"]]) - nrow(graph$get_PtE())

        # Sigma <- Sigma[(nV+1):nrow(Sigma), (nV+1):nrow(Sigma)]
      })

      diag(Sigma) <- diag(Sigma) + sigma_e^2

      loglik_val <- 0

      repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

      if(is.null(repl)){
        u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
      } else{
        u_repl <- unique(repl)
      }


      for(repl_y in 1:length(u_repl)){
          ind_tmp <- (repl_vec %in% u_repl[repl_y])
          y_tmp <- y_graph[ind_tmp]
          na_obs <- is.na(y_tmp)
          Sigma_non_na <- Sigma[!na_obs, !na_obs]
          R <- base::chol(Sigma_non_na)
          v <- y_graph[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y]]

          if(!is.null(X_cov)){
              n_cov <- ncol(X_cov)
              if(n_cov == 0){
                X_cov_repl <- 0
              } else{
                X_cov_repl <- X_cov[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y], ,
                                    drop=FALSE]
                # v <- v - X_cov_repl %*% theta[4:(3+n_cov)]
                v <- v - X_cov_repl %*% theta_covariates
              }
          }

          v <- v[!na_obs]

          loglik_val <- loglik_val + as.double(-sum(log(diag(R))) - 0.5*t(v)%*%solve(Sigma_non_na,v) -
                         length(v)*log(2*pi)/2)
      }

      if(maximize){
        return(loglik_val)
      } else{
        return(-loglik_val)
      }
      }

}


#' Function factory for likelihood evaluation for the graph Laplacian model
#'
#' @param graph metric_graph object
#' @param alpha integer greater or equal to 1. Order of the equation.
#' @param covariates Logical. If `TRUE`, the model will be considered with
#' covariates. It requires `graph` to have covariates included by the method
#' `add_covariates()`.
#' @param log_scale Should the parameters `theta` of the returning function be
#' given in log-scale?
#' @param maximize If `FALSE` the function will return minus the likelihood, so
#' one can directly apply it to the `optim` function.
#' @return The log-likelihood function, which is returned as a function with
#' parameter 'theta'. The parameter `theta` must be supplied as
#' the vector `c(sigma_e, reciprocal_tau, kappa)`.
#'
#' If `covariates` is `TRUE`, then the parameter `theta` must be supplied as the
#' vector `c(sigma_e, reciprocal_tau, kappa, beta[1], ..., beta[p])`,
#' where `beta[1],...,beta[p]` are the coefficients and `p` is the number of
#' covariates.
#' @noRd
likelihood_graph_laplacian <- function(graph, alpha, y_graph, repl,
              X_cov = NULL, maximize = FALSE, parameterization,
              fix_vec = NULL, fix_v_val = NULL) {

  check <- check_graph(graph)

  if(alpha%%1 != 0){
    stop("only integer values of alpha supported")
  }

  if(alpha<= 0){
    stop("alpha must be positive!")
  }

  graph$compute_laplacian(full = FALSE)

  # if(covariates){
  #   if(is.null(graph$covariates)){
  #     stop("If 'covariates' is set to TRUE, the graph must have covariates!")
  #   }
  # }

  loglik <- function(theta){

      if(!is.null(X_cov)){
            n_cov <- ncol(X_cov)
      } else{
            n_cov <- 0
      }

    if(!is.null(fix_v_val)){
      # new_theta <- fix_v_val
      fix_v_val_full <- c(fix_v_val, rep(NA, n_cov))
      fix_vec_full <- c(fix_vec, rep(FALSE, n_cov))
      new_theta <- fix_v_val_full
    }
    if(!is.null(fix_vec)){
      new_theta[!fix_vec_full] <- theta
    } else{
      new_theta <- theta
    }


    repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

    if(is.null(repl)){
      u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
    } else{
      u_repl <- unique(repl)
    }


    sigma_e <- exp(new_theta[1])
    reciprocal_tau <- exp(new_theta[2])
    if(parameterization == "matern"){
      kappa = sqrt(8 * (alpha-0.5)) / exp(new_theta[3])
    } else{
      kappa = exp(new_theta[3])
    }

    y_resp <- y_graph

    l <- 0
    A <- graph$.__enclos_env__$private$A(group = ".all", drop_all_na = FALSE, drop_na = FALSE)

    u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])

    for(repl_y in 1:length(u_repl)){
      K <- kappa^2*Diagonal(graph$nV, 1) + graph$Laplacian[[u_repl[repl_y]]]
      Q <- K
      if (alpha>1) {
        for (k in 2:alpha) {
          Q <- Q %*% K
        }
      }
      Q <- Q / reciprocal_tau^2

      # R <- chol(Q)

      R <- Matrix::Cholesky(Q)

      v <- y_resp[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y]]
      na.obs <- is.na(v)
      A.repl <- A[!na.obs, ]
      v <- v[!na.obs]
      n.o <- length(v)
      Q.p <- Q  + t(A.repl) %*% A.repl/sigma_e^2
      # R.p <- chol(Q.p)
      R.p <- Matrix::Cholesky(Q.p)
      # l <- l + sum(log(diag(R))) - sum(log(diag(R.p))) - n.o*log(sigma_e)
      l <- l + determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus - determinant(R.p, logarithm = TRUE, sqrt=TRUE)$modulus - n.o * log(sigma_e)


      if(!is.null(X_cov)){
          n_cov <- ncol(X_cov)
          if(n_cov == 0){
            X_cov_repl <- 0
          } else{
            X_cov_repl <- X_cov[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y], , drop=FALSE]
            X_cov_repl <- X_cov_repl[!na.obs, , drop = FALSE]
            v <- v - X_cov_repl %*% new_theta[4:(3+n_cov)]
          }
      }

      # mu.p <- solve(Q.p,as.vector(t(A.repl) %*% v / sigma_e^2))
      mu.p <- solve(R.p, as.vector(t(A.repl) %*% v / sigma_e^2), system = "A")
      v <- v - A.repl%*%mu.p
      l <- l - 0.5*(t(mu.p)%*%Q%*%mu.p + t(v)%*%v/sigma_e^2) - 0.5 * n.o*log(2*pi)
    }

    if(maximize){
      return(as.double(l))
    } else{
      return(-as.double(l))
    }
  }
}
