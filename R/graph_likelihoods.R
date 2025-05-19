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
  Tc <- graph$CoB$T[-ind.const, , drop = FALSE]
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

  loglik <- 0

  det_R_count <- NULL
  if(is.null(manual_y)){
    y_resp <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y_resp <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }
  n.o <- 0

  for(repl_y in 1:length(u_repl)){
    loglik <- loglik + det_R
    count <- 0
    Qpmu <- rep(0, 2*nrow(graph$E))

    ind_repl <- graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y]
    y_rep <- y_resp[ind_repl]
    if(!is.null(X_cov)){
      n_cov <- ncol(X_cov)
      if(n_cov == 0){
        X_cov_rep <- 0
      } else{
        X_cov_rep <- X_cov[graph$.__enclos_env__$private$data[[".group"]] == u_repl[repl_y], , drop=FALSE]
      }
    }

    for (e in obs.edges) {

      obs.id <- PtE[,1] == e
      y_i <- y_rep[obs.id]
      idx_na <- is.na(y_i)
      y_i <- y_i[!idx_na]

      if(sum(!idx_na) == 0){
        next
      }

      n.o <- n.o + length(y_i)

      if(!is.null(X_cov)){
        n_cov <- ncol(X_cov)
        if(n_cov == 0){
          X_cov_repl <- 0
        } else{
          X_cov_repl <- X_cov_rep[PtE[,1] == e, ,drop = FALSE]
          X_cov_repl <- X_cov_repl[!idx_na, , drop = FALSE]
          y_i <- y_i - X_cov_repl %*% theta[4:(3+n_cov)]
        }
      }

      l <- graph$edge_lengths[e]

      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]

      # Compute distance matrix
      t <- c(0, l, l*PtE_temp)
      D_matrix <- outer(t, t, `-`)

      S <- r_1(D_matrix, kappa = kappa, tau = 1/reciprocal_tau)

      #covariance update see Art p.17
      E.ind <- c(1:2)
      Obs.ind <- -E.ind

      Bt <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R <- base::chol(Sigma_i)
      Sigma_iB <- backsolve(R, forwardsolve(t(R), t(Bt)))

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

  if(is.null(manual_y)){
    y <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }

  n.o <- 0

  # Precalculate constants
  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]

  # Build Q matrix once
  Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                      alpha = 2, graph = graph, BC=BC)

  R <- Matrix::Cholesky(forceSymmetric(Tc%*%Q%*%t(Tc)),
                        LDL = FALSE, perm = TRUE)

  PtE <- graph$get_PtE()
  obs.edges <- unique(PtE[, 1])

  # Get determinant once
  loglik <- 0
  det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]
  det_R_count <- NULL

  # Cache edge lengths
  edge_lengths <- graph$edge_lengths
  n_edges <- nrow(graph$E)

  for(repl_y in seq_along(u_repl)) {
      loglik <- loglik + det_R

      # Pre-allocate with exact size needed

      Qpmu <- numeric(4 * n_edges)

      # Get data for this replicate only once
      curr_repl <- u_repl[repl_y]
      repl_indices <- (repl_vec == curr_repl)
      y_rep <- y[repl_indices]

      if(!is.null(X_cov)){
        n_cov <- ncol(X_cov)
        if(n_cov > 0){
          X_cov_rep <- X_cov[repl_indices, , drop=FALSE]
        }
      }

      # Pre-allocate with maximum size needed (16 entries per edge)
      max_entries <- 16 * length(obs.edges)
      i_ <- j_ <- x_ <- numeric(max_entries)
      count <- 0

      for (e in obs.edges) {
        obs.id <- PtE[,1] == e
        y_i <- y_rep[obs.id]

        idx_na <- is.na(y_i)
        y_i <- y_i[!idx_na]

        # Skip if no observations
        if(length(y_i) == 0) {
          next
        }

        n.o <- n.o + length(y_i)

        # Handle covariates if present
        if(!is.null(X_cov)){
          n_cov <- ncol(X_cov)
          if(n_cov > 0){
            X_cov_repl <- X_cov_rep[obs.id, , drop=FALSE]
            y_i <- y_i - as.vector(X_cov_repl %*% theta[4:(3+n_cov)])
          }
        }

        # Get edge length once
        l <- edge_lengths[e]

        # Compute distance matrices efficiently
        t <- c(0, l, l*PtE[obs.id, 2])
        D <- outer(t, t, `-`)

        # Pre-allocate matrix
        n_pts <- length(t)
        S <- matrix(0, n_pts + 2, n_pts + 2)

        # Compute all submatrices once
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

        # Cache Cholesky decomposition
        R_i <- base::chol(Sigma_i)

        Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))

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

        # Compute quadratic form directly with Cholesky
        v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))
        quad_form <- sum(y_i * v_i)

        # Update log likelihood
        loglik <- loglik - 0.5 * quad_form - sum(log(diag(R_i)))
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

  return(as.numeric(loglik))
}

#' Precompute data for the alpha2 model
#'
#' @param graph metric_graph object
#' @param data_name name of the response variable
#' @param manual_y manual y values (if data_name is NULL)
#' @param X_cov matrix of covariates
#' @param repl replicates to be considered
#' @return A list with precomputed data that can be passed to likelihood_alpha2_precompute
#' @noRd
precompute_alpha2 <- function(graph, data_name = NULL, manual_y = NULL,
                              X_cov = NULL, repl = NULL) {

  # Ensure we have alpha=2 basis construction
  if(is.null(graph$C)){
    graph$buildC(2)
  } else if(graph$CoB$alpha == 1){
    graph$buildC(2)
  }

  if(is.null(X_cov)){
    n_cov <- 0
  } else{
    n_cov <- ncol(X_cov)
  }

  # Get replication data
  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    u_repl <- unique(repl_vec)
  } else{
    u_repl <- unique(repl)
  }

  # Get response data
  if(is.null(manual_y)){
    y <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }

  # Get observation points
  PtE <- graph$get_PtE()
  obs.edges <- unique(PtE[, 1])

  # Precalculate constants
  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]

  # Cache edge lengths
  edge_lengths <- graph$edge_lengths
  n_edges <- nrow(graph$E)

  # Initialize precomputed data structure
  precomputed <- list()
  precomputed$graph <- graph
  precomputed$u_repl <- u_repl
  precomputed$obs.edges <- obs.edges
  precomputed$n_const <- n_const
  precomputed$ind.const <- ind.const
  precomputed$Tc <- Tc
  precomputed$n_edges <- n_edges
  precomputed$edge_lengths <- edge_lengths
  precomputed$n_cov <- n_cov
  
  # Precompute data for each replicate
  precomputed$y_data <- list()
  precomputed$x_data <- list()
  precomputed$D_data <- list()
  precomputed$t_data <- list()

  for(i in seq_along(u_repl)) {
    curr_repl <- u_repl[i]
    repl_indices <- (repl_vec == curr_repl)
    y_rep <- y[repl_indices]
    
    # Use character names for replicate indices
    repl_name <- paste0("repl_", curr_repl)
    
    precomputed$y_data[[repl_name]] <- list()
    precomputed$x_data[[repl_name]] <- list()
    precomputed$D_data[[repl_name]] <- list()
    precomputed$t_data[[repl_name]] <- list()

    # Process X_cov if provided
    if(!is.null(X_cov)){
      if(n_cov > 0){
        X_cov_rep <- X_cov[repl_indices, , drop=FALSE]
      }
    }

    for(j in seq_along(obs.edges)) {
      e <- obs.edges[j]
      # Use character names for edge indices
      edge_name <- paste0("edge_", e)
      
      obs.id <- PtE[,1] == e
      y_i <- y_rep[obs.id]
      idx_na <- is.na(y_i)
      y_i <- y_i[!idx_na]      

      # Skip if no observations
      if(length(y_i) == 0) {
        precomputed$y_data[[repl_name]][[edge_name]] <- NULL
        precomputed$x_data[[repl_name]][[edge_name]] <- NULL
        precomputed$D_data[[repl_name]][[edge_name]] <- NULL
        precomputed$t_data[[repl_name]][[edge_name]] <- NULL
        next
      }

      # Store y data
      precomputed$y_data[[repl_name]][[edge_name]] <- y_i

      # Store X data if present
      if(!is.null(X_cov) && ncol(X_cov) > 0){
        X_cov_e <- X_cov_rep[obs.id, , drop=FALSE]
        precomputed$x_data[[repl_name]][[edge_name]] <- X_cov_e[!idx_na, , drop=FALSE]
      }

      # Get edge length
      l <- edge_lengths[e]

      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]      

      # Compute and store time points and distance matrix
      t <- c(0, l, l*PtE_temp)
      precomputed$t_data[[repl_name]][[edge_name]] <- t

      D <- outer(t, t, `-`)
      precomputed$D_data[[repl_name]][[edge_name]] <- D
    }
  }

  return(precomputed)
}

#' Log-likelihood calculation for alpha=2 model using precomputed data
#'
#' @param theta parameters (sigma_e, reciprocal_tau, kappa)
#' @param precomputed_data precomputed data from precompute_alpha2
#' @param BC which boundary condition to use (0,1)
#' @param parameterization parameterization to be used
#' @return The log-likelihood
#' @noRd
likelihood_alpha2_precompute <- function(theta, precomputed_data, BC = 1, parameterization = "matern") {
  # Extract parameters
  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  if(parameterization == "matern"){
    kappa = sqrt(8 * 1.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }
  
  # Build Q matrix once
  Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                     alpha = 2, graph = precomputed_data$graph, BC=BC)
  
  R <- Matrix::Cholesky(forceSymmetric(precomputed_data$Tc%*%Q%*%t(precomputed_data$Tc)),
                       LDL = FALSE, perm = TRUE)
  
  # Get determinant once
  loglik <- 0
  det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]
  det_R_count <- NULL
  
  # Creating a large pre-allocated array for all potential entries
  total_max_entries <- 16 * length(precomputed_data$obs.edges)
  all_i <- all_j <- all_x <- numeric(total_max_entries)
  all_count <- 0
  
  n_obs_total <- 0
  
  # Process each replicate
  for(i in seq_along(precomputed_data$u_repl)) {
    curr_repl <- precomputed_data$u_repl[i]
    repl_name <- paste0("repl_", curr_repl)
    
    loglik <- loglik + det_R
    
    # Pre-allocate with exact size needed
    Qpmu <- numeric(4 * precomputed_data$n_edges)
    
    # Process each edge
    for(j in seq_along(precomputed_data$obs.edges)) {
      e <- precomputed_data$obs.edges[j]
      edge_name <- paste0("edge_", e)
      
      # Get data for this edge
      y_i <- precomputed_data$y_data[[repl_name]][[edge_name]]
      
      # Skip if no observations
      if(is.null(y_i) || length(y_i) == 0) {
        next
      }
      
      # Count observations for log determinant term
      n_obs_total <- n_obs_total + length(y_i)
      
      # Handle covariates if present
      if(precomputed_data$n_cov > 0) {
        X_cov_e <- precomputed_data$x_data[[repl_name]][[edge_name]]
        n_cov <- ncol(X_cov_e)
        if(n_cov > 0){
          y_i <- y_i - as.vector(X_cov_e %*% theta[4:(3+n_cov)])
        }
      }
      
      # Get precomputed distance data
      t <- precomputed_data$t_data[[repl_name]][[edge_name]]
      D <- precomputed_data$D_data[[repl_name]][[edge_name]]
      
      # Pre-allocate matrix
      n_pts <- length(t)
      S <- matrix(0, n_pts + 2, n_pts + 2)
      
      # Compute all submatrices
      d.index <- c(1,2)
      S[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                  tau = 1/reciprocal_tau, deriv = 0)
      S[d.index, d.index] <- -r_2(as.matrix(dist(c(0,precomputed_data$edge_lengths[e]))),
                                 kappa = kappa, tau = 1/reciprocal_tau,
                                 deriv = 2)
      S[d.index, -d.index] <- -r_2(D[1:2,], kappa = kappa,
                                 tau = 1/reciprocal_tau, deriv = 1)
      S[-d.index, d.index] <- t(S[d.index, -d.index])
      
      # Covariance updates
      E.ind <- c(1:4)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      
      # Cache Cholesky decomposition
      R_i <- base::chol(Sigma_i)
      
      # Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))
      Sigma_iB <- solve(R_i, forwardsolve(t(R_i), t(Bt)))
      
      BtSinvB <- Bt %*% Sigma_iB
      
      E <- precomputed_data$graph$E[e, ]
      if (E[1] == E[2]) {
        warning("Circle not implemented")
      }
      
      BtSinvB <- BtSinvB[c(3,1,4,2), c(3,1,4,2)]
      Qpmu[4 * (e - 1) + 1:4] <- Qpmu[4 * (e - 1) + 1:4] +
        (t(Sigma_iB) %*% y_i)[c(3, 1, 4, 2)]
      
      # Efficiently add precision matrix entries - use direct indexing instead of loop
      # This is more efficient than using expand.grid in a loop
      idx <- seq(all_count + 1, all_count + 16)
      
      # Lower edge u diagonal
      all_i[idx[1]] <- 4 * (e - 1) + 1
      all_j[idx[1]] <- 4 * (e - 1) + 1
      all_x[idx[1]] <- BtSinvB[1, 1]
      
      # Lower edge u' diagonal  
      all_i[idx[2]] <- 4 * (e - 1) + 2
      all_j[idx[2]] <- 4 * (e - 1) + 2
      all_x[idx[2]] <- BtSinvB[2, 2]
      
      # Upper edge u diagonal
      all_i[idx[3]] <- 4 * (e - 1) + 3
      all_j[idx[3]] <- 4 * (e - 1) + 3
      all_x[idx[3]] <- BtSinvB[3, 3]
      
      # Upper edge u' diagonal
      all_i[idx[4]] <- 4 * (e - 1) + 4
      all_j[idx[4]] <- 4 * (e - 1) + 4
      all_x[idx[4]] <- BtSinvB[4, 4]
      
      # Lower edge (u, u')
      all_i[idx[5]] <- 4 * (e - 1) + 1
      all_j[idx[5]] <- 4 * (e - 1) + 2
      all_x[idx[5]] <- BtSinvB[1, 2]
      
      all_i[idx[6]] <- 4 * (e - 1) + 2
      all_j[idx[6]] <- 4 * (e - 1) + 1
      all_x[idx[6]] <- BtSinvB[1, 2]
      
      # Upper edge (u, u')
      all_i[idx[7]] <- 4 * (e - 1) + 3
      all_j[idx[7]] <- 4 * (e - 1) + 4
      all_x[idx[7]] <- BtSinvB[3, 4]
      
      all_i[idx[8]] <- 4 * (e - 1) + 4
      all_j[idx[8]] <- 4 * (e - 1) + 3
      all_x[idx[8]] <- BtSinvB[3, 4]
      
      # Lower u, upper u
      all_i[idx[9]] <- 4 * (e - 1) + 1
      all_j[idx[9]] <- 4 * (e - 1) + 3
      all_x[idx[9]] <- BtSinvB[1, 3]
      
      all_i[idx[10]] <- 4 * (e - 1) + 3
      all_j[idx[10]] <- 4 * (e - 1) + 1
      all_x[idx[10]] <- BtSinvB[1, 3]
      
      # Lower u, upper u'
      all_i[idx[11]] <- 4 * (e - 1) + 1
      all_j[idx[11]] <- 4 * (e - 1) + 4
      all_x[idx[11]] <- BtSinvB[1, 4]
      
      all_i[idx[12]] <- 4 * (e - 1) + 4
      all_j[idx[12]] <- 4 * (e - 1) + 1
      all_x[idx[12]] <- BtSinvB[1, 4]
      
      # Lower u', upper u
      all_i[idx[13]] <- 4 * (e - 1) + 2
      all_j[idx[13]] <- 4 * (e - 1) + 3
      all_x[idx[13]] <- BtSinvB[2, 3]
      
      all_i[idx[14]] <- 4 * (e - 1) + 3
      all_j[idx[14]] <- 4 * (e - 1) + 2
      all_x[idx[14]] <- BtSinvB[2, 3]
      
      # Lower u', upper u'
      all_i[idx[15]] <- 4 * (e - 1) + 2
      all_j[idx[15]] <- 4 * (e - 1) + 4
      all_x[idx[15]] <- BtSinvB[2, 4]
      
      all_i[idx[16]] <- 4 * (e - 1) + 4
      all_j[idx[16]] <- 4 * (e - 1) + 2
      all_x[idx[16]] <- BtSinvB[2, 4]
      
      all_count <- all_count + 16
      
      # Compute quadratic form directly with Cholesky
      v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))
      quad_form <- sum(y_i * v_i)
      
      # Update log likelihood
      loglik <- loglik - 0.5 * quad_form - sum(log(diag(R_i)))
    }
  }
  
  # Build sparse matrix just once after collecting all entries
  if(all_count > 0) {
    BtSB <- Matrix::sparseMatrix(i = all_i[1:all_count],
                               j = all_j[1:all_count],
                               x = all_x[1:all_count],
                               dims = dim(Q))
    Qp <- Q + BtSB
    Qp <- precomputed_data$Tc %*% Qp %*% t(precomputed_data$Tc)
    R_count <- Matrix::Cholesky(forceSymmetric(Qp), LDL = FALSE, perm = TRUE)
    det_R_count <- Matrix::determinant(R_count, sqrt=TRUE)$modulus[1]
    
    for(i in seq_along(precomputed_data$u_repl)) {
      curr_repl <- precomputed_data$u_repl[i]
      repl_name <- paste0("repl_", curr_repl)
      
      loglik <- loglik - det_R_count
      
      v <- c(as.matrix(Matrix::solve(R_count, Matrix::solve(R_count, precomputed_data$Tc%*%Qpmu, system = 'P'),
                                     system='L')))
      
      loglik <- loglik + 0.5 * t(v) %*% v
    }
  }
  
  # Add constant term
  loglik <- loglik - 0.5 * n_obs_total * log(2 * pi)
  
  return(as.numeric(loglik))
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
        X_cov_tmp <- X_cov_tmp[!na_obs, , drop=FALSE]
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

  loglik <- 0

  det_R_count <- NULL
  if(is.null(manual_y)){
    y_resp <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y_resp <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }
  n.o <- 0

  # Cache some values used in the loop
  nV <- nrow(graph$V)

  for(repl_y in seq_along(u_repl)){
    loglik <- loglik + det_R
    count <- 0
    Qpmu <- numeric(nV)

    # Pre-compute replicate membership only once
    curr_repl <- u_repl[repl_y]
    ind_repl_curr <- (repl_vec == curr_repl)
    y_reply <- y_resp[ind_repl_curr]
    if(!is.null(X_cov)){
      n_cov <- ncol(X_cov)
      if(n_cov == 0){
        X_reply <- 0
      } else{
        X_reply <- X_cov[ind_repl_curr, , drop=FALSE]
      }
    }
    for (e in obs.edges) {
      # Use pre-computed replicate indices
      obs.id <- PtE[,1] == e

      # More efficient indexing
      y_i <- y_reply[obs.id]

      idx_na <- is.na(y_i)
      if(sum(!idx_na) == 0){
        next
      }

      y_i <- y_i[!idx_na]

      n.o <- n.o + length(y_i)      

      if(!is.null(X_cov)){
          n_cov <- ncol(X_cov)
          if(n_cov == 0){
            X_cov_repl <- 0
          } else{
            # Use pre-computed indices
            X_cov_repl <- X_reply[obs.id, , drop=FALSE]
            X_cov_repl <- X_cov_repl[!idx_na, , drop=FALSE]
            # Direct matrix multiplication
            y_i <- y_i - as.vector(X_cov_repl %*% theta[4:(3+n_cov)])
          }
      }

      l <- graph$edge_lengths[e]

      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]

      # Compute distance efficiently
      # Compute distance matrix
      t <- c(0, l, l*PtE_temp)
      D_matrix <- outer(t, t, `-`)

      # Pre-compute S matrix
      S <- r_1(D_matrix, kappa = kappa, tau = 1/reciprocal_tau)

      #covariance update see Art p.17
      E.ind <- c(1:2)
      Obs.ind <- -E.ind

      Bt <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R <- base::chol(Sigma_i)

      Sigma_iB <- backsolve(R, forwardsolve(t(R), t(Bt)))

      BtSinvB <- Bt %*% Sigma_iB

      E <- graph$E[e, ]
      if (E[1] == E[2]) {
        # Pre-compute matrix product
        y_Sigma_iB <- sum(as.vector(t(Sigma_iB) %*% y_i))
        Qpmu[E[1]] <- Qpmu[E[1]] + y_Sigma_iB
        i_[count + 1] <- E[1]
        j_[count + 1] <- E[1]
        x_[count + 1] <- sum(BtSinvB)
        count <- count + 1
      } else {
        # Pre-compute matrix product
        y_prod <- as.vector(t(Sigma_iB) %*% y_i)
        Qpmu[E] <- Qpmu[E] + y_prod

        # More efficient indexing of arrays
        idx <- count + 1:4
        i_[idx] <- c(E[1], E[1], E[2], E[2])
        j_[idx] <- c(E[1], E[2], E[1], E[2])
        x_[idx] <- c(BtSinvB[1, 1], BtSinvB[1, 2],
                           BtSinvB[1, 2], BtSinvB[2, 2])
        count <- count + 4
      }

      # Compute log determinant only once
      log_det <- sum(log(diag(R)))

      # Avoid matrix inversion by using the Cholesky factor directly
      v_i <- backsolve(R, forwardsolve(t(R), y_i))
      quad_form <- sum(y_i * v_i)

      loglik <- loglik - 0.5 * quad_form - log_det
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

#' Precompute Components for the Alpha=1 Likelihood Calculation
#'
#' This function precomputes various components needed to efficiently compute the log-likelihood
#' for an alpha=1 model using a metric graph. The function extracts observed response values, the corresponding
#' covariate matrices (if provided), and constructs distance matrices based on the edge lengths of the graph.
#'
#' @param graph A \code{metric_graph} object that contains the graph structure and observation data.
#' @param data_name Character. The name of the response variable stored in the graph data. This is used
#'   if \code{manual_y} is not provided.
#' @param manual_y Numeric vector. An optional manually provided response vector. If provided, the function uses
#'   this vector instead of the one in the \code{graph} object.
#' @param X_cov Optional matrix of covariates. If provided, the corresponding subset of covariate data is precomputed for
#'   each replicate and edge.
#' @param repl Replicate specification. If \code{NULL}, all replicates present in the graph data are used.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{y}{A nested list with the observed responses split by replicate and edge.}
#'   \item{obs.edges}{A vector of unique edge indices that have observations.}
#'   \item{D_matrix}{A list of distance matrices for each observed edge, computed from the edge lengths and the locations.}
#'   \item{x}{A nested list with subsets of the covariate matrix corresponding to each replicate and edge (if provided).}
#' }
#' @noRd

precompute_alpha1 <- function(graph,data_name = NULL, manual_y = NULL,
                              X_cov = NULL, repl){

  PtE <- graph$get_PtE()
  obs.edges <- unique(PtE[, 1])

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  if(is.null(repl)){
    u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  if(is.null(manual_y)){
    y_resp <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y_resp <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }

  # Cache some values used in the loop
  nV <- nrow(graph$V)

  precomputeddata <- list(y = list(),obs.edges=obs.edges,
                          D_matrix = list(),
                          x = list(),
                          u_repl = u_repl)

  for(j in seq_along(u_repl)){
    curr_repl <- u_repl[j]
    # Use character names for replicate indices
    repl_name <- paste0("repl_", curr_repl)

    # Pre-compute replicate membership only once
    ind_repl_curr <- (repl_vec == curr_repl)
    y_reply <- y_resp[ind_repl_curr]
    if(!is.null(X_cov)){
      n_cov <- ncol(X_cov)
      if(n_cov == 0){
        X_reply <- 0
      } else{
        X_reply <- X_cov[ind_repl_curr, , drop=FALSE]
      }
    }
    precomputeddata$y[[repl_name]] <- list()
    precomputeddata$x[[repl_name]] <- list()
    precomputeddata$D_matrix[[repl_name]] <- list()
    for (i in seq_along(obs.edges)) {
      e <- obs.edges[i]
      # Use character names for edge indices
      edge_name <- paste0("edge_", e)
      
      # Use pre-computed replicate indices
      obs.id <- PtE[,1] == e
      y_i <- y_reply[obs.id]

      idx_na <- is.na(y_i)
      if(sum(!idx_na) == 0){
        precomputeddata$y[[repl_name]][[edge_name]] <- NULL
        precomputeddata$x[[repl_name]][[edge_name]] <- NULL
        precomputeddata$D_matrix[[repl_name]][[edge_name]] <- NULL
        next
      }

      y_i <- y_i[!idx_na]
      precomputeddata$y[[repl_name]][[edge_name]] <- y_i


      if(!is.null(X_cov)){
        n_cov <- ncol(X_cov)
        if(n_cov == 0){
          precomputeddata$x[[repl_name]][[edge_name]] <- 0
        } else{
          X_cov_repl <- X_reply[obs.id, , drop=FALSE]
          precomputeddata$x[[repl_name]][[edge_name]] <- X_cov_repl[!idx_na, , drop=FALSE]
        }
      }

      l <- graph$edge_lengths[e]

      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]

      # Compute and store time points and distance matrix
      t <- c(0, l, l*PtE_temp)
      precomputeddata$D_matrix[[repl_name]][[edge_name]] <- outer(t, t, `-`)
    }
  }
  return(precomputeddata)
}

#' Log-likelihood calculation for alpha=1 model
#' @param theta (sigma_e, reciprocal_tau, kappa)
#' @param graph metric_graph object
#' @param precomputeddata precomputed data
#' @param data_name name of the response variable
#' @param repl replicates
#' @param X_cov matrix of covariates
#' @param BC. - which boundary condition to use (0,1)
#' @noRd
likelihood_alpha1_precompute <- function(theta, graph, precomputeddata ,data_name = NULL, manual_y = NULL,
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
  obs.edges <- precomputeddata$obs.edges

  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))

  if(is.null(repl)){
    u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  loglik <- 0

  det_R_count <- NULL
  n.o <- 0

  # Cache some values used in the loop
  nV <- nrow(graph$V)

  for(j in seq_along(u_repl)){
      curr_repl <- u_repl[j]
      repl_name <- paste0("repl_", curr_repl)
      
      loglik <- loglik + det_R
      count <- 0
      Qpmu <- numeric(nV)
      
    for (i in seq_along(obs.edges)) {
      # Use pre-computed replicate indices
      e <- obs.edges[i]
      edge_name <- paste0("edge_", e)
      
      # Get data for this edge and replicate
      y_i <- precomputeddata$y[[repl_name]][[edge_name]]
      
      # Skip if no observations
      if(is.null(y_i) || length(y_i) == 0){
        next
      }

      n.o <- n.o + length(y_i)

      if(!is.null(X_cov)){
        n_cov <- ncol(X_cov)
        if(n_cov == 0){
          X_cov_repl <- 0
        } else{
          y_i <- y_i - as.vector(precomputeddata$x[[repl_name]][[edge_name]] %*% theta[4:(3+n_cov)])
        }
      }
      D_matrix <- precomputeddata$D_matrix[[repl_name]][[edge_name]]

      # Pre-compute S matrix
      S <- r_1(D_matrix, kappa = kappa, tau = 1/reciprocal_tau)

      #covariance update see Art p.17
      E.ind <- c(1:2)
      Obs.ind <- -E.ind

      Bt <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R <- base::chol(Sigma_i)

      Sigma_iB <- backsolve(R, forwardsolve(t(R), t(Bt)))

      BtSinvB <- Bt %*% Sigma_iB

      E <- graph$E[e, ]
      
      if (E[1] == E[2]) {
        # Pre-compute matrix product
        y_Sigma_iB <- sum(as.vector(t(Sigma_iB) %*% y_i))
        Qpmu[E[1]] <- Qpmu[E[1]] + y_Sigma_iB
        i_[count + 1] <- E[1]
        j_[count + 1] <- E[1]
        x_[count + 1] <- sum(BtSinvB)
        count <- count + 1
      } else {
        # Pre-compute matrix product
        y_prod <- as.vector(t(Sigma_iB) %*% y_i)
        Qpmu[E] <- Qpmu[E] + y_prod

        # More efficient indexing of arrays
        idx <- count + 1:4
        i_[idx] <- c(E[1], E[1], E[2], E[2])
        j_[idx] <- c(E[1], E[2], E[1], E[2])
        x_[idx] <- c(BtSinvB[1, 1], BtSinvB[1, 2],
                     BtSinvB[1, 2], BtSinvB[2, 2])
        count <- count + 4
      }

      # Compute log determinant only once
      log_det <- sum(log(diag(R)))

      # Avoid matrix inversion by using the Cholesky factor directly
      v_i <- backsolve(R, forwardsolve(t(R), y_i))
      quad_form <- sum(y_i * v_i)

      loglik <- loglik - 0.5 * quad_form - log_det
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
#'
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

  # Pre-compute replication data
  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    u_repl <- unique(repl_vec)
  } else{
    u_repl <- unique(repl)
  }

  # Pre-compute covariate dimensions
  n_cov <- if(!is.null(X_cov)) ncol(X_cov) else 0

  loglik <- function(theta){
    # Apply parameter fixes if needed
    if(!is.null(fix_v_val)){
      fix_v_val_full <- c(fix_v_val, rep(NA, n_cov))
      fix_vec_full <- c(fix_vec, rep(FALSE, n_cov))
      new_theta <- fix_v_val_full
      new_theta[!fix_vec_full] <- theta
    } else{
      new_theta <- theta
    }

    # Extract parameters based on model type
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

    # Extract covariate parameters if needed
    theta_covariates <- if(n_cov > 0) {
      new_theta[(length(new_theta)-n_cov+1):length(new_theta)]
    } else {
      NULL
    }

    # Ensure Laplacian is computed if needed
    if(is.null(graph$Laplacian) && (model %in% c("GL1", "GL2"))) {
      graph$compute_laplacian()
    }

    # Build covariance matrix based on model type
    Sigma <- switch(model,
      WM1 = {
        # More efficient precision matrix construction
        Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                          alpha = 1, graph = graph)
        # Use sparse solver when possible
        as.matrix(Matrix::solve(Q))[graph$PtV, graph$PtV]
      },
      WM2 = {
        PtE <- graph$get_PtE()
        n.c <- 1:length(graph$CoB$S)

        # Build precision matrix
        Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau, alpha = 2,
                          graph = graph, BC = 1)

        # Cache matrix products
        TC <- graph$CoB$T
        TCt <- t(TC)
        Qtilde <- TC %*% Q %*% TCt
        Qtilde <- Qtilde[-n.c,-n.c]

        # Use Cholesky for solving when possible
        QChol <- Matrix::Cholesky(Qtilde, LDL = FALSE)
        TC_nc <- TC[-n.c,]

        # Compute using intermediate matrices
        Sigma.overdetermined <- TCt[, -n.c] %*% Matrix::solve(QChol, TC_nc, system = "A")

        # More efficient indexing
        index.obs <- 4 * (PtE[,1] - 1) +
                     1.0 * (abs(PtE[, 2]) < 1e-14) +
                     3.0 * (abs(PtE[, 2]) > 1e-14)

        as.matrix(Sigma.overdetermined[index.obs, index.obs])
      },
      GL1 = {
        # Cache intermediate calculations
        K <- kappa^2 * Matrix::Diagonal(graph$nV, 1)
        L <- graph$Laplacian[[1]]
        Q <- (K + L) / reciprocal_tau^2

        # Use sparse solver
        as.matrix(Matrix::solve(Q))[graph$PtV, graph$PtV]
      },
      GL2 = {
        # Cache intermediate calculations
        K <- kappa^2 * Matrix::Diagonal(graph$nV, 1)
        L <- graph$Laplacian[[1]]
        Q <- (K + L)
        Q <- Q %*% Q / reciprocal_tau^2

        # Use sparse solver
        as.matrix(Matrix::solve(Q))[graph$PtV, graph$PtV]
      },
      isoCov = {
        if(is.null(cov_function)){
          stop("If model is 'isoCov' the covariance function must be supplied!")
        }
        if(!is.function(cov_function)){
          stop("'cov_function' must be a function!")
        }

        # Compute resistance distance if needed
        if(is.null(graph$res_dist)){
          graph$compute_resdist(full = TRUE, check_euclidean = check_euclidean)
        }

        # Compute covariance matrix
        as.matrix(cov_function(as.matrix(graph$res_dist[[".complete"]]), theta_cov))
      }
    )

    # Add measurement error to diagonal
    diag(Sigma) <- diag(Sigma) + sigma_e^2

    # Initialize log-likelihood
    loglik_val <- 0

    # Process each replicate
    for(repl_y in seq_along(u_repl)){
      # Get data for this replicate
      curr_repl <- u_repl[repl_y]
      ind_tmp <- (repl_vec %in% curr_repl)
      y_tmp <- y_graph[ind_tmp]
      na_obs <- is.na(y_tmp)

      # Skip if all observations are NA
      if(all(na_obs)) next

      # Extract observation data
      Sigma_non_na <- Sigma[!na_obs, !na_obs]

      # Compute Cholesky decomposition once
      R <- base::chol(Sigma_non_na)

      # Get data vector
      v <- y_graph[repl_vec == curr_repl]

      # Apply covariate adjustment if needed
      if(!is.null(X_cov) && n_cov > 0){
        X_cov_repl <- X_cov[repl_vec == curr_repl, , drop=FALSE]
        v <- v - X_cov_repl %*% theta_covariates
      } else{
        X_cov_repl <- 0
      }

      # Keep only non-NA observations
      v <- v[!na_obs]

      # Count observations for log determinant term
      n_obs <- length(v)

      # Compute log determinant
      log_det <- sum(log(diag(R)))

      # Solve system efficiently using Cholesky factor
      v_i <- backsolve(R, forwardsolve(t(R), v))
      quad_form <- sum(v * v_i)

      # Update log-likelihood efficiently
      loglik_val <- loglik_val - log_det - 0.5 * quad_form - 0.5 * n_obs * log(2*pi)
    }

    # Return based on maximize flag
    if(maximize){
      return(as.double(loglik_val))
    } else{
      return(-as.double(loglik_val))
    }
  }

  return(loglik)
}

#' Precompute data for the graph covariance model
#'
#' @param graph A `metric_graph` object
#' @param model Type of model: "WM1", "WM2", "GL1", "GL2", or "isoCov"
#' @param y_graph Response vector given in the same order as the internal locations from the graph
#' @param X_cov Matrix with covariates
#' @param repl Vector with the replicates to be considered
#' @param check_euclidean Whether to check if the graph is Euclidean
#' @return A list with precomputed data that can be passed to likelihood_graph_covariance_precompute
#' @noRd
precompute_graph_covariance <- function(graph,
                                        model = "WM1",
                                        y_graph,
                                        X_cov = NULL,
                                        repl = NULL,
                                        check_euclidean = TRUE) {

  # Validate model type
  if(!(model %in% c("WM1", "WM2", "GL1", "GL2", "isoCov"))) {
    stop("The available models are: 'WM1', 'WM2', 'GL1', 'GL2' and 'isoCov'!")
  }

  # Get replication data
  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)) {
    u_repl <- unique(repl_vec)
  } else {
    u_repl <- unique(repl)
  }

  # Prepare data structures
  precomputed <- list()
  precomputed$model <- model
  precomputed$u_repl <- u_repl
  precomputed$y_data <- list()
  precomputed$X_data <- list()
  precomputed$na_indices <- list()

  # Pre-compute covariate dimensions
  precomputed$n_cov <- if(!is.null(X_cov)) ncol(X_cov) else 0

  # Ensure Laplacian is computed if needed for GL1/GL2 models
  if(is.null(graph$Laplacian) && (model %in% c("GL1", "GL2"))) {
    graph$compute_laplacian()
  }

  # For isoCov model, precompute resistance distances
  if(model == "isoCov" && is.null(graph$res_dist)) {
    graph$compute_resdist(full = TRUE, check_euclidean = check_euclidean)
  }

  # Precompute PtV and PtE for WM models
  if(model %in% c("WM1", "WM2")) {
    precomputed$PtV <- graph$PtV
    if(model == "WM2") {
      precomputed$PtE <- graph$get_PtE()
      precomputed$n_constraints <- 1:length(graph$CoB$S)
    }
  }

  # Precompute data for each replicate
  for(i in seq_along(u_repl)) {
    curr_repl <- u_repl[i]
    # Use character names for replicate indices
    repl_name <- paste0("repl_", curr_repl)

    # Get data for this replicate
    ind_tmp <- (repl_vec %in% curr_repl)
    y_tmp <- y_graph[ind_tmp]
    na_obs <- is.na(y_tmp)

    # Store observation mask and non-NA values
    precomputed$na_indices[[repl_name]] <- na_obs
    precomputed$y_data[[repl_name]] <- y_tmp[!na_obs]

    # Store covariate data if present
    if(!is.null(X_cov) && precomputed$n_cov > 0) {
      X_cov_repl <- X_cov[ind_tmp, , drop=FALSE]
      precomputed$X_data[[repl_name]] <- X_cov_repl[!na_obs, , drop=FALSE]
    }
  }

  # Store graph information needed for calculations
  precomputed$nV <- graph$nV

  # Specific model information
  if(model == "GL1" || model == "GL2") {
    precomputed$Laplacian <- graph$Laplacian
  }

  # Store graph components needed for all models
  precomputed$graph <- graph

  return(precomputed)
}

#' Log-likelihood calculation using precomputed data for the graph covariance model
#'
#' @param theta Parameter vector
#' @param precomputed_data Precomputed data from precompute_graph_covariance
#' @param cov_function Covariance function (required for isoCov model)
#' @param log_scale Whether parameters are in log scale
#' @param maximize Whether to return the likelihood (TRUE) or negative likelihood (FALSE)
#' @param fix_vec Vector indicating which parameters are fixed
#' @param fix_v_val Values for fixed parameters
#' @return The log-likelihood
#' @noRd
likelihood_graph_covariance_precompute <- function(theta,
                                                  precomputed_data,
                                                  cov_function = NULL,
                                                  log_scale = TRUE,
                                                  maximize = FALSE) {

  # Extract parameters based on model type
  if(precomputed_data$model == "isoCov"){
    if(log_scale){
      sigma_e <- exp(theta[1])
      theta_cov <- exp(theta[2:(length(theta)-precomputed_data$n_cov)])
    } else{
      sigma_e <- theta[1]
      theta_cov <- theta[2:(length(theta)-precomputed_data$n_cov)]
    }
  } else{
    if(log_scale){
      sigma_e <- exp(theta[1])
      reciprocal_tau <- exp(theta[2])
      kappa <- exp(theta[3])
    } else{
      sigma_e <- theta[1]
      reciprocal_tau <- theta[2]
      kappa <- theta[3]
    }
  }

  # Extract covariate parameters if needed
  theta_covariates <- if(precomputed_data$n_cov > 0) {
    theta[(length(theta)-precomputed_data$n_cov+1):length(theta)]
  } else {
    NULL
  }

  # Build covariance matrix based on model type
  Sigma <- switch(precomputed_data$model,
    WM1 = {
      # More efficient precision matrix construction
      Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                        alpha = 1, graph = precomputed_data$graph)
      # Use sparse solver when possible
      as.matrix(Matrix::solve(Q))[precomputed_data$PtV, precomputed_data$PtV]
    },
    WM2 = {
      PtE <- precomputed_data$PtE
      n.c <- precomputed_data$n_constraints

      # Build precision matrix
      Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau, alpha = 2,
                        graph = precomputed_data$graph, BC = 1)

      # Cache matrix products
      TC <- precomputed_data$graph$CoB$T
      TCt <- t(TC)
      Qtilde <- TC %*% Q %*% TCt
      Qtilde <- Qtilde[-n.c,-n.c]

      # Use Cholesky for solving when possible
      QChol <- Matrix::Cholesky(Qtilde, LDL = FALSE)
      TC_nc <- TC[-n.c,]

      # Compute using intermediate matrices
      Sigma.overdetermined <- TCt[, -n.c] %*% Matrix::solve(QChol, TC_nc, system = "A")

      # More efficient indexing
      index.obs <- 4 * (PtE[,1] - 1) +
                   1.0 * (abs(PtE[, 2]) < 1e-14) +
                   3.0 * (abs(PtE[, 2]) > 1e-14)

      as.matrix(Sigma.overdetermined[index.obs, index.obs])
    },
    GL1 = {
      # Cache intermediate calculations
      K <- kappa^2 * Matrix::Diagonal(precomputed_data$nV, 1)
      L <- precomputed_data$Laplacian[[1]]
      Q <- (K + L) / reciprocal_tau^2

      # Use sparse solver
      as.matrix(Matrix::solve(Q))[precomputed_data$PtV, precomputed_data$PtV]
    },
    GL2 = {
      # Cache intermediate calculations
      K <- kappa^2 * Matrix::Diagonal(precomputed_data$nV, 1)
      L <- precomputed_data$Laplacian[[1]]
      Q <- (K + L)
      Q <- Q %*% Q / reciprocal_tau^2

      # Use sparse solver
      as.matrix(Matrix::solve(Q))[precomputed_data$PtV, precomputed_data$PtV]
    },
    isoCov = {
      if(is.null(cov_function)){
        stop("If model is 'isoCov' the covariance function must be supplied!")
      }
      if(!is.function(cov_function)){
        stop("'cov_function' must be a function!")
      }

      # Compute covariance matrix
      as.matrix(cov_function(as.matrix(precomputed_data$graph$res_dist[[".complete"]]), theta_cov))
    }
  )

  # Add measurement error to diagonal
  diag(Sigma) <- diag(Sigma) + sigma_e^2

  # Initialize log-likelihood
  loglik_val <- 0

  # Process each replicate using precomputed data
  for(i in seq_along(precomputed_data$u_repl)) {
    curr_repl <- precomputed_data$u_repl[i]
    repl_name <- paste0("repl_", curr_repl)
    
    # Get precomputed data for this replicate
    na_obs <- precomputed_data$na_indices[[repl_name]]

    # Skip if all observations are NA
    if(all(na_obs)) next

    # Extract valid observations
    y_i <- precomputed_data$y_data[[repl_name]]

    # Extract observation data - only use the non-NA rows and columns
    Sigma_non_na <- Sigma[!na_obs, !na_obs]

    # Compute Cholesky decomposition once
    R <- base::chol(Sigma_non_na)

    # Get data vector with covariate adjustment if needed
    v <- y_i
    if(!is.null(precomputed_data$X_data) && precomputed_data$n_cov > 0) {
      X_cov_repl <- precomputed_data$X_data[[repl_name]]
      v <- v - X_cov_repl %*% theta_covariates
    }

    # Count observations for log determinant term
    n_obs <- length(v)

    # Compute log determinant
    log_det <- sum(log(diag(R)))

    # Solve system efficiently using Cholesky factor
    v_i <- backsolve(R, forwardsolve(t(R), v))
    quad_form <- sum(v * v_i)

    # Update log-likelihood efficiently
    loglik_val <- loglik_val - log_det - 0.5 * quad_form - 0.5 * n_obs * log(2*pi)
  }

  # Return based on maximize flag
  if(maximize) {
    return(as.double(loglik_val))
  } else {
    return(-as.double(loglik_val))
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

#' Precompute data for the alpha1_directional model
#'
#' @param graph metric_graph object
#' @param data_name name of the response variable
#' @param manual_y manual y values (if data_name is NULL) 
#' @param X_cov matrix of covariates
#' @param repl replicates to be considered
#' @return A list with precomputed data that can be passed to likelihood_alpha1_directional_precompute
#' @noRd
precompute_alpha1_directional <- function(graph, data_name = NULL, manual_y = NULL,
                                         X_cov = NULL, repl = NULL) {
  
  # Ensure directional constraints are built
  if(is.null(graph$C)){
    graph$buildDirectionalConstraints(alpha = 1)
  } else if(graph$CoB$alpha == 2){
    graph$buildDirectionalConstraints(alpha = 1)
  }

  if(is.null(X_cov)){
    n_cov <- 0
  } else{
    n_cov <- ncol(X_cov)
  }
  
  # Get replication data
  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]
  
  if(is.null(repl)){
    u_repl <- unique(repl_vec)
  } else{
    u_repl <- unique(repl)
  }
  
  # Get observation data
  if(is.null(manual_y)){
    y_resp <- graph$.__enclos_env__$private$data[[data_name]]
  } else if(is.null(data_name)){
    y_resp <- manual_y
  } else{
    stop("Either data_name or manual_y must be not NULL")
  }
  
  # Get observation points
  PtE <- graph$get_PtE()
  obs.edges <- unique(PtE[, 1])
  
  # Initialize precomputed data structure
  precomputed <- list()
  precomputed$graph <- graph
  precomputed$u_repl <- u_repl
  precomputed$obs.edges <- obs.edges
  precomputed$n_cov <- n_cov

  # Precalculate constants
  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, , drop = FALSE]
  
  precomputed$n_const <- n_const
  precomputed$ind.const <- ind.const
  precomputed$Tc <- Tc
  
  # Get dimensions
  precomputed$n_edges <- nrow(graph$E)
  precomputed$edge_lengths <- graph$edge_lengths
  
  # Precompute observation data for each replicate
  precomputed$y_data <- list()
  precomputed$x_data <- list()
  precomputed$D_data <- list()
  precomputed$PtE_data <- list()
  precomputed$no_na_indices <- list() # Store indices of non-NA observations
  
  for(i in seq_along(u_repl)) {
    curr_repl <- u_repl[i]
    # Use character names for replicate indices
    repl_name <- paste0("repl_", curr_repl)
    
    ind_repl <- (repl_vec == curr_repl)
    y_rep <- y_resp[ind_repl]
    
    precomputed$y_data[[repl_name]] <- list()
    precomputed$x_data[[repl_name]] <- list()
    precomputed$D_data[[repl_name]] <- list()
    precomputed$PtE_data[[repl_name]] <- list()
    precomputed$no_na_indices[[repl_name]] <- list()
    
    # Process X_cov if provided
    if(!is.null(X_cov)){
      if(n_cov == 0){
        X_reply <- 0
      } else{
        X_reply <- X_cov[ind_repl, , drop=FALSE]
      }
    }
    
    for(j in seq_along(obs.edges)) {
      e <- obs.edges[j]
      # Use character names for edge indices
      edge_name <- paste0("edge_", e)
      
      obs.id <- PtE[,1] == e
      y_i <- y_rep[obs.id]
      
      idx_na <- is.na(y_i)
      
      # Skip if all observations are NA
      if(sum(!idx_na) == 0){
        precomputed$y_data[[repl_name]][[edge_name]] <- NULL
        precomputed$x_data[[repl_name]][[edge_name]] <- NULL
        precomputed$D_data[[repl_name]][[edge_name]] <- NULL
        precomputed$PtE_data[[repl_name]][[edge_name]] <- NULL
        precomputed$no_na_indices[[repl_name]][[edge_name]] <- NULL
        next
      }
      
      # Store non-NA observations
      precomputed$y_data[[repl_name]][[edge_name]] <- y_i[!idx_na]
      precomputed$no_na_indices[[repl_name]][[edge_name]] <- !idx_na
      
      # Store covariate data if present
      if(!is.null(X_cov) && ncol(X_cov) > 0){
        X_cov_e <- X_reply[obs.id, , drop=FALSE]
        precomputed$x_data[[repl_name]][[edge_name]] <- X_cov_e[!idx_na, , drop=FALSE]
      }
      
      # Get observation locations
      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]
      precomputed$PtE_data[[repl_name]][[edge_name]] <- PtE_temp
      
      # Get edge length
      l <- graph$edge_lengths[e]
      
      # Compute distance matrix exactly as in the original function
      t <- c(0, l, l*PtE_temp)
      D <- outer(t, t, `-`)
      precomputed$D_data[[repl_name]][[edge_name]] <- D
    }
  }
  
  return(precomputed)
}

#' Log-likelihood calculation for alpha=1 directional model using precomputed data
#'
#' @param theta parameters (sigma_e, reciprocal_tau, kappa)
#' @param precomputed_data precomputed data from precompute_alpha1_directional
#' @param parameterization parameterization to be used ("matern" or other)
#' @param maximize whether to return the negative log-likelihood
#' @return The log-likelihood
#' @noRd
likelihood_alpha1_directional_precompute <- function(theta,
                                                    precomputed_data,
                                                    parameterization = "matern",
                                                    maximize = FALSE) {
  # Extract parameters
  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  if(parameterization == "matern"){
    kappa = sqrt(8 * 0.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }
  
  # Build Q matrix
  graph <- precomputed_data$graph
  Q.list <- Qalpha1_edges(c(1/reciprocal_tau, kappa),
                          graph,
                          w = 0,
                          BC = 1, 
                          build = FALSE)
  
  Q <- Matrix::sparseMatrix(i = Q.list$i,
                           j = Q.list$j,
                           x = Q.list$x,
                           dims = Q.list$dims)
  
  # Use transformation matrix
  Tc <- precomputed_data$Tc
  Q_transformed <- forceSymmetric(Tc%*%Q%*%t(Tc))
  R <- Matrix::Cholesky(Q_transformed, LDL = FALSE, perm = TRUE)
  det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]
  
  # Initialize log-likelihood
  loglik <- 0
  
  # Pre-allocate arrays for sparse matrix construction
  i_ <- j_ <- x_ <- rep(0, 4 * length(precomputed_data$obs.edges))
  
  # Initialize variables for the precision matrix
  det_R_count <- NULL
  n.o <- 0  # Total observation counter, accumulates across all replicates
  
  # Process each replicate for calculating the likelihood
  for(repl_y in 1:length(precomputed_data$u_repl)){
    curr_repl <- precomputed_data$u_repl[repl_y]
    repl_name <- paste0("repl_", curr_repl)
    
    loglik <- loglik + det_R
    count <- 0  # Reset count for each replicate
    Qpmu <- rep(0, 2*precomputed_data$n_edges)  # Use rep() for consistency
    
    # Process each edge
    for(j in seq_along(precomputed_data$obs.edges)) {
      e <- precomputed_data$obs.edges[j]
      edge_name <- paste0("edge_", e)
      
      # Skip if no observations
      if(is.null(precomputed_data$y_data[[repl_name]][[edge_name]])) {
        next
      }
      
      # Get data for this edge
      y_i <- precomputed_data$y_data[[repl_name]][[edge_name]]
      
      # Count observations
      n.o <- n.o + length(y_i)
      
      # Apply covariate adjustment if needed
      if(precomputed_data$n_cov > 0) {
        X_cov_e <- precomputed_data$x_data[[repl_name]][[edge_name]]
        n_cov <- ncol(X_cov_e)
        if(n_cov > 0){
          y_i <- y_i - X_cov_e %*% theta[4:(3+n_cov)]
        }
      }
      
      # Get precomputed distance matrix
      D_matrix <- precomputed_data$D_data[[repl_name]][[edge_name]]
      
      # Compute covariance function
      S <- r_1(D_matrix, kappa = kappa, tau = 1/reciprocal_tau)
      
      # Covariance updates
      E.ind <- c(1:2)
      Obs.ind <- -E.ind
      
      Bt <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] - 
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R_i <- base::chol(Sigma_i)
      
      Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))
      BtSinvB <- Bt %*% Sigma_iB
      
      E <- graph$E[e, ]
      if(E[1] == E[2]) {
        # Handle self-loop
        Qpmu[2*(e-1)+1] <- Qpmu[2*(e-1)+1] + sum(t(Sigma_iB) %*% y_i)
        i_[count + 1] <- 2*(e-1)+1
        j_[count + 1] <- 2*(e-1)+1
        x_[count + 1] <- sum(Bt %*% Sigma_iB)  # Use Bt %*% Sigma_iB here, not BtSinvB
        count <- count + 1  # Increment count for self-loops
      } else {
        # Handle regular edge
        Qpmu[2*(e-1) + c(1, 2)] <- Qpmu[2*(e-1) + c(1, 2)] + t(Sigma_iB) %*% y_i
        i_[count + (1:4)] <- c(2*(e-1)+1, 2*(e-1)+1, 2*(e-1)+2, 2*(e-1)+2)
        j_[count + (1:4)] <- c(2*(e-1)+1, 2*(e-1)+2, 2*(e-1)+1, 2*(e-1)+2)
        x_[count + (1:4)] <- c(BtSinvB[1, 1], BtSinvB[1, 2],
                               BtSinvB[1, 2], BtSinvB[2, 2])
        count <- count + 4
      }
      
      # Update log likelihood with quadratic term
      loglik <- loglik - 0.5 * t(y_i) %*% solve(Sigma_i, y_i)
      loglik <- loglik - sum(log(diag(R_i)))
    }
    
    # Only compute the Cholesky decomposition once - on the first replicate
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
    
    # Complete the likelihood calculation for this replicate
    loglik <- loglik - det_R_count
    
    v <- c(as.matrix(Matrix::solve(R_count, Matrix::solve(R_count, Tc%*%Qpmu,
                                                        system = "P"),
                                  system = "L")))
    
    loglik <- loglik + 0.5 * t(v) %*% v - 0.5 * n.o * log(2*pi)
  }
  
  result <- loglik[1]  # Get first element
  
  if(maximize) {
    return(result)
  } else {
    return(-result)
  }
}
