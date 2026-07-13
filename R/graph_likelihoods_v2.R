# graph_likelihoods_v2.R
#
# Likelihood evaluation with the fixed effects (beta) profiled out
# analytically, for the bridge/edge representation likelihoods in
# graph_likelihoods.R:
#
#   * alpha = 1               (likelihood_alpha1 / _precompute)
#   * alpha = 2               (likelihood_alpha2 / _precompute)
#   * alpha = 1, directional  (likelihood_alpha1_directional / _precompute)
#
# For fixed covariance parameters theta = (log sigma_e, log 1/tau, log kappa)
# the log-likelihood is quadratic in beta, so beta can be maximized out in
# closed form (GLS). With
#
#   H = sum_e X_e' Sigma_e^-1 X_e - G' Qtilde^-1 G,
#   h = sum_e X_e' Sigma_e^-1 y_e - G' Qtilde^-1 g,
#
# where g/G are the projections of y/X onto the endpoint degrees of freedom,
# the profile estimate is beta_hat(theta) = H^-1 h and
#
#   2 l_p(theta)  = 2 l_0(theta) + h' H^-1 h,
#   2 l_R(theta)  = 2 l_p(theta) - log|H|            (REML, up to a constant),
#
# where l_0 is the corresponding no-covariate likelihood evaluated with the
# raw y. See jonas_local/profile_likelihood_reml.md for the derivation, the
# math-to-code map, and the conventions that are kept identical to the
# originals (per-replicate constants, shared Qtilde factorization, etc.) so
# that l_p(theta) == l(theta, beta_hat(theta)) holds to machine precision.
#
# The existing precompute_alpha1 / precompute_alpha2 /
# precompute_alpha1_directional outputs are consumed as-is; the non-precompute
# variants build the same cache internally on each call.
#
# All functions return the log-likelihood (no sign flip).


#' Check that theta does not contain fixed-effect entries
#' @param theta parameter vector
#' @noRd
check_theta_profile <- function(theta) {
  if (length(theta) != 3) {
    stop("theta must have length 3: (log sigma_e, log reciprocal_tau, log kappa/range). ",
         "beta is profiled out and must not be part of theta.")
  }
}

#' Solve the profiled fixed-effect system via Cholesky
#' @param H p x p information matrix of beta given theta
#' @param h p vector
#' @return list with beta = H^-1 h, vcov = H^-1, quad = h' H^-1 h,
#'   logdetH = log|H|
#' @noRd
profile_beta_solve <- function(H, h) {
  cH <- base::chol(H)
  z <- forwardsolve(t(cH), h)
  return(list(beta = as.vector(backsolve(cH, z)),
              vcov = chol2inv(cH),
              quad = sum(z^2),
              logdetH = 2 * sum(log(diag(cH)))))
}

#' Assemble profile (or REML) log-likelihood from core quantities
#' @param core list with elements base, H, h, n_cov (output of a
#'   profile_lik_core_* function)
#' @param reml if TRUE return the restricted likelihood
#'   2 l_R = 2 l_p - log|H| (no-op when n_cov = 0)
#' @noRd
profile_lik_finish <- function(core, reml) {
  loglik <- as.numeric(core$base)
  if (core$n_cov > 0) {
    bs <- profile_beta_solve(core$H, core$h)
    loglik <- loglik + 0.5 * bs$quad
    if (reml) {
      loglik <- loglik - 0.5 * bs$logdetH
    }
  }
  return(loglik)
}


#' Core quantities for the profiled alpha=1 likelihood
#'
#' Follows likelihood_alpha1_precompute step by step, but evaluates the
#' no-covariate likelihood with the raw y (base) and additionally accumulates
#' H and h for the profiled fixed effects.
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range)
#' @param graph metric_graph object
#' @param precomputeddata output of precompute_alpha1
#' @param BC boundary conditions (0, 1)
#' @param parameterization "matern" or "spde"
#' @return list(base, H, h, n_cov)
#' @noRd
profile_lik_core_alpha1 <- function(theta, graph, precomputeddata, BC,
                                    parameterization) {
  sigma_e <- exp(theta[1])

  if (parameterization == "matern") {
    kappa <- sqrt(8 * 0.5) / exp(theta[3])
  } else {
    kappa <- exp(theta[3])
  }

  reciprocal_tau <- exp(theta[2])

  Q.list <- spde_precision(kappa = kappa, tau = 1 / reciprocal_tau, alpha = 1,
                           graph = graph, build = FALSE, BC = BC)

  Qp <- Matrix::sparseMatrix(i = Q.list$i,
                             j = Q.list$j,
                             x = Q.list$x,
                             dims = Q.list$dims)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  det_R <- Matrix::determinant(R, sqrt = TRUE)$modulus[1]

  obs.edges <- precomputeddata$obs.edges
  u_repl <- precomputeddata$u_repl
  nV <- if (!is.null(precomputeddata$nV)) precomputeddata$nV else nrow(graph$V)
  n_cov <- precomputeddata$n_cov

  base <- 0
  det_R_count <- NULL
  R_count <- NULL
  n.o <- 0

  H <- if (n_cov > 0) matrix(0, n_cov, n_cov) else NULL
  h <- if (n_cov > 0) numeric(n_cov) else NULL

  for (j in seq_along(u_repl)) {
    base <- base + det_R
    count <- 0
    i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))
    Qpmu <- numeric(nV)
    if (n_cov > 0) {
      QpmuX <- matrix(0, nV, n_cov)
      XtSX <- matrix(0, n_cov, n_cov)
      XtSy <- numeric(n_cov)
    }

    edge_cache_j <- precomputeddata$edge_cache[[j]]
    edge_e <- edge_cache_j$e
    edge_E <- edge_cache_j$E
    edge_y <- edge_cache_j$y
    edge_X <- edge_cache_j$X
    edge_D <- edge_cache_j$D

    for (i in seq_along(edge_e)) {
      y_i <- edge_y[[i]]
      n.o <- n.o + length(y_i)
      X_i <- if (n_cov > 0) edge_X[[i]] else NULL
      D_matrix <- edge_D[[i]]

      S <- r_1(D_matrix, kappa = kappa, tau = 1 / reciprocal_tau)

      # covariance update, see likelihood_alpha1
      E.ind <- c(1:2)
      Obs.ind <- -E.ind

      Bt <- solve(S[E.ind, E.ind, drop = FALSE],
                  S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R_i <- base::chol(Sigma_i)

      Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))
      BtSinvB <- Bt %*% Sigma_iB

      E <- edge_E[i, ]

      if (E[1] == E[2]) {
        Qpmu[E[1]] <- Qpmu[E[1]] + sum(as.vector(t(Sigma_iB) %*% y_i))
        if (n_cov > 0) {
          QpmuX[E[1], ] <- QpmuX[E[1], ] + colSums(t(Sigma_iB) %*% X_i)
        }
        i_[count + 1] <- E[1]
        j_[count + 1] <- E[1]
        x_[count + 1] <- sum(BtSinvB)
        count <- count + 1
      } else {
        Qpmu[E] <- Qpmu[E] + as.vector(t(Sigma_iB) %*% y_i)
        if (n_cov > 0) {
          QpmuX[E, ] <- QpmuX[E, , drop = FALSE] + t(Sigma_iB) %*% X_i
        }
        idx <- count + 1:4
        i_[idx] <- c(E[1], E[1], E[2], E[2])
        j_[idx] <- c(E[1], E[2], E[1], E[2])
        x_[idx] <- c(BtSinvB[1, 1], BtSinvB[1, 2],
                     BtSinvB[1, 2], BtSinvB[2, 2])
        count <- count + 4
      }

      v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))
      base <- base - 0.5 * sum(y_i * v_i) - sum(log(diag(R_i)))

      if (n_cov > 0) {
        SinvX <- backsolve(R_i, forwardsolve(t(R_i), X_i))
        XtSX <- XtSX + crossprod(X_i, SinvX)
        XtSy <- XtSy + as.vector(crossprod(X_i, v_i))
      }
    }

    if (is.null(det_R_count)) {
      i_ <- c(Q.list$i, i_[1:count])
      j_ <- c(Q.list$j, j_[1:count])
      x_ <- c(Q.list$x, x_[1:count])

      Qp <- Matrix::sparseMatrix(i = i_,
                                 j = j_,
                                 x = x_,
                                 dims = Q.list$dims)

      R_count <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
      det_R_count <- Matrix::determinant(R_count, sqrt = TRUE)$modulus[1]
    }

    base <- base - det_R_count

    v <- c(as.matrix(Matrix::solve(R_count,
                                   Matrix::solve(R_count, Qpmu, system = "P"),
                                   system = "L")))

    base <- base + 0.5 * sum(v^2) - 0.5 * n.o * log(2 * pi)

    if (n_cov > 0) {
      VX <- as.matrix(Matrix::solve(R_count,
                                    Matrix::solve(R_count, QpmuX, system = "P"),
                                    system = "L"))
      H <- H + XtSX - crossprod(VX)
      h <- h + XtSy - as.vector(crossprod(VX, v))
    }
  }

  return(list(base = base, H = H, h = h, n_cov = n_cov))
}


#' Core quantities for the profiled alpha=2 likelihood
#'
#' Follows likelihood_alpha2_precompute step by step (constrained basis Tc,
#' 4 endpoint dofs per edge, reorder c(3,1,4,2)), evaluating the no-covariate
#' likelihood with the raw y and accumulating H and h.
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range)
#' @param precomputed_data output of precompute_alpha2
#' @param BC boundary conditions (0, 1)
#' @param parameterization "matern" or "spde"
#' @return list(base, H, h, n_cov)
#' @noRd
profile_lik_core_alpha2 <- function(theta, precomputed_data, BC,
                                    parameterization) {
  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  if (parameterization == "matern") {
    kappa <- sqrt(8 * 1.5) / exp(theta[3])
  } else {
    kappa <- exp(theta[3])
  }

  graph <- precomputed_data$graph
  Tc <- precomputed_data$Tc
  n_edges <- precomputed_data$n_edges
  n_cov <- precomputed_data$n_cov
  u_repl <- precomputed_data$u_repl

  Q <- spde_precision(kappa = kappa, tau = 1 / reciprocal_tau,
                      alpha = 2, graph = graph, BC = BC)

  R <- Matrix::Cholesky(forceSymmetric(Tc %*% Q %*% t(Tc)),
                        LDL = FALSE, perm = TRUE)

  det_R <- Matrix::determinant(R, sqrt = TRUE)$modulus[1]
  det_R_count <- NULL
  R_count <- NULL

  base <- 0
  n.o <- 0

  H <- if (n_cov > 0) matrix(0, n_cov, n_cov) else NULL
  h <- if (n_cov > 0) numeric(n_cov) else NULL

  for (i in seq_along(u_repl)) {
    base <- base + det_R

    Qpmu <- numeric(4 * n_edges)
    i_ <- j_ <- x_ <- numeric(16 * length(precomputed_data$obs.edges))
    count <- 0
    if (n_cov > 0) {
      QpmuX <- matrix(0, 4 * n_edges, n_cov)
      XtSX <- matrix(0, n_cov, n_cov)
      XtSy <- numeric(n_cov)
    }

    edge_cache_i <- precomputed_data$edge_cache[[i]]
    edge_e <- edge_cache_i$e
    edge_E <- edge_cache_i$E
    edge_l <- edge_cache_i$l
    edge_y <- edge_cache_i$y
    edge_X <- edge_cache_i$X
    edge_t <- edge_cache_i$t
    edge_D <- edge_cache_i$D

    for (j in seq_along(edge_e)) {
      e <- edge_e[j]
      y_i <- edge_y[[j]]
      n.o <- n.o + length(y_i)
      X_i <- if (n_cov > 0) edge_X[[j]] else NULL

      t_pts <- edge_t[[j]]
      D <- edge_D[[j]]

      n_pts <- length(t_pts)
      S <- matrix(0, n_pts + 2, n_pts + 2)

      d.index <- c(1, 2)
      S[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                   tau = 1 / reciprocal_tau, deriv = 0)
      S[d.index, d.index] <- -r_2(matrix(c(0, -edge_l[j], edge_l[j], 0), 2, 2),
                                  kappa = kappa, tau = 1 / reciprocal_tau,
                                  deriv = 2)
      S[d.index, -d.index] <- -r_2(D[1:2, ], kappa = kappa,
                                   tau = 1 / reciprocal_tau, deriv = 1)
      S[-d.index, d.index] <- t(S[d.index, -d.index])

      E.ind <- c(1:4)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind] - S[Obs.ind, E.ind] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

      R_i <- base::chol(Sigma_i)

      Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))
      BtSinvB <- Bt %*% Sigma_iB

      E <- edge_E[j, ]
      if (E[1] == E[2]) {
        warning("Circle not implemented")
      }

      # reorder endpoint dofs to (u(0), u'(0), u(l), u'(l))
      BtSinvB <- BtSinvB[c(3, 1, 4, 2), c(3, 1, 4, 2)]
      idx4 <- 4 * (e - 1) + 1:4
      Qpmu[idx4] <- Qpmu[idx4] + (t(Sigma_iB) %*% y_i)[c(3, 1, 4, 2)]
      if (n_cov > 0) {
        QpmuX[idx4, ] <- QpmuX[idx4, , drop = FALSE] +
          (t(Sigma_iB) %*% X_i)[c(3, 1, 4, 2), , drop = FALSE]
      }

      # full symmetric 4x4 block of BtSinvB as triplets (equivalent to the
      # 16 explicit assignments in likelihood_alpha2_precompute)
      idx16 <- count + 1:16
      i_[idx16] <- rep(idx4, times = 4)
      j_[idx16] <- rep(idx4, each = 4)
      x_[idx16] <- as.vector(BtSinvB)
      count <- count + 16

      v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))
      base <- base - 0.5 * sum(y_i * v_i) - sum(log(diag(R_i)))

      if (n_cov > 0) {
        SinvX <- backsolve(R_i, forwardsolve(t(R_i), X_i))
        XtSX <- XtSX + crossprod(X_i, SinvX)
        XtSy <- XtSy + as.vector(crossprod(X_i, v_i))
      }
    }

    if (is.null(det_R_count)) {
      idx_count <- seq_len(count)
      BtSB <- Matrix::sparseMatrix(i = i_[idx_count],
                                   j = j_[idx_count],
                                   x = x_[idx_count],
                                   dims = dim(Q))
      Qp <- Q + BtSB
      Qp <- Tc %*% Qp %*% t(Tc)
      R_count <- Matrix::Cholesky(forceSymmetric(Qp), LDL = FALSE, perm = TRUE)
      det_R_count <- Matrix::determinant(R_count, sqrt = TRUE)$modulus[1]
    }

    base <- base - det_R_count

    v <- c(as.matrix(Matrix::solve(R_count,
                                   Matrix::solve(R_count, Tc %*% Qpmu,
                                                 system = "P"),
                                   system = "L")))

    base <- base + 0.5 * sum(v^2) - 0.5 * n.o * log(2 * pi)

    if (n_cov > 0) {
      TcX <- as.matrix(Tc %*% QpmuX)
      VX <- as.matrix(Matrix::solve(R_count,
                                    Matrix::solve(R_count, TcX, system = "P"),
                                    system = "L"))
      H <- H + XtSX - crossprod(VX)
      h <- h + XtSy - as.vector(crossprod(VX, v))
    }
  }

  return(list(base = base, H = H, h = h, n_cov = n_cov))
}


#' Core quantities for the profiled directional alpha=1 likelihood
#'
#' Follows likelihood_alpha1_directional_precompute step by step (edge-wise
#' precision from Qalpha1_edges, directional constraint basis Tc), evaluating
#' the no-covariate likelihood with the raw y and accumulating H and h.
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range)
#' @param precomputed_data output of precompute_alpha1_directional
#' @param parameterization "matern" or "spde"
#' @return list(base, H, h, n_cov)
#' @noRd
profile_lik_core_alpha1_directional <- function(theta, precomputed_data,
                                                parameterization) {
  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  if (parameterization == "matern") {
    kappa <- sqrt(8 * 0.5) / exp(theta[3])
  } else {
    kappa <- exp(theta[3])
  }

  graph <- precomputed_data$graph
  Tc <- precomputed_data$Tc
  n_edges <- precomputed_data$n_edges
  n_cov <- precomputed_data$n_cov
  u_repl <- precomputed_data$u_repl

  Q.list <- Qalpha1_edges(c(1 / reciprocal_tau, kappa),
                          graph,
                          w = 0,
                          BC = 1,
                          build = FALSE)

  Q <- Matrix::sparseMatrix(i = Q.list$i,
                            j = Q.list$j,
                            x = Q.list$x,
                            dims = Q.list$dims)

  R <- Matrix::Cholesky(forceSymmetric(Tc %*% Q %*% t(Tc)),
                        LDL = FALSE, perm = TRUE)
  det_R <- Matrix::determinant(R, sqrt = TRUE)$modulus[1]

  base <- 0
  det_R_count <- NULL
  R_count <- NULL
  n.o <- 0

  H <- if (n_cov > 0) matrix(0, n_cov, n_cov) else NULL
  h <- if (n_cov > 0) numeric(n_cov) else NULL

  for (repl_y in seq_along(u_repl)) {
    curr_repl <- u_repl[repl_y]
    repl_name <- paste0("repl_", curr_repl)

    base <- base + det_R
    count <- 0
    i_ <- j_ <- x_ <- rep(0, 4 * length(precomputed_data$obs.edges))
    Qpmu <- rep(0, 2 * n_edges)
    if (n_cov > 0) {
      QpmuX <- matrix(0, 2 * n_edges, n_cov)
      XtSX <- matrix(0, n_cov, n_cov)
      XtSy <- numeric(n_cov)
    }

    for (j in seq_along(precomputed_data$obs.edges)) {
      e <- precomputed_data$obs.edges[j]
      edge_name <- paste0("edge_", e)

      if (is.null(precomputed_data$y_data[[repl_name]][[edge_name]])) {
        next
      }

      y_i <- precomputed_data$y_data[[repl_name]][[edge_name]]
      n.o <- n.o + length(y_i)
      X_i <- if (n_cov > 0) {
        precomputed_data$x_data[[repl_name]][[edge_name]]
      } else {
        NULL
      }

      D_matrix <- precomputed_data$D_data[[repl_name]][[edge_name]]

      S <- r_1(D_matrix, kappa = kappa, tau = 1 / reciprocal_tau)

      E.ind <- c(1:2)
      Obs.ind <- -E.ind

      Bt <- solve(S[E.ind, E.ind, drop = FALSE],
                  S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt

      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R_i <- base::chol(Sigma_i)

      Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))
      BtSinvB <- Bt %*% Sigma_iB

      E <- graph$E[e, ]
      if (E[1] == E[2]) {
        Qpmu[2 * (e - 1) + 1] <- Qpmu[2 * (e - 1) + 1] +
          sum(as.vector(t(Sigma_iB) %*% y_i))
        if (n_cov > 0) {
          QpmuX[2 * (e - 1) + 1, ] <- QpmuX[2 * (e - 1) + 1, ] +
            colSums(t(Sigma_iB) %*% X_i)
        }
        i_[count + 1] <- 2 * (e - 1) + 1
        j_[count + 1] <- 2 * (e - 1) + 1
        x_[count + 1] <- sum(BtSinvB)
        count <- count + 1
      } else {
        Qpmu[2 * (e - 1) + c(1, 2)] <- Qpmu[2 * (e - 1) + c(1, 2)] +
          as.vector(t(Sigma_iB) %*% y_i)
        if (n_cov > 0) {
          QpmuX[2 * (e - 1) + c(1, 2), ] <-
            QpmuX[2 * (e - 1) + c(1, 2), , drop = FALSE] +
            t(Sigma_iB) %*% X_i
        }
        i_[count + (1:4)] <- c(2 * (e - 1) + 1, 2 * (e - 1) + 1,
                               2 * (e - 1) + 2, 2 * (e - 1) + 2)
        j_[count + (1:4)] <- c(2 * (e - 1) + 1, 2 * (e - 1) + 2,
                               2 * (e - 1) + 1, 2 * (e - 1) + 2)
        x_[count + (1:4)] <- c(BtSinvB[1, 1], BtSinvB[1, 2],
                               BtSinvB[1, 2], BtSinvB[2, 2])
        count <- count + 4
      }

      v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))
      base <- base - 0.5 * sum(y_i * v_i) - sum(log(diag(R_i)))

      if (n_cov > 0) {
        SinvX <- backsolve(R_i, forwardsolve(t(R_i), X_i))
        XtSX <- XtSX + crossprod(X_i, SinvX)
        XtSy <- XtSy + as.vector(crossprod(X_i, v_i))
      }
    }

    if (is.null(det_R_count)) {
      i_ <- c(Q.list$i, i_[1:count])
      j_ <- c(Q.list$j, j_[1:count])
      x_ <- c(Q.list$x, x_[1:count])

      Qp <- Matrix::sparseMatrix(i = i_,
                                 j = j_,
                                 x = x_,
                                 dims = Q.list$dims)

      Qp <- Tc %*% Qp %*% t(Tc)
      R_count <- Matrix::Cholesky(forceSymmetric(Qp), LDL = FALSE, perm = TRUE)
      det_R_count <- Matrix::determinant(R_count, sqrt = TRUE)$modulus[1]
    }

    base <- base - det_R_count

    v <- c(as.matrix(Matrix::solve(R_count,
                                   Matrix::solve(R_count, Tc %*% Qpmu,
                                                 system = "P"),
                                   system = "L")))

    base <- base + 0.5 * sum(v^2) - 0.5 * n.o * log(2 * pi)

    if (n_cov > 0) {
      TcX <- as.matrix(Tc %*% QpmuX)
      VX <- as.matrix(Matrix::solve(R_count,
                                    Matrix::solve(R_count, TcX, system = "P"),
                                    system = "L"))
      H <- H + XtSX - crossprod(VX)
      h <- h + XtSy - as.vector(crossprod(VX, v))
    }
  }

  return(list(base = base, H = H, h = h, n_cov = n_cov))
}


#' Profile log-likelihood for the alpha=1 model (beta maximized out)
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range); length 3,
#'   beta is profiled out analytically.
#' @param graph metric_graph object
#' @param data_name name of the response variable
#' @param manual_y manual y values (if data_name is NULL)
#' @param X_cov matrix of covariates (fixed effects design)
#' @param repl replicates to be considered
#' @param BC which boundary condition to use (0, 1)
#' @param parameterization "matern" or "spde"
#' @param reml if TRUE, return the restricted (REML) log-likelihood
#' @return The profile (or restricted) log-likelihood. Use
#'   profile_beta_estimate() to recover beta_hat(theta).
#' @noRd
likelihood_alpha1_profile <- function(theta, graph, data_name = NULL,
                                      manual_y = NULL, X_cov = NULL,
                                      repl = NULL, BC = 1,
                                      parameterization = "matern",
                                      reml = FALSE) {
  check_theta_profile(theta)
  if (!is.null(X_cov)) {
    X_cov <- as.matrix(X_cov)
  }
  precomputeddata <- precompute_alpha1(graph, data_name = data_name,
                                       manual_y = manual_y, X_cov = X_cov,
                                       repl = repl)
  core <- profile_lik_core_alpha1(theta, graph, precomputeddata, BC,
                                  parameterization)
  return(profile_lik_finish(core, reml))
}

#' Profile log-likelihood for the alpha=1 model using precomputed data
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range); length 3.
#' @param graph metric_graph object
#' @param precomputeddata output of precompute_alpha1 (must have been built
#'   with the X_cov matrix to profile over)
#' @param BC which boundary condition to use (0, 1)
#' @param parameterization "matern" or "spde"
#' @param reml if TRUE, return the restricted (REML) log-likelihood
#' @return The profile (or restricted) log-likelihood.
#' @noRd
likelihood_alpha1_profile_precompute <- function(theta, graph, precomputeddata,
                                                 BC = 1,
                                                 parameterization = "matern",
                                                 reml = FALSE) {
  check_theta_profile(theta)
  core <- profile_lik_core_alpha1(theta, graph, precomputeddata, BC,
                                  parameterization)
  return(profile_lik_finish(core, reml))
}

#' Profile log-likelihood for the alpha=2 model (beta maximized out)
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range); length 3.
#' @param graph metric_graph object
#' @param data_name name of the response variable
#' @param manual_y manual y values (if data_name is NULL)
#' @param X_cov matrix of covariates (fixed effects design)
#' @param repl replicates to be considered
#' @param BC which boundary condition to use (0, 1)
#' @param parameterization "matern" or "spde"
#' @param reml if TRUE, return the restricted (REML) log-likelihood
#' @return The profile (or restricted) log-likelihood.
#' @noRd
likelihood_alpha2_profile <- function(theta, graph, data_name = NULL,
                                      manual_y = NULL, X_cov = NULL,
                                      repl = NULL, BC = 1,
                                      parameterization = "matern",
                                      reml = FALSE) {
  check_theta_profile(theta)
  if (!is.null(X_cov)) {
    X_cov <- as.matrix(X_cov)
  }
  precomputed_data <- precompute_alpha2(graph, data_name = data_name,
                                        manual_y = manual_y, X_cov = X_cov,
                                        repl = repl)
  core <- profile_lik_core_alpha2(theta, precomputed_data, BC,
                                  parameterization)
  return(profile_lik_finish(core, reml))
}

#' Profile log-likelihood for the alpha=2 model using precomputed data
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range); length 3.
#' @param precomputed_data output of precompute_alpha2 (must have been built
#'   with the X_cov matrix to profile over)
#' @param BC which boundary condition to use (0, 1)
#' @param parameterization "matern" or "spde"
#' @param reml if TRUE, return the restricted (REML) log-likelihood
#' @return The profile (or restricted) log-likelihood.
#' @noRd
likelihood_alpha2_profile_precompute <- function(theta, precomputed_data,
                                                 BC = 1,
                                                 parameterization = "matern",
                                                 reml = FALSE) {
  check_theta_profile(theta)
  core <- profile_lik_core_alpha2(theta, precomputed_data, BC,
                                  parameterization)
  return(profile_lik_finish(core, reml))
}

#' Profile log-likelihood for the directional alpha=1 model (beta maximized
#' out)
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range); length 3.
#' @param graph metric_graph object
#' @param data_name name of the response variable
#' @param manual_y manual y values (if data_name is NULL)
#' @param X_cov matrix of covariates (fixed effects design)
#' @param repl replicates to be considered
#' @param parameterization "matern" or "spde"
#' @param reml if TRUE, return the restricted (REML) log-likelihood
#' @return The profile (or restricted) log-likelihood (no sign flip, unlike
#'   likelihood_alpha1_directional_precompute with maximize = FALSE).
#' @noRd
likelihood_alpha1_directional_profile <- function(theta, graph,
                                                  data_name = NULL,
                                                  manual_y = NULL,
                                                  X_cov = NULL, repl = NULL,
                                                  parameterization = "matern",
                                                  reml = FALSE) {
  check_theta_profile(theta)
  if (!is.null(X_cov)) {
    X_cov <- as.matrix(X_cov)
  }
  precomputed_data <- precompute_alpha1_directional(graph,
                                                    data_name = data_name,
                                                    manual_y = manual_y,
                                                    X_cov = X_cov,
                                                    repl = repl)
  core <- profile_lik_core_alpha1_directional(theta, precomputed_data,
                                              parameterization)
  return(profile_lik_finish(core, reml))
}

#' Profile log-likelihood for the directional alpha=1 model using precomputed
#' data
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range); length 3.
#' @param precomputed_data output of precompute_alpha1_directional (must have
#'   been built with the X_cov matrix to profile over)
#' @param parameterization "matern" or "spde"
#' @param reml if TRUE, return the restricted (REML) log-likelihood
#' @return The profile (or restricted) log-likelihood (no sign flip, unlike
#'   likelihood_alpha1_directional_precompute with maximize = FALSE).
#' @noRd
likelihood_alpha1_directional_profile_precompute <- function(theta,
                                                             precomputed_data,
                                                             parameterization = "matern",
                                                             reml = FALSE) {
  check_theta_profile(theta)
  core <- profile_lik_core_alpha1_directional(theta, precomputed_data,
                                              parameterization)
  return(profile_lik_finish(core, reml))
}


#' Profiled fixed-effect estimate beta_hat(theta)
#'
#' Computes the closed-form maximizer of the log-likelihood over the fixed
#' effects for fixed covariance parameters, beta_hat(theta) = H^-1 h, together
#' with its conditional covariance H^-1 (H equals X' V^-1 X, with V the
#' marginal covariance of the observations). Typical use: run the numerical
#' optimization on likelihood_*_profile (or its precompute variant), then call
#' this once at the optimum to recover the fixed effects.
#'
#' @param theta (log sigma_e, log reciprocal_tau, log kappa/range); length 3.
#' @param model one of "alpha1", "alpha2", "alpha1_directional"
#' @param graph metric_graph object (not needed for "alpha2" /
#'   "alpha1_directional" when precomputed_data is supplied)
#' @param precomputed_data optional output of the matching precompute_*
#'   function; when NULL it is built from graph, data_name/manual_y, X_cov and
#'   repl
#' @param data_name name of the response variable (if precomputed_data is NULL)
#' @param manual_y manual y values (if precomputed_data and data_name are NULL)
#' @param X_cov matrix of covariates (if precomputed_data is NULL)
#' @param repl replicates to be considered (if precomputed_data is NULL)
#' @param BC which boundary condition to use (0, 1); ignored for
#'   "alpha1_directional"
#' @param parameterization "matern" or "spde"
#' @return list with elements beta (beta_hat(theta)), vcov (H^-1), H and h.
#'   With zero covariates, beta has length 0.
#' @noRd
profile_beta_estimate <- function(theta,
                                  model = c("alpha1", "alpha2",
                                            "alpha1_directional"),
                                  graph = NULL, precomputed_data = NULL,
                                  data_name = NULL, manual_y = NULL,
                                  X_cov = NULL, repl = NULL, BC = 1,
                                  parameterization = "matern") {
  check_theta_profile(theta)
  model <- match.arg(model)
  if (!is.null(X_cov)) {
    X_cov <- as.matrix(X_cov)
  }

  if (model == "alpha1") {
    if (is.null(precomputed_data)) {
      precomputed_data <- precompute_alpha1(graph, data_name = data_name,
                                            manual_y = manual_y,
                                            X_cov = X_cov, repl = repl)
    }
    core <- profile_lik_core_alpha1(theta, graph, precomputed_data, BC,
                                    parameterization)
  } else if (model == "alpha2") {
    if (is.null(precomputed_data)) {
      precomputed_data <- precompute_alpha2(graph, data_name = data_name,
                                            manual_y = manual_y,
                                            X_cov = X_cov, repl = repl)
    }
    core <- profile_lik_core_alpha2(theta, precomputed_data, BC,
                                    parameterization)
  } else {
    if (is.null(precomputed_data)) {
      precomputed_data <- precompute_alpha1_directional(graph,
                                                        data_name = data_name,
                                                        manual_y = manual_y,
                                                        X_cov = X_cov,
                                                        repl = repl)
    }
    core <- profile_lik_core_alpha1_directional(theta, precomputed_data,
                                                parameterization)
  }

  if (core$n_cov == 0) {
    return(list(beta = numeric(0), vcov = matrix(0, 0, 0),
                H = matrix(0, 0, 0), h = numeric(0)))
  }

  bs <- profile_beta_solve(core$H, core$h)
  return(list(beta = bs$beta, vcov = bs$vcov, H = core$H, h = core$h))
}
