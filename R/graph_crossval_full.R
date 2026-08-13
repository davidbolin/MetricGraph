# Full-graph plug-in leave-one-out cross-validation for alpha=1 models.
#
# This file intentionally contains no subgraph, k-fold, alpha=2, or wrapper
# dispatch. The directional selected-inverse backend is implemented in
# src/crossval_full.cpp.

cv_build_idx_map <- function(precomputed_data) {
  graph     <- precomputed_data$graph
  PtE       <- graph$get_PtE()
  repl_vec  <- graph$.__enclos_env__$private$data[[".group"]]
  u_repl    <- precomputed_data$u_repl
  obs.edges <- precomputed_data$obs.edges

  idx_cache <- integer(0)

  for (repl_y in seq_along(u_repl)) {
    curr_repl  <- u_repl[repl_y]
    repl_name  <- paste0("repl_", curr_repl)

    ind_repl       <- (repl_vec == curr_repl)
    global_in_repl <- which(ind_repl)
    PtE_repl       <- PtE[ind_repl, , drop = FALSE]

    for (j in seq_along(obs.edges)) {
      e         <- obs.edges[j]
      edge_name <- paste0("edge_", e)

      no_na <- precomputed_data$no_na_indices[[repl_name]][[edge_name]]
      if (is.null(no_na)) next

      obs_in_repl <- (PtE_repl[, 1] == e)
      global_on_edge <- global_in_repl[obs_in_repl]
      idx_cache <- c(idx_cache, global_on_edge[no_na])
    }
  }
  idx_cache
}


# ---- Directional CV core ------------------------------------------------

#' Build the full-graph directional alpha=1 precision triplets
#'
#' @param precomputed_data Output from `precompute_alpha1_directional()`.
#' @param kappa Positive SPDE scale parameter.
#' @param reciprocal_tau Positive reciprocal precision parameter.
#' @return A sparse triplet list as returned by `Qalpha1_edges()`.
#' @noRd
.cv_directional_q_list <- function(precomputed_data, kappa, reciprocal_tau) {
  Qalpha1_edges(
    c(1 / reciprocal_tau, kappa),
    precomputed_data$graph,
    w = 0,
    BC = 1,
    build = FALSE
  )
}


#' Full-graph plug-in LOO for the directional alpha=1 model
#'
#' @param theta Length-three vector `(log sigma_e, log reciprocal_tau,
#'   log kappa)`.
#' @param precomputed_data Output from `precompute_alpha1_directional()`.
#' @param parameterization Either `"spde"` or `"matern"`.
#' @param method Either the sparse selected-inverse backend or the dense R
#'   reference implementation.
#' @return A list with `mu`, `var`, `idx`, `beta_hat`, and `H`.
#' @noRd
cv_core_alpha1_directional <- function(
    theta,
    precomputed_data,
    parameterization = "spde",
    method = c("selinv", "reference")) {
  method <- match.arg(method)
  parameterization <- match.arg(parameterization, c("spde", "matern"))

  if (method == "reference") {
    return(.cv_core_alpha1_directional_reference(
      theta,
      precomputed_data,
      parameterization = parameterization
    ))
  }

  sigma_e <- exp(theta[1L])
  reciprocal_tau <- exp(theta[2L])
  kappa <- if (parameterization == "matern") {
    sqrt(8 * 0.5) / exp(theta[3L])
  } else {
    exp(theta[3L])
  }
  Q_list <- .cv_directional_q_list(
    precomputed_data,
    kappa = kappa,
    reciprocal_tau = reciprocal_tau
  )
  cpp_out <- cv_loo_selinv_cpp(
    precomputed_data = precomputed_data,
    edge_endpoints = precomputed_data$graph$E,
    Q_list = Q_list,
    sigma_e = sigma_e,
    reciprocal_tau = reciprocal_tau,
    kappa = kappa
  )

  list(
    mu = cpp_out$mu,
    var = cpp_out$var,
    idx = cv_build_idx_map(precomputed_data),
    beta_hat = as.vector(cpp_out$beta_hat),
    H = if (precomputed_data$n_cov > 0L) cpp_out$H else NULL
  )
}

#' Dense reference plug-in LOO for the directional alpha=1 model
#'
#' Uses a single Cholesky factorization of Q~ and an explicit dense triangular
#' solve. This is intentionally retained only as a small-graph test oracle.
#' Replicates are handled via block-diagonal P across replicates.
#'
#' @param theta Length-three parameter vector containing log measurement error,
#'   log reciprocal precision, and log SPDE scale (or log range).
#' @param precomputed_data output of precompute_alpha1_directional
#' @param parameterization "spde" (default) or "matern"
#' @return list with elements:
#'   mu      - numeric(n_cache): predictive means in cache order
#'   var     - numeric(n_cache): predictive variances (includes sigma_e^2)
#'   idx     - integer(n_cache): data-row index of each cache obs
#'   beta_hat- profiled fixed effects (length n_cov)
#'   H       - n_cov x n_cov information matrix
#' @noRd
.cv_core_alpha1_directional_reference <- function(
    theta,
    precomputed_data,
    parameterization = "spde") {
  parameterization <- match.arg(parameterization, c("spde", "matern"))

  sigma_e        <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  if (parameterization == "matern") {
    kappa <- sqrt(8 * 0.5) / exp(theta[3])
  } else {
    kappa <- exp(theta[3])
  }

  graph     <- precomputed_data$graph
  Tc        <- precomputed_data$Tc
  n_edges   <- precomputed_data$n_edges
  n_cov     <- precomputed_data$n_cov
  u_repl    <- precomputed_data$u_repl
  obs.edges <- precomputed_data$obs.edges
  n_dof     <- 2L * n_edges        # dimension of Q (2*n_edges space)

  # Build Q (once)
  Q.list <- .cv_directional_q_list(precomputed_data, kappa, reciprocal_tau)

  Q <- Matrix::sparseMatrix(i = Q.list$i, j = Q.list$j, x = Q.list$x,
                             dims = Q.list$dims)

  # Global accumulators for H, h (across replicates)
  H <- if (n_cov > 0) matrix(0, n_cov, n_cov) else NULL
  h <- if (n_cov > 0) numeric(n_cov)           else NULL

  # Per-replicate storage (filled in first pass, used in second)
  R_count     <- NULL   # Cholesky of Q~ (built from first replicate)
  vy_list     <- list() # v_y per replicate  (nFree vector)
  VX_list     <- list() # VX per replicate   (nFree x n_cov, only if n_cov > 0)
  Z_list      <- list() # Z per replicate    (nFree x n_cr matrix)

  # Global cache-order arrays (built during first pass)
  all_y           <- numeric(0)
  all_diag_sinv   <- numeric(0)
  all_sinv_r      <- numeric(0)      # Sigma_e^{-1} y, concatenated
  all_sinv_x      <- if (n_cov > 0) matrix(0, 0, n_cov) else NULL
  # Track cache start position of each replicate (1-based)
  repl_cache_starts <- integer(length(u_repl))
  cache_pos <- 0L

  # ---- First pass: build R_count, Z, H, h, per-obs arrays ---------------
  for (repl_y in seq_along(u_repl)) {
    curr_repl  <- u_repl[repl_y]
    repl_name  <- paste0("repl_", curr_repl)
    repl_cache_starts[repl_y] <- cache_pos + 1L

    # BtSinvB accumulation triplets for Q~ (same as likelihood)
    count_btsb <- 0L
    i_b <- j_b <- x_b <- rep(0, 4 * length(obs.edges))

    Qpmu <- numeric(n_dof)
    if (n_cov > 0) {
      QpmuX <- matrix(0, n_dof, n_cov)
      XtSX  <- matrix(0, n_cov, n_cov)
      XtSy  <- numeric(n_cov)
    }

    # C matrix triplets for this replicate
    C_i <- integer(0); C_j <- integer(0); C_x <- numeric(0)

    # Per-obs arrays for this replicate
    diag_sinv_r  <- numeric(0)
    sinv_r_r     <- numeric(0)
    sinv_x_r     <- if (n_cov > 0) matrix(0, 0, n_cov) else NULL
    y_r          <- numeric(0)
    col_local <- 0L   # column within C for this replicate

    for (j in seq_along(obs.edges)) {
      e         <- obs.edges[j]
      edge_name <- paste0("edge_", e)

      y_i <- precomputed_data$y_data[[repl_name]][[edge_name]]
      if (is.null(y_i) || length(y_i) == 0L) next

      n_i <- length(y_i)
      X_i      <- if (n_cov > 0) precomputed_data$x_data[[repl_name]][[edge_name]] else NULL
      D_matrix <- precomputed_data$D_data[[repl_name]][[edge_name]]

      S <- r_1(D_matrix, kappa = kappa, tau = 1 / reciprocal_tau)

      E.ind   <- 1:2
      Obs.ind <- -E.ind

      Bt      <- solve(S[E.ind, E.ind, drop = FALSE],
                       S[E.ind, Obs.ind, drop = FALSE])   # 2 x n_i
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
                 S[Obs.ind, E.ind, drop = FALSE] %*% Bt   # n_i x n_i
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      R_i    <- base::chol(Sigma_i)

      # Sigma_iB = Sigma_i^{-1} B_e^T  (n_i x 2)
      Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))
      BtSinvB  <- Bt %*% Sigma_iB   # 2 x 2

      Sinv_e <- chol2inv(R_i)

      # v_i = Sigma_i^{-1} y_i
      v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))

      # SinvX_i (for H/h and REML PX)
      if (n_cov > 0) {
        SinvX_i <- backsolve(R_i, forwardsolve(t(R_i), X_i))
        XtSX <- XtSX + crossprod(X_i, SinvX_i)
        XtSy <- XtSy + as.vector(crossprod(X_i, v_i))
        sinv_x_r <- rbind(sinv_x_r, SinvX_i)
      }

      # Qpmu and BtSinvB accumulation (mirrors profile_lik_core_alpha1_directional)
      E_e <- graph$E[e, ]
      is_self_loop <- (E_e[1] == E_e[2])

      if (is_self_loop) {
        Qpmu[2 * (e - 1) + 1] <- Qpmu[2 * (e - 1) + 1] +
          sum(as.vector(t(Sigma_iB) %*% y_i))
        if (n_cov > 0)
          QpmuX[2 * (e - 1) + 1, ] <- QpmuX[2 * (e - 1) + 1, ] +
            colSums(t(Sigma_iB) %*% X_i)
        i_b[count_btsb + 1] <- 2 * (e - 1) + 1
        j_b[count_btsb + 1] <- 2 * (e - 1) + 1
        x_b[count_btsb + 1] <- sum(BtSinvB)
        count_btsb <- count_btsb + 1L
      } else {
        dofs <- 2 * (e - 1) + c(1L, 2L)
        Qpmu[dofs] <- Qpmu[dofs] + as.vector(t(Sigma_iB) %*% y_i)
        if (n_cov > 0)
          QpmuX[dofs, ] <- QpmuX[dofs, , drop = FALSE] +
            t(Sigma_iB) %*% X_i
        idx4 <- count_btsb + 1:4
        i_b[idx4] <- c(dofs[1], dofs[1], dofs[2], dofs[2])
        j_b[idx4] <- c(dofs[1], dofs[2], dofs[1], dofs[2])
        x_b[idx4] <- c(BtSinvB[1,1], BtSinvB[1,2], BtSinvB[1,2], BtSinvB[2,2])
        count_btsb <- count_btsb + 4L
      }

      # C matrix columns: one column per obs k, mapping into 2*n_edges dof space
      for (k in seq_len(n_i)) {
        col_k <- col_local + k
        if (is_self_loop) {
          # Both endpoint dofs collapse to position 2*(e-1)+1
          C_i <- c(C_i, 2L * (e - 1L) + 1L)
          C_j <- c(C_j, col_k)
          C_x <- c(C_x, Sigma_iB[k, 1] + Sigma_iB[k, 2])
        } else {
          C_i <- c(C_i, 2L * (e - 1L) + 1L, 2L * (e - 1L) + 2L)
          C_j <- c(C_j, col_k, col_k)
          C_x <- c(C_x, Sigma_iB[k, 1], Sigma_iB[k, 2])
        }
      }

      # Per-obs arrays
      diag_sinv_r  <- c(diag_sinv_r,  diag(Sinv_e))
      sinv_r_r     <- c(sinv_r_r,     as.vector(v_i))
      y_r          <- c(y_r,          y_i)
      col_local    <- col_local + n_i
    }  # end edge loop

    n_cr <- col_local   # number of cache obs in this replicate

    # Build R_count from first replicate (reused for all)
    if (is.null(R_count)) {
      idx_btsb <- seq_len(count_btsb)
      i_full   <- c(Q.list$i, i_b[idx_btsb])
      j_full   <- c(Q.list$j, j_b[idx_btsb])
      x_full   <- c(Q.list$x, x_b[idx_btsb])

      Qp <- Matrix::sparseMatrix(i = i_full, j = j_full, x = x_full,
                                  dims = Q.list$dims)
      Qp <- Matrix::forceSymmetric(Tc %*% Qp %*% t(Tc))
      R_count <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
    }

    # v_y = L^{-1} P Tc Qpmu  (system="P" gives Pb under P A P' = LL')
    Tc_Qpmu <- as.numeric(Tc %*% Qpmu)
    v_y <- c(as.matrix(
      Matrix::solve(R_count,
                    Matrix::solve(R_count, Tc_Qpmu, system = "P"),
                    system = "L")))
    vy_list[[repl_name]] <- v_y

    if (n_cov > 0) {
      TcQpmuX <- as.matrix(Tc %*% QpmuX)
      VX_r <- as.matrix(
        Matrix::solve(R_count,
                      Matrix::solve(R_count, TcQpmuX, system = "P"),
                      system = "L"))
      VX_list[[repl_name]] <- VX_r
      H <- H + XtSX     - crossprod(VX_r)
      h <- h + XtSy - as.vector(crossprod(VX_r, v_y))
    }

    # Z = L^{-1} P Tc C  (system="P" gives Pb)
    if (n_cr > 0L) {
      C_sp  <- Matrix::sparseMatrix(i = C_i, j = C_j, x = C_x,
                                     dims = c(n_dof, n_cr))
      TcC   <- Tc %*% C_sp
      Z_r   <- as.matrix(
        Matrix::solve(R_count,
                      Matrix::solve(R_count, TcC, system = "P"),
                      system = "L"))
    } else {
      Z_r <- matrix(0, nrow(Tc), 0L)
    }
    Z_list[[repl_name]] <- Z_r

    # Accumulate global arrays
    all_y           <- c(all_y,          y_r)
    all_diag_sinv   <- c(all_diag_sinv,  diag_sinv_r)
    all_sinv_r      <- c(all_sinv_r,     sinv_r_r)
    if (n_cov > 0) all_sinv_x <- rbind(all_sinv_x, sinv_x_r)

    cache_pos <- cache_pos + n_cr
  }  # end replicate loop

  n_cache <- cache_pos

  # beta_hat
  if (n_cov > 0) {
    bs       <- profile_beta_solve(H, h)
    beta_hat <- as.vector(bs$beta)   # ensure plain numeric vector
  } else {
    beta_hat <- numeric(0)
  }

  # ---- Second pass: compute d, Pr, PX per replicate ----------------------
  d      <- numeric(n_cache)
  Pr_vec <- numeric(n_cache)

  for (repl_y in seq_along(u_repl)) {
    curr_repl <- u_repl[repl_y]
    repl_name <- paste0("repl_", curr_repl)
    Z_r       <- Z_list[[repl_name]]
    n_cr      <- ncol(Z_r)
    if (n_cr == 0L) next

    start_r <- repl_cache_starts[repl_y]
    idx_r   <- start_r:(start_r + n_cr - 1L)

    v_y_r <- vy_list[[repl_name]]
    if (n_cov > 0) {
      VX_r <- VX_list[[repl_name]]
      v_r  <- v_y_r - as.vector(VX_r %*% beta_hat)
    } else {
      v_r <- v_y_r
    }

    # Sinv_r for the residual r = y - X beta_hat
    sinv_r_adj <- all_sinv_r[idx_r]
    if (n_cov > 0) {
      sinv_r_adj <- sinv_r_adj -
        as.vector(all_sinv_x[idx_r, , drop = FALSE] %*% beta_hat)
    }

    d[idx_r]      <- all_diag_sinv[idx_r] - colSums(Z_r^2)
    Pr_vec[idx_r] <- sinv_r_adj           - as.vector(t(Z_r) %*% v_r)

  }

  mu_cache <- all_y - Pr_vec / d
  var_cache <- 1 / d

  # Build data-row index map
  idx <- cv_build_idx_map(precomputed_data)

  list(mu       = mu_cache,
       var      = var_cache,
       idx      = idx,
       beta_hat = beta_hat,
       H        = H)
}


# ---- idx-map helpers ---------------------------------------------------

#' Build cache-order → data-row index for alpha=1 (non-directional)
#'
#' Replays the precompute_alpha1 replicate × edge loop and records which
#' global data row each cache position corresponds to.
#'
#' @param graph metric_graph object (same as used for precomputation)
#' @param precomputed_data output of precompute_alpha1
#' @param y_resp_full full y vector (same length as nrow(graph$get_PtE()));
#'   NAs are used to identify non-observed positions.  If NULL an attempt is
#'   made to infer it from the graph's private data; if that fails a warning
#'   is issued and 1:n_cache is returned as a fallback.
#' @return integer vector of length n_cache
#' @noRd
cv_build_idx_map_alpha1 <- function(graph, precomputed_data,
                                    y_resp_full = NULL) {
  PtE       <- graph$get_PtE()
  repl_vec  <- graph$.__enclos_env__$private$data[[".group"]]
  u_repl    <- precomputed_data$u_repl
  if (is.null(y_resp_full)) {
    data_names <- names(graph$.__enclos_env__$private$data)
    data_names <- data_names[!startsWith(data_names, ".")]
    if (length(data_names) > 0L) {
      y_resp_full <- graph$.__enclos_env__$private$data[[data_names[1L]]]
    } else {
      n_cache <- sum(vapply(
        precomputed_data$edge_cache,
        function(edge_cache) sum(lengths(edge_cache$y)),
        integer(1L)
      ))
      warning("cv_build_idx_map_alpha1: cannot determine y_resp_full; ",
              "returning 1:n_cache as fallback")
      return(seq_len(n_cache))
    }
  }

  idx_cache <- integer(0)

  for (repl_y in seq_along(u_repl)) {
    ind_repl       <- (repl_vec == u_repl[repl_y])
    global_in_repl <- which(ind_repl)
    PtE_repl       <- PtE[ind_repl, , drop = FALSE]
    y_repl         <- y_resp_full[ind_repl]
    edge_cache     <- precomputed_data$edge_cache[[repl_y]]

    for (e in edge_cache$e) {
      obs_in_repl <- (PtE_repl[, 1L] == e)
      no_na       <- !is.na(y_repl[obs_in_repl])
      idx_cache   <- c(idx_cache, global_in_repl[obs_in_repl][no_na])
    }
  }
  idx_cache
}


#' Full-graph plug-in LOO for the non-directional alpha=1 model
#'
#' @param theta Length-three vector `(log sigma_e, log reciprocal_tau,
#'   log kappa)`.
#' @param graph The full `metric_graph` used for precomputation.
#' @param precomputed_data Output from `precompute_alpha1()`.
#' @param y_resp_full Full response vector, used to restore graph-data order.
#' @param parameterization Either `"spde"` or `"matern"`.
#' @param BC Boundary condition used by the WM1 precision.
#' @return A list with `mu`, `var`, `idx`, `beta_hat`, and `H`.
#' @noRd
cv_core_alpha1 <- function(
    theta,
    graph,
    precomputed_data,
    y_resp_full = NULL,
    parameterization = "spde",
    BC = 1L) {
  parameterization <- match.arg(parameterization, c("spde", "matern"))

  sigma_e        <- exp(theta[1L])
  reciprocal_tau <- exp(theta[2L])
  if (parameterization == "matern") {
    kappa <- sqrt(8L * 0.5) / exp(theta[3L])
  } else {
    kappa <- exp(theta[3L])
  }

  nV        <- nrow(graph$V)
  n_dof     <- nV
  u_repl    <- precomputed_data$u_repl
  obs.edges <- precomputed_data$obs.edges
  n_cov     <- precomputed_data$n_cov

  # Build Q once (vertex space, nV x nV)
  Q.list <- spde_precision(kappa = kappa, tau = 1L / reciprocal_tau, alpha = 1L,
                           graph = graph, build = FALSE, BC = BC)
  Q <- Matrix::sparseMatrix(i = Q.list$i, j = Q.list$j, x = Q.list$x,
                             dims = Q.list$dims)

  # Global H, h accumulators (across replicates)
  H <- if (n_cov > 0L) matrix(0, n_cov, n_cov) else NULL
  h <- if (n_cov > 0L) numeric(n_cov)           else NULL

  # Per-replicate storage
  R_count   <- NULL
  vy_list   <- list()
  VX_list   <- list()
  Z_list    <- list()

  # Global cache-order arrays
  all_y           <- numeric(0L)
  all_diag_sinv   <- numeric(0L)
  all_sinv_r      <- numeric(0L)
  all_sinv_x      <- if (n_cov > 0L) matrix(0, 0L, n_cov) else NULL

  repl_cache_starts <- integer(length(u_repl))
  cache_pos         <- 0L

  # ---- First pass: build R_count, Z, H, h, per-obs arrays ---------------
  for (repl_y in seq_along(u_repl)) {
    curr_repl  <- u_repl[repl_y]
    repl_name  <- paste0("repl_", curr_repl)
    repl_cache_starts[repl_y] <- cache_pos + 1L

    count_btsb <- 0L
    i_b <- j_b <- x_b <- rep(0, 4L * length(obs.edges))

    Qpmu <- numeric(n_dof)
    if (n_cov > 0L) {
      QpmuX <- matrix(0, n_dof, n_cov)
      XtSX  <- matrix(0, n_cov, n_cov)
      XtSy  <- numeric(n_cov)
    }

    C_i <- integer(0L); C_j <- integer(0L); C_x <- numeric(0L)

    diag_sinv_r <- numeric(0L)
    sinv_r_r    <- numeric(0L)
    sinv_x_r    <- if (n_cov > 0L) matrix(0, 0L, n_cov) else NULL
    y_r         <- numeric(0L)

    col_local <- 0L

    edge_cache <- precomputed_data$edge_cache[[repl_y]]
    for (j in seq_along(edge_cache$e)) {
      e <- edge_cache$e[j]
      y_i <- edge_cache$y[[j]]
      n_i      <- length(y_i)
      X_i      <- if (n_cov > 0L) edge_cache$X[[j]] else NULL
      D_matrix <- edge_cache$D[[j]]

      S <- r_1(D_matrix, kappa = kappa, tau = 1L / reciprocal_tau)

      E.ind   <- 1L:2L
      Obs.ind <- -E.ind

      Bt      <- solve(S[E.ind, E.ind, drop = FALSE],
                       S[E.ind, Obs.ind, drop = FALSE])   # 2 x n_i
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
                 S[Obs.ind, E.ind,   drop = FALSE] %*% Bt  # n_i x n_i
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2L
      R_i      <- base::chol(Sigma_i)

      Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))  # n_i x 2
      BtSinvB  <- Bt %*% Sigma_iB                              # 2 x 2

      Sinv_e <- chol2inv(R_i)

      v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))

      if (n_cov > 0L) {
        SinvX_i <- backsolve(R_i, forwardsolve(t(R_i), X_i))
        XtSX     <- XtSX + crossprod(X_i, SinvX_i)
        XtSy     <- XtSy + as.vector(crossprod(X_i, v_i))
        sinv_x_r <- rbind(sinv_x_r, SinvX_i)
      }

      # Accumulate into vertex space (E[e,] are 1-based vertex indices)
      E_e          <- edge_cache$E[j, ]
      is_self_loop <- (E_e[1L] == E_e[2L])

      if (is_self_loop) {
        Qpmu[E_e[1L]] <- Qpmu[E_e[1L]] + sum(as.vector(t(Sigma_iB) %*% y_i))
        if (n_cov > 0L)
          QpmuX[E_e[1L], ] <- QpmuX[E_e[1L], ] + colSums(t(Sigma_iB) %*% X_i)
        i_b[count_btsb + 1L] <- E_e[1L]
        j_b[count_btsb + 1L] <- E_e[1L]
        x_b[count_btsb + 1L] <- sum(BtSinvB)
        count_btsb <- count_btsb + 1L
      } else {
        Qpmu[E_e] <- Qpmu[E_e] + as.vector(t(Sigma_iB) %*% y_i)
        if (n_cov > 0L)
          QpmuX[E_e, ] <- QpmuX[E_e, , drop = FALSE] + t(Sigma_iB) %*% X_i
        idx4         <- count_btsb + 1L:4L
        i_b[idx4]   <- c(E_e[1L], E_e[1L], E_e[2L], E_e[2L])
        j_b[idx4]   <- c(E_e[1L], E_e[2L], E_e[1L], E_e[2L])
        x_b[idx4]   <- c(BtSinvB[1L, 1L], BtSinvB[1L, 2L],
                          BtSinvB[1L, 2L], BtSinvB[2L, 2L])
        count_btsb  <- count_btsb + 4L
      }

      # C matrix columns: one column per observation
      for (k in seq_len(n_i)) {
        col_k <- col_local + k
        if (is_self_loop) {
          C_i <- c(C_i, E_e[1L])
          C_j <- c(C_j, col_k)
          C_x <- c(C_x, Sigma_iB[k, 1L] + Sigma_iB[k, 2L])
        } else {
          C_i <- c(C_i, E_e[1L], E_e[2L])
          C_j <- c(C_j, col_k,   col_k)
          C_x <- c(C_x, Sigma_iB[k, 1L], Sigma_iB[k, 2L])
        }
      }

      diag_sinv_r <- c(diag_sinv_r, diag(Sinv_e))
      sinv_r_r    <- c(sinv_r_r,    as.vector(v_i))
      y_r         <- c(y_r,         y_i)
      col_local   <- col_local + n_i
    }  # end edge loop

    n_cr <- col_local

    # Build R_count from the first replicate (reused for all replicates)
    if (is.null(R_count)) {
      idx_btsb <- seq_len(count_btsb)
      i_full   <- c(Q.list$i, i_b[idx_btsb])
      j_full   <- c(Q.list$j, j_b[idx_btsb])
      x_full   <- c(Q.list$x, x_b[idx_btsb])

      Qp <- Matrix::sparseMatrix(i = i_full, j = j_full, x = x_full,
                                  dims = Q.list$dims)
      # alpha=1: Q_tilde is already in vertex space — no Tc
      R_count <- Matrix::Cholesky(Matrix::forceSymmetric(Qp),
                                   LDL = FALSE, perm = TRUE)
    }

    # v_y = L^{-1} P Qpmu  (system="P" gives Pb; no Tc for alpha=1)
    v_y <- c(as.matrix(
      Matrix::solve(R_count,
                    Matrix::solve(R_count, Qpmu, system = "P"),
                    system = "L")))
    vy_list[[repl_name]] <- v_y

    if (n_cov > 0L) {
      VX_r <- as.matrix(
        Matrix::solve(R_count,
                      Matrix::solve(R_count, QpmuX, system = "P"),
                      system = "L"))
      VX_list[[repl_name]] <- VX_r
      H <- H + XtSX - crossprod(VX_r)
      h <- h + XtSy - as.vector(crossprod(VX_r, v_y))
    }

    # Z = L^{-1} P C  (system="P" gives Pb; no Tc for alpha=1)
    if (n_cr > 0L) {
      C_sp <- Matrix::sparseMatrix(i = C_i, j = C_j, x = C_x,
                                    dims = c(n_dof, n_cr))
      Z_r  <- as.matrix(
        Matrix::solve(R_count,
                      Matrix::solve(R_count, C_sp, system = "P"),
                      system = "L"))
    } else {
      Z_r <- matrix(0, n_dof, 0L)
    }
    Z_list[[repl_name]] <- Z_r

    all_y           <- c(all_y,          y_r)
    all_diag_sinv   <- c(all_diag_sinv,  diag_sinv_r)
    all_sinv_r      <- c(all_sinv_r,     sinv_r_r)
    if (n_cov > 0L) all_sinv_x <- rbind(all_sinv_x, sinv_x_r)

    cache_pos <- cache_pos + n_cr
  }  # end replicate loop

  n_cache <- cache_pos

  if (n_cov > 0L) {
    bs       <- profile_beta_solve(H, h)
    beta_hat <- as.vector(bs$beta)
  } else {
    beta_hat <- numeric(0L)
  }

  # ---- Second pass: compute the diagonal precision and residual score ----
  d      <- numeric(n_cache)
  Pr_vec <- numeric(n_cache)

  for (repl_y in seq_along(u_repl)) {
    curr_repl <- u_repl[repl_y]
    repl_name <- paste0("repl_", curr_repl)
    Z_r       <- Z_list[[repl_name]]
    n_cr      <- ncol(Z_r)
    if (n_cr == 0L) next

    start_r <- repl_cache_starts[repl_y]
    idx_r   <- start_r:(start_r + n_cr - 1L)

    v_y_r <- vy_list[[repl_name]]
    if (n_cov > 0L) {
      VX_r <- VX_list[[repl_name]]
      v_r  <- v_y_r - as.vector(VX_r %*% beta_hat)
    } else {
      v_r <- v_y_r
    }

    sinv_r_adj <- all_sinv_r[idx_r]
    if (n_cov > 0L) {
      sinv_r_adj <- sinv_r_adj -
        as.vector(all_sinv_x[idx_r, , drop = FALSE] %*% beta_hat)
    }

    d[idx_r]      <- all_diag_sinv[idx_r] - colSums(Z_r^2L)
    Pr_vec[idx_r] <- sinv_r_adj           - as.vector(t(Z_r) %*% v_r)

  }

  mu_cache <- all_y - Pr_vec / d
  var_cache <- 1 / d

  idx <- cv_build_idx_map_alpha1(graph, precomputed_data, y_resp_full)

  list(mu       = mu_cache,
       var      = var_cache,
       idx      = idx,
       beta_hat = beta_hat,
       H        = H)
}
