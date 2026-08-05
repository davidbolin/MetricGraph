# Unconditional prior simulation of Whittle-Matérn fields on metric graphs.
#
#   Method A ("direct")  – O(m^3) per edge via dense Cholesky of Sigma*_e.
#   Method B ("kriging") – O(m)   per edge via Markov recursion + correction.
#
# Public methods:
#   simulate.metric_graph()  – serial, or parallel when parallel=TRUE
#   simulate_parallel()      – thin wrapper around simulate.metric_graph



# ============================================================
# Bridge predictor and covariance
# ============================================================

#' @noRd
bridge_alpha1 <- function(t_abs, l_e, kappa, tau) {
  t_ends   <- c(0, l_e)
  D_xb     <- outer(t_abs, t_ends, `-`)
  D_bb     <- outer(t_ends, t_ends, `-`)
  D_xx     <- outer(t_abs, t_abs,  `-`)
  Sigma_bb <- r_1(D_bb, kappa, tau)
  Sigma_xb <- r_1(D_xb, kappa, tau)
  Sigma_xx <- r_1(D_xx, kappa, tau)
  S        <- Sigma_xb %*% solve(Sigma_bb)
  list(S = S, Sigma_star = Sigma_xx - S %*% t(Sigma_xb))
}

#' @noRd
bridge_alpha2 <- function(t_abs, l_e, kappa, tau) {
  m      <- length(t_abs)
  t_ends <- c(0, l_e)

  Sigma_bb <- matrix(0, 4, 4)
  D_vv     <- outer(t_ends, t_ends, `-`)
  Sigma_bb[c(1,3), c(1,3)] <- r_2(D_vv, kappa, tau, 0)
  D_dd     <- as.matrix(dist(t_ends))
  Sigma_bb[c(2,4), c(2,4)] <- -r_2(D_dd, kappa, tau, 2)
  D_dv     <- outer(t_ends, t_ends, `-`)
  Sigma_bb[c(2,4), c(1,3)] <- r_2(D_dv, kappa, tau, 1)
  Sigma_bb[c(1,3), c(2,4)] <- t(Sigma_bb[c(2,4), c(1,3)])

  D_xv     <- outer(t_abs, t_ends, `-`)
  Sigma_xb <- matrix(0, m, 4)
  Sigma_xb[, c(1,3)] <-  r_2(D_xv,  kappa, tau, 0)
  Sigma_xb[, c(2,4)] <-  r_2(-D_xv, kappa, tau, 1)

  D_xx     <- outer(t_abs, t_abs, `-`)
  Sigma_xx <- r_2(D_xx, kappa, tau, 0)
  S        <- Sigma_xb %*% solve(Sigma_bb)
  list(S = S, Sigma_star = Sigma_xx - S %*% t(Sigma_xb))
}


# ============================================================
# Markov chain simulation (Method B internals)
# ============================================================

#' @noRd
sim_markov_alpha1 <- function(t_aug, kappa, tau) {
  n      <- length(t_aug)
  r0     <- 1 / (2 * kappa * tau^2)
  h_vals <- diff(t_aug)
  phi    <- exp(-kappa * h_vals)
  sd_vec <- sqrt(r0 * (1 - phi * phi))
  x      <- numeric(n)
  x[1]   <- rnorm(1, 0, sqrt(r0))
  for (i in seq_len(n - 1L))
    x[i + 1L] <- phi[i] * x[i] + sd_vec[i] * rnorm(1)
  x
}

#' @noRd
sim_markov_alpha2 <- function(t_aug, kappa, tau) {
  n      <- length(t_aug)
  h_vals <- diff(t_aug)
  kappa2 <- kappa * kappa
  c_val  <- r_2(0, kappa, tau, 0)
  r0_du  <- -r_2(0, kappa, tau, 2)
  R0_diag <- diag(c(c_val, r0_du))
  R0     <- chol(R0_diag)
  h_uniq <- unique(h_vals)
  cache  <- lapply(h_uniq, function(h) {
    phi <- exp(-kappa * h); kh <- kappa * h
    A   <- phi * matrix(c(1 + kh, h, -kappa2 * h, 1 - kh), 2, 2, byrow = TRUE)
    Ch  <- matrix(c( r_2(h, kappa, tau, 0), -r_2(h, kappa, tau, 1),
                     r_2(h, kappa, tau, 1), -r_2(h, kappa, tau, 2)),
                  2, 2, byrow = TRUE)
    list(A = A, R_om = chol(R0_diag - A %*% t(Ch)))
  })
  idx     <- match(h_vals, h_uniq)
  X       <- matrix(0, 2, n)
  X[, 1]  <- t(R0) %*% rnorm(2)
  for (i in seq_len(n - 1L)) {
    mk       <- cache[[idx[i]]]
    X[, i+1] <- mk$A %*% X[, i] + t(mk$R_om) %*% rnorm(2)
  }
  X
}


# ============================================================
# Per-edge draws
# ============================================================

#' @noRd
draw_edge_direct <- function(kappa, tau, b_e, l_e, t_abs, alpha) {
  m <- length(t_abs)
  if (m == 0L) return(numeric(0))
  br <- if (alpha == 1L) bridge_alpha1(t_abs, l_e, kappa, tau)
        else             bridge_alpha2(t_abs, l_e, kappa, tau)
  mu <- as.vector(br$S %*% b_e)
  R  <- tryCatch(chol(br$Sigma_star), error = function(e) NULL)
  if (!is.null(R)) {
    return(mu + as.vector(t(R) %*% rnorm(m)))
  }
  # Cholesky failed: Sigma* is near-singular due to floating-point error.
  # Fall back to eigendecomposition, clamping small negative eigenvalues to zero.
  Sigma_star <- 0.5 * (br$Sigma_star + t(br$Sigma_star))
  eig <- eigen(Sigma_star, symmetric = TRUE)
  tol <- max(1e-14, 1e-6 * max(abs(eig$values)))
  if (min(eig$values) < -tol) {
    warning("Bridge covariance has materially negative eigenvalues")
  }
  lambda <- sqrt(pmax(eig$values, 0))
  mu + as.vector(eig$vectors %*% (lambda * rnorm(m)))
}

#' @noRd
bridge_s_alpha1 <- function(t_abs, l_e, kappa, tau) {
  t_ends   <- c(0, l_e)
  Sigma_bb <- r_1(outer(t_ends, t_ends, `-`), kappa, tau)
  Sigma_xb <- r_1(outer(t_abs,  t_ends, `-`), kappa, tau)
  Sigma_xb %*% solve(Sigma_bb)
}

bridge_s_alpha2 <- function(t_abs, l_e, kappa, tau) {
  m      <- length(t_abs)
  t_ends <- c(0, l_e)
  Sigma_bb <- matrix(0, 4, 4)
  D_vv     <- outer(t_ends, t_ends, `-`)
  Sigma_bb[c(1,3), c(1,3)] <-  r_2(D_vv,  kappa, tau, 0)
  Sigma_bb[c(2,4), c(2,4)] <- -r_2(abs(outer(t_ends, t_ends, `-`)), kappa, tau, 2)
  Sigma_bb[c(2,4), c(1,3)] <-  r_2(outer(t_ends, t_ends, `-`), kappa, tau, 1)
  Sigma_bb[c(1,3), c(2,4)] <-  t(Sigma_bb[c(2,4), c(1,3)])
  Sigma_xb <- matrix(0, m, 4)
  D_xv     <- outer(t_abs, t_ends, `-`)
  Sigma_xb[, c(1,3)] <-  r_2(D_xv,  kappa, tau, 0)
  Sigma_xb[, c(2,4)] <-  r_2(-D_xv, kappa, tau, 1)
  Sigma_xb %*% solve(Sigma_bb)
}

#' @noRd
draw_edge_kriging <- function(kappa, tau, b_e, l_e, t_abs, alpha) {
  m <- length(t_abs)
  if (m == 0L) return(numeric(0))
  t_aug  <- c(0, sort(t_abs), l_e)
  n_aug  <- length(t_aug)
  idx_int <- seq(2L, m + 1L)
  if (alpha == 1L) {
    x_aug  <- sim_markov_alpha1(t_aug, kappa, tau)
    x_int  <- x_aug[idx_int]
    x_ends <- x_aug[c(1L, n_aug)]
    S      <- bridge_s_alpha1(t_abs, l_e, kappa, tau)
    x_int  + as.vector(S %*% (b_e - x_ends))
  } else {
    X_aug  <- sim_markov_alpha2(t_aug, kappa, tau)
    x_int  <- X_aug[1, idx_int]
    x_ends_state <- c(X_aug[1, 1], X_aug[2, 1],
                      X_aug[1, n_aug], X_aug[2, n_aug])
    S      <- bridge_s_alpha2(t_abs, l_e, kappa, tau)
    x_int  + as.vector(S %*% (b_e - x_ends_state))
  }
}


# ============================================================
# Extended method: add PtE locations as vertices, single sparse Cholesky draw
# ============================================================

#' @noRd
.simulate_extended_wm <- function(graph, kappa, tau, alpha, BC, PtE, nsim = 1L) {
  order_PtE <- order(PtE[, 1], PtE[, 2])
  n_pts     <- nrow(PtE)

  # Build extended graph with PtE promoted to vertices (once, amortized over nsim)
  g_ext   <- graph$get_initial_graph()
  nV_orig <- g_ext$nV
  df_obs  <- data.frame(
    y                = rep(NA_real_, n_pts),
    edge_number      = PtE[, 1],
    distance_on_edge = PtE[, 2]
  )
  g_ext$add_observations(data = df_obs, normalized = TRUE,
                         suppress_warnings = TRUE, verbose = 0)
  g_ext$observation_to_vertex()
  nV_ext  <- g_ext$nV
  idx_out <- seq.int(nV_orig + 1L, nV_ext)

  na_out <- function() if (nsim == 1L) rep(NA_real_, n_pts) else matrix(NA_real_, n_pts, nsim)

  if (alpha == 1L) {
    Q_ext <- spde_precision(kappa = kappa, tau = tau, alpha = 1L,
                            graph = g_ext, BC = BC)
    R_ext <- tryCatch(
      Matrix::Cholesky(Q_ext, LDL = FALSE, perm = TRUE),
      error = function(e) NULL
    )
    if (is.null(R_ext)) return(na_out())
    draw_one <- function() {
      u_ext <- as.vector(Matrix::solve(R_ext,
                         Matrix::solve(R_ext, rnorm(nV_ext), system = "Lt"),
                         system = "Pt"))
      u_ext[idx_out]
    }
  } else {
    # Pre-compute CoB transformation and Cholesky for alpha=2
    if (is.null(g_ext$CoB) || g_ext$CoB$alpha != 2) g_ext$buildC(2)
    n_con  <- nrow(g_ext$CoB$U)
    Q_ext  <- spde_precision(kappa = kappa, tau = tau, alpha = 2L,
                             graph = g_ext, BC = BC)
    Qtile  <- (g_ext$CoB$T) %*% Q_ext %*% t(g_ext$CoB$T)
    if (n_con > 0L) Qtile <- Qtile[-seq_len(n_con), -seq_len(n_con)]
    R_ext  <- tryCatch(
      Matrix::Cholesky(Matrix::forceSymmetric(Qtile), LDL = FALSE, perm = TRUE),
      error = function(e) NULL
    )
    if (is.null(R_ext)) return(na_out())
    n_draw <- 4L * g_ext$nE - n_con
    CoB_T  <- g_ext$CoB$T
    VtE    <- g_ext$VtEfirst()
    draw_one <- function() {
      V0    <- as.vector(Matrix::solve(R_ext,
                         Matrix::solve(R_ext, rnorm(n_draw), system = "Lt"),
                         system = "Pt"))
      u_ext <- as.vector(t(CoB_T) %*% c(rep(0, n_con), V0))
      vapply(idx_out, function(v) {
        e   <- VtE[v, 1L]
        pos <- VtE[v, 2L]
        u_ext[4L * (e - 1L) + if (pos == 0L) 1L else 3L]
      }, numeric(1L))
    }
  }

  if (nsim == 1L) {
    u_out            <- numeric(n_pts)
    u_out[order_PtE] <- draw_one()
    return(u_out)
  }

  # Multiple samples: graph and Cholesky built once above, only re-solve per draw
  mat <- matrix(0.0, n_pts, nsim)
  for (s in seq_len(nsim)) {
    u_out            <- numeric(n_pts)
    u_out[order_PtE] <- draw_one()
    mat[, s]         <- u_out
  }
  mat
}


# ============================================================
# Global vertex state
# ============================================================

#' @noRd
.draw_vertex_state_wm <- function(graph, kappa, tau, alpha, BC) {
  if (alpha == 1L) {
    Q <- spde_precision(kappa = kappa, tau = tau, alpha = 1,
                        graph = graph, BC = BC)
    R <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)
    return(as.vector(Matrix::solve(R,
                     Matrix::solve(R, rnorm(graph$nV), system = "Lt"),
                     system = "Pt")))
  }
  if (is.null(graph$CoB) || graph$CoB$alpha != 2) graph$buildC(2)
  n_con <- nrow(graph$CoB$U)
  Q     <- spde_precision(kappa = kappa, tau = tau, alpha = 2,
                          graph = graph, BC = BC)
  Qtile <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
  if (n_con > 0L) Qtile <- Qtile[-seq_len(n_con), -seq_len(n_con)]
  R <- Matrix::Cholesky(Matrix::forceSymmetric(Qtile), LDL = FALSE, perm = TRUE)
  V0 <- as.vector(Matrix::solve(R,
                  Matrix::solve(R, rnorm(4L * graph$nE - n_con), system = "Lt"),
                  system = "Pt"))
  as.vector(t(graph$CoB$T) %*% c(rep(0, n_con), V0))
}


# ============================================================
# Public: simulate.metric_graph
# ============================================================

#' Simulate a Whittle-Matérn field on a metric graph
#'
#' Draws unconditional (prior) samples of a Whittle-Matérn field on a metric
#' graph using one of two exact algorithms from Section 6.6 of the paper.
#'
#' @details
#' **Algorithm A ("direct")**: interior values on each edge are drawn from the
#' bridge conditional \eqn{\mathcal{N}(S_e b_e, \Sigma^*_e)} via a dense
#' Cholesky of the \eqn{m \times m} bridge covariance.  Cost: \eqn{O(m^3)}
#' per edge.  Preferred when \eqn{m} is small (roughly \eqn{m \lesssim 100}).
#'
#' **Algorithm B ("kriging")**: an unconditioned Matérn process is simulated on
#' each edge by a linear Markov recursion, then the bridge boundary data are
#' imposed by a kriging correction (Theorem 2 / Proposition 4).  Cost:
#' \eqn{O(m)} per edge.  Preferred for large \eqn{m}.
#'
#' Both algorithms produce draws from exactly the same finite-dimensional
#' distribution; they differ only in computational cost.
#'
#' The two-step procedure is:
#' \enumerate{
#'   \item Draw the vertex (boundary) state from the vertex precision via a
#'     sparse Cholesky.
#'   \item For each edge (independently given the vertex state), draw the
#'     interior values using the chosen algorithm.
#' }
#'
#' When `parallel = TRUE` the edge loop is distributed over a
#' `doParallel`/`foreach` cluster.  Reproducibility requires passing a `seed`
#' so that per-edge seeds are derived deterministically.
#'
#' @param object A `metric_graph` object.
#' @param nsim Number of samples.  Returns a vector (nsim = 1) or matrix
#'   (nsim > 1) of sampled values.
#' @param seed Integer seed for `set.seed()` before the vertex draw.  Pass
#'   `NULL` (default) to use the current RNG state.
#' @param alpha Smoothness parameter (1 or 2).
#' @param method Simulation algorithm: `"direct"` (Method A, \eqn{O(m^3)}),
#'   `"kriging"` (Method B, \eqn{O(m)}), or `"extended"` (add simulation
#'   locations as graph vertices, single sparse Cholesky draw).
#' @param impl Implementation: `"cpp"` (default, Rcpp/Eigen) or `"R"` (pure
#'   R reference, for testing/comparison).  Ignored for `method = "extended"`.
#' @param kappa Range parameter.
#' @param tau Precision parameter.
#' @param range Practical correlation range (alternative to `kappa`/`tau`).
#' @param sigma Marginal standard deviation (alternative to `kappa`/`tau`).
#' @param PtE Matrix with columns `(edge_number, normalised_position)`.
#'   Required when `type = "manual"` (the default).
#' @param type Location specification: `"manual"` (use `PtE`), `"mesh"` (mesh
#'   nodes), or `"obs"` (observation locations).
#' @param BC Boundary condition for degree-1 vertices: 1 = stationary
#'   (default), 0 = Neumann.
#' @param parallel Logical.  Use a parallel cluster for the edge loop?
#' @param n_cores Number of parallel workers.  Ignored if `cluster` is
#'   supplied.  Defaults to `detectCores() - 1`.
#' @param cluster An already-registered `doParallel` cluster.  If `NULL` and
#'   `parallel = TRUE`, a cluster is created and stopped automatically.
#' @param ... Ignored.
#' @return Numeric vector (nsim = 1) or matrix with nsim columns of field
#'   values at the requested locations, ordered to match `PtE`.
#' @seealso [sample_spde()]
#' @examples
#' # Small square graph, alpha = 1
#' V <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
#' E <- rbind(c(1,2), c(2,3), c(3,4), c(4,1))
#' g <- metric_graph$new(V = V, E = E)
#' t_norm <- seq(0.1, 0.9, by = 0.2)
#' PtE    <- do.call(rbind, lapply(1:4, function(e) cbind(e, t_norm)))
#' u <- simulate(g, alpha = 1, method = "direct", kappa = 1, tau = 1, PtE = PtE)
#' \donttest{
#' # Method B on a larger graph, parallel
#' u2 <- simulate(g, alpha = 2, method = "kriging", range = 1, sigma = 1,
#'                PtE = PtE, parallel = TRUE, n_cores = 2)
#' }
#' @method simulate metric_graph
#' @export
simulate.metric_graph <- function(object, nsim = 1, seed = NULL,
                                   alpha = 1,
                                   method = c("direct", "kriging", "extended"),
                                   impl   = c("cpp", "R"),
                                   kappa, tau, range, sigma,
                                   PtE = NULL, type = "manual",
                                   BC = 1,
                                   parallel = FALSE,
                                   n_cores = NULL,
                                   cluster = NULL,
                                   ...) {
  method <- match.arg(method)
  impl   <- match.arg(impl)
  alpha  <- as.integer(alpha)
  if (!alpha %in% c(1L, 2L)) stop("alpha must be 1 or 2")
  if (!(type %in% c("manual", "mesh", "obs")))
    stop("type must be 'manual', 'mesh', or 'obs'")

  # --- Parameter conversion ---
  has_kappa <- !missing(kappa); has_tau  <- !missing(tau)
  has_range <- !missing(range); has_sigma <- !missing(sigma)
  if ((!has_kappa || !has_tau) && (!has_range || !has_sigma))
    stop("Supply either kappa and tau, or range and sigma.")
  if (has_range && has_sigma) {
    nu    <- alpha - 0.5
    kappa <- sqrt(8 * nu) / range
    tau   <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
                               (4 * pi)^0.5 * gamma(nu + 0.5)))
  }

  # --- Locations ---
  if (type == "manual") {
    if (is.null(PtE)) stop("PtE must be supplied when type = 'manual'")
    PtE <- as.matrix(PtE)
  } else if (type == "mesh") {
    if (is.null(object$mesh)) stop("Graph has no mesh; build one first.")
    PtE <- object$mesh$PtE
  } else {
    if (is.null(object$PtE)) stop("Graph has no observations.")
    PtE <- object$PtE
  }

  # --- Multiple samples ---
  if (nsim > 1L) {
    if (method == "extended") {
      # Build extended graph and Cholesky once, amortize over all samples
      if (!is.null(seed)) set.seed(seed)
      return(.simulate_extended_wm(object, kappa, tau, alpha, BC, PtE, nsim = nsim))
    }
    samps <- lapply(seq_len(nsim), function(s) {
      sd_s <- if (!is.null(seed)) seed + s - 1L else NULL
      simulate.metric_graph(object, nsim = 1L, seed = sd_s,
                             alpha = alpha, method = method, impl = impl,
                             kappa = kappa, tau = tau,
                             PtE = PtE, type = "manual",
                             BC = BC, parallel = FALSE)
    })
    return(do.call(cbind, samps))
  }

  if (!is.null(seed)) set.seed(seed)

  # --- Extended method: bypass per-edge loop ---
  if (method == "extended") {
    return(.simulate_extended_wm(object, kappa, tau, alpha, BC, PtE))
  }

  # --- Cluster setup ---
  cl_created <- FALSE
  if (parallel) {
    if (is.null(cluster)) {
      if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
      cluster    <- parallel::makeCluster(n_cores)
      cl_created <- TRUE
    }
    doParallel::registerDoParallel(cluster)
    on.exit(if (cl_created) parallel::stopCluster(cluster), add = TRUE)
    # Export package functions once per cluster lifetime; skip on subsequent calls
    if (!isTRUE(attr(cluster, "MetricGraph_exported"))) {
      pkg_env  <- environment()
      base_fns <- c("draw_edge_direct", "draw_edge_kriging",
                    "bridge_alpha1", "bridge_alpha2",
                    "sim_markov_alpha1", "sim_markov_alpha2",
                    "r_1", "r_2")
      cpp_fns  <- c("draw_edge_direct_cpp", "draw_edge_kriging_cpp")
      parallel::clusterExport(cluster,
        varlist = c(base_fns, cpp_fns),
        envir = pkg_env)
      attr(cluster, "MetricGraph_exported") <- TRUE
    }
  }

  # --- 1. Vertex draw ---
  b_all <- .draw_vertex_state_wm(object, kappa, tau, alpha, BC)

  # --- 2. Edge-loop setup ---
  order_PtE   <- order(PtE[, 1], PtE[, 2])
  ordered_PtE <- PtE[order_PtE, , drop = FALSE]
  inds_e      <- unique(ordered_PtE[, 1])
  t_by_edge   <- split(ordered_PtE[, 2], ordered_PtE[, 1])

  if (!is.null(seed)) {
    set.seed(seed + 1L)
    edge_seeds <- sample.int(.Machine$integer.max, object$nE)
  } else {
    edge_seeds <- NULL
  }

  edge_data <- lapply(seq_along(inds_e), function(k) {
    i <- inds_e[k]
    t_k <- t_by_edge[[as.character(i)]]
    list(
      t_abs = t_k * object$edge_lengths[i],
      l_e   = object$edge_lengths[i],
      b_e   = if (alpha == 1L) b_all[object$E[i, ]]
              else b_all[4L * (i - 1L) + 1:4],
      seed  = if (!is.null(edge_seeds)) edge_seeds[i] else NULL
    )
  })

  # --- 3. Edge loop ---
  .draw_edge <- function(d) {
    if (!is.null(d$seed)) set.seed(d$seed)
    if (impl == "cpp") {
      if (method == "direct")
        draw_edge_direct_cpp(kappa, tau, d$b_e, d$l_e, d$t_abs, alpha)
      else
        draw_edge_kriging_cpp(kappa, tau, d$b_e, d$l_e, d$t_abs, alpha)
    } else {
      if (method == "direct")
        draw_edge_direct(kappa, tau, d$b_e, d$l_e, d$t_abs, alpha)
      else
        draw_edge_kriging(kappa, tau, d$b_e, d$l_e, d$t_abs, alpha)
    }
  }

  if (parallel) {
    # Chunk edges into one batch per worker to eliminate per-task dispatch overhead
    n_par  <- length(cluster)
    n_par  <- min(n_par, length(edge_data))
    chunks <- split(seq_along(edge_data),
                    ceiling(seq_along(edge_data) / ceiling(length(edge_data) / n_par)))
    u_chunks <- foreach::foreach(
      ed_ch = lapply(chunks, function(i) edge_data[i])
    ) %dopar% {
      lapply(ed_ch, .draw_edge)
    }
    u_list <- unlist(u_chunks, recursive = FALSE)
  } else {
    u_list <- lapply(edge_data, .draw_edge)
  }

  # --- 4. Reassemble in original order ---
  u_out              <- numeric(nrow(PtE))
  u_out[order_PtE]  <- unlist(u_list, use.names = FALSE)
  u_out
}


# ============================================================
# Public: simulate_parallel (thin wrapper)
# ============================================================

#' Parallel simulation of a Whittle-Matérn field on a metric graph
#'
#' Calls `simulate.metric_graph()` with `parallel = TRUE`.  The edge loop is
#' distributed over a `doParallel`/`foreach` cluster.
#'
#' @param graph A `metric_graph` object.
#' @param ... Further arguments passed to `simulate.metric_graph()`, including
#'   `alpha`, `method`, `kappa`/`tau` or `range`/`sigma`, `PtE`, `type`,
#'   `BC`, `nsim`, and `seed`.
#' @param n_cores Number of parallel workers.
#' @param cluster An already-registered `doParallel` cluster.
#' @return Numeric vector (nsim = 1) or matrix with nsim columns.
#' @seealso [sample_spde()]
#' @examples
#' \donttest{
#' V <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
#' E <- rbind(c(1,2), c(2,3), c(3,4), c(4,1))
#' g <- metric_graph$new(V = V, E = E)
#' t_norm <- seq(0.1, 0.9, by = 0.2)
#' PtE    <- do.call(rbind, lapply(1:4, function(e) cbind(e, t_norm)))
#' u <- simulate_parallel(g, alpha = 1, method = "kriging",
#'                        kappa = 1, tau = 1, PtE = PtE, n_cores = 2)
#' }
#' @export
simulate_parallel <- function(graph, ..., n_cores = NULL, cluster = NULL) {
  simulate.metric_graph(graph, ...,
                        parallel = TRUE,
                        n_cores  = n_cores,
                        cluster  = cluster)
}
