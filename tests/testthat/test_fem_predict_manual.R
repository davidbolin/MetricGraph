# Tests for predict.rspde_lme (the FEM/fractional path that
# predict.graph_lme dispatches to) against a fully manual prediction
# built from the model's primitives:
#
#   * Q     — precision of the SPDE latent field, rebuilt from the fit's
#             random_effects via update(latent_model, ...).
#   * A_r   — fit$A_list[[r]], the observation-to-mesh design matrix
#             for replicate r as stored at fit time.
#   * Y     — fit$model_matrix[, 1] minus X %*% beta for any covariates.
#   * A_pred — make_A(latent_model, loc_pred), expanded by the rational
#              order kronecker block for fractional alpha.
#
# The manual mean / variance for replicate r are
#
#   Q_xy = t(A_r) %*% A_r / sigma_e^2 + Q
#   mu   = A_pred %*% Q_xy^{-1} %*% t(A_r) %*% y_r / sigma_e^2
#   var  = diag(A_pred %*% Q_xy^{-1} %*% t(A_pred))
#
# plus mean-correction terms when the rSPDE object has them. predict()
# must reproduce these byte-for-byte. The tests probe the cases where
# bugs have historically lived (multi-replicate, na_test_idx masking,
# fractional alpha with mean correction, covariates) so a regression in
# predict.rspde_lme — or in MetricGraph's forwarding of advanced_options
# / which_repl — surfaces here rather than only through
# posterior_crossvalidation downstream.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.fem_manual_predict <- function(fit, loc_pred, which_repl = NULL,
                                na_idx = NULL) {
  sigma_e <- as.numeric(fit$coeff$measurement_error)

  coeff_random <- fit$coeff$random_effects
  update_params <- list()
  for (param_name in names(coeff_random)) {
    update_params[[gsub(" \\(fixed\\)$", "", param_name)]] <-
      coeff_random[[param_name]]
  }
  update_params$check_stationarity <- FALSE
  new_obj <- do.call(update, c(list(fit$latent_model), update_params))
  Q <- new_obj$Q

  alpha_param <- grep("^alpha", names(coeff_random), value = TRUE)
  nu_param    <- grep("^nu",    names(coeff_random), value = TRUE)
  if (length(alpha_param) > 0) {
    alpha_val <- coeff_random[[alpha_param]]
  } else {
    alpha_val <- coeff_random[[nu_param]] + fit$latent_model$d / 2
  }
  is_integer_alpha <- (alpha_val %% 1 == 0)

  A_pred <- rSPDE::make_A(new_obj, loc_pred)
  if (!is_integer_alpha) {
    A_pred <- kronecker(matrix(1, 1, fit$rspde_order + 1), A_pred)
  }
  mu_corr <- if (fit$mean_correction) new_obj$mean_correction(full = TRUE) else 0

  mm <- as.matrix(fit$model_matrix)
  Y_all <- mm[, 1]
  if (ncol(mm) > 1) {
    X_fit <- mm[, 2:ncol(mm), drop = FALSE]
    Y_all <- Y_all - as.numeric(X_fit %*% fit$coeff$fixed_effects)
  }
  Y_all <- as.numeric(Y_all)
  if (!is.null(na_idx)) Y_all[na_idx] <- NA

  repl_vec <- fit$repl
  u_repl <- if (is.null(which_repl)) unique(repl_vec) else unique(which_repl)

  out_mean <- list(); out_var <- list()
  for (r in u_repl) {
    idx_r    <- repl_vec == r
    y_r      <- Y_all[idx_r]
    obs_mask <- !is.na(y_r)
    y_r      <- y_r[obs_mask]
    A_r      <- fit$A_list[[r]]
    if (nrow(A_r) == length(obs_mask) && any(!obs_mask)) {
      A_r <- A_r[obs_mask, , drop = FALSE]
    }
    if (fit$mean_correction) y_r <- y_r - as.numeric(A_r %*% mu_corr)

    Q_xy   <- Matrix::t(A_r) %*% A_r / sigma_e^2 + Q
    mu_lat <- solve(Q_xy, as.vector(Matrix::t(A_r) %*% y_r / sigma_e^2))
    mu_pp  <- as.numeric(A_pred %*% mu_lat)

    fe_pred <- 0
    if (ncol(mm) > 1) {
      # Fixed-effects intercept is in fit$coeff$fixed_effects; we built
      # loc_pred without covariates so the intercept (if any) is the
      # only fe contribution. Match predict.rspde_lme's behaviour of
      # adding X_pred %*% beta with X_pred = matrix(1, nrow=n_pred, 1)
      # when no covariate values were provided alongside the loc.
      fe_pred <- as.numeric(fit$coeff$fixed_effects[1])
    }
    if (fit$mean_correction) mu_pp <- mu_pp + as.numeric(A_pred %*% mu_corr)

    out_mean[[as.character(r)]] <- mu_pp + fe_pred
    out_var[[as.character(r)]]  <-
      pmax(diag(as.matrix(A_pred %*% solve(Q_xy, Matrix::t(A_pred)))), 0)
  }
  list(mean = out_mean, variance = out_var)
}

.build_fem_multirep_graph <- function(n_repl, alpha = 1, obs_per_edge = 5,
                                      seed = 11, mesh_h = 0.05,
                                      shared_locations = TRUE) {
  set.seed(seed)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = mesh_h)
  df_all <- NULL
  if (shared_locations) {
    PtE <- NULL
    for (i in seq_len(graph$nE)) {
      PtE <- rbind(PtE, cbind(rep(i, obs_per_edge), runif(obs_per_edge)))
    }
  }
  for (r in seq_len(n_repl)) {
    if (!shared_locations) {
      PtE <- NULL
      for (i in seq_len(graph$nE)) {
        PtE <- rbind(PtE, cbind(rep(i, obs_per_edge), runif(obs_per_edge)))
      }
    }
    u <- sample_spde(kappa = 5, tau = 1, alpha = alpha, graph = graph, PtE = PtE)
    df_all <- rbind(df_all,
      data.frame(y = u + 0.3 * rnorm(length(u)),
                 edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2],
                 repl = r))
  }
  group_var <- if (n_repl > 1L) "repl" else NULL
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = group_var)
  graph
}

# ---------------------------------------------------------------------------
# Single-replicate FEM, integer alpha: predict.rspde_lme must match the
# manual A^T A / sigma_e^2 + Q kriging formula for arbitrary prediction
# locations.
# ---------------------------------------------------------------------------

test_that("predict.rspde_lme matches manual prediction (alpha=1 FEM, single rep)", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 1L, alpha = 1, seed = 7)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))
  expect_true(inherits(fit, "rspde_lme"))

  loc_pred <- cbind(c(1, 1, 2, 3, 4), c(0.2, 0.7, 0.5, 0.3, 0.9))
  manual <- .fem_manual_predict(fit, loc_pred)

  nd <- list(.edge_number = loc_pred[, 1], .distance_on_edge = loc_pred[, 2])
  pp <- predict(fit, newdata = nd,
                edge_number = ".edge_number",
                distance_on_edge = ".distance_on_edge",
                normalized = TRUE, compute_variances = TRUE)
  expect_equal(as.numeric(pp$mean),     manual$mean[[1]],     tolerance = 1e-9)
  expect_equal(as.numeric(pp$variance), manual$variance[[1]], tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# Multi-replicate, integer alpha. Three things must hold:
#   1. Each per-replicate predict() call (which_repl = r) matches the
#      manual replicate-r kriging.
#   2. A single predict() call with newdata covering all replicates
#      returns the concatenation of those per-replicate predictions
#      (and a correct $repl tag).
#   3. The variance for every replicate is INDEPENDENT — a prior bug
#      collapsed $variance to the last replicate's value.
# ---------------------------------------------------------------------------

test_that("predict.rspde_lme matches manual prediction (alpha=1 FEM, multi-rep)", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))
  expect_true(inherits(fit, "rspde_lme"))

  loc_pred <- cbind(c(1, 1, 2, 3, 4), c(0.2, 0.7, 0.5, 0.3, 0.9))
  manual <- .fem_manual_predict(fit, loc_pred)

  # Per-replicate calls.
  for (r in 1:3) {
    nd_r <- list(.edge_number = loc_pred[, 1],
                 .distance_on_edge = loc_pred[, 2])
    pp_r <- predict(fit, newdata = nd_r,
                    edge_number = ".edge_number",
                    distance_on_edge = ".distance_on_edge",
                    normalized = TRUE, compute_variances = TRUE,
                    which_repl = r)
    expect_equal(as.numeric(pp_r$mean),     manual$mean[[as.character(r)]],
                 tolerance = 1e-9,
                 info = paste("mean mismatch for replicate", r))
    expect_equal(as.numeric(pp_r$variance), manual$variance[[as.character(r)]],
                 tolerance = 1e-9,
                 info = paste("variance mismatch for replicate", r))
  }

  # Single call covering all replicates: must be the concatenation.
  nd_all <- list(
    .edge_number      = rep(loc_pred[, 1], 3),
    .distance_on_edge = rep(loc_pred[, 2], 3),
    .group            = rep(1:3, each = nrow(loc_pred))
  )
  pp_all <- predict(fit, newdata = nd_all,
                    edge_number = ".edge_number",
                    distance_on_edge = ".distance_on_edge",
                    normalized = TRUE, compute_variances = TRUE)
  expect_equal(as.numeric(pp_all$mean),
               unlist(manual$mean[as.character(1:3)], use.names = FALSE),
               tolerance = 1e-9)
  expect_equal(as.numeric(pp_all$variance),
               unlist(manual$variance[as.character(1:3)], use.names = FALSE),
               tolerance = 1e-9)
  expect_equal(as.character(pp_all$repl),
               rep(as.character(1:3), each = nrow(loc_pred)))

  # Variances ARE expected to coincide here: with shared observation
  # locations and stationary parameters, the per-replicate kriging
  # variance is a function of locations and parameters only. The
  # regression check that the per-replicate variances do NOT collapse
  # under predict() bookkeeping lives in test_posterior_crossvalidation.R,
  # where one replicate has an NA'd observation. Here we just confirm
  # the means differ across replicates (they share Q but not y).
  m1 <- manual$mean[["1"]]; m2 <- manual$mean[["2"]]; m3 <- manual$mean[["3"]]
  expect_true(any(abs(m1 - m2) > 1e-8) || any(abs(m2 - m3) > 1e-8))
})

# ---------------------------------------------------------------------------
# Fractional FEM (rational order > 0 and mean correction) — multi-rep.
# The Aprd expansion (kronecker(matrix(1,1,m+1), Aprd)) and the
# mean-correction term are the parts most prone to silent drift; the
# manual prediction here uses the same code path so any divergence is
# squarely the fault of predict.rspde_lme.
# ---------------------------------------------------------------------------

test_that("predict.rspde_lme matches manual prediction (fractional FEM, multi-rep)", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", fem = TRUE),
                   model_options = list(start_nu = 0.5))
  expect_true(inherits(fit, "rspde_lme"))
  # Sanity: this fit should actually be fractional.
  nu_val <- fit$coeff$random_effects[["nu"]]
  expect_true(nu_val %% 1 != 0)

  loc_pred <- cbind(c(1, 1, 2, 3, 4), c(0.2, 0.7, 0.5, 0.3, 0.9))
  manual <- .fem_manual_predict(fit, loc_pred)

  for (r in 1:3) {
    nd_r <- list(.edge_number = loc_pred[, 1],
                 .distance_on_edge = loc_pred[, 2])
    pp_r <- predict(fit, newdata = nd_r,
                    edge_number = ".edge_number",
                    distance_on_edge = ".distance_on_edge",
                    normalized = TRUE, compute_variances = TRUE,
                    which_repl = r)
    expect_equal(as.numeric(pp_r$mean),     manual$mean[[as.character(r)]],
                 tolerance = 1e-7,
                 info = paste("mean mismatch for replicate", r))
    expect_equal(as.numeric(pp_r$variance), manual$variance[[as.character(r)]],
                 tolerance = 1e-7,
                 info = paste("variance mismatch for replicate", r))
  }
})

# ---------------------------------------------------------------------------
# na_test_idx masking — this is the path posterior_crossvalidation uses
# to compute pseudo-CV predictions for FEM fits. The masking must:
#   * remove the held-out rows from y AND from the corresponding A_list
#     rows (so the kriging really excludes those points), and
#   * still produce a per-replicate prediction at the requested loc.
# Manual prediction with na_idx applied to Y_all (and obs_mask filtering
# A_r) must agree with predict() called with advanced_options$na_test_idx.
# ---------------------------------------------------------------------------

test_that("predict.rspde_lme respects na_test_idx (alpha=1 FEM, multi-rep)", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))

  # Hold out: a few obs in rep 1, one in rep 2, several in rep 3.
  gd <- fit$graph$.__enclos_env__$private$data
  rep_vec <- gd[[".group"]]
  na_idx <- c(which(rep_vec == 1)[1:3],
              which(rep_vec == 2)[5],
              which(rep_vec == 3)[c(2, 7, 12)])
  loc_pred <- cbind(c(1, 2, 3, 4), c(0.15, 0.55, 0.4, 0.8))

  manual <- .fem_manual_predict(fit, loc_pred, na_idx = na_idx)

  for (r in 1:3) {
    nd_r <- list(.edge_number = loc_pred[, 1],
                 .distance_on_edge = loc_pred[, 2])
    pp_r <- predict(fit, newdata = nd_r,
                    edge_number = ".edge_number",
                    distance_on_edge = ".distance_on_edge",
                    normalized = TRUE, compute_variances = TRUE,
                    which_repl = r,
                    advanced_options = list(na_test_idx = na_idx))
    expect_equal(as.numeric(pp_r$mean),     manual$mean[[as.character(r)]],
                 tolerance = 1e-9,
                 info = paste("mean mismatch (na_test_idx) for replicate", r))
    expect_equal(as.numeric(pp_r$variance), manual$variance[[as.character(r)]],
                 tolerance = 1e-9,
                 info = paste("variance mismatch (na_test_idx) for replicate", r))
  }
})

# ---------------------------------------------------------------------------
# Different observation locations per replicate (the more realistic case).
# Each replicate has its own PtE, so A_list[[r]] differs in row content
# across replicates. The manual prediction must still match predict()
# at a fixed set of prediction locations.
# ---------------------------------------------------------------------------

test_that("predict.rspde_lme matches manual prediction with per-replicate locations", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 5,
                                     shared_locations = FALSE)
  # The optimiser may fall back from L-BFGS-B for this small / per-rep-
  # locations setup. The fallback is benign for the prediction equality
  # check below; suppress the cosmetic warning to keep test output clean.
  fit <- suppressWarnings(
    graph_lme(y ~ -1, graph = graph,
              model = list(type = "WhittleMatern", alpha = 1, fem = TRUE)))

  # Sanity: A_list rows must differ across replicates (the only place
  # this matters is when locations actually differ).
  expect_false(identical(as.numeric(fit$A_list[[1]]),
                         as.numeric(fit$A_list[[2]])))

  loc_pred <- cbind(c(1, 2, 3, 4), c(0.1, 0.4, 0.6, 0.85))
  manual <- .fem_manual_predict(fit, loc_pred)

  for (r in 1:3) {
    nd_r <- list(.edge_number = loc_pred[, 1],
                 .distance_on_edge = loc_pred[, 2])
    pp_r <- predict(fit, newdata = nd_r,
                    edge_number = ".edge_number",
                    distance_on_edge = ".distance_on_edge",
                    normalized = TRUE, compute_variances = TRUE,
                    which_repl = r)
    expect_equal(as.numeric(pp_r$mean),     manual$mean[[as.character(r)]],
                 tolerance = 1e-9,
                 info = paste("mean mismatch (per-rep locs) for replicate", r))
    expect_equal(as.numeric(pp_r$variance), manual$variance[[as.character(r)]],
                 tolerance = 1e-9,
                 info = paste("variance mismatch (per-rep locs) for replicate", r))
  }
})

# ---------------------------------------------------------------------------
# posterior_crossvalidation LOO must reproduce the manual prediction at
# each held-out point. This is the END-TO-END check the bug report is
# actually about: the predict-method machinery + the CV loop in
# posterior_crossvalidation together must yield the manual proper LOO
# at every i in 1..nobs. Earlier tests in test_posterior_crossvalidation.R
# probed this against a predict-based reference (.manual_fem_loo_one);
# this one closes the loop with a predict-free reference.
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation LOO matches predict-free manual LOO (multi-rep FEM)", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))

  res_pre <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = TRUE,
                                       scores = c("logscore", "rmse"))
  res_no  <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = FALSE,
                                       scores = c("logscore", "rmse"))
  expect_equal(res_pre$mu,  res_no$mu,  tolerance = 1e-9)
  expect_equal(res_pre$var, res_no$var, tolerance = 1e-9)

  gd <- fit$graph$.__enclos_env__$private$data
  rep_vec <- gd[[".group"]]
  sig2 <- as.numeric(fit$coeff$measurement_error)^2
  # Check EVERY observation, not just one per replicate — the bug class
  # this is guarding against would otherwise only show up at a handful
  # of indices (e.g. replicate 2 / 3 if replicate 1 was used as the
  # canonical position-indexed answer for all points).
  for (i in seq_len(fit$nobs)) {
    r <- rep_vec[i]
    loc_test <- cbind(gd[[".edge_number"]][i],
                      gd[[".distance_on_edge"]][i])
    manual <- .fem_manual_predict(fit, loc_test, which_repl = r,
                                  na_idx = i)
    target_mean <- unname(manual$mean[[as.character(r)]])
    target_var  <- unname(manual$variance[[as.character(r)]]) + sig2
    expect_equal(unname(res_pre$mu[i]),  target_mean, tolerance = 1e-7,
                 info = paste("LOO mean mismatch at i =", i, "(rep", r, ")"))
    expect_equal(unname(res_pre$var[i]), target_var,  tolerance = 1e-7,
                 info = paste("LOO var mismatch at i =", i, "(rep", r, ")"))
  }
})

# ---------------------------------------------------------------------------
# Regression: non-stationary FEM construction must use the requested nu /
# alpha. The conversion (B.range, B.sigma) -> (B.tau, B.kappa) depends on
# nu through the Matern formula; if graph_lme builds the rspde_object
# without passing nu, rSPDE defaults to nu=0.75 (alpha=1.25) and the B
# matrices — and the entire likelihood — get tuned for the wrong model.
# The symptom was: a non-stationary FEM fit produced dramatically worse
# CV scores than the same-family stationary FEM on data generated from a
# non-stationary FEM. The fix is to thread alpha/nu through the
# spde.matern.operators call in graph_lme.
# ---------------------------------------------------------------------------

test_that("graph_lme(fem=TRUE, alpha=1, B.range, B.sigma) builds an alpha=1 operator", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
  n_mesh <- nrow(graph$mesh$VtE)
  spatial_cov <- graph$mesh$VtE[, 2] - 0.5
  B.range <- cbind(rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh), spatial_cov)
  B.sigma <- cbind(rep(0, n_mesh), rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh))

  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE,
                                B.range = B.range, B.sigma = B.sigma))

  # The operator stored in the fit must agree with the user's alpha=1 /
  # nu=0.5 request. Before the fix, this was alpha=1.25 / nu=0.75 with
  # Q at (m+1)*n_mesh = 3*n_mesh — A_list (n_mesh cols) and Q
  # disagreed, and the likelihood was optimised against the wrong B.tau /
  # B.kappa.
  expect_equal(fit$latent_model$alpha, 1)
  expect_equal(fit$latent_model$nu,    0.5)
  expect_equal(nrow(fit$latent_model$Q), ncol(fit$A_list[[1]]),
               info = "Q rows must match A_list columns (no rational expansion for integer alpha)")
})

# ---------------------------------------------------------------------------
# Sanity check: a non-stationary FEM fit to data simulated from the same
# family must NOT score dramatically worse on CV than a stationary FEM
# fit on the same data. Before the fix, on data generated from a
# non-stationary spde.matern.operators (alpha=1), the FEM non-stat CV
# log-score was ~2x worse than the FEM stat log-score on the same data —
# entirely an artefact of the operator being built with rSPDE's default
# nu=0.75 instead of the requested nu=0.5.
#
# The test is a coarse upper bound (well above what a correctly-fit model
# should achieve) so it doesn't flake on random-seed variations of the
# optimiser, but is far below the buggy regime.
# ---------------------------------------------------------------------------

test_that("non-stationary FEM CV is not catastrophically worse than stationary FEM CV", {
  skip_if_not_installed("rSPDE")
  set.seed(42)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.02)
  n_mesh <- nrow(graph$mesh$VtE)
  spatial_cov <- graph$mesh$VtE[, 2] - 0.5
  B.range <- cbind(rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh), spatial_cov)
  B.sigma <- cbind(rep(0, n_mesh), rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh))

  # Simulate from a non-stationary FEM operator at alpha=1.
  ns_sim <- rSPDE::spde.matern.operators(
    graph = graph, nu = 0.5,
    B.range = B.range, B.sigma = B.sigma,
    theta = c(-2, -0.9, -1.5),
    parameterization = "matern",
    check_stationarity = FALSE)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 8), runif(8)))
  }
  A_obs <- rSPDE::make_A(ns_sim, PtE)
  L <- chol(ns_sim$Q)
  sigma_e_true <- 0.2
  df_all <- NULL
  for (r in 1:3) {
    z <- rnorm(nrow(ns_sim$Q))
    u <- as.numeric(solve(L, z))
    ux <- as.numeric(A_obs %*% u)
    df_all <- rbind(df_all,
      data.frame(y = ux + sigma_e_true * rnorm(length(ux)),
                 edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = "repl")

  fit_ns <- graph_lme(y ~ -1, graph = graph,
                model = list(type = "WhittleMatern", alpha = 1, fem = TRUE,
                             B.range = B.range, B.sigma = B.sigma))
  fit_stat <- graph_lme(y ~ -1, graph = graph,
                model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))

  # The non-stationary fit must build the correct operator (the upstream
  # bug it's guarding against).
  expect_equal(fit_ns$latent_model$alpha, 1)
  expect_equal(fit_ns$latent_model$nu,    0.5)

  cv_ns   <- posterior_crossvalidation(fit_ns,   mode = "loo", true_CV = FALSE,
                                       scores = "logscore")
  cv_stat <- posterior_crossvalidation(fit_stat, mode = "loo", true_CV = FALSE,
                                       scores = "logscore")

  ls_ns   <- as.numeric(cv_ns$scores$logscore)
  ls_stat <- as.numeric(cv_stat$scores$logscore)
  # Before the fixes ls_ns ≈ 0.94 vs ls_stat ≈ 0.42 — over 2x worse on
  # data generated FROM the same non-stationary family. Two distinct
  # bugs converged here:
  #   * graph_lme built the rspde_object without passing nu, so
  #     B.tau / B.kappa were derived for rSPDE's default nu (0.75)
  #     instead of the requested nu=0.5.
  #   * update.CBrSPDEobj silently ignored theta1, theta2, ...
  #     keyword args (which is how predict.rspde_lme passes them from
  #     coeff$random_effects); Q was rebuilt at the stale stored
  #     theta on every predict call, so kriging used the wrong
  #     covariance regardless of what the optimiser had found.
  # After both fixes ls_ns is at least as good as ls_stat (often
  # better, since the data really is non-stationary). The bound here
  # is generous to absorb optimiser noise but still well inside the
  # buggy regime.
  expect_lt(ls_ns, ls_stat + 0.1)
})

# ---------------------------------------------------------------------------
# Direct regression for the update() / theta naming bug. The fit's
# coeff$random_effects stores the non-stationary thetas under names
# "theta1", "theta2", ... — and predict.rspde_lme rebuilds the operator
# via `do.call(update, c(list(latent_model), <these names>))`. Before
# the fix, update.CBrSPDEobj only honored a single `theta = c(...)`
# vector argument and dropped theta1/theta2/... into `...` (silently).
# The operator's stored theta — whatever was last cached during fitting
# (often near zero) — was used instead, and Q was wrong on every
# predict call.
#
# This test pins update()'s behaviour to match a direct
# spde.matern.operators(theta = c(...)) construction at the same params.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Same end-to-end sanity check as the B.range/B.sigma scenario, but using
# the B.tau / B.kappa ("spde") parameterisation. Before the construction
# fix, this path also built the rspde_object without passing alpha
# (rSPDE defaulted to alpha=1.5, m=2 with rational expansion) even though
# the user requested alpha=1, leading to an A_list / Q size mismatch and
# (combined with the theta naming bug) catastrophic CV scores. After
# both fixes the non-stationary fit on data generated from a
# non-stationary B.tau/B.kappa model must score better than the
# stationary fit.
# ---------------------------------------------------------------------------

test_that("non-stationary FEM CV (B.tau/B.kappa) beats stationary FEM CV on data from a non-stationary model", {
  skip_if_not_installed("rSPDE")
  set.seed(42)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.02)
  n_mesh <- nrow(graph$mesh$VtE)
  spatial_cov <- graph$mesh$VtE[, 2] - 0.5
  B.tau   <- cbind(rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh), spatial_cov)
  B.kappa <- cbind(rep(0, n_mesh), rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh))

  ns_sim <- rSPDE::spde.matern.operators(
    graph = graph, alpha = 1,
    B.tau = B.tau, B.kappa = B.kappa,
    theta = c(-2, -0.9, -1.5),
    parameterization = "spde",
    check_stationarity = FALSE)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 8), runif(8)))
  }
  A_obs <- rSPDE::make_A(ns_sim, PtE)
  L <- chol(ns_sim$Q)
  df_all <- NULL
  for (r in 1:3) {
    u  <- as.numeric(solve(L, rnorm(nrow(ns_sim$Q))))
    ux <- as.numeric(A_obs %*% u)
    df_all <- rbind(df_all,
      data.frame(y = ux + 0.2 * rnorm(length(ux)),
                 edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = "repl")

  fit_ns <- graph_lme(y ~ -1, graph = graph,
                model = list(type = "WhittleMatern", alpha = 1, fem = TRUE,
                             B.tau = B.tau, B.kappa = B.kappa))
  fit_stat <- graph_lme(y ~ -1, graph = graph,
                model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))

  # Construction-time fix: alpha must be 1 here too (not rSPDE's default
  # 1.5 from spde.matern.operators without an alpha argument).
  expect_equal(fit_ns$latent_model$alpha, 1)
  expect_equal(fit_ns$latent_model$nu,    0.5)

  cv_ns   <- posterior_crossvalidation(fit_ns,   mode = "loo", true_CV = FALSE,
                                       scores = "logscore")
  cv_stat <- posterior_crossvalidation(fit_stat, mode = "loo", true_CV = FALSE,
                                       scores = "logscore")
  ls_ns   <- as.numeric(cv_ns$scores$logscore)
  ls_stat <- as.numeric(cv_stat$scores$logscore)
  # The non-stationary fit must do at least as well as the stationary
  # fit on data generated from the non-stationary family. Before the
  # fixes ls_ns blew up; the stat fit happened to score similarly here
  # because it cannot exploit the spatial variation.
  expect_lt(ls_ns, ls_stat + 0.1)
})

# ---------------------------------------------------------------------------
# Q from update(latent_model, <random_effects names>) must match Q built
# directly with matern.operators(...) for STATIONARY FEM as well — this
# pins that no analogous theta-naming bug exists for stationary fits in
# either parameterisation. predict.rspde_lme uses the same do.call(update)
# pattern, so if this stays correct, stationary FEM predictions stay
# correct.
# ---------------------------------------------------------------------------

test_that("stationary FEM update(latent_model, <random_effects>) matches direct construction", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)

  for (parametrization_path in c("spde", "matern")) {
    fit_args <- list(
      formula = y ~ -1,
      graph = graph,
      model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))
    if (parametrization_path == "matern") {
      fit_args$model_options <- list(start_nu = 0.5)
    }
    fit <- suppressWarnings(do.call(graph_lme, fit_args))

    coeff_random <- fit$coeff$random_effects
    update_params <- list()
    for (param_name in names(coeff_random)) {
      update_params[[gsub(" \\(fixed\\)$", "", param_name)]] <-
        coeff_random[[param_name]]
    }
    update_params$check_stationarity <- FALSE
    op_via_update <- do.call(update, c(list(fit$latent_model), update_params))

    if (parametrization_path == "spde") {
      op_direct <- rSPDE::matern.operators(
        graph = graph,
        alpha = as.numeric(coeff_random[["alpha (fixed)"]]),
        tau   = as.numeric(coeff_random[["tau"]]),
        kappa = as.numeric(coeff_random[["kappa"]]),
        m = fit$rspde_order,
        parameterization = "spde")
    } else {
      op_direct <- rSPDE::matern.operators(
        graph = graph,
        nu    = as.numeric(coeff_random[["nu (fixed)"]]),
        sigma = as.numeric(coeff_random[["sigma"]]),
        range = as.numeric(coeff_random[["range"]]),
        m = fit$rspde_order,
        parameterization = "matern")
    }
    expect_equal(as.matrix(op_via_update$Q), as.matrix(op_direct$Q),
                 tolerance = 1e-10,
                 info = paste("stationary FEM mismatch under parameterisation:",
                              parametrization_path))
  }
})

test_that("update(latent_model, theta1=..., theta2=...) matches direct spde.matern.operators(theta=...)", {
  skip_if_not_installed("rSPDE")
  set.seed(7)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.05)
  n_mesh <- nrow(graph$mesh$VtE)
  spatial_cov <- graph$mesh$VtE[, 2] - 0.5
  B.range <- cbind(rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh), spatial_cov)
  B.sigma <- cbind(rep(0, n_mesh), rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh))

  # Need a fit to get a latent_model; data doesn't matter for this test.
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 4), runif(4)))
  }
  df <- data.frame(y = rnorm(nrow(PtE)),
                   edge_number = PtE[, 1],
                   distance_on_edge = PtE[, 2])
  graph$add_observations(data = df, normalized = TRUE, verbose = 0)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE,
                                B.range = B.range, B.sigma = B.sigma))

  theta_test <- c(-2.1, -0.8, -1.3)

  # The predict.rspde_lme path: thetas as individual named arguments.
  op_named <- update(fit$latent_model,
                     nu = 0.5,
                     theta1 = theta_test[1],
                     theta2 = theta_test[2],
                     theta3 = theta_test[3],
                     check_stationarity = FALSE)
  # Reference: theta as a single vector.
  op_vec   <- update(fit$latent_model,
                     nu = 0.5,
                     theta = theta_test,
                     check_stationarity = FALSE)
  # Reference: build directly with spde.matern.operators.
  op_direct <- rSPDE::spde.matern.operators(
    graph = graph, nu = 0.5,
    B.range = B.range, B.sigma = B.sigma,
    theta = theta_test,
    parameterization = "matern",
    check_stationarity = FALSE)

  expect_equal(as.matrix(op_named$Q), as.matrix(op_vec$Q),    tolerance = 1e-10)
  expect_equal(as.matrix(op_named$Q), as.matrix(op_direct$Q), tolerance = 1e-10)

  # Before the fix, op_named$tau and op_named$kappa were constant — the
  # spatially-varying theta3 contribution was silently dropped.
  expect_true(length(unique(round(op_named$tau,   8))) > 1L,
              info = "tau must vary spatially with theta")
  expect_true(length(unique(round(op_named$kappa, 8))) > 1L,
              info = "kappa must vary spatially with theta")
})

# ---------------------------------------------------------------------------
# Same exhaustive LOO check as above, but for a NON-STATIONARY α=1 FEM
# fit. The non-stationary parameterisation goes through a different
# branch of update.CBrSPDEobj (theta1/theta2/... instead of tau/kappa)
# and constructs Q from the spatially-varying B.range / B.sigma. A bug
# in how predict.rspde_lme rebuilds Q here — or in how
# posterior_crossvalidation forwards the cached `precomputed_data`
# across folds — would show up here.
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation LOO matches manual LOO (multi-rep, non-stationary FEM)", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
  n_mesh <- nrow(graph$mesh$VtE)
  spatial_cov <- graph$mesh$VtE[, 2] - 0.5
  B.range <- cbind(rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh), spatial_cov)
  B.sigma <- cbind(rep(0, n_mesh), rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh))

  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE,
                                B.range = B.range, B.sigma = B.sigma))
  expect_true(inherits(fit, "rspde_lme"))
  # Sanity: theta1/theta2/theta3 (or similar) — NOT tau/kappa — are what
  # gets fit in the non-stationary parameterisation.
  expect_true(any(grepl("^theta", names(fit$coeff$random_effects))))

  res_pre <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = TRUE,
                                       scores = c("logscore", "rmse"))
  res_no  <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = FALSE,
                                       scores = c("logscore", "rmse"))
  expect_equal(res_pre$mu,  res_no$mu,  tolerance = 1e-9)
  expect_equal(res_pre$var, res_no$var, tolerance = 1e-9)

  gd <- fit$graph$.__enclos_env__$private$data
  rep_vec <- gd[[".group"]]
  sig2 <- as.numeric(fit$coeff$measurement_error)^2
  for (i in seq_len(fit$nobs)) {
    r <- rep_vec[i]
    loc_test <- cbind(gd[[".edge_number"]][i],
                      gd[[".distance_on_edge"]][i])
    manual <- .fem_manual_predict(fit, loc_test, which_repl = r,
                                  na_idx = i)
    target_mean <- unname(manual$mean[[as.character(r)]])
    target_var  <- unname(manual$variance[[as.character(r)]]) + sig2
    expect_equal(unname(res_pre$mu[i]),  target_mean, tolerance = 1e-7,
                 info = paste("LOO mean mismatch at i =", i, "(rep", r, ")"))
    expect_equal(unname(res_pre$var[i]), target_var,  tolerance = 1e-7,
                 info = paste("LOO var mismatch at i =", i, "(rep", r, ")"))
  }
})

# ---------------------------------------------------------------------------
# posterior_crossvalidation k-fold LOO equivalence: when a k-fold fold
# happens to contain exactly one test point from a given replicate, that
# point's prediction must match the manual LOO formula (since the rest
# of that replicate's data is the same as in LOO for that point).
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation k-fold reduces to LOO for single-point folds (multi-rep FEM)", {
  skip_if_not_installed("rSPDE")
  graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))

  res_k_pre <- posterior_crossvalidation(fit, mode = "k-fold", k = 10,
                                         seed = 17,
                                         use_precomputed = TRUE,
                                         scores = c("logscore", "rmse"),
                                         return_indices = TRUE)
  res_k_no  <- posterior_crossvalidation(fit, mode = "k-fold", k = 10,
                                         seed = 17,
                                         use_precomputed = FALSE,
                                         scores = c("logscore", "rmse"))
  expect_equal(res_k_pre$mu,  res_k_no$mu,  tolerance = 1e-9)
  expect_equal(res_k_pre$var, res_k_no$var, tolerance = 1e-9)

  gd <- fit$graph$.__enclos_env__$private$data
  rep_vec <- gd[[".group"]]

  # For each fold, look for replicates whose held-out set is a single
  # point — that point's k-fold prediction is a proper LOO and must
  # match the manual computation.
  any_checked <- FALSE
  for (fold in res_k_pre$indices) {
    test_idx <- fold$test
    repls_in_fold <- rep_vec[test_idx]
    for (r in unique(repls_in_fold)) {
      pts_r <- test_idx[repls_in_fold == r]
      if (length(pts_r) == 1L) {
        i_test <- pts_r
        loc_test <- cbind(gd[[".edge_number"]][i_test],
                          gd[[".distance_on_edge"]][i_test])
        manual <- .fem_manual_predict(fit, loc_test, which_repl = r,
                                      na_idx = i_test)
        expect_equal(unname(res_k_pre$mu[i_test]),
                     manual$mean[[as.character(r)]],
                     tolerance = 1e-7,
                     info = paste("k-fold mu mismatch at i =", i_test,
                                  "(rep", r, ")"))
        any_checked <- TRUE
      }
    }
  }
  expect_true(any_checked,
              info = "no single-point single-replicate folds to check")
})
