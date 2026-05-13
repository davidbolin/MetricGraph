# Tests for posterior_crossvalidation() comparing against a hand-rolled
# leave-one-out reference. These pin two classes of bugs that were live
# at the time these tests were written:
#
#   (a) precomputation returned NULL for rspde_lme objects, which then
#       sent process_fold down a fallback path that referenced an
#       undefined `new_data` variable;
#   (b) prediction output for multi-replicate models was read by linear
#       position into pred$mean, so the LOO of an observation in
#       replicate r was scored against the prediction for replicate 1.
#
# The reference implementation builds the joint covariance of all
# observations from the fitted parameters and does LOO with the standard
# Schur-complement formulas. For replicates we LOO within each replicate
# since replicates are conditionally independent.

# ---------------------------------------------------------------------------
# Reference LOO using the model's joint covariance at observation locations.
# ---------------------------------------------------------------------------

.loo_reference <- function(fit) {
  graph <- fit$graph$clone()
  response_name <- as.character(fit$response_var)
  data <- graph$.__enclos_env__$private$data
  y_all   <- data[[response_name]]
  grp_all <- data[[".group"]]

  sigma_e <- fit$coeff$measurement_error
  tau     <- fit$coeff$random_effects[1]
  kappa   <- fit$coeff$random_effects[2]

  # Build the full covariance matrix at the (within-replicate) observation
  # locations. The graph is the same for every replicate, so we build it
  # once and re-use it per replicate.
  alpha <- fit$latent_model$alpha
  g_v <- graph$clone()
  g_v$observation_to_vertex()
  if (alpha == 1) {
    Q <- MetricGraph:::Qalpha1(c(tau, kappa), g_v, BC = fit$BC)
    Sigma <- as.matrix(solve(Q))[g_v$PtV, g_v$PtV]
  } else if (alpha == 2) {
    g_v$buildC(2, fit$BC == 0)
    n.c <- 1:length(g_v$CoB$S)
    Q <- MetricGraph:::Qalpha2(c(tau, kappa), g_v, BC = fit$BC)
    Qtilde <- (g_v$CoB$T) %*% Q %*% t(g_v$CoB$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.full <- t(g_v$CoB$T[-n.c, ]) %*% solve(Qtilde) %*% (g_v$CoB$T[-n.c, ])
    PtE <- g_v$get_PtE()
    index.obs <- 4 * (PtE[, 1] - 1) +
      (1 * (PtE[, 2] == 0)) + (3 * (PtE[, 2] != 0))
    Sigma <- as.matrix(Sigma.full[index.obs, index.obs])
  } else {
    stop("alpha must be 1 or 2 for the reference")
  }
  Sigma_o <- Sigma
  diag(Sigma_o) <- diag(Sigma_o) + sigma_e^2

  # Per-replicate LOO using the same Sigma (graph is shared) but only the
  # rows/columns of observations in that replicate.
  n_total <- length(y_all)
  mu  <- rep(NA_real_, n_total)
  var <- rep(NA_real_, n_total)
  for (r in unique(grp_all)) {
    pos_r <- which(grp_all == r)
    # Within-replicate observations sit at the first length(pos_r) entries
    # of Sigma — every replicate uses the same set of observation
    # locations on the (shared) graph.
    nr <- length(pos_r)
    S_r   <- Sigma[seq_len(nr), seq_len(nr), drop = FALSE]
    So_r  <- Sigma_o[seq_len(nr), seq_len(nr), drop = FALSE]
    y_r   <- y_all[pos_r]
    for (i in seq_len(nr)) {
      mu[pos_r[i]]  <- S_r[i, -i, drop = FALSE] %*%
                         solve(So_r[-i, -i, drop = FALSE], y_r[-i])
      var[pos_r[i]] <- So_r[i, i] - S_r[i, -i, drop = FALSE] %*%
                         solve(So_r[-i, -i, drop = FALSE],
                               t(S_r[i, -i, drop = FALSE]))
    }
  }
  list(mu = mu, var = var)
}

.fit_small <- function(alpha = 1, n_repl = 1, seed = 1L,
                       obs_per_edge = 4L) {
  set.seed(seed)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)

  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, obs_per_edge), runif(obs_per_edge)))
  }

  kappa <- 5; tau <- 1; sigma_e <- 0.3
  df_all <- NULL
  for (r in seq_len(n_repl)) {
    u <- sample_spde(kappa = kappa, tau = tau, alpha = alpha,
                     graph = graph, PtE = PtE)
    y <- u + sigma_e * rnorm(length(u))
    df_all <- rbind(df_all,
      data.frame(y = y, edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
  group_var <- if (n_repl > 1) "repl" else NULL
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = group_var)

  # No intercept so that the reference (which treats y as zero-mean)
  # matches posterior_crossvalidation byte-for-byte.
  graph_lme(y ~ -1, graph = graph,
            model = list(type = "WhittleMatern", alpha = alpha))
}

# ---------------------------------------------------------------------------
# Single-replicate, alpha = 1: LOO mean and variance match the reference,
# and the use_precomputed = TRUE / FALSE paths agree with each other.
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation LOO matches manual LOO (alpha=1, no replicates)", {
  fit <- .fit_small(alpha = 1, n_repl = 1)
  ref <- .loo_reference(fit)

  res_pre <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = TRUE,
                                       scores = c("logscore", "crps", "mae", "rmse"))
  res_no  <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = FALSE,
                                       scores = c("logscore", "crps", "mae", "rmse"))

  expect_equal(as.numeric(res_pre$mu),  ref$mu,  tolerance = 1e-8)
  expect_equal(as.numeric(res_pre$var), ref$var, tolerance = 1e-8)
  expect_equal(as.numeric(res_no$mu),   ref$mu,  tolerance = 1e-8)
  expect_equal(as.numeric(res_no$var),  ref$var, tolerance = 1e-8)

  # The two paths must produce the same result.
  expect_equal(res_pre$mu,  res_no$mu,  tolerance = 1e-10)
  expect_equal(res_pre$var, res_no$var, tolerance = 1e-10)
  expect_equal(as.numeric(res_pre$scores$logscore),
               as.numeric(res_no$scores$logscore),  tolerance = 1e-10)
  expect_equal(as.numeric(res_pre$scores$crps),
               as.numeric(res_no$scores$crps),      tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Single-replicate, alpha = 2.
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation LOO matches manual LOO (alpha=2, no replicates)", {
  fit <- .fit_small(alpha = 2, n_repl = 1)
  ref <- .loo_reference(fit)

  # logscore is requested so that variances are actually computed —
  # posterior_crossvalidation only fills in $var when some
  # variance-dependent score is asked for.
  res <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                   scores = c("logscore", "mae", "rmse"))
  expect_equal(as.numeric(res$mu),  ref$mu,  tolerance = 1e-7)
  expect_equal(as.numeric(res$var), ref$var, tolerance = 1e-7)
})

# ---------------------------------------------------------------------------
# Multi-replicate, alpha = 1: two prior bugs converged here.
#   1. process_fold indexed pred$mean by position, so the LOO of an
#      observation in replicate r was scored against replicate 1's
#      prediction.
#   2. predict.graph_lme overwrote $variance with the LAST replicate's
#      kriging variance, so every replicate reported the same wrong
#      variance and the logscore blew up.
# Both fixes are exercised together: $mu AND $var must match the
# per-replicate reference computation.
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation LOO matches manual LOO (alpha=1, with replicates)", {
  fit <- .fit_small(alpha = 1, n_repl = 3)
  ref <- .loo_reference(fit)

  # Both the precomputed (fast) and non-precomputed (fallback) paths
  # must agree with the manual per-replicate LOO. For LOO every fold
  # has a single test point in a single replicate, so the precomputed
  # path is always exercised.
  res_pre <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = TRUE,
                                       scores = c("logscore", "crps", "mae", "rmse"))
  res_no  <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = FALSE,
                                       scores = c("logscore", "crps", "mae", "rmse"))

  expect_equal(as.numeric(res_pre$mu),  ref$mu,  tolerance = 1e-7)
  expect_equal(as.numeric(res_pre$var), ref$var, tolerance = 1e-7)
  expect_equal(as.numeric(res_no$mu),   ref$mu,  tolerance = 1e-7)
  expect_equal(as.numeric(res_no$var),  ref$var, tolerance = 1e-7)

  # The two paths must give identical numerics — the precomputed path is
  # an optimisation, not a different algorithm.
  expect_equal(res_pre$mu,  res_no$mu,  tolerance = 1e-10)
  expect_equal(res_pre$var, res_no$var, tolerance = 1e-10)

  # RMSE is a function of the residuals only, so it should match the
  # reference exactly up to numerical noise.
  ref_rmse <- sqrt(mean((ref$mu - fit$graph$.__enclos_env__$private$data[["y"]])^2))
  expect_equal(as.numeric(res_pre$scores$rmse), ref_rmse, tolerance = 1e-7)

  # logscore must be finite — before the variance fix it was ~1e6 because
  # var.p collapsed to ~0 for most folds.
  expect_true(is.finite(as.numeric(res_pre$scores$logscore)))
  expect_lt(as.numeric(res_pre$scores$logscore), 10)
})

# ---------------------------------------------------------------------------
# Focused regression test for the predict.graph_lme variance bug.
# Construct two replicates that share observation locations, NA out one
# observation in replicate 1, fit, and predict at that exact location.
# The kriging variance must be sizeable for replicate 1 (the NA'd one)
# and essentially zero for replicate 2 (still observed there). Before the
# fix, both replicates reported the same value — whichever replicate the
# per-replicate loop visited LAST.
# ---------------------------------------------------------------------------

test_that("predict.graph_lme reports per-replicate kriging variances independently", {
  set.seed(2)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)

  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 4), runif(4)))
  }

  df_all <- NULL
  for (r in 1:2) {
    u <- sample_spde(kappa = 5, tau = 1, alpha = 1, graph = graph, PtE = PtE)
    y <- u + 0.3 * rnorm(length(u))
    df_all <- rbind(df_all,
      data.frame(y = y, edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = "repl")

  # NA out a single observation in replicate 1 only.
  i_na <- 5
  graph$.__enclos_env__$private$data$y[i_na] <- NA

  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1))

  gd <- graph$.__enclos_env__$private$data
  nd <- list(
    .edge_number      = gd[[".edge_number"]][i_na],
    .distance_on_edge = gd[[".distance_on_edge"]][i_na],
    .group            = "1",
    y                 = NA_real_
  )
  pred <- predict(fit, newdata = nd,
                  edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge",
                  normalized = TRUE, compute_variances = TRUE)

  v1 <- pred$variance[pred$repl == "1"]
  v2 <- pred$variance[pred$repl == "2"]
  expect_length(v1, 1L)
  expect_length(v2, 1L)

  # Sanity floor: replicate 1's obs is NA'd, so its kriging variance must
  # not be near zero.
  expect_gt(v1, 1e-3)

  # Core regression check: the two replicates must NOT collapse to the
  # same value. Before the fix, $variance was overwritten with the last
  # replicate's value, so v1 / v2 were bit-identical. After the fix the
  # NA'd replicate's variance is substantially larger.
  expect_gt(v1 / v2, 5)
})

# ---------------------------------------------------------------------------
# rspde_lme (FEM-based, fractional nu) — the original bug report.
# Before the fix this raised: object 'new_data' not found, after warning
# "Precomputation for pseudo-CV failed, falling back to standard prediction".
# The precomputation cannot succeed for rspde_lme because predict.graph_lme
# dispatches to predict.rspde_lme which does not produce a precomputed_data
# field. posterior_crossvalidation must skip precomputation and just
# produce finite predictions.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# posterior_crossvalidation_loo is the older covariance-algebra LOO. It has
# two analogous bugs to the ones fixed in posterior_crossvalidation /
# predict.graph_lme:
#
#   (A) rspde_lme objects fall through the model-type dispatch and get
#       silently re-classified as a graph-Laplacian model with the wrong
#       parameters, returning finite but meaningless scores.
#   (B) for multi-replicate fits, mu.p[i] is overwritten on every
#       replicate iteration of the inner loop, so the returned $mu / $var
#       have length n_obs_per_replicate instead of nobs and only contain
#       the LAST replicate's predictions.
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation_loo matches posterior_crossvalidation for single-replicate alpha=1", {
  fit <- .fit_small(alpha = 1, n_repl = 1)
  res_loo <- posterior_crossvalidation_loo(fit)
  res_pc  <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       scores = c("logscore", "crps", "mae", "rmse"))
  expect_equal(as.numeric(res_loo$mu),  as.numeric(res_pc$mu),  tolerance = 1e-8)
  expect_equal(as.numeric(res_loo$var), as.numeric(res_pc$var), tolerance = 1e-8)
})

test_that("posterior_crossvalidation_loo returns per-observation mu/var for all replicates", {
  fit <- .fit_small(alpha = 1, n_repl = 3)
  ref <- .loo_reference(fit)

  res <- posterior_crossvalidation_loo(fit)

  # Bug B regression: $mu used to have length 16 (one replicate) for a
  # 48-observation fit.
  expect_equal(length(res$mu),  fit$nobs)
  expect_equal(length(res$var), fit$nobs)
  expect_equal(as.numeric(res$mu),  ref$mu,  tolerance = 1e-7)
  expect_equal(as.numeric(res$var), ref$var, tolerance = 1e-7)
})

test_that("posterior_crossvalidation_loo refuses to silently mis-classify rspde_lme fits", {
  skip_if_not_installed("rSPDE")
  set.seed(11)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.05)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 5), runif(5)))
  }
  u <- sample_spde(kappa = 5, tau = 1, alpha = 1, graph = graph, PtE = PtE)
  y <- u + 0.3 * rnorm(length(u))
  graph$add_observations(
    data = data.frame(y = y,
                      edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))
  expect_true(inherits(fit, "rspde_lme"))

  # Should error explicitly rather than return scores from a wrong model.
  # Before the fix, it silently produced wildly off scores (logscore ~ 0.42
  # vs the correct ~ 0.02 from posterior_crossvalidation).
  expect_error(
    posterior_crossvalidation_loo(fit),
    regexp = "rspde|fem|fractional|posterior_crossvalidation",
    ignore.case = TRUE
  )
})

test_that("posterior_crossvalidation LOO runs cleanly on an rspde_lme (FEM) fit", {
  skip_if_not_installed("rSPDE")
  set.seed(11)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.05)

  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 5), runif(5)))
  }
  u <- sample_spde(kappa = 5, tau = 1, alpha = 1, graph = graph, PtE = PtE)
  y <- u + 0.3 * rnorm(length(u))
  graph$add_observations(
    data = data.frame(y = y,
                      edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0)

  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))
  expect_true(inherits(fit, "rspde_lme"))

  # Must not warn about precomputation failing or error with the old
  # "object 'new_data' not found".
  expect_no_warning(
    res <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                     scores = c("mae", "rmse"))
  )
  expect_equal(length(res$mu), fit$nobs)
  expect_true(all(is.finite(res$mu)))
  expect_true(is.finite(as.numeric(res$scores$rmse)))
})

# ---------------------------------------------------------------------------
# Once rSPDE::predict.rspde_lme grew an `na_test_idx` option (so the
# response and A_list are correctly masked for held-out points) and a
# `precompute_data` option (so the parameter-dependent Q + mean
# correction can be reused across folds), posterior_crossvalidation can
# enable precomputation for rspde_lme fits. The two paths
# (use_precomputed TRUE / FALSE) must produce the same predictions, and
# the predictions must be a *true* LOO (no leak of the test value into
# training).
# ---------------------------------------------------------------------------

test_that("posterior_crossvalidation LOO is a true LOO on an rspde_lme (FEM) fit", {
  skip_if_not_installed("rSPDE")
  set.seed(11)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.05)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 5), runif(5)))
  }
  u <- sample_spde(kappa = 5, tau = 1, alpha = 1, graph = graph, PtE = PtE)
  y <- u + 0.3 * rnorm(length(u))
  graph$add_observations(
    data = data.frame(y = y, edge_number = PtE[, 1], distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))

  # Manual proper LOO for a handful of folds: mask the i-th observation
  # in both model_matrix and A_list, then call predict directly.
  gd <- fit$graph$.__enclos_env__$private$data
  manual_loo <- function(i_test) {
    nd <- lapply(gd, function(x) x[i_test])
    cv_model <- fit
    mm <- as.matrix(cv_model$model_matrix)
    mm[i_test, 1] <- NA
    cv_model$model_matrix <- mm
    cv_model$A_list <- fit$A_list
    cv_model$A_list[[1]] <- cv_model$A_list[[1]][-i_test, , drop = FALSE]
    p <- predict(cv_model,
                 newdata = nd,
                 edge_number = ".edge_number",
                 distance_on_edge = ".distance_on_edge",
                 normalized = TRUE,
                 compute_variances = TRUE)
    list(mean = as.numeric(p$mean), variance = as.numeric(p$variance))
  }

  res_pre <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = TRUE,
                                       scores = c("logscore", "crps", "mae", "rmse"))
  res_no  <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = FALSE,
                                       scores = c("logscore", "crps", "mae", "rmse"))

  # The two paths must agree numerically — precompute is an optimisation.
  expect_equal(res_pre$mu,  res_no$mu,  tolerance = 1e-9)
  expect_equal(res_pre$var, res_no$var, tolerance = 1e-9)

  # Each path's per-observation predictions must match the manual proper LOO.
  for(i in c(1L, 5L, 13L, fit$nobs)) {
    m <- manual_loo(i)
    # measurement-error variance is added by posterior_crossvalidation
    # but is NOT part of predict's kriging variance.
    target_var <- m$variance + fit$coeff$measurement_error^2
    expect_equal(unname(res_pre$mu[i]),  m$mean,     tolerance = 1e-7,
                 info = paste("mean mismatch at i =", i))
    expect_equal(unname(res_pre$var[i]), unname(target_var), tolerance = 1e-7,
                 info = paste("variance mismatch at i =", i))
  }
})

test_that("predict.rspde_lme reports per-replicate kriging variances independently", {
  skip_if_not_installed("rSPDE")
  set.seed(2)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = 0.1)
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 4), runif(4)))
  }
  df_all <- NULL
  for (r in 1:2) {
    u <- sample_spde(kappa = 5, tau = 1, alpha = 1, graph = graph, PtE = PtE)
    y <- u + 0.3 * rnorm(length(u))
    df_all <- rbind(df_all,
      data.frame(y = y, edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = "repl")
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))
  expect_true(inherits(fit, "rspde_lme"))

  # NA one observation in replicate 1 only via the cv helper.
  i_na <- 5
  cv_model <- MetricGraph:::update_graph_lme_with_na(fit, i_na)
  gd <- fit$graph$.__enclos_env__$private$data
  nd <- list(
    .edge_number      = gd[[".edge_number"]][i_na],
    .distance_on_edge = gd[[".distance_on_edge"]][i_na],
    .group            = "1",
    y                 = NA_real_
  )
  pred <- predict(cv_model,
                  newdata = nd,
                  edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge",
                  normalized = TRUE, compute_variances = TRUE)

  v1 <- pred$variance[pred$repl == "1"]
  v2 <- pred$variance[pred$repl == "2"]
  expect_length(v1, 1L)
  expect_length(v2, 1L)
  # Replicate 1 holds out its observation -> variance should be sizeable.
  expect_gt(v1, 1e-3)
  # Core regression check: the two replicates must NOT collapse to the
  # same value. Before the rspde fix, $variance was overwritten with
  # the last replicate's value.
  expect_gt(v1 / v2, 2)
})
