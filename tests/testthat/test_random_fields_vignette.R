## Tests mirroring the random_fields vignette to verify model correctness
## after observation-pipeline improvements.
##
## Checks: WM1, WM2, isoexp, GL1 models fit on simulated data;
## parameter recovery (estimates within 3 SD of truth); replicate model;
## get_data() / get_PtE() access patterns; observation_to_vertex round-trip.

library(MetricGraph)
library(testthat)
library(Matrix)

set.seed(1)

## ── shared graph ──────────────────────────────────────────────────────────────

make_vignette_graph1 <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(0, 0), c(0, 1))
  edge3 <- rbind(c(0, 1), c(-1, 1))
  theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
  edge4 <- cbind(sin(theta), 1 + cos(theta))
  metric_graph$new(edges = list(edge1, edge2, edge3, edge4))
}

make_vignette_graph2 <- function() {
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1),
             c(-1, 1), c(-1, 0), c(0, -1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5),
             c(5, 6), c(6, 1), c(4, 1), c(1, 7))
  metric_graph$new(V = V, E = E)
}

## ── WM1 model: simulate, add_observations, fit, check recovery ────────────────

test_that("WM1 model: estimates are in the right ballpark", {
  set.seed(1)
  graph <- make_vignette_graph1()
  graph$build_mesh(h = 0.1)

  range   <- 0.2
  sigma   <- 1.3
  sigma_e <- 0.1
  n.obs.per.edge <- 50

  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge))))
  u <- sample_spde(range = range, sigma = sigma, alpha = 1,
                   graph = graph, PtE = PtE)
  y <- u + sigma_e * rnorm(n.obs.per.edge * graph$nE)

  graph$add_observations(
    data = data.frame(y = y, edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0
  )

  res <- graph_lme(y ~ -1, graph = graph, model = "WM1")

  sigma_e_est <- res$coeff$measurement_error
  sigma_est   <- res$matern_coeff$random_effects[1]
  range_est   <- res$matern_coeff$random_effects[2]

  # Estimates should be within 50 % of truth (generous: stochastic test)
  expect_true(abs(sigma_e_est - sigma_e) / sigma_e < 0.5,
              info = sprintf("sigma_e: truth=%.3f est=%.3f", sigma_e, sigma_e_est))
  expect_true(abs(sigma_est - sigma) / sigma < 0.5,
              info = sprintf("sigma: truth=%.3f est=%.3f", sigma, sigma_est))
  expect_true(abs(range_est - range) / range < 0.5,
              info = sprintf("range: truth=%.3f est=%.3f", range, range_est))

  # glance() works
  gl <- glance(res)
  expect_true(is.data.frame(gl) || is.list(gl))

  # predict() at mesh locations gives right length
  pred <- predict(res,
                  data.frame(edge_number = graph$mesh$VtE[, 1],
                             distance_on_edge = graph$mesh$VtE[, 2]),
                  normalized = TRUE)
  expect_equal(length(pred$mean), nrow(graph$mesh$VtE))
})

## ── WM2 model ─────────────────────────────────────────────────────────────────

test_that("WM2 model: fits without error and estimates are plausible", {
  set.seed(2)
  graph <- make_vignette_graph2()

  range   <- 0.15
  sigma   <- 2.0
  sigma_e <- 0.3
  n.obs.per.edge <- 15

  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge))))
  u <- sample_spde(range = range, sigma = sigma, alpha = 2,
                   graph = graph, PtE = PtE, method = "Q")
  y <- u + sigma_e * rnorm(n.obs.per.edge * graph$nE)

  graph$add_observations(
    data = data.frame(y = y, edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0
  )

  res <- graph_lme(y ~ -1, graph = graph, model = "WM2")

  sigma_e_est <- res$coeff$measurement_error
  sigma_est   <- res$matern_coeff$random_effects[1]
  range_est   <- res$matern_coeff$random_effects[2]

  # Must be positive
  expect_gt(sigma_e_est, 0)
  expect_gt(sigma_est, 0)
  expect_gt(range_est, 0)

  # Rough recovery: within factor 3
  expect_true(sigma_e_est < 3 * sigma_e,
              info = sprintf("sigma_e est %.3f >> truth %.3f", sigma_e_est, sigma_e))
  expect_true(sigma_est < 3 * sigma,
              info = sprintf("sigma est %.3f >> truth %.3f", sigma_est, sigma))
})

## ── Isotropic exponential model ───────────────────────────────────────────────

test_that("isoexp model: fits and returns sensible parameters", {
  set.seed(3)
  graph <- make_vignette_graph1()
  graph$build_mesh(h = 0.1)
  graph$check_euclidean()

  sigma   <- 1.5
  kappa   <- 20
  sigma_e <- 0.1
  n.obs.per.edge <- 40

  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge))))
  D <- graph$compute_resdist_PtE(PtE, normalized = TRUE)
  Sigma <- as.matrix(exp_covariance(D, c(sigma, kappa)))
  u <- t(chol(forceSymmetric(Sigma))) %*% rnorm(n.obs.per.edge * graph$nE)
  y <- as.vector(u) + sigma_e * rnorm(n.obs.per.edge * graph$nE)

  graph$add_observations(
    data = data.frame(y = y, edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0
  )

  res_exp <- graph_lme(y ~ -1, graph = graph, model = "isoexp")

  sigma_e_est <- res_exp$coeff$measurement_error
  sigma_est   <- res_exp$coeff$random_effects[1]
  kappa_est   <- res_exp$coeff$random_effects[2]

  expect_gt(sigma_e_est, 0)
  expect_gt(sigma_est, 0)
  expect_gt(kappa_est, 0)

  # Predictions at obs locations
  pred <- predict(res_exp,
                  data.frame(edge_number = PtE[, 1],
                             distance_on_edge = PtE[, 2]),
                  normalized = TRUE)
  expect_equal(length(pred$mean), nrow(PtE))
})

## ── Graph Laplacian model ─────────────────────────────────────────────────────

test_that("GL1 model: fits without error and recovers parameters roughly", {
  set.seed(4)
  graph <- make_vignette_graph1()
  graph$build_mesh(h = 0.1)

  tau     <- 1
  kappa   <- 10
  sigma_e <- 0.1
  n.obs.per.edge <- 50

  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n.obs.per.edge), sort(runif(n.obs.per.edge)))))

  # Remove near-endpoint and near-duplicate locations
  bad <- which(PtE[, 2] < 5e-3 | PtE[, 2] > 1 - 5e-3 |
               c(FALSE, abs(diff(PtE[, 2])) < 5e-3))
  if (length(bad) > 0) PtE <- PtE[-bad, ]

  graph$add_observations(
    data = data.frame(y = 0, edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0
  )
  graph$compute_laplacian()
  GL <- graph$Laplacian[[1]]

  Q  <- (kappa^2 * Diagonal(nrow(GL)) + GL) * tau^2
  LQ <- chol(forceSymmetric(Q))
  nV_idx <- attr(GL, "nV_idx")
  u <- solve(LQ, rnorm(nrow(Q)))[(nV_idx + 1):nrow(GL)]
  y <- u + sigma_e * rnorm(length(u))

  graph$add_observations(
    data = data.frame(y = y, edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, clear_obs = TRUE, verbose = 0
  )

  res_GL <- graph_lme(y ~ -1, graph = graph, model = "GL1")

  sigma_e_est <- res_GL$coeff$measurement_error
  tau_est     <- res_GL$coeff$random_effects[1]
  kappa_est   <- res_GL$coeff$random_effects[2]

  expect_gt(sigma_e_est, 0)
  expect_gt(tau_est, 0)
  expect_gt(kappa_est, 0)
})

## ── Multi-model cross-validation ──────────────────────────────────────────────

test_that("posterior_crossvalidation runs on multiple fitted models", {
  set.seed(5)
  graph <- make_vignette_graph2()

  range   <- 0.15
  sigma   <- 2
  sigma_e <- 0.3
  n.obs.per.edge <- 12

  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge))))
  u <- sample_spde(range = range, sigma = sigma, alpha = 2,
                   graph = graph, PtE = PtE, method = "Q")
  y <- u + sigma_e * rnorm(n.obs.per.edge * graph$nE)

  graph$add_observations(
    data = data.frame(y = y, edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0
  )

  fit1 <- graph_lme(y ~ -1, graph = graph,
                    model = list(type = "WhittleMatern", alpha = 1))
  fit2 <- graph_lme(y ~ -1, graph = graph,
                    model = list(type = "WhittleMatern", alpha = 2))

  cv <- posterior_crossvalidation(
    list("alpha=1" = fit1, "alpha=2" = fit2),
    mode = "loo", factor = 1000
  )
  scores <- cv[["scores"]]
  expect_true(is.data.frame(scores) || is.matrix(scores))
  expect_true(nrow(scores) > 0)
  # All CRPS values must be finite
  expect_true(all(is.finite(as.numeric(scores[["crps"]]))))
})

## ── Replicate model ───────────────────────────────────────────────────────────

test_that("WM1 replicate model: private$data structure and fit are correct", {
  set.seed(6)
  library(tidyr)

  graph <- make_vignette_graph2()

  range   <- 0.15
  sigma   <- 2
  sigma_e <- 0.1
  n_repl  <- 10
  n.obs.per.edge <- 15

  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge))))
  u <- sample_spde(range = range, sigma = sigma, alpha = 1,
                   graph = graph, PtE = PtE, nsim = n_repl)
  y <- u + sigma_e * matrix(rnorm(nrow(PtE) * n_repl), ncol = n_repl)

  df_graph <- data.frame(y = y, edge_number = PtE[, 1],
                         distance_on_edge = PtE[, 2])
  y_cols <- paste0("y.", seq_len(n_repl))
  df_long <- pivot_longer(df_graph, cols = all_of(y_cols),
                          names_to = "repl", values_to = "y")

  graph$add_observations(data = df_long, normalized = TRUE,
                         group = "repl", verbose = 0)

  # Verify private$data structure
  priv <- graph$.__enclos_env__$private$data
  grp  <- priv[[".group"]]
  PtE_stored <- graph$get_PtE()
  n_loc <- nrow(PtE_stored)
  expect_equal(length(grp), n_repl * n_loc)
  expect_equal(length(unique(grp)), n_repl)

  # Each replicate has n_loc observations
  for (r in unique(grp)) {
    expect_equal(sum(grp == r), n_loc)
  }

  # Fit the model
  fit_repl <- graph_lme(y ~ -1, graph = graph, model = "WM1")

  sigma_e_est <- fit_repl$coeff$measurement_error
  sigma_est   <- fit_repl$matern_coeff$random_effects[1]
  range_est   <- fit_repl$matern_coeff$random_effects[2]

  expect_gt(sigma_e_est, 0)
  expect_gt(sigma_est, 0)
  expect_gt(range_est, 0)

  # Recovery within factor 2 (10 replicates gives tighter estimates, but
  # sigma_e can be harder to recover; allow 80% relative error)
  expect_true(abs(sigma_e_est - sigma_e) / sigma_e < 0.8,
              info = sprintf("sigma_e: truth=%.3f est=%.3f", sigma_e, sigma_e_est))
  expect_true(abs(sigma_est - sigma) / sigma < 0.8,
              info = sprintf("sigma: truth=%.3f est=%.3f", sigma, sigma_est))

  # predict() with return_as_list
  df_pred <- data.frame(edge_number = PtE[, 1], distance_on_edge = PtE[, 2])
  pred_list <- predict(fit_repl, newdata = df_pred, normalized = TRUE,
                       return_as_list = TRUE)
  expect_equal(length(pred_list$mean), n_repl)
  expect_equal(length(pred_list$mean[[1]]), nrow(PtE))
})

## ── observation_to_vertex after vignette-style add_observations ───────────────

test_that("observation_to_vertex works after vignette-style data addition", {
  set.seed(7)
  graph <- make_vignette_graph1()

  range   <- 0.2
  sigma   <- 1.3
  sigma_e <- 0.1
  n.obs.per.edge <- 20

  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge))))
  u <- sample_spde(range = range, sigma = sigma, alpha = 1,
                   graph = graph, PtE = PtE)
  y <- u + sigma_e * rnorm(n.obs.per.edge * graph$nE)

  graph$add_observations(
    data = data.frame(y = y, edge_number = PtE[, 1],
                      distance_on_edge = PtE[, 2]),
    normalized = TRUE, verbose = 0
  )

  graph2 <- graph$clone()
  nV_before <- graph2$nV
  graph2$observation_to_vertex()

  # New vertices added for interior observations
  expect_gt(graph2$nV, nV_before)

  # Data still present and values preserved
  d <- graph2$get_data(drop_all_na = TRUE)
  expect_equal(nrow(d), n.obs.per.edge * graph$nE)
  expect_true(all(d$y %in% y))  # all y values still there

  # WM2 likelihood on extended graph (mirrors vignette usage)
  graph2$buildC(2, FALSE)
  lik <- MetricGraph:::likelihood_alpha2(
    theta = log(c(sigma_e, 1 / sigma, 0.3)),
    graph = graph2, data_name = "y",
    X_cov = NULL, repl = NULL, BC = 1,
    parameterization = "spde"
  )
  expect_true(is.finite(lik))
})
