## Regression tests for removing Kirchhoff-constraint rows when a graph has no
## constraints (e.g. a single edge). The idiom `T[-c(1:n_const), ]` silently
## dropped the first row when n_const == 0 (since 1:0 is c(1, 0)), giving
## dimension errors or wrong numbers. Each quantity on a single edge is compared
## with the same interval split at a degree-2 vertex, which has constraints and
## so takes the other code path.

line_graph <- function(split = FALSE) {
  if (split) {
    edges <- list(rbind(c(0, 0), c(0.5, 0)), rbind(c(0.5, 0), c(1, 0)))
  } else {
    edges <- list(rbind(c(0, 0), c(1, 0)))
  }
  metric_graph$new(edges = edges, verbose = 0)
}

# (edge_number, normalized distance_on_edge) for positions x in [0, 1]
line_PtE <- function(x, split = FALSE) {
  if (!split) {
    return(data.frame(edge_number = 1, distance_on_edge = x))
  }
  data.frame(edge_number = ifelse(x < 0.5, 1, 2),
             distance_on_edge = ifelse(x < 0.5, 2 * x, 2 * x - 1))
}

line_graph_obs <- function(x, y, split = FALSE) {
  g <- line_graph(split)
  df <- line_PtE(x, split)
  df$y <- y
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  g
}

mesh_values_by_x <- function(g, v) {
  as.vector(v)[order(g$mesh$V[, 1])]
}

obs_x <- c(0.1, 0.3, 0.45, 0.7, 0.9)
obs_y <- c(-0.96, -0.29, 0.26, -1.15, 0.2)
theta_log <- log(c(0.3, 0.8, 2)) # (sigma_e, reciprocal_tau, kappa), log scale
theta_nat <- c(0.3, 0.8, 2)      # (sigma_e, tau, kappa)

test_that("buildDirectionalConstraints works on a single edge", {
  g <- line_graph()
  expect_no_error(g$buildDirectionalConstraints(alpha = 1))
  expect_length(g$CoB$S, 0)
  expect_equal(dim(g$CoB$T), c(2L, 2L))
})

test_that("alpha = 2 likelihoods on a single edge match the split edge", {
  g1 <- line_graph_obs(obs_x, obs_y)
  g1$buildC(2)
  g2 <- line_graph_obs(obs_x, obs_y, split = TRUE)
  g2$buildC(2)
  expect_length(g1$CoB$S, 0)

  lik <- function(g) {
    likelihood_alpha2(theta_log, g, data_name = "y", repl = NULL, BC = 1,
                      parameterization = "spde")
  }
  expect_equal(lik(g1), lik(g2), tolerance = 1e-8)

  lik_pre <- function(g) {
    likelihood_alpha2_precompute(theta_log, precompute_alpha2(g, data_name = "y"),
                                 BC = 1, parameterization = "spde")
  }
  expect_equal(lik_pre(g1), lik_pre(g2), tolerance = 1e-8)
})

test_that("directional alpha = 1 likelihoods on a single edge match the split edge", {
  g1 <- line_graph_obs(obs_x, obs_y)
  g1$buildDirectionalConstraints(alpha = 1)
  g2 <- line_graph_obs(obs_x, obs_y, split = TRUE)
  g2$buildDirectionalConstraints(alpha = 1)

  lik <- function(g) {
    likelihood_alpha1_directional(theta_log, g, data_name = "y",
                                  parameterization = "spde")
  }
  expect_equal(lik(g1), lik(g2), tolerance = 1e-8)

  lik_pre <- function(g) {
    likelihood_alpha1_directional_precompute(
      theta_log, precompute_alpha1_directional(g, data_name = "y"),
      parameterization = "spde")
  }
  expect_equal(lik_pre(g1), lik_pre(g2), tolerance = 1e-8)
})

test_that("spde_covariance on a single edge matches the split edge", {
  kappa <- 2
  tau <- 1
  for (directional in c(FALSE, TRUE)) {
    alpha <- if (directional) 1 else 2
    m1 <- line_graph()
    m1$build_mesh(h = 0.1)
    m2 <- line_graph(split = TRUE)
    m2$build_mesh(h = 0.1)
    c1 <- spde_covariance(c(1, 0.2), kappa = kappa, tau = tau, alpha = alpha,
                          graph = m1, directional = directional)
    c2 <- spde_covariance(c(1, 0.4), kappa = kappa, tau = tau, alpha = alpha,
                          graph = m2, directional = directional)
    expect_equal(mesh_values_by_x(m1, c1), mesh_values_by_x(m2, c2),
                 tolerance = 1e-8)

    # Marginal variance at P: 1/(4 kappa^3 tau^2) for alpha = 2 and
    # 1/(2 kappa tau^2) for the directional (OU) alpha = 1 model.
    r0 <- if (alpha == 2) 1 / (4 * kappa^3 * tau^2) else 1 / (2 * kappa * tau^2)
    expect_equal(as.vector(c1)[which.min(abs(m1$mesh$V[, 1] - 0.2))], r0,
                 tolerance = 1e-8)
  }
})

test_that("spde_variance (directional) on a single edge matches the split edge", {
  m1 <- line_graph()
  m1$build_mesh(h = 0.1)
  m2 <- line_graph(split = TRUE)
  m2$build_mesh(h = 0.1)
  v1 <- spde_variance(kappa = 2, tau = 1, alpha = 1, graph = m1,
                      directional = TRUE)
  v2 <- spde_variance(kappa = 2, tau = 1, alpha = 1, graph = m2,
                      directional = TRUE)
  x1 <- round(m1$mesh$PtE[, 2], 8)
  x2 <- round((m2$mesh$PtE[, 1] - 1) / 2 + m2$mesh$PtE[, 2] / 2, 8)
  common <- intersect(x1, x2)
  expect_gt(length(common), 0)
  expect_equal(v1[match(common, x1)], v2[match(common, x2)], tolerance = 1e-8)
})

test_that("sample_spde (directional, alpha = 1) works on a single edge", {
  u <- sample_spde(kappa = 1, tau = 1, alpha = 1, directional = TRUE,
                   graph = line_graph(), PtE = cbind(1, c(0.2, 0.5)),
                   type = "manual")
  expect_length(u, 2)
  expect_true(all(is.finite(u)))
})

test_that("posterior means on a single edge match the split edge", {
  g1 <- line_graph_obs(obs_x, obs_y)
  g1$buildC(2)
  g2 <- line_graph_obs(obs_x, obs_y, split = TRUE)
  g2$buildC(2)
  pm2 <- function(g) {
    posterior_mean_obs_alpha2(theta_nat, g, resp = g$get_data()$y,
                              PtE_resp = g$get_PtE(), type = "obs")
  }
  expect_equal(pm2(g1), pm2(g2), tolerance = 1e-8)

  d1 <- line_graph_obs(obs_x, obs_y)
  d1$buildDirectionalConstraints(alpha = 1)
  d2 <- line_graph_obs(obs_x, obs_y, split = TRUE)
  d2$buildDirectionalConstraints(alpha = 1)
  pm1 <- function(g) {
    posterior_mean_obs_alpha1(theta_nat, g, resp = g$get_data()$y,
                              PtE_resp = g$get_PtE(), type = "obs",
                              directional = TRUE)
  }
  expect_equal(pm1(d1), pm1(d2), tolerance = 1e-8)
})

test_that("covariance-based alpha = 2 functions on a single edge match the split edge", {
  # These functions assume observations at vertices; observing only the two
  # end points keeps the single edge unsplit by observation_to_vertex().
  make <- function(split) {
    g <- line_graph_obs(c(0, 1), c(0.4, -0.7), split = split)
    g$observation_to_vertex()
    g$buildC(2)
    g
  }
  g1 <- make(FALSE)
  g2 <- make(TRUE)
  expect_equal(g1$nE, 1)
  expect_length(g1$CoB$S, 0)

  lik <- function(g) {
    likelihood_graph_covariance(g, model = "WM2", y_graph = g$get_data()$y,
                                repl = NULL, check_euclidean = FALSE)(theta_log)
  }
  expect_equal(lik(g1), lik(g2), tolerance = 1e-8)

  lik_pre <- function(g) {
    pre <- precompute_graph_covariance(g, model = "WM2",
                                       y_graph = g$get_data()$y,
                                       check_euclidean = FALSE)
    likelihood_graph_covariance_precompute(theta_log, pre)
  }
  expect_equal(lik_pre(g1), lik_pre(g2), tolerance = 1e-8)

  for (cv_fun in list(posterior_crossvalidation_manual,
                      posterior_crossvalidation_covariance_manual)) {
    expect_equal(cv_fun(theta_nat, make(FALSE), data_name = "y", model = "alpha2")$mu,
                 cv_fun(theta_nat, make(TRUE), data_name = "y", model = "alpha2")$mu,
                 tolerance = 1e-8)
  }
})

test_that("graph_lme predict and LOO CV on a single edge match the split edge", {
  skip_on_cran()
  # Observations and predictions only at the end points (replicated), so that
  # observation_to_vertex() inside predict/CV keeps the graph a single edge.
  n_rep <- 6
  y_rep <- c(-0.84, 1.38, -1.26, 0.07, 1.71, -0.6, -0.47, -0.64, -0.29, 0.14,
             0.78, 0.03)
  make <- function(split) {
    g <- line_graph(split)
    df <- line_PtE(rep(c(0, 1), n_rep), split)
    df$y <- y_rep
    df$repl <- rep(seq_len(n_rep), each = 2)
    g$add_observations(data = df, normalized = TRUE, group = "repl", verbose = 0)
    g
  }
  newdata <- function(split) cbind(line_PtE(c(0, 1), split), repl = 1)

  for (model in c("WM2", "WMD1")) {
    fit1 <- suppressWarnings(graph_lme(y ~ -1, graph = make(FALSE),
                                       model = model, parallel = FALSE))
    fit2 <- suppressWarnings(graph_lme(y ~ -1, graph = make(TRUE),
                                       model = model, parallel = FALSE))
    expect_equal(fit1$graph$nE, 1)
    # Use identical parameters so the two fits are directly comparable.
    fit2$coeff <- fit1$coeff

    p1 <- predict(fit1, newdata = newdata(FALSE), normalized = TRUE,
                  compute_variances = TRUE)
    p2 <- predict(fit2, newdata = newdata(TRUE), normalized = TRUE,
                  compute_variances = TRUE)
    expect_equal(p1$mean, p2$mean, tolerance = 1e-8)
    expect_equal(p1$variance, p2$variance, tolerance = 1e-8)

    cv1 <- posterior_crossvalidation_loo(fit1, tibble = FALSE)
    cv2 <- posterior_crossvalidation_loo(fit2, tibble = FALSE)
    expect_equal(cv1$mu, cv2$mu, tolerance = 1e-8)
    expect_equal(cv1$var, cv2$var, tolerance = 1e-8)

    if (model == "WM2") {
      expect_equal(as.matrix(get_covariance_precision(fit1)$prec_cov),
                   as.matrix(get_covariance_precision(fit2)$prec_cov),
                   tolerance = 1e-8)
    }
  }
})
