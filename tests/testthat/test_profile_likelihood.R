library(testthat)
library(MetricGraph)

# Tests for the profiled-beta likelihoods in R/graph_likelihoods_v2.R:
#   likelihood_alpha1_profile / _precompute
#   likelihood_alpha2_profile / _precompute
#   likelihood_alpha1_directional_profile / _precompute
#   profile_beta_estimate
#
# theta is always length 3, log scale: (log sigma_e, log reciprocal_tau,
# log kappa), used with parameterization = "spde" so kappa = exp(theta[3]).
# All functions under test are internal (MetricGraph:::).

# ---- graph helpers ---------------------------------------------------

# 4-edge graph (3 straight edges + 1 arc), as in test_precomputed_likelihood.R
make_graph_alpha <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(0, 0), c(0, 1))
  edge3 <- rbind(c(0, 1), c(-1, 1))
  theta_seq <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
  edge4 <- cbind(sin(theta_seq), 1 + cos(theta_seq))
  edges <- list(edge1, edge2, edge3, edge4)
  metric_graph$new(edges = edges)
}

# 3-edge graph oriented into a confluence, as in test_directional.R
make_graph_directional <- function() {
  edge1 <- rbind(c(1, 0), c(0, 0))
  edge2 <- rbind(c(1, 1), c(1, 0))
  edge3 <- rbind(c(1, -1), c(1, 0))
  edges <- list(edge1, edge2, edge3)
  metric_graph$new(edges = edges)
}

# Adds n_per_edge observations per edge (strictly inside (0,1)) to graph
# and returns the y / X_cov vectors used (aligned with add_observations'
# internal ordering since positions are already sorted per edge).
add_test_observations <- function(graph, n_per_edge = 10, n_cov = TRUE) {
  nE <- graph$nE
  y_list <- vector("list", nE)
  x_list <- vector("list", nE)
  edge_list <- vector("list", nE)
  dist_list <- vector("list", nE)

  for (e in seq_len(nE)) {
    pos <- sort(runif(n_per_edge, 0.02, 0.98))
    y_e <- rnorm(n_per_edge)
    edge_list[[e]] <- rep(e, n_per_edge)
    dist_list[[e]] <- pos
    y_list[[e]] <- y_e
  }

  y <- unlist(y_list)
  edge_number <- unlist(edge_list)
  distance_on_edge <- unlist(dist_list)
  n <- length(y)
  X_cov <- cbind(1, runif(n))

  df_data <- data.frame(y = y, edge_number = edge_number,
                        distance_on_edge = distance_on_edge)

  graph$clear_observations()
  graph$add_observations(data = df_data, normalized = TRUE)

  list(y = y, X_cov = X_cov, n = n)
}

theta3 <- c(log(0.2), log(0.8), log(1.3))

# ---- 1. profile equals joint at beta_hat ------------------------------

test_that("profile equals joint at beta_hat", {
  set.seed(123)

  ## alpha1
  graph <- make_graph_alpha()
  obs <- add_test_observations(graph)
  y <- obs$y; X <- obs$X_cov

  bhat <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha1",
                                              graph = graph, manual_y = y,
                                              X_cov = X, repl = NULL, BC = 1,
                                              parameterization = "spde")$beta

  prof <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                  manual_y = y, X_cov = X,
                                                  repl = NULL, BC = 1,
                                                  parameterization = "spde")
  joint <- MetricGraph:::likelihood_alpha1(c(theta3, bhat), graph = graph,
                                           manual_y = y, X_cov = X,
                                           repl = NULL, BC = 1,
                                           parameterization = "spde")
  expect_equal(prof, joint, tolerance = 1e-8)

  ## alpha2
  graph2 <- make_graph_alpha()
  obs2 <- add_test_observations(graph2)
  y2 <- obs2$y; X2 <- obs2$X_cov

  bhat2 <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha2",
                                               graph = graph2, manual_y = y2,
                                               X_cov = X2, repl = NULL, BC = 1,
                                               parameterization = "spde")$beta

  prof2 <- MetricGraph:::likelihood_alpha2_profile(theta3, graph = graph2,
                                                   manual_y = y2, X_cov = X2,
                                                   repl = NULL, BC = 1,
                                                   parameterization = "spde")
  joint2 <- MetricGraph:::likelihood_alpha2(c(theta3, bhat2), graph = graph2,
                                            manual_y = y2, X_cov = X2,
                                            repl = NULL, BC = 1,
                                            parameterization = "spde")
  expect_equal(prof2, joint2, tolerance = 1e-8)

  ## alpha1 directional
  graphd <- make_graph_directional()
  obsd <- add_test_observations(graphd)
  yd <- obsd$y; Xd <- obsd$X_cov

  bhatd <- MetricGraph:::profile_beta_estimate(theta3,
                                               model = "alpha1_directional",
                                               graph = graphd, manual_y = yd,
                                               X_cov = Xd, repl = NULL,
                                               parameterization = "spde")$beta

  profd <- MetricGraph:::likelihood_alpha1_directional_profile(
    theta3, graph = graphd, manual_y = yd, X_cov = Xd, repl = NULL,
    parameterization = "spde")
  jointd <- MetricGraph:::likelihood_alpha1_directional(
    c(theta3, bhatd), graph = graphd, manual_y = yd, X_cov = Xd, repl = NULL,
    parameterization = "spde")
  expect_equal(profd, jointd, tolerance = 1e-8)
})

# ---- 2. profile dominates joint in beta -------------------------------

test_that("profile dominates joint in beta", {
  set.seed(123)

  ## alpha1
  graph <- make_graph_alpha()
  obs <- add_test_observations(graph)
  y <- obs$y; X <- obs$X_cov
  bhat <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha1",
                                              graph = graph, manual_y = y,
                                              X_cov = X, repl = NULL, BC = 1,
                                              parameterization = "spde")$beta
  prof <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                  manual_y = y, X_cov = X,
                                                  repl = NULL, BC = 1,
                                                  parameterization = "spde")
  for (k in 1:3) {
    delta <- 0.5 * rnorm(length(bhat))
    joint <- MetricGraph:::likelihood_alpha1(c(theta3, bhat + delta),
                                             graph = graph, manual_y = y,
                                             X_cov = X, repl = NULL, BC = 1,
                                             parameterization = "spde")
    expect_true(prof >= joint - 1e-10)
  }

  ## alpha2
  graph2 <- make_graph_alpha()
  obs2 <- add_test_observations(graph2)
  y2 <- obs2$y; X2 <- obs2$X_cov
  bhat2 <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha2",
                                               graph = graph2, manual_y = y2,
                                               X_cov = X2, repl = NULL, BC = 1,
                                               parameterization = "spde")$beta
  prof2 <- MetricGraph:::likelihood_alpha2_profile(theta3, graph = graph2,
                                                   manual_y = y2, X_cov = X2,
                                                   repl = NULL, BC = 1,
                                                   parameterization = "spde")
  for (k in 1:3) {
    delta <- 0.5 * rnorm(length(bhat2))
    joint2 <- MetricGraph:::likelihood_alpha2(c(theta3, bhat2 + delta),
                                              graph = graph2, manual_y = y2,
                                              X_cov = X2, repl = NULL, BC = 1,
                                              parameterization = "spde")
    expect_true(prof2 >= joint2 - 1e-10)
  }

  ## alpha1 directional
  graphd <- make_graph_directional()
  obsd <- add_test_observations(graphd)
  yd <- obsd$y; Xd <- obsd$X_cov
  bhatd <- MetricGraph:::profile_beta_estimate(theta3,
                                               model = "alpha1_directional",
                                               graph = graphd, manual_y = yd,
                                               X_cov = Xd, repl = NULL,
                                               parameterization = "spde")$beta
  profd <- MetricGraph:::likelihood_alpha1_directional_profile(
    theta3, graph = graphd, manual_y = yd, X_cov = Xd, repl = NULL,
    parameterization = "spde")
  for (k in 1:3) {
    delta <- 0.5 * rnorm(length(bhatd))
    jointd <- MetricGraph:::likelihood_alpha1_directional(
      c(theta3, bhatd + delta), graph = graphd, manual_y = yd, X_cov = Xd,
      repl = NULL, parameterization = "spde")
    expect_true(profd >= jointd - 1e-10)
  }
})

# ---- 3. beta_hat matches numerical optimum (alpha1 only) --------------

test_that("beta_hat matches numerical optimum", {
  set.seed(123)
  graph <- make_graph_alpha()
  obs <- add_test_observations(graph)
  y <- obs$y; X <- obs$X_cov

  bhat <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha1",
                                              graph = graph, manual_y = y,
                                              X_cov = X, repl = NULL, BC = 1,
                                              parameterization = "spde")$beta
  prof <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                  manual_y = y, X_cov = X,
                                                  repl = NULL, BC = 1,
                                                  parameterization = "spde")

  fn <- function(b) {
    -MetricGraph:::likelihood_alpha1(c(theta3, b), graph = graph,
                                     manual_y = y, X_cov = X, repl = NULL,
                                     BC = 1, parameterization = "spde")
  }
  opt <- optim(c(0, 0), fn, method = "BFGS",
              control = list(reltol = 1e-14, maxit = 500))

  expect_equal(opt$par, bhat, tolerance = 1e-4)
  expect_equal(-opt$value, prof, tolerance = 1e-6)
})

# ---- 4. precompute equals non-precompute ------------------------------

test_that("precompute equals non-precompute", {
  set.seed(123)

  ## alpha1
  graph <- make_graph_alpha()
  obs <- add_test_observations(graph)
  y <- obs$y; X <- obs$X_cov

  pcomp1 <- MetricGraph:::precompute_alpha1(graph, manual_y = y, X_cov = X,
                                            repl = NULL)
  for (reml in c(FALSE, TRUE)) {
    prof <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                    manual_y = y, X_cov = X,
                                                    repl = NULL, BC = 1,
                                                    parameterization = "spde",
                                                    reml = reml)
    prof_pre <- MetricGraph:::likelihood_alpha1_profile_precompute(
      theta3, graph = graph, precomputeddata = pcomp1, BC = 1,
      parameterization = "spde", reml = reml)
    expect_equal(prof, prof_pre, tolerance = 1e-10)
  }

  ## alpha2
  graph2 <- make_graph_alpha()
  obs2 <- add_test_observations(graph2)
  y2 <- obs2$y; X2 <- obs2$X_cov

  pcomp2 <- MetricGraph:::precompute_alpha2(graph2, manual_y = y2,
                                            X_cov = X2, repl = NULL)
  for (reml in c(FALSE, TRUE)) {
    prof2 <- MetricGraph:::likelihood_alpha2_profile(theta3, graph = graph2,
                                                     manual_y = y2, X_cov = X2,
                                                     repl = NULL, BC = 1,
                                                     parameterization = "spde",
                                                     reml = reml)
    prof2_pre <- MetricGraph:::likelihood_alpha2_profile_precompute(
      theta3, precomputed_data = pcomp2, BC = 1, parameterization = "spde",
      reml = reml)
    expect_equal(prof2, prof2_pre, tolerance = 1e-10)
  }

  ## alpha1 directional
  graphd <- make_graph_directional()
  obsd <- add_test_observations(graphd)
  yd <- obsd$y; Xd <- obsd$X_cov

  pcompd <- MetricGraph:::precompute_alpha1_directional(graphd,
                                                        manual_y = yd,
                                                        X_cov = Xd,
                                                        repl = NULL)
  for (reml in c(FALSE, TRUE)) {
    profd <- MetricGraph:::likelihood_alpha1_directional_profile(
      theta3, graph = graphd, manual_y = yd, X_cov = Xd, repl = NULL,
      parameterization = "spde", reml = reml)
    profd_pre <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
      theta3, precomputed_data = pcompd, parameterization = "spde",
      reml = reml)
    expect_equal(profd, profd_pre, tolerance = 1e-10)
  }
})

# ---- 5. REML equals profile minus half logdet H (finite differences) --

test_that("REML equals profile minus half logdet H (finite differences)", {
  set.seed(123)
  graph <- make_graph_alpha()
  obs <- add_test_observations(graph)
  y <- obs$y; X <- obs$X_cov

  bhat <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha1",
                                              graph = graph, manual_y = y,
                                              X_cov = X, repl = NULL, BC = 1,
                                              parameterization = "spde")$beta

  # The joint log-likelihood in beta is exactly quadratic, so exact finite
  # differences (any step d) recover the Hessian H = -d^2 l / d beta^2.
  lfun <- function(b) {
    MetricGraph:::likelihood_alpha1(c(theta3, b), graph = graph,
                                    manual_y = y, X_cov = X, repl = NULL,
                                    BC = 1, parameterization = "spde")
  }

  p <- length(bhat)
  d <- 0.1
  l0 <- lfun(bhat)
  H_fd <- matrix(0, p, p)
  for (i in 1:p) {
    ei <- rep(0, p); ei[i] <- d
    for (j in 1:p) {
      ej <- rep(0, p); ej[j] <- d
      H_fd[i, j] <- -(lfun(bhat + ei + ej) - lfun(bhat + ei) -
                        lfun(bhat + ej) + l0) / d^2
    }
  }

  prof <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                  manual_y = y, X_cov = X,
                                                  repl = NULL, BC = 1,
                                                  parameterization = "spde",
                                                  reml = FALSE)
  reml_val <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                      manual_y = y, X_cov = X,
                                                      repl = NULL, BC = 1,
                                                      parameterization = "spde",
                                                      reml = TRUE)

  expected_reml <- prof - 0.5 * determinant(H_fd)$modulus
  expect_equal(as.numeric(reml_val), as.numeric(expected_reml),
              tolerance = 1e-6)
})

# ---- 6. two replicates (alpha1) ---------------------------------------

test_that("two replicates", {
  set.seed(123)
  # add_observations(group = <colname>) stacks a ".group" column from the
  # named column of `data` (confirmed in R/metric_graph.R roxygen for
  # add_observations and used e.g. in tests/testthat/test_obs_improvements.R
  # and test_posterior_crossvalidation.R as
  # `graph$add_observations(data = df, group = "repl", normalized = TRUE)`).
  # We build a long-format data frame with the same locations repeated for
  # group values 1 and 2, each with an independent y draw.
  graph <- make_graph_alpha()
  nE <- graph$nE
  n_per_edge <- 10

  edge_number <- NULL
  distance_on_edge <- NULL
  for (e in seq_len(nE)) {
    pos <- sort(runif(n_per_edge, 0.02, 0.98))
    edge_number <- c(edge_number, rep(e, n_per_edge))
    distance_on_edge <- c(distance_on_edge, pos)
  }
  n_loc <- length(edge_number)

  y1 <- rnorm(n_loc)
  y2 <- rnorm(n_loc)

  df_data <- data.frame(
    y = c(y1, y2),
    edge_number = rep(edge_number, 2),
    distance_on_edge = rep(distance_on_edge, 2),
    repl = rep(c(1, 2), each = n_loc)
  )

  graph$clear_observations()
  graph$add_observations(data = df_data, normalized = TRUE, group = "repl")

  y_full <- c(y1, y2)
  n <- length(y_full)
  X_full <- cbind(1, runif(n))

  bhat <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha1",
                                              graph = graph,
                                              manual_y = y_full,
                                              X_cov = X_full, repl = NULL,
                                              BC = 1,
                                              parameterization = "spde")$beta

  prof <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                  manual_y = y_full,
                                                  X_cov = X_full, repl = NULL,
                                                  BC = 1,
                                                  parameterization = "spde")
  joint <- MetricGraph:::likelihood_alpha1(c(theta3, bhat), graph = graph,
                                           manual_y = y_full, X_cov = X_full,
                                           repl = NULL, BC = 1,
                                           parameterization = "spde")
  expect_equal(prof, joint, tolerance = 1e-8)
})

# ---- 7. no covariates (alpha1) -----------------------------------------

test_that("no covariates", {
  set.seed(123)
  graph <- make_graph_alpha()
  obs <- add_test_observations(graph)
  y <- obs$y

  prof <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                  manual_y = y, X_cov = NULL,
                                                  repl = NULL, BC = 1,
                                                  parameterization = "spde")
  joint <- MetricGraph:::likelihood_alpha1(theta3, graph = graph,
                                           manual_y = y, X_cov = NULL,
                                           repl = NULL, BC = 1,
                                           parameterization = "spde")
  expect_equal(prof, joint, tolerance = 1e-8)

  prof_reml <- MetricGraph:::likelihood_alpha1_profile(theta3, graph = graph,
                                                       manual_y = y,
                                                       X_cov = NULL,
                                                       repl = NULL, BC = 1,
                                                       parameterization = "spde",
                                                       reml = TRUE)
  expect_equal(prof_reml, prof, tolerance = 1e-10)

  bs <- MetricGraph:::profile_beta_estimate(theta3, model = "alpha1",
                                            graph = graph, manual_y = y,
                                            X_cov = NULL, repl = NULL, BC = 1,
                                            parameterization = "spde")
  expect_length(bs$beta, 0)
})

# ---- theta length validation --------------------------------------------

test_that("theta length validation", {
  set.seed(123)
  graph <- make_graph_alpha()
  obs <- add_test_observations(graph)
  y <- obs$y; X <- obs$X_cov

  expect_error(
    MetricGraph:::likelihood_alpha1_profile(c(theta3, 0.5), graph = graph,
                                            manual_y = y, X_cov = X,
                                            repl = NULL, BC = 1,
                                            parameterization = "spde")
  )
})
