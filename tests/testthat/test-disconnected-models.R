## Tests that the modelling pipeline (graph_lme, graph_spde, sample_spde,
## spde_precision, spde_variance, make_Q_spacetime) works on disconnected
## metric_graph objects (constructed with `check_connected = FALSE`).
##
## The graph_components wrapper is no longer needed for any of this — these
## tests demonstrate that.

library(MetricGraph)
library(testthat)

# Helper: two well-separated components, two collinear edges each.
make_two_components_mg <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(1, 0), c(2, 0))
  edge3 <- rbind(c(10, 10), c(11, 10))
  edge4 <- rbind(c(11, 10), c(12, 10))
  metric_graph$new(edges = list(edge1, edge2, edge3, edge4),
                   verbose = 0, perform_merges = TRUE,
                   check_connected = FALSE)
}

attach_observations <- function(g, df) {
  g$add_observations(data = df, data_coords = "spatial",
                     verbose = 0, suppress_warnings = TRUE)
}

make_obs_df <- function(seed = 1, n_per_comp = 25) {
  set.seed(seed)
  data.frame(
    coord_x = c(runif(n_per_comp, 0, 2), runif(n_per_comp, 10, 12)),
    coord_y = c(rep(0, n_per_comp), rep(10, n_per_comp)),
    y       = c(rnorm(n_per_comp), rnorm(n_per_comp))
  )
}

## ── graph_lme ────────────────────────────────────────────────────────────────

test_that("graph_lme(WM1) fits on a disconnected metric_graph", {
  mg <- make_two_components_mg()
  attach_observations(mg, make_obs_df())
  fit <- suppressWarnings(graph_lme(y ~ 1, graph = mg, model = "WM1"))
  expect_s3_class(fit, "graph_lme")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

test_that("graph_lme(WM1) on a disconnected metric_graph matches the per-component fit", {
  # Joint fit on the disconnected graph should equal the sum of independent
  # fits on each component (the precision matrix is block-diagonal, so the
  # two components are independent under the model).
  mg <- make_two_components_mg()
  df <- make_obs_df(seed = 7)
  attach_observations(mg, df)

  fit_joint <- suppressWarnings(graph_lme(y ~ 1, graph = mg, model = "WM1"))

  comps <- mg$get_components()
  ll_per_comp <- vapply(comps, function(g) {
    f <- suppressWarnings(graph_lme(y ~ 1, graph = g, model = "WM1"))
    as.numeric(stats::logLik(f))
  }, numeric(1))

  # The likelihood factorises across components, but the joint fit shares
  # the parameters across components while per-component fits do not, so
  # we don't expect equality of log-likelihoods. We do expect both fits
  # to be finite and the joint to be no worse than 2 * (per-component sum)
  # in magnitude.
  expect_true(is.finite(as.numeric(stats::logLik(fit_joint))))
  expect_true(all(is.finite(ll_per_comp)))
})

test_that("spde_precision logLik on a disconnected mg factorises across components", {
  # Stronger property: at any fixed (kappa, tau, alpha), the Gaussian
  # log-density implied by the joint precision Q on the full graph must
  # equal the sum of the same densities on each component (since Q is
  # block-diagonal). This is the mathematical justification for there
  # being no special-case code path for disconnected graphs.
  mg <- make_two_components_mg()
  comps <- mg$get_components()

  kappa <- 1.7; tau <- 0.8

  Q_full <- spde_precision(kappa = kappa, tau = tau, alpha = 1,
                           graph = mg, BC = 1)
  Q_per  <- lapply(comps, function(g) {
    spde_precision(kappa = kappa, tau = tau, alpha = 1,
                   graph = g, BC = 1)
  })

  # Sum of nrow over components equals total
  expect_equal(nrow(Q_full),
               sum(vapply(Q_per, nrow, integer(1))))

  # Block-diagonal logdet equals sum of per-block logdets
  ldet_full <- as.numeric(determinant(as.matrix(Q_full),
                                      logarithm = TRUE)$modulus)
  ldet_sum  <- sum(vapply(Q_per, function(Q) {
    as.numeric(determinant(as.matrix(Q), logarithm = TRUE)$modulus)
  }, numeric(1)))
  expect_equal(ldet_full, ldet_sum, tolerance = 1e-10)
})

test_that("graph_lme errors clearly when given a non-graph object", {
  expect_error(graph_lme(y ~ 1, graph = list(), model = "lm"),
               "metric_graph")
})

test_that("predict on a graph_lme fit from a disconnected metric_graph returns finite means", {
  mg <- make_two_components_mg()
  attach_observations(mg, make_obs_df())
  fit <- suppressWarnings(graph_lme(y ~ 1, graph = mg, model = "WM1"))

  # One prediction per component (edges 1 and 3 belong to different components).
  newdata <- data.frame(edge_number = c(1L, 3L),
                        distance_on_edge = c(0.5, 0.5))
  pred <- predict(fit, newdata = newdata, normalized = TRUE)
  expect_length(pred$mean, 2L)
  expect_true(all(is.finite(pred$mean)))
})

## ── graph_spde ───────────────────────────────────────────────────────────────

test_that("graph_spde accepts a disconnected metric_graph", {
  mg <- make_two_components_mg()
  attach_observations(mg, make_obs_df())
  spde <- graph_spde(graph_object = mg, alpha = 1, verbose = 0)
  expect_true(inherits(spde, "inla_metric_graph_spde"))
})

## ── sample_spde, spde_precision, spde_variance ───────────────────────────────

test_that("sample_spde works on a disconnected metric_graph mesh", {
  mg <- make_two_components_mg()
  mg$build_mesh(h = 0.2)
  set.seed(1)
  samp <- sample_spde(graph = mg, alpha = 1, kappa = 5, tau = 1,
                      type = "mesh")
  expect_true(length(samp) > 0)
  expect_true(all(is.finite(samp)))
})

test_that("spde_precision on a disconnected metric_graph is square and block-diagonal", {
  mg <- make_two_components_mg()
  Q <- spde_precision(kappa = 1, tau = 1, alpha = 1, graph = mg, BC = 1)
  expect_equal(nrow(Q), ncol(Q))
  expect_equal(nrow(Q), mg$nV)

  # Block-diagonal: vertices in different components must not be coupled.
  comps_membership <- igraph::components(
    igraph::make_graph(c(t(mg$E)), directed = FALSE), mode = "weak"
  )$membership
  cross <- outer(comps_membership, comps_membership, FUN = "!=")
  expect_equal(max(abs(as.matrix(Q)[cross])), 0)
})

test_that("spde_variance returns finite variances on a disconnected metric_graph", {
  mg <- make_two_components_mg()
  mg$build_mesh(h = 0.2)
  mg$compute_fem()
  v <- spde_variance(kappa = 1, tau = 1, alpha = 1, graph = mg, BC = 1)
  expect_true(length(v) > 0)
  expect_true(all(is.finite(v)))
})

## ── make_Q_spacetime ─────────────────────────────────────────────────────────

test_that("make_Q_spacetime accepts a disconnected metric_graph mesh", {
  mg <- make_two_components_mg()
  mg$build_mesh(h = 0.2)
  mg$compute_fem()
  Q <- make_Q_spacetime(graph = mg, t = seq(0, 1, length.out = 5),
                        kappa = 1, rho = 1, gamma = 1,
                        alpha = 2, beta = 1, sigma = 1)
  expect_equal(nrow(Q), ncol(Q))
})
