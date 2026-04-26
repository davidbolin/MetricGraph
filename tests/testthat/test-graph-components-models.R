## Unit tests for graph_components compatibility with the modelling pipeline:
##   as_metric_graph(), graph_lme(), graph_spde(), sample_spde(),
##   spde_precision(), spde_variance(), make_Q_spacetime().

library(MetricGraph)
library(testthat)

# Helper: two well-separated components, two collinear edges each.
make_two_components <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(1, 0), c(2, 0))
  edge3 <- rbind(c(10, 10), c(11, 10))
  edge4 <- rbind(c(11, 10), c(12, 10))
  graph_components$new(edges = list(edge1, edge2, edge3, edge4),
                       verbose = 0, perform_merges = TRUE)
}

# Helper: matched single metric_graph that should give numerically identical
# fits to the components-based version (same edges in the same order, plus
# `check_connected = FALSE`).
make_matched_metric_graph <- function() {
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

## ── as_metric_graph ──────────────────────────────────────────────────────────

test_that("as_metric_graph builds a single graph with combined edges and data", {
  gc <- make_two_components()
  df <- make_obs_df()
  attach_observations(gc, df)

  mg <- gc$as_metric_graph()
  expect_s3_class(mg, "metric_graph")
  expect_equal(mg$nE, 4L)
  expect_equal(mg$nV, 6L)

  obs <- mg$get_data()
  expect_equal(nrow(obs), nrow(df))
  # Component 2 observations should now live on edges 3 and 4 of the combined
  # graph.
  expect_true(all(obs$.edge_number[obs$.coord_x > 5] %in% c(3L, 4L)))
  expect_true(all(obs$.edge_number[obs$.coord_x < 5] %in% c(1L, 2L)))
})

test_that("as_metric_graph block-stacks mesh and FEM matrices when present", {
  gc <- make_two_components()
  attach_observations(gc, make_obs_df())
  gc$build_mesh(h = 0.25)
  gc$compute_fem()

  mg <- gc$as_metric_graph()
  expect_false(is.null(mg$mesh))
  expect_false(is.null(mg$mesh$C))
  expect_false(is.null(mg$mesh$G))
  expect_equal(nrow(mg$mesh$V), nrow(mg$mesh$C))
  # FEM C should be the bdiag of per-component mesh$C — sum of per-component
  # nrow equals combined nrow.
  expect_equal(nrow(mg$mesh$C),
               sum(vapply(gc$graphs, function(g) nrow(g$mesh$C), integer(1))))
})

test_that("as_metric_graph returns no mesh when components have no mesh", {
  gc <- make_two_components()
  mg <- gc$as_metric_graph()
  expect_null(mg$mesh)
})

test_that("as_metric_graph marks the result as disconnected", {
  gc <- make_two_components()
  mg <- gc$as_metric_graph()
  expect_true(mg$is_disconnected())

  # Regular metric_graph instances are not disconnected
  g <- metric_graph$new(edges = list(rbind(c(0, 0), c(1, 0))),
                         check_connected = FALSE, verbose = 0)
  expect_false(g$is_disconnected())
})

test_that("add_observations is rejected on a disconnected metric_graph", {
  gc <- make_two_components()
  attach_observations(gc, make_obs_df())
  mg <- gc$as_metric_graph()

  expect_error(
    mg$add_observations(
      data = data.frame(edge_number = 1, distance_on_edge = 0.5, y = 1),
      verbose = 0
    ),
    "assembled from a 'graph_components'"
  )

  # Internal callers can bypass via the .allow_disconnected argument
  expect_silent(
    mg$add_observations(
      data = data.frame(edge_number = 1, distance_on_edge = 0.5, y = 1),
      verbose = 0, suppress_warnings = TRUE,
      .allow_disconnected = TRUE
    )
  )
})

test_that("as_metric_graph yields fields equivalent to the constructor route", {
  gc <- make_two_components()
  attach_observations(gc, make_obs_df(seed = 11))
  gc$build_mesh(h = 0.25)
  gc$compute_fem()

  mg_fast <- gc$as_metric_graph()

  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(2, 0)),
    rbind(c(10, 10), c(11, 10)),
    rbind(c(11, 10), c(12, 10))
  )
  mg_slow <- metric_graph$new(edges = edges, perform_merges = FALSE,
                              check_connected = FALSE, verbose = 0)
  mg_slow$add_observations(data = make_obs_df(seed = 11),
                           data_coords = "spatial", verbose = 0)
  mg_slow$build_mesh(h = 0.25)
  mg_slow$compute_fem()

  expect_equal(mg_fast$nE, mg_slow$nE)
  expect_equal(mg_fast$nV, mg_slow$nV)
  expect_equal(unname(mg_fast$V), unname(mg_slow$V))
  expect_equal(unname(mg_fast$E), unname(mg_slow$E))
  expect_equal(as.numeric(mg_fast$edge_lengths),
               as.numeric(mg_slow$edge_lengths))
  expect_equal(dim(mg_fast$mesh$V), dim(mg_slow$mesh$V))
  expect_equal(dim(mg_fast$mesh$C), dim(mg_slow$mesh$C))
})

test_that("as_metric_graph preserves edge weights", {
  edges <- list(rbind(c(0, 0), c(1, 0)),
                rbind(c(1, 0), c(2, 0)),
                rbind(c(10, 10), c(11, 10)))
  ew <- c(0.5, 1.5, 2.5)
  gc <- graph_components$new(edges = edges, edge_weights = ew, verbose = 0)
  mg <- gc$as_metric_graph()
  expect_equal(mg$nE, 3L)
  expect_setequal(unname(mg$get_edge_weights()$.weights), ew)
})

## ── graph_lme ────────────────────────────────────────────────────────────────

test_that("graph_lme accepts a graph_components object (WM1)", {
  gc <- make_two_components()
  attach_observations(gc, make_obs_df())
  fit <- suppressWarnings(graph_lme(y ~ 1, graph = gc, model = "WM1"))
  expect_s3_class(fit, "graph_lme")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

test_that("graph_lme on graph_components matches identical metric_graph", {
  gc <- make_two_components()
  mg <- make_matched_metric_graph()

  df <- make_obs_df(seed = 7)
  attach_observations(gc, df)
  attach_observations(mg, df)

  fit_gc <- suppressWarnings(graph_lme(y ~ 1, graph = gc, model = "WM1"))
  fit_mg <- suppressWarnings(graph_lme(y ~ 1, graph = mg, model = "WM1"))

  expect_equal(as.numeric(stats::logLik(fit_gc)),
               as.numeric(stats::logLik(fit_mg)),
               tolerance = 1e-4)
  expect_equal(unname(fit_gc$coeff$fixed_effects),
               unname(fit_mg$coeff$fixed_effects),
               tolerance = 1e-4)
  expect_equal(unname(fit_gc$coeff$random_effects),
               unname(fit_mg$coeff$random_effects),
               tolerance = 1e-3)
})

test_that("graph_lme errors clearly when given a non-graph object", {
  expect_error(graph_lme(y ~ 1, graph = list(), model = "lm"),
               "metric_graph")
})

test_that("predict on a graph_lme fit from graph_components returns finite means", {
  gc <- make_two_components()
  attach_observations(gc, make_obs_df())
  fit <- suppressWarnings(graph_lme(y ~ 1, graph = gc, model = "WM1"))

  # Edge numbering after as_metric_graph: comp1 edges 1-2, comp2 edges 3-4.
  newdata <- data.frame(edge_number = c(1L, 3L),
                        distance_on_edge = c(0.5, 0.5))
  pred <- predict(fit, newdata = newdata, normalized = TRUE)
  expect_length(pred$mean, 2L)
  expect_true(all(is.finite(pred$mean)))
})

## ── graph_spde ───────────────────────────────────────────────────────────────

test_that("graph_spde accepts a graph_components object", {
  gc <- make_two_components()
  attach_observations(gc, make_obs_df())
  spde <- graph_spde(graph_object = gc, alpha = 1, verbose = 0)
  expect_true(inherits(spde, "inla_metric_graph_spde"))
})

## ── sample_spde, spde_precision, spde_variance ───────────────────────────────

test_that("sample_spde works on a graph_components mesh", {
  gc <- make_two_components()
  gc$build_mesh(h = 0.2)
  set.seed(1)
  samp <- sample_spde(graph = gc, alpha = 1, kappa = 5, tau = 1,
                      type = "mesh")
  expect_true(length(samp) > 0)
  expect_true(all(is.finite(samp)))
})

test_that("spde_precision returns a square matrix on graph_components", {
  gc <- make_two_components()
  Q <- spde_precision(kappa = 1, tau = 1, alpha = 1, graph = gc, BC = 1)
  expect_equal(nrow(Q), ncol(Q))
  expect_equal(nrow(Q), gc$graphs[[1]]$nV + gc$graphs[[2]]$nV)
})

test_that("spde_variance returns finite variances on graph_components", {
  gc <- make_two_components()
  gc$build_mesh(h = 0.2)
  gc$compute_fem()
  v <- spde_variance(kappa = 1, tau = 1, alpha = 1, graph = gc, BC = 1)
  expect_true(length(v) > 0)
  expect_true(all(is.finite(v)))
})

## ── make_Q_spacetime ─────────────────────────────────────────────────────────

test_that("make_Q_spacetime accepts a graph_components mesh", {
  gc <- make_two_components()
  gc$build_mesh(h = 0.2)
  gc$compute_fem()
  Q <- make_Q_spacetime(graph = gc, t = seq(0, 1, length.out = 5),
                        kappa = 1, rho = 1, gamma = 1,
                        alpha = 2, beta = 1, sigma = 1)
  expect_equal(nrow(Q), ncol(Q))
})
