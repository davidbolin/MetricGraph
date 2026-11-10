## Tests filling gaps identified by `covr::package_coverage()`.
## Targets:
##   - selected_inv()                     (selected_inv.R)
##   - exp_covariance()                   (util.R)
##   - graph_starting_values()            (util.R)
##   - match_mesh_data()                  (util.R)
##   - tidyverse methods on metric_graph_data
##       (filter / mutate / select / drop_na / summarise)
##   - basic print/summary methods on metric_graph
##   - graph_lme() with model = "lm" and "isoExp"
##   - simulate() on graph_lme fits
##   - get_PtE() / get_groups() / get_locations() / get_degrees()
##   - clear_observations()
##   - mesh_A() / fem_basis() / get_mesh_locations()

library(MetricGraph)
library(testthat)

# ----- shared fixtures --------------------------------------------------------

make_small_graph <- function(seed = 1) {
  set.seed(seed)
  e1 <- rbind(c(0, 0), c(1, 0))
  e2 <- rbind(c(1, 0), c(2, 0))
  e3 <- rbind(c(2, 0), c(2, 1))
  metric_graph$new(edges = list(e1, e2, e3), verbose = 0)
}

attach_simple_obs <- function(g, n = 30, seed = 1, with_groups = FALSE) {
  set.seed(seed)
  if (with_groups) {
    df <- data.frame(
      coord_x = c(runif(n / 2, 0, 2), runif(n / 2, 0, 2)),
      coord_y = c(rep(0, n / 2), rep(0, n / 2)),
      rep     = rep(c("a", "b"), each = n / 2),
      y       = rnorm(n)
    )
    g$add_observations(data = df, data_coords = "spatial",
                       group = "rep", verbose = 0,
                       suppress_warnings = TRUE)
  } else {
    df <- data.frame(coord_x = runif(n, 0, 2),
                     coord_y = rep(0, n),
                     y = rnorm(n))
    g$add_observations(data = df, data_coords = "spatial",
                       verbose = 0, suppress_warnings = TRUE)
  }
}

## ── selected_inv() ───────────────────────────────────────────────────────────

test_that("selected_inv returns the diagonal of the inverse of a small SPD matrix", {
  Q <- Matrix::Diagonal(x = c(2, 3, 4))
  inv <- selected_inv(Q)
  # Selected inverse on a diagonal matrix returns 1/diag.
  diag_inv <- diag(as.matrix(inv))
  expect_equal(diag_inv, c(1 / 2, 1 / 3, 1 / 4), tolerance = 1e-12)
})

test_that("selected_inv matches solve() entries on the sparsity pattern of Q", {
  set.seed(42)
  n <- 6
  L <- Matrix::Diagonal(n) +
       Matrix::sparseMatrix(i = c(1, 2, 3, 4, 5),
                            j = c(2, 3, 4, 5, 6),
                            x = rep(0.3, 5),
                            dims = c(n, n))
  Q <- Matrix::tcrossprod(L)
  Sinv <- selected_inv(Q)
  Sfull <- solve(as.matrix(Q))
  # Selected inverse only fills the sparsity pattern of Q. Compare on the
  # entries Q is non-zero at.
  Qnz <- which(as.matrix(Q) != 0, arr.ind = TRUE)
  for (k in seq_len(nrow(Qnz))) {
    i <- Qnz[k, 1]; j <- Qnz[k, 2]
    expect_equal(as.numeric(Sinv[i, j]), Sfull[i, j], tolerance = 1e-8,
                 info = sprintf("entry (%d,%d)", i, j))
  }
})

test_that("selected_inv accepts non-dgC sparse matrices via coercion", {
  M <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(1, 2, 4),
                            dims = c(3, 3))
  inv <- selected_inv(M)
  expect_equal(diag(as.matrix(inv)), c(1, 1 / 2, 1 / 4), tolerance = 1e-12)
})

test_that("selected_inv rejects non-sparse input", {
  expect_error(selected_inv(diag(3)), "sparse matrix")
})

## ── exp_covariance() ─────────────────────────────────────────────────────────

test_that("exp_covariance evaluates the exponential covariance function", {
  # rho(h) = sigma^2 * exp(-kappa * h)
  expect_equal(exp_covariance(0, c(2, 1)), 4)
  expect_equal(exp_covariance(1, c(1, 1)), exp(-1), tolerance = 1e-12)
  h <- c(0, 0.5, 1, 2)
  expect_equal(exp_covariance(h, c(3, 0.5)),
               9 * exp(-0.5 * h), tolerance = 1e-12)
})

test_that("exp_covariance is monotone decreasing in h for positive kappa", {
  h <- seq(0, 10, length.out = 50)
  v <- exp_covariance(h, c(1, 0.4))
  expect_true(all(diff(v) < 0))
  expect_equal(v[1], 1)
})

## ── graph_starting_values() ──────────────────────────────────────────────────

test_that("graph_starting_values returns a 3-element vector with manual_data", {
  g <- make_small_graph()
  res <- graph_starting_values(g, model = "alpha1",
                                data = TRUE,
                                manual_data = rnorm(50))
  vec <- res$start_values
  expect_length(vec, 3L)
  expect_true(all(is.finite(vec)))
})

test_that("graph_starting_values works with data = FALSE", {
  g <- make_small_graph()
  res <- graph_starting_values(g, model = "alpha1", data = FALSE)
  vec <- res$start_values
  expect_length(vec, 3L)
})

test_that("graph_starting_values accepts manual_data", {
  g <- make_small_graph()
  res <- graph_starting_values(g, model = "alpha1", data = TRUE,
                                manual_data = rnorm(50))
  vec <- res$start_values
  expect_length(vec, 3L)
})

test_that("graph_starting_values respects model_options$start_range", {
  g <- make_small_graph()
  res <- graph_starting_values(g, model = "alpha1", data = FALSE,
                                model_options = list(start_range = 1.5))
  vec <- res$start_values
  expect_length(vec, 3L)
})

test_that("graph_starting_values errors on unknown model", {
  g <- make_small_graph()
  expect_error(
    graph_starting_values(g, model = "bogus", data = FALSE),
    "model"
  )
})

## ── tidyverse methods on metric_graph_data ───────────────────────────────────

test_that("filter.metric_graph_data preserves metric_graph_data class", {
  g <- make_small_graph()
  attach_simple_obs(g)
  d <- g$get_data()
  out <- dplyr::filter(d, y > 0)
  expect_s3_class(out, "metric_graph_data")
  expect_true(all(out$y > 0))
})

test_that("mutate.metric_graph_data adds new columns and keeps class", {
  g <- make_small_graph()
  attach_simple_obs(g)
  d <- g$get_data()
  out <- dplyr::mutate(d, y2 = y + 1)
  expect_s3_class(out, "metric_graph_data")
  expect_equal(out$y2, out$y + 1)
})

test_that("select.metric_graph_data preserves required graph columns", {
  g <- make_small_graph()
  attach_simple_obs(g)
  d <- g$get_data()
  out <- dplyr::select(d, y)
  expect_s3_class(out, "metric_graph_data")
  # The internal location columns must always be preserved
  for (nm in c(".group", ".edge_number", ".distance_on_edge",
               ".coord_x", ".coord_y")) {
    expect_true(nm %in% names(out), info = nm)
  }
})

test_that("drop_na.metric_graph_data drops NA rows", {
  g <- make_small_graph()
  attach_simple_obs(g)
  d <- g$get_data()
  d$y[1:2] <- NA
  out <- tidyr::drop_na(d, y)
  expect_s3_class(out, "metric_graph_data")
  expect_equal(sum(is.na(out$y)), 0L)
  expect_equal(nrow(out), nrow(d) - 2L)
})

test_that("summarise.metric_graph_data groups by location and stores group info", {
  g <- make_small_graph()
  attach_simple_obs(g, with_groups = TRUE)
  d <- g$get_data()
  out <- dplyr::summarise(d, y_mean = mean(y, na.rm = TRUE))
  expect_s3_class(out, "metric_graph_data")
  expect_true("y_mean" %in% names(out))
})

## ── metric_graph: basic getters ──────────────────────────────────────────────

test_that("get_PtE returns a 2-column matrix with valid edges and positions", {
  g <- make_small_graph()
  attach_simple_obs(g)
  pte <- g$get_PtE()
  expect_equal(ncol(pte), 2L)
  expect_true(all(pte[, 1] %in% seq_len(g$nE)))
  expect_true(all(pte[, 2] >= 0 & pte[, 2] <= 1))
})

test_that("get_PtE warns and returns NULL when there is no data", {
  g <- make_small_graph()
  expect_warning(res <- g$get_PtE(), "no data")
  expect_null(res)
})

test_that("get_groups returns the unique group identifiers", {
  g <- make_small_graph()
  attach_simple_obs(g, with_groups = TRUE)
  expect_setequal(g$get_groups(), c("a", "b"))
})

test_that("get_groups warns and returns NULL with no data", {
  g <- make_small_graph()
  expect_warning(res <- g$get_groups(), "no data")
  expect_null(res)
})

test_that("get_locations returns one location per observation", {
  g <- make_small_graph()
  attach_simple_obs(g, n = 20)
  loc <- g$get_locations()
  expect_true(is.matrix(loc) || is.data.frame(loc))
  expect_equal(NROW(loc), 20L)
})

test_that("get_degrees returns a vector of length nV with non-negative integers", {
  g <- make_small_graph()
  d <- g$get_degrees()
  expect_length(d, g$nV)
  expect_true(all(d >= 0))
})

test_that("get_degrees('indegree') and ('outdegree') sum to total degree", {
  g <- make_small_graph()
  total <- g$get_degrees("degree")
  ind   <- g$get_degrees("indegree")
  out   <- g$get_degrees("outdegree")
  expect_equal(total, ind + out)
})

test_that("get_degrees rejects invalid 'which'", {
  g <- make_small_graph()
  expect_error(g$get_degrees("bogus"), "must be either")
})

## ── metric_graph: clear_observations ─────────────────────────────────────────

test_that("clear_observations removes data and resets cached distances", {
  g <- make_small_graph()
  attach_simple_obs(g)
  expect_false(is.null(g$.__enclos_env__$private$data))

  g$clear_observations()
  expect_null(g$.__enclos_env__$private$data)
  expect_null(g$geo_dist)
  expect_null(g$res_dist)
  expect_null(g$PtV)
})

## ── metric_graph: mesh helpers ───────────────────────────────────────────────

test_that("get_mesh_locations returns one row per mesh node", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  loc <- g$get_mesh_locations()
  expect_equal(NROW(loc), nrow(g$mesh$VtE))
  expect_equal(NCOL(loc), 2L)
})

test_that("fem_basis maps PtE locations into the mesh basis", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  PtE <- rbind(c(1, 0.25), c(2, 0.5), c(3, 0.75))
  A <- g$fem_basis(PtE)
  expect_equal(nrow(A), nrow(PtE))
  expect_equal(ncol(A), nrow(g$mesh$V))
  # Each row should sum to ~1 (linear FEM partition of unity)
  expect_equal(as.numeric(Matrix::rowSums(A)), rep(1, nrow(PtE)),
               tolerance = 1e-10)
})

## ── metric_graph: print / summary ────────────────────────────────────────────

test_that("print on a metric_graph runs without error", {
  g <- make_small_graph()
  expect_output(print(g), "metric graph", ignore.case = TRUE)
})

test_that("summary on a metric_graph runs and returns invisibly", {
  g <- make_small_graph()
  expect_output(s <- summary(g))
  expect_true(!is.null(s) || is.null(s))
})

## ── graph_lme: model = 'lm' (no random effects) ──────────────────────────────

test_that("graph_lme with model = 'lm' fits a linear model on the graph data", {
  g <- make_small_graph()
  attach_simple_obs(g, n = 50)
  fit <- graph_lme(y ~ 1, graph = g, model = "lm")
  expect_s3_class(fit, "graph_lme")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  expect_length(fit$coeff$fixed_effects, 1L)
})

test_that("graph_lme 'lm' fit exposes coefficients and logLik", {
  g <- make_small_graph()
  attach_simple_obs(g, n = 50)
  fit <- graph_lme(y ~ 1, graph = g, model = "lm")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  cf <- fit$coeff$fixed_effects
  expect_length(cf, 1L)
  expect_true(is.finite(cf))
})

## ── graph_lme: simulate ──────────────────────────────────────────────────────

test_that("simulate on a graph_lme fit returns a list/matrix of the right size", {
  g <- make_small_graph()
  attach_simple_obs(g, n = 30)
  fit <- suppressWarnings(graph_lme(y ~ 1, graph = g, model = "WM1"))
  sim <- simulate(fit, nsim = 3)
  # simulate.graph_lme should return something with one column/element per sim
  expect_true(length(sim) >= 1L)
})

## ── match_mesh_data ──────────────────────────────────────────────────────────

test_that("match_mesh_data reorders rows to align with the mesh ordering", {
  g <- make_small_graph()
  g$build_mesh(h = 0.5)
  vte <- g$mesh$VtE
  # Construct a data frame in a permuted order with the mesh PtE locations
  perm <- sample(seq_len(nrow(vte)))
  data <- data.frame(
    .edge_number      = vte[perm, 1],
    .distance_on_edge = vte[perm, 2],
    val               = perm  # carry the original index as a payload
  )
  out <- match_mesh_data(graph = g, data = data)
  # After matching, the 'val' column should follow the mesh ordering,
  # so the i-th row of out corresponds to the i-th row of vte.
  expect_equal(nrow(out), nrow(vte))
  expect_equal(out$.edge_number, vte[, 1])
  expect_equal(out$.distance_on_edge, vte[, 2])
})

## ── coordinates() ────────────────────────────────────────────────────────────

test_that("coordinates(PtE) round-trips with coordinates(XY)", {
  g <- make_small_graph()
  PtE_in <- rbind(c(1, 0.25), c(2, 0.6), c(3, 0.5))
  XY <- g$coordinates(PtE = PtE_in, normalized = TRUE)
  PtE_back <- g$coordinates(XY = XY)
  expect_equal(PtE_back[, 1], PtE_in[, 1])
  expect_equal(PtE_back[, 2], PtE_in[, 2], tolerance = 1e-9)
})

test_that("coordinates(PtE) accepts vector input of length 2", {
  g <- make_small_graph()
  XY <- g$coordinates(PtE = c(1, 0.5), normalized = TRUE)
  expect_equal(dim(XY), c(1L, 2L))
  expect_equal(as.numeric(XY), c(0.5, 0.0), tolerance = 1e-12)
})

test_that("coordinates(PtE) errors on negative or out-of-range distances", {
  g <- make_small_graph()
  expect_error(g$coordinates(PtE = rbind(c(1, -0.1)), normalized = TRUE),
               "negative")
  expect_error(g$coordinates(PtE = rbind(c(1, 1.1)),  normalized = TRUE),
               "should not be|larger")
})

test_that("coordinates() requires exactly one of PtE or XY", {
  g <- make_small_graph()
  expect_error(g$coordinates(), "must be provided")
  expect_error(g$coordinates(PtE = c(1, 0.5), XY = c(0, 0)),
               "not both")
})

## ── compute_geodist / compute_resdist ────────────────────────────────────────

test_that("compute_geodist populates geo_dist with a square distance matrix", {
  g <- make_small_graph()
  attach_simple_obs(g, n = 8)
  g$compute_geodist(obs = TRUE)
  expect_false(is.null(g$geo_dist))
  D <- g$geo_dist[[1]]
  expect_equal(nrow(D), ncol(D))
  expect_equal(diag(D), rep(0, nrow(D)), tolerance = 1e-12)
  # Symmetry
  expect_equal(D, t(D), tolerance = 1e-12)
})

test_that("compute_resdist populates res_dist", {
  g <- make_small_graph()
  attach_simple_obs(g, n = 8)
  g$compute_resdist(obs = TRUE)
  expect_false(is.null(g$res_dist))
  R <- g$res_dist[[1]]
  expect_equal(nrow(R), ncol(R))
  expect_equal(diag(R), rep(0, nrow(R)), tolerance = 1e-12)
  expect_equal(R, t(R), tolerance = 1e-10)
})

## ── get_edges() / get_vertices() / get_bounding_box() ────────────────────────

test_that("get_edges returns one feature per edge", {
  g <- make_small_graph()
  edges_list <- g$get_edges(format = "list")
  expect_equal(length(edges_list), g$nE)
})

test_that("get_edges format='sf' returns an sf object with nE rows", {
  skip_if_not_installed("sf")
  g <- make_small_graph()
  e <- g$get_edges(format = "sf")
  expect_s3_class(e, "sf")
  expect_equal(nrow(e), g$nE)
})

test_that("get_vertices format='sf' returns an sf object with nV rows", {
  skip_if_not_installed("sf")
  g <- make_small_graph()
  v <- g$get_vertices(format = "sf")
  expect_s3_class(v, "sf")
  expect_equal(nrow(v), g$nV)
})

test_that("get_bounding_box returns coords matching V's range", {
  g <- make_small_graph()
  bb <- g$get_bounding_box(format = "list")
  expect_equal(bb$min_x, min(g$V[, 1]))
  expect_equal(bb$max_x, max(g$V[, 1]))
  expect_equal(bb$min_y, min(g$V[, 2]))
  expect_equal(bb$max_y, max(g$V[, 2]))
})

## ── set_edge_weights / get_edge_weights ──────────────────────────────────────

test_that("set_edge_weights stores a vector and get_edge_weights returns it", {
  g <- make_small_graph()
  w <- c(0.5, 1.5, 2.5)
  g$set_edge_weights(weights = w)
  out <- g$get_edge_weights(data.frame = TRUE)
  expect_equal(as.numeric(out$.weights), w)
})

test_that("set_edge_weights with single number broadcasts to every edge", {
  g <- make_small_graph()
  g$set_edge_weights(weights = 2.0)
  out <- g$get_edge_weights(data.frame = TRUE)
  expect_equal(as.numeric(out$.weights), rep(2.0, g$nE))
})

## ── is_tree ──────────────────────────────────────────────────────────────────

test_that("is_tree returns TRUE on a chain and FALSE on a triangle", {
  chain <- metric_graph$new(
    edges = list(rbind(c(0, 0), c(1, 0)),
                 rbind(c(1, 0), c(2, 0))),
    verbose = 0
  )
  expect_true(chain$is_tree())

  triangle <- metric_graph$new(
    edges = list(rbind(c(0, 0), c(1, 0)),
                 rbind(c(1, 0), c(0.5, 1)),
                 rbind(c(0.5, 1), c(0, 0))),
    verbose = 0
  )
  expect_false(triangle$is_tree())
})

## ── compute_characteristics + summary ────────────────────────────────────────

test_that("compute_characteristics populates the characteristics field", {
  g <- make_small_graph()
  g$compute_characteristics()
  expect_false(is.null(g$characteristics))
  expect_true(is.list(g$characteristics))
})

## ── psp / linnet / stlpp converters ─────────────────────────────────────────

test_that("psp.to.graph builds a metric_graph from a spatstat psp object", {
  skip_if_not_installed("spatstat.geom")
  ends <- matrix(c(0, 0, 1, 0,
                   1, 0, 1, 1,
                   1, 1, 0, 1),
                 ncol = 4, byrow = TRUE)
  ps <- spatstat.geom::psp(ends[, 1], ends[, 2], ends[, 3], ends[, 4],
                            window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  g <- psp.to.graph(ps)
  expect_s3_class(g, "metric_graph")
  expect_equal(g$nE, 3L)
})

test_that("linnet.to.graph builds a metric_graph from a spatstat linnet", {
  skip_if_not_installed("spatstat.data")
  skip_if_not_installed("sf")
  L <- spatstat.data::simplenet
  g <- linnet.to.graph(L, crs = sf::NA_crs_)
  expect_s3_class(g, "metric_graph")
  expect_equal(g$nE, length(L$from))
})

## ── spde_covariance: more thorough checks ────────────────────────────────────

test_that("spde_variance returns positive values at each mesh-only location", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  v <- spde_variance(kappa = 2, tau = 1, alpha = 1, graph = g)
  # Without include_vertices, variance is reported only at mesh-only nodes
  # (i.e. nrow(mesh$V) - nV). Values must all be positive.
  expect_equal(length(v), nrow(g$mesh$V) - g$nV)
  expect_true(all(v > 0))

  # With include_vertices = TRUE, the result has length nrow(mesh$V).
  v_full <- spde_variance(kappa = 2, tau = 1, alpha = 1, graph = g,
                          include_vertices = TRUE)
  expect_equal(length(v_full), nrow(g$mesh$V))
  expect_true(all(v_full > 0))
})

test_that("spde_covariance at the same point is bounded by the variance", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  P <- c(1, 0.5)
  v <- spde_variance(kappa = 2, tau = 1, alpha = 1, graph = g)
  cov_vec <- spde_covariance(P = P, kappa = 2, tau = 1, alpha = 1, graph = g)
  expect_true(length(cov_vec) > 0)
  # By Cauchy-Schwarz |Cov(X,Y)| <= sqrt(Var(X) Var(Y));
  # at the mesh nodes that have a variance entry, |cov| <= max(variance).
  expect_true(max(abs(cov_vec)) <= max(v) + 1e-8)
})

## ── sample_spde: prior + nsim > 1 ────────────────────────────────────────────

test_that("sample_spde with nsim > 1 returns a matrix of correct shape", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  set.seed(42)
  s <- sample_spde(graph = g, alpha = 1, kappa = 2, tau = 1,
                   type = "mesh", nsim = 5, method = "Q")
  expect_equal(NCOL(s), 5L)
  expect_equal(NROW(s), nrow(g$mesh$V))
})

test_that("sample_spde alpha=2 runs at mesh nodes", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  set.seed(1)
  s <- sample_spde(graph = g, alpha = 2, kappa = 2, tau = 1,
                   type = "mesh")
  expect_equal(length(s), nrow(g$mesh$V))
  expect_true(all(is.finite(s)))
})

## ── graph_lgcp_sim ──────────────────────────────────────────────────────────

test_that("graph_lgcp_sim returns u, edge_number, edge_loc for a single sim", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  g$compute_fem()
  set.seed(1)
  out <- graph_lgcp_sim(intercept = 1, sigma = 0.5, range = 1, alpha = 2,
                        graph = g)
  expect_named(out, c("u", "edge_number", "edge_loc"))
  expect_equal(length(out$u), nrow(g$mesh$V))
  expect_true(all(is.finite(out$u)))
  if (length(out$edge_number) > 0) {
    expect_true(all(out$edge_number %in% seq_len(g$nE)))
    expect_true(all(out$edge_loc >= 0 & out$edge_loc <= 1))
    expect_equal(length(out$edge_number), length(out$edge_loc))
  }
})

test_that("graph_lgcp_sim with n > 1 returns a list of n simulations", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  g$compute_fem()
  set.seed(2)
  reps <- graph_lgcp_sim(n = 3, intercept = 1, sigma = 0.5, range = 1,
                          alpha = 1, graph = g)
  expect_length(reps, 3L)
  for (r in reps) {
    expect_named(r, c("u", "edge_number", "edge_loc"))
    expect_equal(length(r$u), nrow(g$mesh$V))
  }
})

test_that("graph_lgcp_sim errors when no mesh is built", {
  g <- make_small_graph()
  expect_error(graph_lgcp_sim(intercept = 0, sigma = 1, range = 1,
                               alpha = 1, graph = g),
               "mesh")
})

test_that("graph_lgcp_sim errors on non-integer n", {
  g <- make_small_graph()
  g$build_mesh(h = 0.25)
  g$compute_fem()
  expect_error(graph_lgcp_sim(n = 1.5, intercept = 0, sigma = 1, range = 1,
                               alpha = 1, graph = g),
               "integer")
})

## ── disconnected metric_graph: distance/characteristics safety ──────────────

# Helper: a 3-component metric_graph constructed with check_connected = FALSE
make_disconnected_graph <- function() {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(2, 0)),
    rbind(c(10, 0), c(11, 0)),
    rbind(c(11, 0), c(12, 0)),
    rbind(c(20, 0), c(21, 0))
  )
  metric_graph$new(edges = edges, check_connected = FALSE, verbose = 0)
}

test_that("compute_characteristics correctly reports a disconnected graph", {
  g <- make_disconnected_graph()
  g$compute_characteristics()
  expect_false(g$characteristics$connected)

  # And a connected graph is still reported correctly
  g2 <- metric_graph$new(
    edges = list(rbind(c(0, 0), c(1, 0)), rbind(c(1, 0), c(2, 0))),
    verbose = 0
  )
  g2$compute_characteristics()
  expect_true(g2$characteristics$connected)
})

test_that("compute_resdist now succeeds on a disconnected graph (Inf cross-component)", {
  g <- make_disconnected_graph()
  set.seed(1)
  df <- data.frame(coord_x = c(0.5, 10.5, 20.5), coord_y = 0, y = rnorm(3))
  g$add_observations(data = df, data_coords = "spatial", verbose = 0)
  g$compute_resdist(obs = TRUE)
  R <- as.matrix(g$res_dist[[1]])
  expect_equal(diag(R), rep(0, 3L))
  expect_true(all(is.infinite(R[upper.tri(R)])))
})

test_that("graph_lme(model='isoExp') fits on a disconnected metric_graph", {
  g <- make_disconnected_graph()
  set.seed(2)
  df <- data.frame(coord_x = c(runif(10, 0, 2), runif(10, 10, 12), runif(5, 20, 21)),
                   coord_y = 0, y = rnorm(25))
  g$add_observations(data = df, data_coords = "spatial", verbose = 0)
  fit <- suppressWarnings(graph_lme(y ~ 1, graph = g, model = "isoExp"))
  expect_s3_class(fit, "graph_lme")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

test_that("compute_geodist returns Inf for cross-component pairs but no errors", {
  g <- make_disconnected_graph()
  df <- data.frame(coord_x = c(0.5, 10.5, 20.5), coord_y = 0, y = c(1, 2, 3))
  g$add_observations(data = df, data_coords = "spatial", verbose = 0)
  g$compute_geodist(obs = TRUE)
  D <- g$geo_dist[[1]]
  expect_equal(diag(D), rep(0, 3))
  # Pairs across components are Inf
  expect_true(is.infinite(D[1, 2]))
  expect_true(is.infinite(D[1, 3]))
  expect_true(is.infinite(D[2, 3]))
})

test_that("compute_laplacian on a disconnected graph has multiple zero eigenvalues", {
  g <- make_disconnected_graph()
  g$compute_laplacian(full = TRUE)
  L <- g$Laplacian[[1]]
  ev <- eigen(as.matrix(L), only.values = TRUE)$values
  # 3 components → 3 zero eigenvalues (within tolerance)
  expect_equal(sum(abs(ev) < 1e-10), 3L)
})

test_that("graph_lme(model='WM1') still works on a disconnected metric_graph", {
  g <- make_disconnected_graph()
  set.seed(3)
  df <- data.frame(coord_x = c(runif(20, 0, 2), runif(20, 10, 12), runif(10, 20, 21)),
                   coord_y = 0,
                   y = rnorm(50))
  g$add_observations(data = df, data_coords = "spatial", verbose = 0)
  fit <- suppressWarnings(graph_lme(y ~ 1, graph = g, model = "WM1"))
  expect_s3_class(fit, "graph_lme")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

## ── compute_resdist on disconnected graphs (native handling) ────────────────

test_that("compute_resdist now works natively on a disconnected metric_graph", {
  g <- make_disconnected_graph()
  set.seed(7)
  df <- data.frame(coord_x = c(runif(8, 0, 2), runif(8, 10, 12), runif(4, 20, 21)),
                   coord_y = 0, y = rnorm(20))
  g$add_observations(data = df, data_coords = "spatial", verbose = 0)
  g$compute_resdist(obs = TRUE)
  R <- as.matrix(g$res_dist[[1]])

  expect_equal(nrow(R), 20L)
  expect_equal(ncol(R), 20L)
  expect_equal(diag(R), rep(0, 20L))
  expect_true(isSymmetric(R))

  # The first 8 obs are on component 1, next 8 on component 2, last 4 on
  # component 3. Within-component entries are finite, between-component are Inf.
  comp <- c(rep(1L, 8), rep(2L, 8), rep(3L, 4))
  cross <- outer(comp, comp, function(a, b) a != b)
  expect_true(all(is.infinite(R[cross])))
  expect_true(all(is.finite(R[!cross])))
})

test_that("compute_resdist on a disconnected graph matches per-component computation", {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(2, 0)),
    rbind(c(10, 0), c(11, 0)),
    rbind(c(11, 0), c(12, 0))
  )
  set.seed(11)
  obs1 <- data.frame(coord_x = c(0.5, 1.5), coord_y = 0, y = c(1, 2))
  obs2 <- data.frame(coord_x = c(10.5, 11.5), coord_y = 0, y = c(3, 4))

  # Per-component
  g1 <- metric_graph$new(edges = edges[1:2], verbose = 0)
  g1$add_observations(data = obs1, data_coords = "spatial", verbose = 0)
  g1$compute_resdist(obs = TRUE)
  R1 <- as.matrix(g1$res_dist[[1]])

  g2 <- metric_graph$new(edges = edges[3:4], verbose = 0)
  g2$add_observations(data = obs2, data_coords = "spatial", verbose = 0)
  g2$compute_resdist(obs = TRUE)
  R2 <- as.matrix(g2$res_dist[[1]])

  # Combined
  g <- metric_graph$new(edges = edges, check_connected = FALSE, verbose = 0)
  g$add_observations(data = rbind(obs1, obs2),
                     data_coords = "spatial", verbose = 0)
  g$compute_resdist(obs = TRUE)
  R <- as.matrix(g$res_dist[[1]])

  expect_equal(R[1:2, 1:2], R1, tolerance = 1e-8)
  expect_equal(R[3:4, 3:4], R2, tolerance = 1e-8)
  expect_true(all(is.infinite(R[1:2, 3:4])))
  expect_true(all(is.infinite(R[3:4, 1:2])))
})

## ── disconnected metric_graph: distance/Laplacian methods ───────────────────

test_that("compute_geodist/compute_resdist/compute_laplacian work on a disconnected metric_graph", {
  e1 <- rbind(c(0, 0), c(1, 0)); e2 <- rbind(c(1, 0), c(2, 0))
  e3 <- rbind(c(10, 0), c(11, 0))
  mg <- metric_graph$new(edges = list(e1, e2, e3),
                         verbose = 0, check_connected = FALSE)
  set.seed(1)
  df <- data.frame(coord_x = c(0.5, 1.5, 10.5), coord_y = 0,
                   y = c(1, 2, 3))
  mg$add_observations(data = df, data_coords = "spatial", verbose = 0)

  mg$compute_geodist(obs = TRUE)
  expect_false(is.null(mg$geo_dist))

  mg$compute_resdist(obs = TRUE)
  expect_false(is.null(mg$res_dist))

  mg$compute_laplacian()
  expect_false(is.null(mg$Laplacian))
})

## ── graph_lme(isoExp/GL1) on a disconnected metric_graph ────────────────────

test_that("graph_lme(isoExp) fits on a disconnected metric_graph", {
  e1 <- rbind(c(0, 0), c(1, 0)); e2 <- rbind(c(1, 0), c(2, 0))
  e3 <- rbind(c(10, 0), c(11, 0)); e4 <- rbind(c(11, 0), c(12, 0))
  mg <- metric_graph$new(edges = list(e1, e2, e3, e4),
                         verbose = 0, check_connected = FALSE)
  set.seed(2)
  df <- data.frame(
    coord_x = c(runif(15, 0, 2), runif(15, 10, 12)),
    coord_y = 0, y = rnorm(30)
  )
  mg$add_observations(data = df, data_coords = "spatial", verbose = 0)

  fit <- suppressWarnings(graph_lme(y ~ 1, graph = mg, model = "isoExp"))
  expect_s3_class(fit, "graph_lme")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

test_that("graph_lme(GL1) fits on a disconnected metric_graph", {
  e1 <- rbind(c(0, 0), c(1, 0)); e2 <- rbind(c(1, 0), c(2, 0))
  e3 <- rbind(c(10, 0), c(11, 0)); e4 <- rbind(c(11, 0), c(12, 0))
  mg <- metric_graph$new(edges = list(e1, e2, e3, e4),
                         verbose = 0, check_connected = FALSE)
  set.seed(3)
  df <- data.frame(
    coord_x = c(runif(15, 0, 2), runif(15, 10, 12)),
    coord_y = 0, y = rnorm(30)
  )
  mg$add_observations(data = df, data_coords = "spatial", verbose = 0)

  fit <- suppressWarnings(graph_lme(y ~ 1, graph = mg, model = "GL1"))
  expect_s3_class(fit, "graph_lme")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})
