test_that("single-edge covariance matches the OU formula", {
  graph <- metric_graph$new(
    edges = list(rbind(c(0, 0), c(1, 0))), verbose = 0
  )
  graph$set_edge_weights(
    weights = data.frame(w = 1), directional_weights = "w"
  )
  kappa <- 0.7
  tau <- 1.3
  points <- rbind(c(1, 0.2), c(1, 0.7))
  covariance <- directional_ou_covariance(
    graph, kappa, tau, PtE = points
  )
  expected <- exp(-kappa * 0.5) / (2 * kappa * tau^2)
  expect_equal(covariance[1, 2], expected, tolerance = 1e-12)
})

test_that("public covariance C++ and R paths agree", {
  graph <- make_directional_test_graph()
  points <- rbind(
    c(1, 0.2), c(2, 0.6), c(3, 0.4), c(4, 0.3), c(7, 0.8)
  )
  points_two <- points[c(5, 2, 4), , drop = FALSE]

  covariance_cpp <- directional_ou_covariance(
    graph, 0.8, 1.2, PtE = points, PtE2 = points_two, cpp = TRUE
  )
  covariance_r <- directional_ou_covariance(
    graph, 0.8, 1.2, PtE = points, PtE2 = points_two, cpp = FALSE
  )
  expect_equal(covariance_cpp, covariance_r, tolerance = 1e-10)

  points_absolute <- points
  points_absolute[, 2] <- points[, 2] * graph$edge_lengths[points[, 1]]
  expect_equal(
    directional_ou_covariance(
      graph, 0.8, 1.2, PtE = points, normalized = TRUE
    ),
    directional_ou_covariance(
      graph, 0.8, 1.2, PtE = points_absolute, normalized = FALSE
    ),
    tolerance = 1e-12
  )

  via_method <- graph$compute_directional_covariance(
    0.8, 1.2, PtE = points, cpp = TRUE
  )
  expect_equal(
    via_method,
    directional_ou_covariance(graph, 0.8, 1.2, PtE = points),
    tolerance = 1e-12
  )
})

test_that("C++ numeric covariance setup matches the R recursion", {
  graph <- make_directional_test_graph()
  structure <- MetricGraph:::directional_ou_setup_structure(graph)
  fields <- c(
    "var_tail", "var_head", "logG_head", "logG_tail",
    "signG_head", "signG_tail"
  )
  for (parameters in list(c(0.5, 1), c(1.4, 0.7))) {
    kappa <- parameters[1]
    tau <- parameters[2]
    source_values <- seq(
      0.8, by = 0.2,
      length.out = sum(structure$V_indegree == 0)
    )
    cpp <- MetricGraph:::directional_ou_setup_numeric_dispatch(
      structure, kappa, tau, source_values, cpp = TRUE
    )
    reference <- MetricGraph:::directional_ou_setup_numeric_dispatch(
      structure, kappa, tau, source_values, cpp = FALSE
    )
    for (field in fields) {
      expect_equal(cpp[[field]], reference[[field]], tolerance = 1e-11)
    }
  }
})

test_that("covariance is symmetric positive semidefinite", {
  graph <- make_directional_test_graph()
  set.seed(1)
  points <- cbind(sample.int(graph$nE, 30, replace = TRUE), runif(30))
  covariance <- directional_ou_covariance(
    graph, 1.1, 0.9, PtE = points
  )
  expect_equal(covariance, t(covariance), tolerance = 1e-12)
  expect_gte(
    min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
    -1e-9
  )
})

test_that("K2 has stationary variance while K1 need not", {
  graph <- make_directional_test_graph()
  points <- cbind(rep(seq_len(graph$nE), each = 3L),
                  rep(c(0.1, 0.5, 0.9), graph$nE))
  kappa <- 0.6
  tau <- 1.4

  graph$setDirectionalWeightFunction(
    f_in = function(weight) sqrt(weight / sum(weight))
  )
  variance_k2 <- directional_ou_variance(
    graph, kappa, tau, PtE = points, cpp = TRUE
  )
  expect_equal(
    variance_k2, rep(1 / (2 * kappa * tau^2), nrow(points)),
    tolerance = 1e-11
  )
  expect_equal(
    variance_k2,
    directional_ou_variance(
      graph, kappa, tau, PtE = points, cpp = FALSE
    ),
    tolerance = 1e-11
  )

  graph$setDirectionalWeightFunction()
  variance_k1 <- directional_ou_variance(graph, kappa, tau, PtE = points)
  expect_gt(diff(range(variance_k1)), 1e-6)
})

test_that("closed-form endpoint covariance matches constrained precision", {
  for (weight_rule in c("K1", "K2")) {
    graph <- make_directional_test_graph()
    if (identical(weight_rule, "K2")) {
      graph$setDirectionalWeightFunction(
        f_in = function(weight) sqrt(weight / sum(weight))
      )
    }
    kappa <- 0.9
    tau <- 1.1
    precision_covariance <- directional_endpoint_covariance_from_precision(
      graph, kappa, tau
    )
    closed_form <- directional_ou_covariance(
      graph, kappa, tau,
      PtE = directional_endpoint_points(graph), normalized = FALSE
    )
    expect_equal(
      closed_form, precision_covariance, tolerance = 2e-9,
      label = weight_rule
    )
  }
})

test_that("source variances are validated and applied", {
  graph <- make_directional_test_graph()
  structure <- MetricGraph:::directional_ou_setup_structure(graph)
  sources <- which(structure$V_indegree == 0)
  expect_error(
    MetricGraph:::directional_ou_resolve_source_var(
      structure, sigma_source = 2, sigma_stationary = 1
    ),
    "length"
  )
  expect_error(
    MetricGraph:::directional_ou_resolve_source_var(
      structure, sigma_source = c(not_a_vertex = 2),
      sigma_stationary = 1
    ),
    "not_a_vertex"
  )

  source_variances <- stats::setNames(
    seq(1.2, by = 0.3, length.out = length(sources)),
    as.character(sources)
  )
  source_edges <- match(sources, graph$E[, 1])
  variances <- directional_ou_variance(
    graph, 0.8, 1, PtE = cbind(source_edges, 0),
    sigma_source = source_variances
  )
  expect_equal(variances, unname(source_variances), tolerance = 1e-12)
})

test_that("non-dendritic covariance falls back to R", {
  graph <- make_nondendritic_directional_graph()
  setup <- MetricGraph:::directional_ou_setup(graph, 0.7, 1.1, NULL)
  expect_false(setup$is_dendritic)
  points <- rbind(c(1, 0.3), c(2, 0.4), c(2, 0.6), c(3, 0.7))
  expect_identical(
    MetricGraph:::directional_ou_covariance_from_setup_cpp(
      setup, points
    ),
    MetricGraph:::directional_ou_covariance_from_setup(setup, points)
  )
})

test_that("a branching-source LCA gives a clear unsupported-case error", {
  graph <- metric_graph$new(
    edges = list(
      rbind(c(0, 0), c(1, 0)),
      rbind(c(0, 0), c(-1, 0))
    ),
    verbose = 0
  )
  graph$set_edge_weights(
    weights = data.frame(w = c(1, 1)), directional_weights = "w"
  )
  expect_error(
    directional_ou_covariance(
      graph, 0.8, 1, PtE = rbind(c(1, 0.4), c(2, 0.6))
    ),
    "branching source vertex"
  )
})

test_that("invalid topology, deep trees, and dense-size guard are handled", {
  cycle_edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(0.5, 1)),
    rbind(c(0.5, 1), c(0, 0))
  )
  cycle_graph <- metric_graph$new(
    edges = cycle_edges, verbose = 0, check_connected = FALSE
  )
  cycle_graph$set_edge_weights(
    weights = data.frame(w = rep(1, cycle_graph$nE)),
    directional_weights = "w"
  )
  expect_error(
    directional_ou_covariance(cycle_graph, 1, 1, PtE = cbind(1, 0.5)),
    "cycle|tree"
  )

  n_edges <- 5000L
  E <- cbind(seq_len(n_edges), seq_len(n_edges) + 1L)
  outdegree <- c(rep(1L, n_edges), 0L)
  in_edges <- split(seq_len(n_edges), E[, 2])
  labels <- MetricGraph:::directional_ou_ancestor_labels(
    E, outdegree, in_edges, is_dendritic = TRUE
  )
  expect_length(labels$enter, n_edges)
  expect_true(all(labels$enter < labels$exit))

  graph <- make_directional_test_graph()
  setup <- MetricGraph:::directional_ou_setup(graph, 1, 1, NULL)
  too_many_points <- cbind(rep(1L, 10001L), rep(0.5, 10001L))
  expect_error(
    MetricGraph:::directional_ou_covariance_from_setup_cpp(
      setup, too_many_points
    ),
    "must be <= 10000"
  )
})
