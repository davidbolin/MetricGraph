test_that("covariance likelihood C++ and R paths agree", {
  graph <- make_directional_test_graph(add_observations = TRUE)
  response <- graph$.__enclos_env__$private$data[["y"]]
  design <- cbind(1, seq_along(response) / length(response))
  theta_values <- list(
    c(log(0.15), log(0.8), log(1.1)),
    c(log(0.3), log(1.3), log(0.6))
  )

  precomputed <- MetricGraph:::precompute_directional_ou_covariance(
    graph, data_name = "y", X_cov = design
  )
  for (theta in theta_values) {
    for (reml in c(FALSE, TRUE)) {
      value_cpp <- MetricGraph:::directional_ou_covariance_loglik_precompute(
        theta, precomputed, reml = reml, cpp = TRUE
      )
      value_r <- MetricGraph:::directional_ou_covariance_loglik_precompute(
        theta, precomputed, reml = reml, cpp = FALSE
      )
      expect_equal(value_cpp, value_r, tolerance = 1e-9)
    }

    direct <- MetricGraph:::directional_ou_covariance_loglik(
      theta, graph, manual_y = response, X_cov = design, cpp = TRUE
    )
    cached <- MetricGraph:::directional_ou_covariance_loglik_precompute(
      theta, precomputed, cpp = TRUE
    )
    expect_equal(direct, cached, tolerance = 1e-10)
  }
})

test_that("covariance likelihood agrees with sparse precision likelihood", {
  graph <- make_directional_test_graph(add_observations = TRUE)
  graph$setDirectionalWeightFunction(
    f_in = function(weight) sqrt(weight / sum(weight))
  )
  response <- graph$.__enclos_env__$private$data[["y"]]
  design <- cbind(1, sin(seq_along(response)))
  theta <- c(log(0.2), log(0.9), log(1.2))

  covariance_precomputed <- MetricGraph:::precompute_directional_ou_covariance(
    graph, data_name = "y", X_cov = design
  )
  precision_precomputed <- MetricGraph:::precompute_alpha1_directional(
    graph, data_name = "y", X_cov = design
  )
  covariance_value <- MetricGraph:::directional_ou_covariance_loglik_precompute(
    theta, covariance_precomputed, cpp = TRUE
  )
  precision_value <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
    theta, precision_precomputed, parameterization = "spde", cpp = TRUE
  )
  expect_equal(covariance_value, precision_value, tolerance = 1e-6)
})

test_that("sparsity-aware and forced-dense GLS calculations agree", {
  block <- matrix(c(1.5, 0.2, 0.2, 1.2), 2, 2)
  covariance <- as.matrix(Matrix::bdiag(replicate(10, block, simplify = FALSE)))
  response <- seq_len(nrow(covariance)) / 10
  design <- cbind(1, cos(response))

  routed <- MetricGraph:::directional_ou_gls_core_from_sigma(
    covariance, response, design, n_cov = 2L, force_dense = FALSE
  )
  dense <- MetricGraph:::directional_ou_gls_core_from_sigma(
    covariance, response, design, n_cov = 2L, force_dense = TRUE
  )
  expect_equal(routed, dense, tolerance = 1e-10)
})

test_that("covariance likelihood C++ and R paths agree on out-trees", {
  graph <- make_reversed_directional_test_graph(add_observations = TRUE)
  precomputed <- MetricGraph:::precompute_directional_ou_covariance(
    graph, data_name = "y"
  )
  expect_identical(precomputed$structure$tree_orientation, "out")
  theta <- c(log(0.2), log(1.1), log(0.8))
  expect_equal(
    MetricGraph:::directional_ou_covariance_loglik_precompute(
      theta, precomputed, cpp = TRUE
    ),
    MetricGraph:::directional_ou_covariance_loglik_precompute(
      theta, precomputed, cpp = FALSE
    ),
    tolerance = 1e-9
  )
})

test_that("small reversed Columbia covariance smoke test is finite", {
  skip_on_cran()
  helper <- system.file(
    "examples/directional/columbia_full_graph_helpers.R",
    package = "MetricGraph"
  )
  skip_if(!nzchar(helper), "Installed Columbia helper is unavailable")

  helper_environment <- new.env(parent = globalenv())
  sys.source(helper, envir = helper_environment)
  utils::data(
    "columbia_main_component", package = "MetricGraph",
    envir = environment()
  )
  graph <- helper_environment$columbia_make_graph(
    columbia_main_component, reversed = TRUE
  )
  graph$setDirectionalWeightFunction()
  structure <- MetricGraph:::directional_ou_setup_structure(graph)
  expect_identical(structure$tree_orientation, "out")
  expect_gt(
    sum(structure$V_outdegree > 1L),
    0L
  )
  set.seed(10)
  points <- cbind(sample.int(graph$nE, 12L), runif(12L, 0.1, 0.9))
  covariance <- directional_ou_covariance(
    graph, kappa = 0.7, tau = 1, PtE = points, cpp = TRUE
  )
  expect_true(all(is.finite(covariance)))
  expect_equal(covariance, t(covariance), tolerance = 1e-10)
})
