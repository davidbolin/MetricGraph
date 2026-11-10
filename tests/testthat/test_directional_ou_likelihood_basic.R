# Small-graph equivalence checks for the closed-form directional OU covariance
# and its dense likelihood against the independent sparse precision route.

theta_test <- c(log(0.3), log(0.8), log(1.2))

test_that("outward flux fixture is a genuine out-tree", {
  graph <- subgraph_equivalence_cases()[["outward_flux"]]
  structure <- MetricGraph:::directional_ou_setup_structure(graph)

  expect_identical(structure$tree_orientation, "out")
  expect_true(any(structure$V_outdegree > 1L))
  expect_true(all(structure$V_indegree <= 1L))
})

test_that("closed-form endpoint covariance equals constrained precision inverse on all small graphs", {
  cases <- subgraph_equivalence_cases()
  parameter_grid <- list(
    c(kappa = 0.5, tau = 0.8),
    c(kappa = 2, tau = 1.5)
  )

  for (case_name in names(cases)) {
    graph <- cases[[case_name]]
    expect_lte(
      graph$nE, 10L,
      label = paste(case_name, "is within the dense-inverse size limit")
    )
    for (parameters in parameter_grid) {
      kappa <- parameters[["kappa"]]
      tau <- parameters[["tau"]]
      precision_inverse <- directional_endpoint_covariance_from_precision(
        graph, kappa, tau
      )
      closed_form <- directional_ou_covariance(
        graph, kappa, tau,
        PtE = directional_endpoint_points(graph), normalized = FALSE
      )
      expect_equal(
        closed_form, precision_inverse, tolerance = 2e-9,
        label = paste(case_name, "kappa", kappa, "tau", tau)
      )
    }
  }
})

test_that("dense directional OU likelihood matches sparse likelihood without covariates", {
  cases <- subgraph_equivalence_cases()
  case_names <- c(
    "no_pruning", "two_source_merge", "two_source_k2",
    "nested_three_hop"
  )
  for (case_name in case_names) {
    graph <- cases[[case_name]]
    response <- graph$.__enclos_env__$private$data[["temp"]]
    covariance_value <- MetricGraph:::directional_ou_covariance_loglik(
      theta_test, graph, manual_y = response, X_cov = NULL,
      sigma_source = NULL
    )
    precision_precomputed <- MetricGraph:::precompute_alpha1_directional(
      graph, data_name = "temp"
    )
    precision_value <-
      MetricGraph:::likelihood_alpha1_directional_profile_precompute(
        theta_test, precomputed_data = precision_precomputed,
        parameterization = "spde", cpp = TRUE
      )
    expect_equal(
      covariance_value, precision_value, tolerance = 1e-6,
      label = paste("dense versus sparse log likelihood", case_name)
    )
  }
})

test_that("dense directional OU likelihood matches sparse likelihood with covariates", {
  graph <- subgraph_equivalence_cases()[["no_pruning"]]
  response <- graph$.__enclos_env__$private$data[["temp"]]
  set.seed(1)
  design <- cbind(1, rnorm(length(response)))

  covariance_value <- MetricGraph:::directional_ou_covariance_loglik(
    theta_test, graph, manual_y = response, X_cov = design,
    sigma_source = NULL
  )
  precision_precomputed <- MetricGraph:::precompute_alpha1_directional(
    graph, data_name = "temp", X_cov = design
  )
  precision_value <-
    MetricGraph:::likelihood_alpha1_directional_profile_precompute(
      theta_test, precomputed_data = precision_precomputed,
      parameterization = "spde", cpp = TRUE
    )
  expect_equal(
    covariance_value, precision_value, tolerance = 1e-6,
    label = "dense versus sparse log likelihood with two covariates"
  )
})

test_that("directional OU covariance likelihood rejects wrong-length theta", {
  graph <- subgraph_equivalence_cases()[["no_pruning"]]
  response <- graph$.__enclos_env__$private$data[["temp"]]
  expect_error(
    MetricGraph:::directional_ou_covariance_loglik(
      c(log(0.3), log(0.8)), graph, manual_y = response,
      X_cov = NULL, sigma_source = NULL
    ),
    "theta must have length 3"
  )
})

test_that("REML likelihood adjustment matches the determinant formula", {
  graph <- subgraph_equivalence_cases()[["no_pruning"]]
  response <- graph$.__enclos_env__$private$data[["temp"]]
  set.seed(2)
  design <- cbind(1, rnorm(length(response)))

  ml_value <- MetricGraph:::directional_ou_covariance_loglik(
    theta_test, graph, manual_y = response, X_cov = design,
    sigma_source = NULL, reml = FALSE
  )
  reml_value <- MetricGraph:::directional_ou_covariance_loglik(
    theta_test, graph, manual_y = response, X_cov = design,
    sigma_source = NULL, reml = TRUE
  )
  core <- MetricGraph:::directional_ou_covariance_loglik_core(
    theta_test, graph, manual_y = response, X_cov = design,
    sigma_source = NULL
  )
  expected_adjustment <- -0.5 * as.numeric(
    determinant(core$H, logarithm = TRUE)$modulus
  )
  expect_equal(
    reml_value - ml_value, expected_adjustment, tolerance = 1e-8
  )
})
