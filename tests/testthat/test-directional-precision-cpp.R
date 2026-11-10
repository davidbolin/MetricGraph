test_that("C++ edge precision matches the R reference", {
  graph <- make_directional_test_graph()
  source_vertices <- which(graph$get_degrees("indegree") == 0)
  stationary_cases <- list("all", "none", source_vertices[1:2])

  for (theta in list(c(1, 0.4), c(0.7, 1.3))) {
    for (weight in c(0, 0.25, 1)) {
      for (stationary_points in stationary_cases) {
        cpp_triplets <- MetricGraph:::Qalpha1_edges(
          theta, graph, weight, stationary_points = stationary_points,
          build = FALSE, cpp = TRUE
        )
        r_triplets <- MetricGraph:::Qalpha1_edges(
          theta, graph, weight, stationary_points = stationary_points,
          build = FALSE, cpp = FALSE
        )
        expect_equal(cpp_triplets, r_triplets, tolerance = 1e-13)

        cpp_matrix <- MetricGraph:::Qalpha1_edges(
          theta, graph, weight, stationary_points = stationary_points,
          cpp = TRUE
        )
        r_matrix <- MetricGraph:::Qalpha1_edges(
          theta, graph, weight, stationary_points = stationary_points,
          cpp = FALSE
        )
        expect_equal(cpp_matrix, r_matrix, tolerance = 1e-13)
      }
    }
  }
})

test_that("edge precision validates stationary sources", {
  graph <- make_directional_test_graph()
  non_source <- which(graph$get_degrees("indegree") > 0)[1]
  expect_error(
    MetricGraph:::Qalpha1_edges(
      c(1, 1), graph, 0, stationary_points = non_source
    ),
    "inward degree zero"
  )
  expect_error(
    MetricGraph:::Qalpha1_edges(
      c(1, 1), graph, 0, stationary_points = "invalid"
    ),
    "either 'all' or 'none'"
  )
})

test_that("C++ edge precision supports graph_components", {
  graph_one <- make_directional_test_graph()
  graph_two <- make_directional_test_graph()
  components <- structure(
    list(graphs = list(graph_one, graph_two)),
    class = "graph_components"
  )

  cpp <- MetricGraph:::Qalpha1_edges(c(0.8, 1.1), components, w = 0,
                                     cpp = TRUE)
  reference <- MetricGraph:::Qalpha1_edges(c(0.8, 1.1), components, w = 0,
                                           cpp = FALSE)
  expect_equal(cpp, reference, tolerance = 1e-13)
  expect_error(
    MetricGraph:::Qalpha1_edges(
      c(0.8, 1.1), components, w = 0, build = FALSE
    ),
    "not supported"
  )
})

test_that("active directional constraint builder matches its R reference", {
  graph <- make_directional_test_graph()
  weights <- graph$get_edge_weights()[["w"]]
  indegree <- graph$get_degrees("indegree")
  outdegree <- graph$get_degrees("outdegree")

  reference <- MetricGraph:::construct_directional_constraint_matrix(
    E = graph$E, nV = graph$nV, nE = graph$nE, alpha = 1L,
    V_indegree = indegree, V_outdegree = outdegree, weight = weights,
    DirectionalWeightFunction_out = graph$DirectionalWeightFunction_out,
    DirectionalWeightFunction_in = graph$DirectionalWeightFunction_in
  )
  graph$buildDirectionalConstraints(alpha = 1)
  expect_equal(graph$C, reference, tolerance = 1e-13)
})

test_that("directional profile likelihood defaults to C++ and retains R path", {
  graph <- make_directional_test_graph(add_observations = TRUE)
  design <- cbind(1, seq_len(10) / 10)
  precomputed <- MetricGraph:::precompute_alpha1_directional(
    graph, data_name = "y", X_cov = design
  )
  theta <- c(log(0.2), log(0.9), log(1.2))

  value_default <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
    theta, precomputed, parameterization = "spde"
  )
  value_cpp <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
    theta, precomputed, parameterization = "spde", cpp = TRUE
  )
  value_r <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
    theta, precomputed, parameterization = "spde", cpp = FALSE
  )
  expect_identical(value_default, value_cpp)
  expect_equal(value_cpp, value_r, tolerance = 1e-10)

  fit <- optim(
    theta,
    function(parameters) {
      -MetricGraph:::likelihood_alpha1_directional_profile_precompute(
        parameters, precomputed, parameterization = "spde", cpp = TRUE
      )
    },
    method = "Nelder-Mead",
    control = list(maxit = 6)
  )
  expect_true(is.finite(fit$value))
  expect_true(all(is.finite(fit$par)))
  expect_equal(
    MetricGraph:::likelihood_alpha1_directional_profile_precompute(
      fit$par, precomputed, parameterization = "spde", cpp = TRUE
    ),
    MetricGraph:::likelihood_alpha1_directional_profile_precompute(
      fit$par, precomputed, parameterization = "spde", cpp = FALSE
    ),
    tolerance = 1e-9
  )
})
