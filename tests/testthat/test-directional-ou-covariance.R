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

test_that("public covariance C++ and R paths agree on oriented trees", {
  graphs_and_points <- directional_covariance_tree_cases()

  for (case_name in names(graphs_and_points)) {
    case <- graphs_and_points[[case_name]]
    points_two <- case$points[c(5, 2, 4), , drop = FALSE]
    covariance_cpp <- directional_ou_covariance(
      case$graph, 0.8, 1.2, PtE = case$points,
      PtE2 = points_two, cpp = TRUE
    )
    covariance_r <- directional_ou_covariance(
      case$graph, 0.8, 1.2, PtE = case$points,
      PtE2 = points_two, cpp = FALSE
    )
    expect_equal(
      covariance_cpp, covariance_r, tolerance = 1e-10,
      label = case_name
    )
  }
})

test_that("normalized and absolute point coordinates agree", {
  graphs_and_points <- directional_covariance_tree_cases()

  for (case_name in names(graphs_and_points)) {
    case <- graphs_and_points[[case_name]]
    absolute_points <- case$points
    absolute_points[, 2] <- case$points[, 2] *
      case$graph$edge_lengths[case$points[, 1]]
    expect_equal(
      directional_ou_covariance(
        case$graph, 0.8, 1.2, PtE = case$points, normalized = TRUE
      ),
      directional_ou_covariance(
        case$graph, 0.8, 1.2, PtE = absolute_points,
        normalized = FALSE
      ),
      tolerance = 1e-12,
      label = case_name
    )
  }
})

test_that("metric graph covariance method matches the public function", {
  graph <- make_directional_test_graph()
  points <- rbind(
    c(1, 0.2), c(2, 0.6), c(3, 0.4), c(4, 0.3), c(7, 0.8)
  )

  expect_equal(
    graph$compute_directional_covariance(
      0.8, 1.2, PtE = points, cpp = TRUE
    ),
    directional_ou_covariance(graph, 0.8, 1.2, PtE = points),
    tolerance = 1e-12
  )
})

test_that("C++ numeric covariance setup matches the R recursion", {
  fields <- c(
    "var_tail", "var_head", "logG_head", "logG_tail",
    "signG_head", "signG_tail"
  )
  graphs <- list(
    in_tree = make_directional_test_graph(),
    out_tree = make_reversed_directional_test_graph()
  )
  for (graph_name in names(graphs)) {
    structure <- MetricGraph:::directional_ou_setup_structure(
      graphs[[graph_name]]
    )
    expected_orientation <- if (graph_name == "in_tree") "in" else "out"
    expect_identical(structure$tree_orientation, expected_orientation)
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
        expect_equal(
          cpp[[field]], reference[[field]], tolerance = 1e-11,
          label = paste(graph_name, field)
        )
      }
    }
  }
})

test_that("covariance is symmetric positive semidefinite", {
  graph <- make_directional_test_graph()
  points <- cbind(
    rep(seq_len(graph$nE), each = 3L),
    rep(c(0.15, 0.5, 0.85), graph$nE)
  )
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

test_that("reversed K1 selects the out-tree fast path", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  kappa <- 0.6
  tau <- 1.4
  setup <- MetricGraph:::directional_ou_setup(
    graph, kappa, tau, NULL, cpp = TRUE
  )

  expect_identical(setup$tree_orientation, "out")
  expect_false(setup$is_dendritic)
  expect_length(setup$depth, graph$nE)
  expect_length(setup$parent_edge, graph$nE)
  expect_named(
    setup$out_tree_lca_index,
    c("first", "root_edge", "log2_floor", "rmq")
  )
})

test_that("out-tree parent edges match their incoming vertex edges", {
  E <- rbind(
    c(1L, 2L),
    c(2L, 3L),
    c(2L, 4L),
    c(4L, 5L),
    c(6L, 7L)
  )
  n_edges <- nrow(E)
  indegree <- tabulate(E[, 2], nbins = 7L)
  outdegree <- tabulate(E[, 1], nbins = 7L)
  labels <- MetricGraph:::directional_ou_ancestor_labels(
    E,
    indegree,
    outdegree,
    split(seq_len(n_edges), E[, 2]),
    split(seq_len(n_edges), E[, 1]),
    "out"
  )

  expect_identical(
    labels$parent_edge,
    c(NA_integer_, 1L, 1L, 3L, NA_integer_)
  )
})

test_that("out-tree C++ covariance bypasses the R pairwise fallback", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  setup <- MetricGraph:::directional_ou_setup(
    graph, 0.6, 1.4, NULL, cpp = TRUE
  )
  points <- rbind(
    c(2, 0.2), c(1, 0.6), c(3, 0.4), c(4, 0.3), c(7, 0.8)
  )
  expected <- MetricGraph:::directional_ou_covariance_from_setup(
    setup, points
  )

  testthat::local_mocked_bindings(
    directional_ou_covariance_from_setup = function(...) {
      stop("unexpected R covariance fallback")
    },
    .package = "MetricGraph"
  )
  expect_equal(
    MetricGraph:::directional_ou_covariance_from_setup_cpp(setup, points),
    expected,
    tolerance = 1e-11
  )
})

test_that("reversed K1 retains stationary variance", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  kappa <- 0.6
  tau <- 1.4

  points <- cbind(
    rep(seq_len(graph$nE), each = 3L),
    rep(c(0.1, 0.5, 0.9), graph$nE)
  )
  expect_equal(
    directional_ou_variance(graph, kappa, tau, PtE = points),
    rep(1 / (2 * kappa * tau^2), nrow(points)),
    tolerance = 1e-11
  )
})

test_that("out-tree two-pointer LCA agrees with the generic oracle", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  setup <- MetricGraph:::directional_ou_setup(
    graph, 0.6, 1.4, NULL, cpp = TRUE
  )

  for (edge_s in seq_len(graph$nE - 1L)) {
    for (edge_t in seq.int(edge_s + 1L, graph$nE)) {
      s_contains_t <- MetricGraph:::directional_edge_downstream_of(
        setup, edge_s, edge_t
      )
      t_contains_s <- MetricGraph:::directional_edge_downstream_of(
        setup, edge_t, edge_s
      )
      if (!s_contains_t && !t_contains_s) {
        expect_identical(
          MetricGraph:::directional_lca_edge_out_tree(
            setup, edge_s, edge_t
          ),
          MetricGraph:::directional_lca_edge_nondendritic(
            setup, edge_s, edge_t
          ),
          label = paste("edges", edge_s, "and", edge_t)
        )
      }
    }
  }
})

test_that("out-tree sibling covariance uses a third-edge LCA", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  kappa <- 0.6
  tau <- 1.4
  setup <- MetricGraph:::directional_ou_setup(
    graph, kappa, tau, NULL, cpp = TRUE
  )

  lca <- MetricGraph:::directional_lca_edge_out_tree(setup, 4L, 5L)
  expect_identical(lca$edge, 1L)

  covariance <- directional_ou_covariance(
    graph, kappa, tau,
    PtE = rbind(c(4, 0.3), c(5, 0.3)), normalized = FALSE
  )
  expected <- exp(-2 * kappa * 0.3) / (2 * kappa * tau^2)
  expect_equal(covariance[1, 2], expected, tolerance = 1e-11)
})

test_that("out-tree fast transfer agrees with the explicit walk", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  setup <- MetricGraph:::directional_ou_setup(
    graph, 0.6, 1.4, NULL, cpp = TRUE
  )

  for (from_edge in seq_len(graph$nE)) {
    for (to_edge in seq_len(graph$nE)) {
      if (from_edge != to_edge &&
          MetricGraph:::directional_edge_downstream_of(
            setup, from_edge, to_edge
          )) {
        from <- c(from_edge, 0.2 * graph$edge_lengths[from_edge])
        to <- c(to_edge, 0.7 * graph$edge_lengths[to_edge])
        fast <- MetricGraph:::directional_log_transfer(setup, from, to)
        walk <- MetricGraph:::directional_transfer_walk(setup, from, to)
        expect_equal(
          fast, walk, tolerance = 1e-11,
          label = paste("edges", from_edge, "to", to_edge)
        )
      }
    }
  }
})

test_that("out-tree transfer is multiplicative along chains", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  setup <- MetricGraph:::directional_ou_setup(
    graph, 0.6, 1.4, NULL, cpp = TRUE
  )

  from <- c(2, 0.2)
  middle <- c(1, 0.5)
  to <- c(4, 0.3)
  transfer <- function(first, second) {
    value <- MetricGraph:::directional_log_transfer(setup, first, second)
    value$sign * exp(value$log_abs)
  }
  expect_equal(
    transfer(from, to),
    transfer(from, middle) * transfer(middle, to),
    tolerance = 1e-11
  )
})

test_that("out-tree Euler reachability agrees with forward search", {
  graph <- make_reversed_directional_test_graph()
  graph$setDirectionalWeightFunction()
  setup <- MetricGraph:::directional_ou_setup(
    graph, 0.6, 1.4, NULL, cpp = TRUE
  )

  for (from_edge in seq_len(graph$nE)) {
    for (to_edge in setdiff(seq_len(graph$nE), from_edge)) {
      expect_identical(
        MetricGraph:::directional_edge_downstream_of(
          setup, from_edge, to_edge
        ),
        MetricGraph:::directional_edge_reachable_forward(
          graph$E, setup$out_edges_by_vertex, from_edge, to_edge
        ),
        label = paste("edges", from_edge, "to", to_edge)
      )
    }
  }
})

test_that("closed-form endpoint covariance matches constrained precision", {
  graph_makers <- list(
    "in-tree" = make_directional_test_graph,
    reversed = make_reversed_directional_test_graph
  )
  parameter_pairs <- list(c(0.5, 0.8), c(2, 1.5))

  for (graph_name in names(graph_makers)) {
    for (weight_rule in c("K1", "K2")) {
      graph <- graph_makers[[graph_name]]()
      if (weight_rule == "K2") {
        graph$setDirectionalWeightFunction(
          f_in = function(weight) sqrt(weight / sum(weight))
        )
      }
      for (parameters in parameter_pairs) {
        kappa <- parameters[1]
        tau <- parameters[2]
        precision_covariance <- directional_endpoint_covariance_from_precision(
          graph, kappa, tau
        )
        closed_form <- directional_ou_covariance(
          graph, kappa, tau,
          PtE = directional_endpoint_points(graph), normalized = FALSE
        )
        expect_equal(
          closed_form, precision_covariance, tolerance = 2e-9,
          label = paste(graph_name, weight_rule, kappa, tau)
        )
      }
    }
  }
})

test_that("source variance input is validated", {
  graph <- make_directional_test_graph()
  structure <- MetricGraph:::directional_ou_setup_structure(graph)
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
})

test_that("named source variances are applied at source endpoints", {
  graph <- make_directional_test_graph()
  structure <- MetricGraph:::directional_ou_setup_structure(graph)
  sources <- which(structure$V_indegree == 0)

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

test_that("irregular covariance retains the generic R fallback", {
  graph <- make_irregular_directional_graph()
  setup <- MetricGraph:::directional_ou_setup(graph, 0.7, 1.1, NULL)
  expect_false(setup$is_dendritic)
  expect_identical(setup$tree_orientation, "irregular")
  expect_null(setup$enter)
  expect_null(setup$parent_edge)
  points <- rbind(c(1, 0.3), c(2, 0.4), c(3, 0.6), c(5, 0.7))
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

test_that("disconnected out-tree components have zero cross-covariance", {
  graph <- metric_graph$new(
    edges = list(
      rbind(c(0, 0), c(1, 0)),
      rbind(c(1, 0), c(2, -1)),
      rbind(c(1, 0), c(2, 1)),
      rbind(c(10, 0), c(11, 0)),
      rbind(c(11, 0), c(12, -1)),
      rbind(c(11, 0), c(12, 1))
    ),
    verbose = 0,
    check_connected = FALSE
  )
  graph$set_edge_weights(
    weights = data.frame(w = rep(1, 6)), directional_weights = "w"
  )
  setup <- MetricGraph:::directional_ou_setup(graph, 0.8, 1, NULL)
  expect_identical(setup$tree_orientation, "out")
  expect_null(MetricGraph:::directional_lca_edge_out_tree(setup, 2L, 5L))
  covariance <- directional_ou_covariance(
    graph, 0.8, 1, PtE = rbind(c(2, 0.4), c(5, 0.6))
  )
  expect_equal(covariance[1, 2], 0)
})

test_that("directional covariance rejects cyclic topology", {
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
})

test_that("ancestor labeling handles deep in-trees and out-trees iteratively", {
  n_edges <- 5000L
  E <- cbind(seq_len(n_edges), seq_len(n_edges) + 1L)
  indegree <- c(0L, rep(1L, n_edges))
  outdegree <- c(rep(1L, n_edges), 0L)
  in_edges <- split(seq_len(n_edges), E[, 2])
  out_edges <- split(seq_len(n_edges), E[, 1])
  in_labels <- MetricGraph:::directional_ou_ancestor_labels(
    E, indegree, outdegree, in_edges, out_edges, "in"
  )
  expect_length(in_labels$enter, n_edges)
  expect_true(all(in_labels$enter < in_labels$exit))

  out_labels <- MetricGraph:::directional_ou_ancestor_labels(
    E, indegree, outdegree, in_edges, out_edges, "out"
  )
  expect_length(out_labels$enter, n_edges)
  expect_identical(out_labels$depth, seq.int(0L, n_edges - 1L))
  expect_identical(out_labels$parent_edge[-1L], seq_len(n_edges - 1L))
})

test_that("C++ dense covariance enforces its point-count guard", {
  old_limit <- getOption("DIRECTIONAL_OU_MAX_POINTS")
  on.exit(options(DIRECTIONAL_OU_MAX_POINTS = old_limit), add = TRUE)
  options(DIRECTIONAL_OU_MAX_POINTS = 10000L)

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

test_that("dense covariance point-count guard honors its R option", {
  old_limit <- getOption("DIRECTIONAL_OU_MAX_POINTS")
  on.exit(options(DIRECTIONAL_OU_MAX_POINTS = old_limit), add = TRUE)

  graph <- make_directional_test_graph()
  setup <- MetricGraph:::directional_ou_setup(graph, 1, 1, NULL)
  points <- cbind(rep(1L, 3L), c(0.2, 0.5, 0.8))

  options(DIRECTIONAL_OU_MAX_POINTS = 2L)
  expect_error(
    MetricGraph:::directional_ou_covariance_from_setup_cpp(setup, points),
    "must be <= 2"
  )

  options(DIRECTIONAL_OU_MAX_POINTS = 3L)
  expect_no_error({
    covariance <- MetricGraph:::directional_ou_covariance_from_setup_cpp(
      setup, points
    )
  })
  expect_identical(dim(covariance), c(3L, 3L))
})
