make_directional_test_graph <- function(add_observations = FALSE) {
  edges <- list(
    rbind(c(0, 1), c(1, 1)),
    rbind(c(1, 1), c(2, 1)),
    rbind(c(1, 3), c(1, 1)),
    rbind(c(-1, 0), c(0, 1)),
    rbind(c(-1, 2), c(0, 1)),
    rbind(c(0, 3), c(1, 3)),
    rbind(c(-2, 0), c(-1, 0))
  )
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_edge_weights(
    weights = data.frame(w = c(1, 1, 1, 2, 1, 3, 1)),
    directional_weights = "w"
  )
  if (add_observations) {
    graph$add_observations(
      data = data.frame(
        y = seq_len(10) / 7,
        edge_number = rep(1:5, each = 2L),
        distance_on_edge = rep(c(0.23, 0.71), 5L)
      ),
      normalized = TRUE,
      verbose = 0
    )
  }
  graph
}

make_nondendritic_directional_graph <- function(add_observations = FALSE) {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(2, 0)),
    rbind(c(1, 0), c(1, 1))
  )
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_edge_weights(
    weights = data.frame(w = c(1, 1, 1)),
    directional_weights = "w"
  )
  if (add_observations) {
    graph$add_observations(
      data = data.frame(
        y = c(-0.2, 0.4, 0.7, -0.1),
        edge_number = c(1, 2, 2, 3),
        distance_on_edge = c(0.2, 0.3, 0.6, 0.7)
      ),
      normalized = TRUE,
      verbose = 0
    )
  }
  graph
}

make_reversed_directional_test_graph <- function(add_observations = FALSE) {
  edges <- list(
    rbind(c(1, 1), c(0, 1)),
    rbind(c(2, 1), c(1, 1)),
    rbind(c(1, 1), c(1, 3)),
    rbind(c(0, 1), c(-1, 0)),
    rbind(c(0, 1), c(-1, 2)),
    rbind(c(1, 3), c(0, 3)),
    rbind(c(-1, 0), c(-2, 0))
  )
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_edge_weights(
    weights = data.frame(w = c(1, 1, 1, 2, 1, 3, 1)),
    directional_weights = "w"
  )
  if (add_observations) {
    graph$add_observations(
      data = data.frame(
        y = seq_len(10) / 7,
        edge_number = rep(1:5, each = 2L),
        distance_on_edge = 1 - rep(c(0.23, 0.71), 5L)
      ),
      normalized = TRUE,
      verbose = 0
    )
  }
  graph
}

make_irregular_directional_graph <- function(add_observations = FALSE) {
  edges <- list(
    rbind(c(0, 0), c(1, 1)),
    rbind(c(0, 2), c(1, 1)),
    rbind(c(1, 1), c(2, 1)),
    rbind(c(2, 1), c(3, 0)),
    rbind(c(2, 1), c(3, 2))
  )
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_edge_weights(
    weights = data.frame(w = c(1, 2, 3, 2, 1)),
    directional_weights = "w"
  )
  if (add_observations) {
    graph$add_observations(
      data = data.frame(
        y = c(-0.2, 0.4, 0.7, -0.1, 0.3),
        edge_number = seq_len(5),
        distance_on_edge = c(0.2, 0.3, 0.6, 0.7, 0.4)
      ),
      normalized = TRUE,
      verbose = 0
    )
  }
  graph
}

directional_covariance_tree_cases <- function() {
  list(
    in_tree = list(
      graph = make_directional_test_graph(),
      points = rbind(
        c(1, 0.2), c(2, 0.6), c(3, 0.4), c(4, 0.3), c(7, 0.8)
      )
    ),
    out_tree = list(
      graph = make_reversed_directional_test_graph(),
      points = rbind(
        c(2, 0.2), c(1, 0.6), c(3, 0.4), c(4, 0.3), c(7, 0.8)
      )
    )
  )
}

directional_endpoint_covariance_from_precision <- function(graph, kappa,
                                                           tau, cpp = TRUE) {
  graph$buildDirectionalConstraints(alpha = 1)
  n_constraints <- length(graph$CoB$S)
  free_basis <- graph$CoB$T[-seq_len(n_constraints), , drop = FALSE]
  edge_precision <- MetricGraph:::Qalpha1_edges(
    c(tau, kappa), graph, w = 0, BC = 1, build = TRUE, cpp = cpp
  )
  free_precision <- as.matrix(
    free_basis %*% edge_precision %*% Matrix::t(free_basis)
  )
  as.matrix(Matrix::t(free_basis)) %*%
    solve(free_precision, as.matrix(free_basis))
}

directional_endpoint_points <- function(graph) {
  points <- matrix(0, nrow = 2L * graph$nE, ncol = 2L)
  for (edge in seq_len(graph$nE)) {
    points[2L * edge - 1L, ] <- c(edge, 0)
    points[2L * edge, ] <- c(edge, graph$edge_lengths[edge])
  }
  points
}
