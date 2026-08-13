# Shared small directed-tree fixtures for covariance/precision equivalence
# tests. testthat loads helper-*.R files automatically before test files.

make_equivalence_graph <- function(edges, observed_edges, weights = NULL) {
  graph <- metric_graph$new(edges = edges, verbose = 0)
  if (is.null(weights)) {
    weights <- 1 + (seq_len(graph$nE) %% 3)
  }
  graph$set_edge_weights(
    weights = data.frame(w = weights),
    directional_weights = "w"
  )
  graph$setDirectionalWeightFunction()
  n_observations <- 2L * length(observed_edges)
  graph$add_observations(
    data = data.frame(
      temp = seq_len(n_observations) / 7,
      obs_id = sprintf("id-%02d", seq_len(n_observations)),
      edge_number = rep(observed_edges, each = 2L),
      distance_on_edge = rep(c(0.23, 0.71), length(observed_edges))
    ),
    normalized = TRUE,
    verbose = 0
  )
  graph
}

make_equivalence_tree_graph <- function(extra_upstream_source = FALSE) {
  edges <- list(
    rbind(c(0, 1), c(1, 1)),
    rbind(c(1, 1), c(2, 1)),
    rbind(c(1, 3), c(1, 1)),
    rbind(c(-1, 0), c(0, 1)),
    rbind(c(-1, 2), c(0, 1)),
    rbind(c(0, 3), c(1, 3)),
    rbind(c(-2, 0), c(-1, 0))
  )
  weights <- c(1, 1, 1, 2, 1, 3, 1)
  if (extra_upstream_source) {
    edges[[8]] <- rbind(c(-2, 2), c(-1, 0))
    weights <- c(weights, 2)
  }
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_edge_weights(
    weights = data.frame(w = weights), directional_weights = "w"
  )
  graph$setDirectionalWeightFunction()
  graph$add_observations(
    data = data.frame(
      temp = seq_len(6L) / 7,
      obs_id = sprintf("tree-%02d", seq_len(6L)),
      edge_number = rep(1:3, each = 2L),
      distance_on_edge = rep(c(0.23, 0.71), 3L)
    ),
    normalized = TRUE,
    verbose = 0
  )
  graph
}

subgraph_equivalence_cases <- function() {
  chain_edges <- lapply(0:3, function(index) {
    rbind(c(index, 0), c(index + 1, 0))
  })
  two_source_edges <- list(
    rbind(c(-1, 1), c(0, 0)),
    rbind(c(-1, -1), c(0, 0)),
    rbind(c(0, 0), c(1, 0))
  )
  nested_edges <- list(
    rbind(c(1, 0), c(2, 0)),
    rbind(c(2, 0), c(3, 0)),
    rbind(c(0, 0), c(1, 0)),
    rbind(c(-1, 0), c(0, 0)),
    rbind(c(-2, 1), c(-1, 0)),
    rbind(c(-2, -1), c(-1, 0)),
    rbind(c(0, -1), c(0, 0))
  )
  confluence_edges <- list(
    rbind(c(-1, 1), c(0, 0)),
    rbind(c(-1, -1), c(0, 0)),
    rbind(c(0, 0), c(1, 0))
  )
  outward_flux_edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(2, 1)),
    rbind(c(1, 0), c(2, -1))
  )

  list(
    no_pruning = make_equivalence_graph(chain_edges, 1:4),
    one_soft_boundary = make_equivalence_graph(chain_edges, 4),
    two_source_merge = make_equivalence_graph(
      two_source_edges, 3, weights = c(1, 2, 1)
    ),
    two_source_k2 = {
      graph <- make_equivalence_graph(
        two_source_edges, 3, weights = c(1, 4, 1)
      )
      graph$setDirectionalWeightFunction(
        f_in = function(weight) sqrt(weight / sum(weight))
      )
      graph
    },
    two_source_custom = {
      graph <- make_equivalence_graph(
        two_source_edges, 3, weights = c(1, 4, 1)
      )
      graph$setDirectionalWeightFunction(
        f_in = function(weight) 0.6 * sqrt(weight / sum(weight)),
        f_out = function(weight) rep(-0.8, length(weight))
      )
      graph
    },
    nested_three_hop = make_equivalence_graph(nested_edges, 1:2),
    inward_confluence = make_equivalence_graph(confluence_edges, c(1, 3)),
    outward_flux = make_equivalence_graph(
      outward_flux_edges, 1:3, weights = c(1, 1, 1)
    ),
    seven_edge = make_equivalence_tree_graph(),
    shared_suffix_eight_edge = make_equivalence_tree_graph(
      extra_upstream_source = TRUE
    )
  )
}
