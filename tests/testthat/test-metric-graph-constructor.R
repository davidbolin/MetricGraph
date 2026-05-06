test_that("constructor builds planar grid with expected dimensions", {
  edges <- make_grid_edges(5)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  # 5x5 grid has 25 unique vertices and 2 * 5 * 4 = 40 edges
  expect_equal(g$nV, 25L)
  expect_equal(g$nE, 40L)
})

test_that("constructor merges near-duplicate vertices when perform_merges = TRUE", {
  # Two edges sharing an endpoint that differ by eps < tolerance
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1 + 1e-6, 0), c(2, 0))
  )
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         tolerance = list(vertex_vertex = 1e-4,
                                           vertex_edge = 0,
                                           edge_edge = 0),
                         check_connected = FALSE, verbose = 0)
  # Should merge: 3 unique vertices, 2 edges
  expect_equal(g$nV, 3L)
  expect_equal(g$nE, 2L)
})

test_that("constructor does NOT merge vertices when perform_merges = FALSE", {
  # Two edges with endpoints separated by a distance that should NOT be merged
  # (eps must be larger than any epsilon the constructor uses for deduping
  # coincident points — 1e-2 is safely above any numerical-coincidence threshold).
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1 + 1e-2, 0), c(2, 0))
  )
  g <- metric_graph$new(edges = edges, perform_merges = FALSE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nV, 4L)
  expect_equal(g$nE, 2L)
})

test_that("constructor accepts sf LINESTRING input and sets longlat = TRUE", {
  skip_if_not_installed("sf")
  edges_sf <- make_longlat_grid(5)
  g <- metric_graph$new(edges = edges_sf, perform_merges = TRUE,
                         tolerance = list(vertex_vertex = 1e-5,
                                           vertex_edge = 0,
                                           edge_edge = 0),
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nV, 25L)
  expect_equal(g$nE, 40L)
  # longlat flag should propagate to edge attribute after set_edge_weights
  expect_true(isTRUE(attr(g$edges[[1L]], "longlat")))
})

test_that("constructor computes edge lengths matching polyline Euclidean lengths (planar)", {
  edges <- make_grid_edges(4)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  for (i in seq_len(g$nE)) {
    expected <- polyline_length(g$edges[[i]])
    expect_equal(g$edge_lengths[i], expected, tolerance = 1e-12,
                 info = sprintf("edge %d length mismatch", i))
  }
})

test_that("constructor produces edges with valid PtE attributes", {
  edges <- make_grid_edges(4)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  ok <- all_pte_valid(g)
  expect_true(all(ok),
              info = sprintf("edges with invalid PtE: %s",
                             paste(which(!ok), collapse = ", ")))
})

test_that("constructor produces valid PtE on polylines with multiple interior points", {
  # A single edge with 4 polyline points
  edges <- list(
    cbind(c(0, 1, 3, 5), c(0, 0, 1, 1))
  )
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nE, 1L)
  pte <- attr(g$edges[[1L]], "PtE")
  expect_equal(pte[1], 0)
  expect_equal(pte[length(pte)], 1)
  expect_true(all(diff(pte) >= 0))
  # Interior PtE values should be cumulative arc lengths normalized
  expected_pte <- cumsum(c(0, sqrt(c(1^2, 2^2 + 1^2, 2^2))))
  expected_pte <- expected_pte / expected_pte[length(expected_pte)]
  expect_equal(pte, expected_pte, tolerance = 1e-12)
})

test_that("constructor removes degenerate (zero-length) edges", {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(1, 0)),  # zero-length: both points identical
    rbind(c(1, 0), c(2, 0))
  )
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nE, 2L)
  expect_true(all(g$edge_lengths > 0))
})

test_that("constructor deduplicates consecutive duplicate points in polylines", {
  # A polyline with a consecutive duplicate point
  edges <- list(
    rbind(c(0, 0), c(1, 0), c(1, 0), c(2, 0))
  )
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nE, 1L)
  # The duplicate point should be removed
  expect_true(nrow(g$edges[[1L]]) < 4L)
})

test_that("constructor produces consistent E matrix linking to edges", {
  edges <- make_grid_edges(4)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  # E[i,] = (start_vertex, end_vertex) for edge i; both must be in [1, nV]
  expect_equal(nrow(g$E), g$nE)
  expect_true(all(g$E >= 1L & g$E <= g$nV))
  # Each edge's polyline start/end should match the V coordinates at E[i,1] and E[i,2]
  for (i in seq_len(g$nE)) {
    expect_equal(unname(g$V[g$E[i, 1L], ]),
                 unname(g$edges[[i]][1L, ]),
                 tolerance = 1e-10)
    expect_equal(unname(g$V[g$E[i, 2L], ]),
                 unname(g$edges[[i]][nrow(g$edges[[i]]), ]),
                 tolerance = 1e-10)
  }
})

test_that("constructor is deterministic across repeated calls with same input", {
  edges <- make_grid_edges(5)
  g1 <- metric_graph$new(edges = edges, perform_merges = TRUE,
                          check_connected = FALSE, verbose = 0)
  g2 <- metric_graph$new(edges = edges, perform_merges = TRUE,
                          check_connected = FALSE, verbose = 0)
  expect_equal(g1$nV, g2$nV)
  expect_equal(g1$nE, g2$nE)
  expect_equal(sort(g1$edge_lengths), sort(g2$edge_lengths))
})
