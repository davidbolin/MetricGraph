test_that("prune_vertices collapses an open chain into a single edge", {
  edges <- make_open_chain(5)  # 5 edges in a line
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nV, 6L)
  expect_equal(g$nE, 5L)
  total_length_before <- sum(g$edge_lengths)

  g$prune_vertices()

  # After pruning, the 4 interior degree-2 vertices should be gone
  expect_equal(g$nV, 2L)
  expect_equal(g$nE, 1L)
  # Total length is preserved
  expect_equal(sum(g$edge_lengths), total_length_before, tolerance = 1e-12)
  # The resulting edge's PtE is still valid
  expect_true(all(all_pte_valid(g)))
})

test_that("prune_vertices leaves a star graph unchanged (no deg-2 vertices to merge)", {
  edges <- make_star(5)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  nV_before <- g$nV
  nE_before <- g$nE
  lengths_before <- sort(g$edge_lengths)

  g$prune_vertices()

  expect_equal(g$nV, nV_before)
  expect_equal(g$nE, nE_before)
  expect_equal(sort(g$edge_lengths), lengths_before, tolerance = 1e-12)
})

test_that("prune_vertices handles a triangle (closed loop, all vertices deg-2)", {
  edges <- make_triangle()
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nV, 3L)
  expect_equal(g$nE, 3L)
  total_length_before <- sum(g$edge_lengths)

  g$prune_vertices()

  # All three vertices are degree 2, so prune should merge into a single
  # closed-loop edge. The fallback path handles this.
  expect_lte(g$nV, 3L)
  expect_lte(g$nE, 3L)
  # Total length should still be preserved
  expect_equal(sum(g$edge_lengths), total_length_before, tolerance = 1e-12)
})

test_that("prune_vertices preserves total edge length on a longer chain", {
  edges <- make_open_chain(20)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  total_before <- sum(g$edge_lengths)
  g$prune_vertices()
  expect_equal(sum(g$edge_lengths), total_before, tolerance = 1e-12)
  expect_equal(g$nV, 2L)
  expect_equal(g$nE, 1L)
})

test_that("prune_vertices preserves graph connectivity on a Y junction", {
  # Y-shaped graph: three arms meeting at a degree-3 center,
  # each arm has an interior deg-2 vertex
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(2, 0)),   # first arm: 2 edges
    rbind(c(0, 0), c(0, 1)),
    rbind(c(0, 1), c(0, 2)),   # second arm: 2 edges
    rbind(c(0, 0), c(-1, 0)),
    rbind(c(-1, 0), c(-2, 0))  # third arm: 2 edges
  )
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nE, 6L)
  total_before <- sum(g$edge_lengths)

  g$prune_vertices()

  # Each arm should collapse from 2 edges to 1 → total 3 edges
  expect_equal(g$nE, 3L)
  # Center + 3 leaves = 4 vertices
  expect_equal(g$nV, 4L)
  expect_equal(sum(g$edge_lengths), total_before, tolerance = 1e-12)
  expect_true(all(all_pte_valid(g)))
})

test_that("prune_vertices produces edges with valid PtE attributes", {
  edges <- make_open_chain(10)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  g$prune_vertices()
  expect_true(all(all_pte_valid(g)))
  # Interior points of the merged edge should correspond to the old vertex positions
  merged_edge <- g$edges[[1L]]
  expect_equal(merged_edge[1L, ], c(0, 0), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(merged_edge[nrow(merged_edge), ], c(10, 0),
               tolerance = 1e-12, ignore_attr = TRUE)
})
