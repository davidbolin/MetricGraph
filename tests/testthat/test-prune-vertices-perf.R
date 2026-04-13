test_that("benchmark: large open chain (2000 edges)", {
  skip_on_cran()
  edges <- make_open_chain(2000)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  expect_equal(g$nV, 2001L)
  expect_equal(g$nE, 2000L)
  total_before <- sum(g$edge_lengths)

  t <- system.time(g$prune_vertices())

  expect_equal(g$nV, 2L)
  expect_equal(g$nE, 1L)
  expect_equal(sum(g$edge_lengths), total_before, tolerance = 1e-12)
  expect_true(all(all_pte_valid(g)))
  message(sprintf("  chain-2000: %.3f s", t[["elapsed"]]))
})

test_that("benchmark: large grid (30x30)", {
  skip_on_cran()
  edges <- make_grid_edges(30)
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  nV_before <- g$nV
  nE_before <- g$nE
  total_before <- sum(g$edge_lengths)

  t <- system.time(g$prune_vertices())

  expect_lte(g$nV, nV_before)
  expect_lte(g$nE, nE_before)
  expect_equal(sum(g$edge_lengths), total_before, tolerance = 1e-12)
  expect_true(all(all_pte_valid(g)))
  message(sprintf("  grid-30x30: %.3f s  (V: %d->%d, E: %d->%d)",
                  t[["elapsed"]], nV_before, g$nV, nE_before, g$nE))
})

test_that("benchmark: mixed problematic and unproblematic vertices", {
  skip_on_cran()
  # Build a graph with both "problematic" (directional) and normal degree-2 vertices.
  # Problematic vertices arise when both incident edges point the same direction
  # (both outgoing or both incoming), giving indegree==0 or outdegree==0 at a
  # degree-2 vertex.
  #
  # Layout: multiple star hubs connected by chains, with some chains having
  # edges that both point toward a shared vertex (creating problematic vertices).
  edges <- list()
  k <- 1L

  # Create 10 hub-and-spoke clusters connected by chains

  for (hub in 0:9) {
    cx <- hub * 50
    # Hub at (cx, 0) with 3 spokes (degree-3, never pruned)
    for (angle_idx in 1:3) {
      ang <- 2 * pi * angle_idx / 3
      edges[[k]] <- rbind(c(cx, 0), c(cx + 5 * cos(ang), 5 * sin(ang)))
      k <- k + 1L
    }

    if (hub < 9) {
      # Chain of 20 normal edges connecting hub to next hub
      chain_start <- cx + 5
      chain_end   <- (hub + 1) * 50 - 5
      xs <- seq(chain_start, chain_end, length.out = 21)
      for (i in seq_len(20)) {
        edges[[k]] <- rbind(c(xs[i], 0), c(xs[i + 1], 0))
        k <- k + 1L
      }

      # Add a pair of edges that both point TOWARD a shared vertex
      # (creating a problematic degree-2 vertex at (cx + 25, 10))
      mid_x <- cx + 25
      edges[[k]] <- rbind(c(mid_x - 3, 10), c(mid_x, 10))  # points right
      k <- k + 1L
      edges[[k]] <- rbind(c(mid_x + 3, 10), c(mid_x, 10))  # points left
      k <- k + 1L
    }
  }

  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  nV_before <- g$nV
  nE_before <- g$nE
  total_before <- sum(g$edge_lengths)

  # Count problematic vertices
  n_problematic <- sum(vapply(g$vertices,
    function(v) isTRUE(attr(v, "problematic")), logical(1)))

  t <- system.time(g$prune_vertices())

  expect_lte(g$nV, nV_before)
  expect_lte(g$nE, nE_before)
  expect_equal(sum(g$edge_lengths), total_before, tolerance = 1e-12)
  expect_true(all(all_pte_valid(g)))
  message(sprintf("  mixed: %.3f s  (V: %d->%d, E: %d->%d, problematic: %d)",
                  t[["elapsed"]], nV_before, g$nV, nE_before, g$nE, n_problematic))
})

test_that("benchmark: large cycle (500 edges)", {
  skip_on_cran()
  # Closed loop with 500 edges — all vertices degree-2, exercises fallback path
  n <- 500
  angles <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  pts <- cbind(cos(angles), sin(angles))
  edges <- lapply(seq_len(n), function(i) {
    j <- if (i < n) i + 1L else 1L
    rbind(pts[i, ], pts[j, ])
  })
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  total_before <- sum(g$edge_lengths)

  t <- system.time(g$prune_vertices())

  expect_equal(sum(g$edge_lengths), total_before, tolerance = 1e-10)
  expect_true(all(all_pte_valid(g)))
  message(sprintf("  cycle-500: %.3f s  (V: %d, E: %d)",
                  t[["elapsed"]], g$nV, g$nE))
})
