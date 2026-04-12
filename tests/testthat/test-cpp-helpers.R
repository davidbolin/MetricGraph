# Direct tests of the C++ helpers.

test_that("compute_edge_lengths_cpp matches polyline Euclidean length (planar)", {
  edges <- list(
    rbind(c(0, 0), c(3, 4)),               # length 5
    cbind(c(0, 1, 3, 5), c(0, 0, 0, 0)),   # length 5
    rbind(c(0, 0), c(1, 1), c(2, 0))       # length 2*sqrt(2)
  )
  lens <- MetricGraph:::compute_edge_lengths_cpp(edges, FALSE)
  expect_equal(lens[1], 5, tolerance = 1e-12)
  expect_equal(lens[2], 5, tolerance = 1e-12)
  expect_equal(lens[3], 2 * sqrt(2), tolerance = 1e-12)
})

test_that("compute_edge_lengths_cpp returns non-negative lengths on lon/lat", {
  # A 1-degree edge at the equator should be ~111 km
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(0, 0), c(0, 1))
  )
  lens <- MetricGraph:::compute_edge_lengths_cpp(edges, TRUE)
  # Haversine on sphere of radius 6371008.8 m gives ~111.195 km per degree
  expect_equal(lens[1] / 1000, 111.195, tolerance = 0.5)
  expect_equal(lens[2] / 1000, 111.195, tolerance = 0.5)
})

test_that("compute_PtE_edges_cpp produces valid PtE attributes anchored at [0, 1]", {
  edges <- list(
    rbind(c(0, 0), c(1, 0), c(3, 0), c(5, 0)),
    rbind(c(0, 0), c(3, 4))
  )
  pte_list <- MetricGraph:::compute_PtE_edges_cpp(edges, FALSE)
  for (i in seq_along(pte_list)) {
    pte <- pte_list[[i]]
    expect_equal(length(pte), nrow(edges[[i]]))
    expect_equal(pte[1], 0, tolerance = 1e-12)
    expect_equal(pte[length(pte)], 1, tolerance = 1e-12)
    expect_true(all(diff(pte) >= -1e-12))
  }
  # First edge: cumulative lengths 0, 1, 3, 5 → PtE = 0, 0.2, 0.6, 1
  expect_equal(pte_list[[1]], c(0, 0.2, 0.6, 1), tolerance = 1e-12)
})

test_that("aeqd_project_cpp projects the origin to (0, 0)", {
  pts <- matrix(c(2.35, 48.86), nrow = 1)
  proj <- MetricGraph:::aeqd_project_cpp(pts, 2.35, 48.86)
  expect_equal(proj[1, 1], 0, tolerance = 1e-6)
  expect_equal(proj[1, 2], 0, tolerance = 1e-6)
})

test_that("aeqd_project_cpp projects eastward and northward points with expected signs", {
  # Point 0.001 degrees east of origin should have positive x, near-zero y
  pts <- matrix(c(2.35 + 0.001, 48.86), nrow = 1)
  proj <- MetricGraph:::aeqd_project_cpp(pts, 2.35, 48.86)
  expect_true(proj[1, 1] > 0)
  expect_true(abs(proj[1, 2]) < 1)
  # Point 0.001 degrees north of origin should have positive y, near-zero x
  pts <- matrix(c(2.35, 48.86 + 0.001), nrow = 1)
  proj <- MetricGraph:::aeqd_project_cpp(pts, 2.35, 48.86)
  expect_true(proj[1, 2] > 0)
  expect_true(abs(proj[1, 1]) < 1)
})

test_that("postprocess_edges_cpp removes consecutive duplicate interior points", {
  edges <- list(
    rbind(c(0, 0), c(1, 0), c(1, 0), c(2, 0)),        # duplicate interior
    rbind(c(0, 0), c(1, 1), c(2, 2))                  # no duplicates
  )
  out <- MetricGraph:::postprocess_edges_cpp(edges)
  expect_lte(nrow(out[[1]]), 3L)
  expect_equal(nrow(out[[2]]), 3L)
  # No two consecutive rows of any output should be identical
  for (i in seq_along(out)) {
    e <- out[[i]]
    if (nrow(e) >= 2L) {
      diffs <- rowSums((e[-1, , drop = FALSE] - e[-nrow(e), , drop = FALSE])^2)
      expect_true(all(diffs > 0),
                  info = sprintf("edge %d still has consecutive duplicates", i))
    }
  }
})

# -----------------------------------------------------------------------------
# split_one_edge_cpp — the per-edge helper used by the new split_edge_batch
# -----------------------------------------------------------------------------

test_that("split_one_edge_cpp with a single split point creates two segments of correct length", {
  edge <- rbind(c(0, 0), c(10, 0))
  PtE  <- c(0, 1)
  res  <- MetricGraph:::split_one_edge_cpp(edge, PtE, t_values = 0.3,
                                             edge_len = 10, E_row = c(1L, 2L),
                                             first_new_v = 5L)
  expect_equal(res$segment_lengths, c(3, 7), tolerance = 1e-12)
  # coords_list1 should end at the split point (3, 0)
  expect_equal(res$coords_list1[nrow(res$coords_list1), ], c(3, 0),
               tolerance = 1e-12, ignore_attr = TRUE)
  # coords_list2[[1]] should start at (3, 0) and end at (10, 0)
  cl2 <- res$coords_list2[[1]]
  expect_equal(cl2[1, ], c(3, 0), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(cl2[nrow(cl2), ], c(10, 0), tolerance = 1e-12, ignore_attr = TRUE)
  # aux_matrix
  expect_equal(res$aux_matrix, matrix(c(1L, 5L, 5L, 2L), ncol = 2L, byrow = TRUE))
})

test_that("split_one_edge_cpp with three splits produces four segments, lengths summing to edge_len", {
  edge <- rbind(c(0, 0), c(10, 0))
  PtE  <- c(0, 1)
  t_values <- c(0.2, 0.5, 0.8)
  res <- MetricGraph:::split_one_edge_cpp(edge, PtE, t_values,
                                             edge_len = 10, E_row = c(3L, 7L),
                                             first_new_v = 20L)
  expect_equal(sum(res$segment_lengths), 10, tolerance = 1e-12)
  expect_length(res$segment_lengths, 4L)
  expect_length(res$coords_list2, 3L)
  # aux_matrix should have 4 rows chaining the new vertices
  expect_equal(dim(res$aux_matrix), c(4L, 2L))
  expect_equal(res$aux_matrix[1, ], c(3L, 20L))
  expect_equal(res$aux_matrix[2, ], c(20L, 21L))
  expect_equal(res$aux_matrix[3, ], c(21L, 22L))
  expect_equal(res$aux_matrix[4, ], c(22L, 7L))
})

test_that("split_one_edge_cpp PtE attributes on output segments are anchored at [0, 1]", {
  edge <- cbind(c(0, 1, 3, 5), c(0, 0, 0, 0))
  PtE  <- c(0, 0.2, 0.6, 1)
  res <- MetricGraph:::split_one_edge_cpp(edge, PtE, t_values = c(0.3, 0.7),
                                            edge_len = 5, E_row = c(1L, 2L),
                                            first_new_v = 10L)
  pte1 <- attr(res$coords_list1, "PtE")
  expect_equal(pte1[1], 0, tolerance = 1e-12)
  expect_equal(pte1[length(pte1)], 1, tolerance = 1e-12)
  expect_true(all(diff(pte1) >= -1e-12))
  for (j in seq_along(res$coords_list2)) {
    p <- attr(res$coords_list2[[j]], "PtE")
    expect_equal(p[1], 0, tolerance = 1e-12)
    expect_equal(p[length(p)], 1, tolerance = 1e-12)
    expect_true(all(diff(p) >= -1e-12))
  }
})

test_that("split_one_edge_cpp handles a multi-point polyline with one split inside segment 2", {
  edge <- cbind(c(0, 1, 3, 5), c(0, 0, 0, 0))
  PtE  <- c(0, 0.2, 0.6, 1)
  res <- MetricGraph:::split_one_edge_cpp(edge, PtE, t_values = 0.4,
                                             edge_len = 5, E_row = c(1L, 2L),
                                             first_new_v = 10L)
  # Split at t = 0.4 (which is halfway through segment 2: between PtE=0.2 and PtE=0.6)
  expect_equal(sum(res$segment_lengths), 5, tolerance = 1e-12)
  expect_equal(res$segment_lengths, c(2, 3), tolerance = 1e-12)
  # coords_list1 should include the original first polyline point (0,0), the
  # interior vertex (1,0) and the split point (2,0)
  cl1 <- res$coords_list1
  expect_equal(cl1[nrow(cl1), ], c(2, 0), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("split_one_edge_cpp with split exactly on an interior vertex returns a valid structure", {
  edge <- cbind(c(0, 1, 3, 5), c(0, 0, 0, 0))
  PtE  <- c(0, 0.2, 0.6, 1)
  # t = 0.6 is exactly on the third polyline point
  res <- MetricGraph:::split_one_edge_cpp(edge, PtE, t_values = 0.6,
                                            edge_len = 5, E_row = c(1L, 2L),
                                            first_new_v = 10L)
  expect_equal(sum(res$segment_lengths), 5, tolerance = 1e-12)
  # No negative lengths or NaN PtE
  expect_true(all(res$segment_lengths >= 0))
  expect_true(all(!is.na(attr(res$coords_list1, "PtE"))))
  expect_true(all(!is.na(attr(res$coords_list2[[1]], "PtE"))))
})

test_that("split_edges_batch_cpp matches split_one_edge_cpp on each group", {
  edges_full <- list(
    rbind(c(0, 0), c(10, 0)),
    cbind(c(0, 1, 3, 5), c(0, 0, 0, 0)),
    rbind(c(0, 0), c(1, 1), c(2, 0))
  )
  PtE_full <- list(
    c(0, 1),
    c(0, 0.2, 0.6, 1),
    c(0, 0.5, 1)
  )
  E_full <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 2L, byrow = TRUE)
  edge_lens <- c(10, 5, 2 * sqrt(2))
  edge_ids  <- c(1L, 2L, 3L)
  t_values_list <- list(
    c(0.3, 0.7),
    0.5,
    c(0.25, 0.75)
  )
  first_new_vs <- c(100L, 200L, 300L)

  batch_res <- MetricGraph:::split_edges_batch_cpp(
    edges_full, PtE_full, E_full, edge_lens, edge_ids, t_values_list, first_new_vs
  )
  expect_length(batch_res, 3L)
  # Each group result matches a direct split_one_edge_cpp call
  for (g in seq_len(3L)) {
    one_res <- MetricGraph:::split_one_edge_cpp(
      edges_full[[g]], PtE_full[[g]], t_values_list[[g]],
      edge_lens[g], E_full[g, ], first_new_vs[g]
    )
    expect_equal(batch_res[[g]]$val_lines, one_res$val_lines, tolerance = 1e-12)
    expect_equal(batch_res[[g]]$segment_lengths, one_res$segment_lengths, tolerance = 1e-12)
    expect_equal(batch_res[[g]]$aux_matrix, one_res$aux_matrix)
    expect_equal(unname(batch_res[[g]]$coords_list1),
                 unname(one_res$coords_list1), tolerance = 1e-12)
    expect_equal(attr(batch_res[[g]]$coords_list1, "PtE"),
                 attr(one_res$coords_list1, "PtE"), tolerance = 1e-12)
    for (j in seq_along(batch_res[[g]]$coords_list2)) {
      expect_equal(unname(batch_res[[g]]$coords_list2[[j]]),
                   unname(one_res$coords_list2[[j]]), tolerance = 1e-12)
      expect_equal(attr(batch_res[[g]]$coords_list2[[j]], "PtE"),
                   attr(one_res$coords_list2[[j]], "PtE"), tolerance = 1e-12)
    }
  }
})

test_that("compute_edge_lengths_cpp handles a polyline with many points", {
  # Zigzag with 10 points
  coords <- cbind(0:9, rep(c(0, 1), times = 5))
  edges <- list(coords)
  expected <- sum(sqrt(rowSums(diff(coords)^2)))
  lens <- MetricGraph:::compute_edge_lengths_cpp(edges, FALSE)
  expect_equal(lens[1], expected, tolerance = 1e-12)
})

test_that("compute_PtE_edges_cpp with a 2-point edge returns c(0, 1)", {
  edges <- list(rbind(c(0, 0), c(5, 0)))
  pte <- MetricGraph:::compute_PtE_edges_cpp(edges, FALSE)
  expect_equal(pte[[1]], c(0, 1), tolerance = 1e-12)
})

# -----------------------------------------------------------------------------
# split_one_edge_cpp edge-case tests
# -----------------------------------------------------------------------------

test_that("split_one_edge_cpp produces val_lines matching interpolation", {
  edge <- cbind(c(0, 2, 5, 10), c(0, 0, 0, 0))
  PtE  <- c(0, 0.2, 0.5, 1)
  t_vals <- c(0.1, 0.35, 0.75)
  res <- MetricGraph:::split_one_edge_cpp(edge, PtE, t_vals,
                                            edge_len = 10,
                                            E_row = c(1L, 2L),
                                            first_new_v = 10L)
  # Each val_line should lie on the original polyline and its x-coordinate
  # should be t_vals[i] * 10 (for this purely horizontal edge)
  for (i in seq_along(t_vals)) {
    expect_equal(res$val_lines[i, 1], t_vals[i] * 10, tolerance = 1e-12)
    expect_equal(res$val_lines[i, 2], 0, tolerance = 1e-12)
  }
})

test_that("split_one_edge_cpp aux_matrix connects start to first new, new to new, last new to end", {
  edge <- rbind(c(0, 0), c(10, 0))
  res <- MetricGraph:::split_one_edge_cpp(edge, c(0, 1), c(0.25, 0.5, 0.75),
                                            edge_len = 10,
                                            E_row = c(42L, 99L),
                                            first_new_v = 500L)
  aux <- res$aux_matrix
  # Row 1: start_v -> first_new_v
  expect_equal(aux[1, ], c(42L, 500L))
  # Middle rows: chain of new vertices
  expect_equal(aux[2, ], c(500L, 501L))
  expect_equal(aux[3, ], c(501L, 502L))
  # Last row: last new vertex -> end_v
  expect_equal(aux[4, ], c(502L, 99L))
})

test_that("split_edges_batch_cpp preserves per-group ordering in the output list", {
  edges_full <- list(
    rbind(c(0, 0), c(10, 0)),
    rbind(c(0, 0), c(20, 0))
  )
  PtE_full <- list(c(0, 1), c(0, 1))
  E_full <- matrix(c(1L, 2L, 3L, 4L), ncol = 2L, byrow = TRUE)
  # Process the edges in a non-sorted edge_ids order to check output ordering
  res <- MetricGraph:::split_edges_batch_cpp(
    edges_full, PtE_full, E_full, c(10, 20),
    edge_ids = c(2L, 1L),  # reversed
    t_values_list = list(0.5, 0.3),
    first_new_vs = c(100L, 200L)
  )
  # First group (edge_ids = 2, t = 0.5, edge len = 20) -> split point at x = 10
  expect_equal(res[[1]]$val_lines[1, 1], 10, tolerance = 1e-12)
  expect_equal(res[[1]]$segment_lengths, c(10, 10), tolerance = 1e-12)
  # Second group (edge_ids = 1, t = 0.3, edge len = 10) -> split point at x = 3
  expect_equal(res[[2]]$val_lines[1, 1], 3, tolerance = 1e-12)
  expect_equal(res[[2]]$segment_lengths, c(3, 7), tolerance = 1e-12)
})

test_that("split_edges_batch_cpp handles a mix of 1-split and many-split groups", {
  edge_mat <- rbind(c(0, 0), c(10, 0))
  edges_full <- list(edge_mat, edge_mat, edge_mat)
  PtE_full <- list(c(0, 1), c(0, 1), c(0, 1))
  E_full <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 2L, byrow = TRUE)
  t_values_list <- list(
    0.5,                    # 1 split
    seq(0.1, 0.9, by = 0.1), # 9 splits
    c(0.25, 0.5, 0.75)      # 3 splits
  )
  first_new_vs <- c(100L, 110L, 130L)
  res <- MetricGraph:::split_edges_batch_cpp(
    edges_full, PtE_full, E_full, rep(10, 3), c(1L, 2L, 3L),
    t_values_list, first_new_vs
  )
  expect_length(res[[1]]$coords_list2, 1L)
  expect_length(res[[2]]$coords_list2, 9L)
  expect_length(res[[3]]$coords_list2, 3L)
  for (g in seq_along(res)) {
    expect_equal(sum(res[[g]]$segment_lengths), 10, tolerance = 1e-12,
                 info = sprintf("group %d total length mismatch", g))
  }
})
