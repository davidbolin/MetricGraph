## Tests for the get_components / which_component / plot(components=...)
## methods added to metric_graph (so users can work with disconnected
## graphs without the graph_components wrapper).

library(MetricGraph)
library(testthat)

make_disconnected_mg <- function() {
  e1 <- rbind(c(0, 0), c(1, 0))
  e2 <- rbind(c(1, 0), c(2, 0))
  e3 <- rbind(c(10, 10), c(11, 10))
  e4 <- rbind(c(11, 10), c(12, 10))
  metric_graph$new(edges = list(e1, e2, e3, e4),
                   verbose = 0, perform_merges = TRUE,
                   check_connected = FALSE)
}

make_connected_mg <- function() {
  e1 <- rbind(c(0, 0), c(1, 0))
  e2 <- rbind(c(1, 0), c(2, 0))
  metric_graph$new(edges = list(e1, e2), verbose = 0,
                   perform_merges = TRUE)
}

test_that("get_components returns one entry per connected component", {
  mg <- make_disconnected_mg()
  comps <- mg$get_components()
  expect_length(comps, 2L)
  expect_true(all(vapply(comps, inherits, logical(1), what = "metric_graph")))
  expect_equal(sort(vapply(comps, function(g) g$nE, integer(1))),
               c(2L, 2L))
})

test_that("get_components is identity for a connected graph", {
  mg <- make_connected_mg()
  comps <- mg$get_components()
  expect_length(comps, 1L)
  expect_identical(comps[[1]], mg)
})

test_that("get_components sorts by total edge length, descending", {
  e1 <- rbind(c(0, 0), c(5, 0))   # length 5
  e2 <- rbind(c(5, 0), c(10, 0))  # length 5  -> comp A, total 10
  e3 <- rbind(c(20, 0), c(23, 0)) # comp B, total 3
  e4 <- rbind(c(40, 0), c(41, 0)) # comp C, total 1
  mg <- metric_graph$new(edges = list(e1, e2, e3, e4),
                         verbose = 0, perform_merges = TRUE,
                         check_connected = FALSE)
  comps <- mg$get_components()
  expect_length(comps, 3L)
  totals <- vapply(comps,
                   function(g) sum(as.numeric(g$edge_lengths)),
                   numeric(1))
  expect_equal(totals, c(10, 3, 1))
})

test_that("get_components routes existing observations to the right component", {
  mg <- make_disconnected_mg()
  df <- data.frame(
    coord_x = c(0.5, 10.5, 1.5, 11.5),
    coord_y = c(0,   10,   0,   10),
    y = c(1, 2, 3, 4)
  )
  mg$add_observations(data = df, data_coords = "spatial",
                      verbose = 0, suppress_warnings = TRUE)
  comps <- mg$get_components()
  d1 <- comps[[1]]$get_data(format = "list")
  d2 <- comps[[2]]$get_data(format = "list")
  expect_equal(sort(d1$y), c(1, 3))
  expect_equal(sort(d2$y), c(2, 4))
})

## ── add_observations: spatial / PtE / sf on a disconnected graph ────────────

test_that("add_observations (spatial) does not warn or error on a disconnected graph", {
  mg <- make_disconnected_mg()
  df <- data.frame(coord_x = c(0.5, 10.5),
                   coord_y = c(0, 10), y = c(1, 2))
  expect_no_warning(
    mg$add_observations(data = df, data_coords = "spatial",
                        verbose = 0, suppress_warnings = TRUE)
  )
  expect_equal(nrow(as.data.frame(mg$get_data())), 2L)
})

test_that("add_observations (PtE) lands on the correct edges of a disconnected graph", {
  # Edge indexing is global. With `make_disconnected_mg()`, edges 1-2 are on
  # component 1 (rows 1,0)-(2,0), edges 3-4 are on component 2 (rows
  # 10,10)-(12,10). Using PtE input with edge_number = c(1, 4) should
  # snap obs to component 1 and component 2 respectively.
  mg <- make_disconnected_mg()
  df <- data.frame(
    edge_number = c(1L, 2L, 3L, 4L),
    distance_on_edge = c(0.25, 0.5, 0.5, 0.75),
    y = c(10, 20, 30, 40)
  )
  mg$add_observations(data = df, data_coords = "PtE",
                      normalized = TRUE, verbose = 0,
                      suppress_warnings = TRUE)

  out <- mg$get_data(format = "list")
  expect_equal(length(out$y), 4L)
  expect_setequal(out$.edge_number, c(1L, 2L, 3L, 4L))

  # Round-trip via get_components: edges 1-2 land in component 1, edges 3-4
  # in component 2 (per-component edge numbers are 1 and 2 within each).
  comps <- mg$get_components()
  d1 <- comps[[1]]$get_data(format = "list")
  d2 <- comps[[2]]$get_data(format = "list")
  expect_setequal(d1$y, c(10, 20))
  expect_setequal(d2$y, c(30, 40))
})

test_that("add_observations errors on out-of-range PtE edge_number", {
  mg <- make_disconnected_mg()
  df <- data.frame(edge_number = c(1L, 999L),
                   distance_on_edge = c(0.5, 0.5), y = c(1, 2))
  expect_error(
    mg$add_observations(data = df, data_coords = "PtE",
                        normalized = TRUE, verbose = 0,
                        suppress_warnings = TRUE)
  )
})

test_that("add_observations (sf) on a disconnected graph routes to the right component", {
  skip_if_not_installed("sf")
  mg <- make_disconnected_mg()
  pts <- sf::st_as_sf(
    data.frame(x = c(0.5, 10.5), y = c(0, 10), val = c(7, 8)),
    coords = c("x", "y")
  )
  mg$add_observations(data = pts, verbose = 0, suppress_warnings = TRUE)
  comps <- mg$get_components()
  d1 <- comps[[1]]$get_data(format = "list")
  d2 <- comps[[2]]$get_data(format = "list")
  expect_equal(d1$val, 7)
  expect_equal(d2$val, 8)
})

## ── components_cache invalidation ───────────────────────────────────────────

test_that("get_components reflects new observations added after a previous call", {
  # The lazy components_cache must be invalidated by add_observations(); a
  # stale cache would return per-component graphs missing the new data.
  mg <- make_disconnected_mg()
  mg$add_observations(
    data = data.frame(coord_x = 0.5, coord_y = 0, y = 1),
    data_coords = "spatial", verbose = 0, suppress_warnings = TRUE
  )
  comps_first <- mg$get_components()
  expect_equal(length(comps_first[[1]]$get_data(format = "list")$y), 1L)

  mg$add_observations(
    data = data.frame(coord_x = 10.5, coord_y = 10, y = 2),
    data_coords = "spatial", verbose = 0, suppress_warnings = TRUE
  )
  comps_second <- mg$get_components()
  # Component 2 now has the second observation
  d2 <- comps_second[[2]]$get_data(format = "list")
  expect_equal(d2$y, 2)
})

test_that("get_components reflects clear_observations", {
  mg <- make_disconnected_mg()
  mg$add_observations(
    data = data.frame(coord_x = c(0.5, 10.5), coord_y = c(0, 10),
                      y = c(1, 2)),
    data_coords = "spatial", verbose = 0, suppress_warnings = TRUE
  )
  comps_before <- mg$get_components()
  expect_false(is.null(comps_before[[1]]$.__enclos_env__$private$data))

  mg$clear_observations()
  comps_after <- mg$get_components()
  for (k in seq_along(comps_after)) {
    expect_true(is.null(comps_after[[k]]$.__enclos_env__$private$data))
  }
})

## ── is_disconnected ─────────────────────────────────────────────────────────

test_that("is_disconnected reflects the cached component count", {
  expect_true(make_disconnected_mg()$is_disconnected())
  expect_false(make_connected_mg()$is_disconnected())
})

## ── cache invalidation under structural mutators ────────────────────────────

test_that("prune_vertices refreshes the component cache (membership length matches nV)", {
  edges <- list(
    rbind(c(0, 0),   c(0.5, 0)),  rbind(c(0.5, 0),  c(1, 0)),
    rbind(c(1, 0),   c(1.5, 0)),
    rbind(c(10, 0),  c(10.5, 0)), rbind(c(10.5, 0), c(11, 0))
  )
  mg <- metric_graph$new(edges = edges, check_connected = FALSE,
                         perform_merges = TRUE, verbose = 0)
  expect_true(mg$is_disconnected())
  nV_before <- mg$nV
  expect_equal(length(mg$.__enclos_env__$private$component_membership),
               nV_before)

  mg$prune_vertices(verbose = FALSE)

  # Membership vector must match the new nV (otherwise which_component()
  # would index past the end of the new vertex list).
  expect_equal(length(mg$.__enclos_env__$private$component_membership),
               mg$nV)
  expect_true(mg$is_disconnected())
  # which_component must still route correctly to the surviving components.
  expect_equal(mg$which_component(rbind(c(0.5, 0), c(10.5, 0))),
               c(1L, 2L))
})

test_that("remove_small_circles refreshes the component cache", {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(2, 0)),
    rbind(c(10, 0), c(11, 0))
  )
  mg <- metric_graph$new(edges = edges, check_connected = FALSE,
                         verbose = 0)
  # No-op call (no circular edges below the tolerance) — the cache should
  # still be valid afterwards, and the membership vector size should match
  # nV.
  mg$remove_small_circles(tolerance = 1e-12, verbose = 0)
  expect_equal(length(mg$.__enclos_env__$private$component_membership),
               mg$nV)
  expect_true(mg$is_disconnected())
})

test_that("set_manual_edge_lengths refreshes the component ordering", {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),     # comp A — total length 1
    rbind(c(10, 0), c(15, 0))    # comp B — total length 5
  )
  mg <- metric_graph$new(edges = edges, check_connected = FALSE,
                         verbose = 0)
  # By default comp B (length 5) is component 1 (largest); a point near it
  # should map to component 1.
  expect_equal(mg$which_component(c(12, 0)), 1L)

  # Flip the lengths so that comp A is now larger.
  mg$set_manual_edge_lengths(edge_lengths = c(10, 1))
  # Component ordering should be re-sorted — comp A (now length 10) is
  # component 1; the point at (12, 0) is still spatially closest to comp B,
  # which is now component 2.
  expect_equal(mg$which_component(c(12, 0)), 2L)
  expect_equal(mg$which_component(c(0.5, 0)), 1L)
})

test_that("which_component routes spatial points to the closest component", {
  mg <- make_disconnected_mg()
  XY <- rbind(c(0.5, 0.0), c(10.5, 10.0),
              c(1.7, 0.1), c(11.9, 9.9))
  expect_equal(mg$which_component(XY), c(1L, 2L, 1L, 2L))
})

test_that("which_component handles a vector of length 2 and returns 1 for connected graph", {
  mg <- make_disconnected_mg()
  expect_equal(mg$which_component(c(11, 10)), 2L)
  expect_equal(mg$which_component(c(0.5, 0)), 1L)

  mg_c <- make_connected_mg()
  expect_equal(mg_c$which_component(rbind(c(100, 100), c(-50, 0))),
               c(1L, 1L))
})

test_that("which_component validates input shape", {
  mg <- make_disconnected_mg()
  expect_error(mg$which_component(c(1, 2, 3)), "length 2")
  expect_error(mg$which_component(matrix(1:6, nrow = 2)), "two columns")
})

test_that("plot(components = TRUE) returns a ggplot for disconnected graphs", {
  mg <- make_disconnected_mg()
  p <- mg$plot(components = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot(components = TRUE) on a connected graph falls back to default", {
  mg <- make_connected_mg()
  p <- mg$plot(components = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot(components = matrix) accepts user-supplied colors", {
  mg <- make_disconnected_mg()
  cols <- rbind(c(1, 0, 0), c(0, 0, 1))
  p <- mg$plot(components = cols)
  expect_s3_class(p, "ggplot")
})

test_that("plot(components = matrix) errors on wrong shape", {
  mg <- make_disconnected_mg()
  expect_error(mg$plot(components = rbind(c(1, 0, 0))),
               "n x 3 matrix")
})
