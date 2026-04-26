## Unit tests for graph_components methods that mirror metric_graph behaviour.
## Covers which_component, add_observations (PtE and spatial), get_data,
## clear_observations, get_PtE, get_groups, coordinates, build_mesh,
## compute_fem, plot, plot_function.

library(MetricGraph)
library(testthat)

# Helper: two well-separated components, each consisting of two collinear edges.
make_two_components <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(1, 0), c(2, 0))
  edge3 <- rbind(c(10, 10), c(11, 10))
  edge4 <- rbind(c(11, 10), c(12, 10))
  graph_components$new(edges = list(edge1, edge2, edge3, edge4),
                       verbose = 0, perform_merges = TRUE)
}

# Helper: three components of different sizes (sorted by length).
make_three_components <- function() {
  e1 <- rbind(c(0, 0), c(5, 0))   # comp A, length 5
  e2 <- rbind(c(5, 0), c(10, 0))  # comp A, length 5
  e3 <- rbind(c(20, 0), c(23, 0)) # comp B, length 3
  e4 <- rbind(c(40, 0), c(41, 0)) # comp C, length 1
  graph_components$new(edges = list(e1, e2, e3, e4),
                       verbose = 0, perform_merges = TRUE)
}

## ── which_component ──────────────────────────────────────────────────────────

test_that("which_component routes points to the closest component", {
  gc <- make_two_components()
  XY <- rbind(c(0.5, 0.0),   # comp 1
              c(10.5, 10.0), # comp 2
              c(1.7, 0.1),   # comp 1
              c(11.9, 9.9))  # comp 2
  wc <- gc$which_component(XY)
  expect_equal(wc, c(1L, 2L, 1L, 2L))
})

test_that("which_component handles a vector of length 2", {
  gc <- make_two_components()
  expect_equal(gc$which_component(c(11, 10)), 2L)
  expect_equal(gc$which_component(c(0.5, 0)), 1L)
})

test_that("which_component returns 1 for single-component graph", {
  edges <- list(rbind(c(0, 0), c(1, 0)), rbind(c(1, 0), c(2, 0)))
  gc <- graph_components$new(edges = edges, verbose = 0)
  expect_equal(gc$n, 1L)
  XY <- rbind(c(100, 100), c(-50, 0))
  expect_equal(gc$which_component(XY), c(1L, 1L))
})

test_that("which_component picks midpoint between components consistently", {
  gc <- make_two_components()
  # Point exactly at vertex of one component
  expect_equal(gc$which_component(c(2, 0)), 1L)
  expect_equal(gc$which_component(c(10, 10)), 2L)
})

test_that("which_component validates input shape", {
  gc <- make_two_components()
  expect_error(gc$which_component(c(1, 2, 3)), "length 2")
  expect_error(gc$which_component(matrix(1:6, nrow = 2)), "two columns")
})

## ── add_observations: spatial routing ────────────────────────────────────────

test_that("add_observations (spatial) routes each point to the nearest component", {
  gc <- make_two_components()
  df <- data.frame(
    coord_x = c(0.5, 10.5, 1.5, 11.5),
    coord_y = c(0,   10,   0,   10),
    y = c(1, 2, 3, 4)
  )
  gc$add_observations(data = df, data_coords = "spatial", verbose = 0)

  d1 <- gc$graphs[[1]]$get_data(format = "list")
  d2 <- gc$graphs[[2]]$get_data(format = "list")
  expect_equal(sort(d1$y), c(1, 3))
  expect_equal(sort(d2$y), c(2, 4))
})

test_that("add_observations (spatial) handles all points to a single component", {
  gc <- make_two_components()
  df <- data.frame(coord_x = c(0.2, 0.8, 1.5),
                   coord_y = c(0, 0, 0),
                   y = c(1, 2, 3))
  gc$add_observations(data = df, data_coords = "spatial", verbose = 0)
  expect_false(is.null(gc$graphs[[1]]$.__enclos_env__$private$data))
  expect_true(is.null(gc$graphs[[2]]$.__enclos_env__$private$data))
})

test_that("add_observations (spatial) supports sf input", {
  skip_if_not_installed("sf")
  gc <- make_two_components()
  pts <- sf::st_as_sf(
    data.frame(x = c(0.5, 10.5), y = c(0, 10), val = c(7, 8)),
    coords = c("x", "y"))
  gc$add_observations(data = pts, verbose = 0)
  d1 <- gc$graphs[[1]]$get_data(format = "list")
  d2 <- gc$graphs[[2]]$get_data(format = "list")
  expect_equal(d1$val, 7)
  expect_equal(d2$val, 8)
})

## ── add_observations: PtE convention ─────────────────────────────────────────

test_that("add_observations (PtE) splits data by component column", {
  gc <- make_two_components()
  df <- data.frame(
    component = c(1, 1, 2, 2),
    edge_number = c(1, 2, 1, 2),
    distance_on_edge = c(0.2, 0.7, 0.4, 0.5),
    y = c(10, 20, 30, 40)
  )
  gc$add_observations(data = df, data_coords = "PtE",
                      normalized = TRUE, verbose = 0)

  d1 <- gc$graphs[[1]]$get_data(format = "list")
  d2 <- gc$graphs[[2]]$get_data(format = "list")
  expect_equal(sort(d1$y), c(10, 20))
  expect_equal(sort(d2$y), c(30, 40))
  # Edge numbering on each component should be preserved
  expect_setequal(d1$.edge_number, c(1L, 2L))
  expect_setequal(d2$.edge_number, c(1L, 2L))
})

test_that("add_observations (PtE) errors on missing component column", {
  gc <- make_two_components()
  df <- data.frame(edge_number = 1, distance_on_edge = 0.5, y = 1)
  expect_error(
    gc$add_observations(data = df, data_coords = "PtE",
                        normalized = TRUE, verbose = 0),
    "component"
  )
})

test_that("add_observations (PtE) errors on out-of-range component", {
  gc <- make_two_components()
  df <- data.frame(component = c(1, 999), edge_number = c(1, 1),
                   distance_on_edge = c(0.5, 0.5), y = c(1, 2))
  expect_error(
    gc$add_observations(data = df, data_coords = "PtE",
                        normalized = TRUE, verbose = 0),
    "between 1 and"
  )
})

test_that("add_observations clear_obs wipes existing data first", {
  gc <- make_two_components()
  gc$add_observations(
    data = data.frame(coord_x = c(0.5, 10.5), coord_y = c(0, 10),
                      y = c(1, 2)),
    data_coords = "spatial", verbose = 0)
  gc$add_observations(
    data = data.frame(coord_x = 1.5, coord_y = 0, y = 99),
    data_coords = "spatial", clear_obs = TRUE, verbose = 0)
  d1 <- gc$graphs[[1]]$get_data(format = "list")
  expect_equal(d1$y, 99)
  expect_true(is.null(gc$graphs[[2]]$.__enclos_env__$private$data))
})

## ── get_data / get_PtE / get_groups / clear_observations ─────────────────────

test_that("get_data combines per-component data and adds a .component column", {
  gc <- make_two_components()
  df <- data.frame(
    coord_x = c(0.5, 10.5, 1.5, 11.5),
    coord_y = c(0, 10, 0, 10),
    y = c(1, 2, 3, 4)
  )
  gc$add_observations(data = df, data_coords = "spatial", verbose = 0)
  out <- gc$get_data(format = "tibble")
  expect_true(".component" %in% names(out))
  expect_equal(nrow(out), 4L)
  expect_setequal(out$.component, c(1L, 2L))
  # All comp 1 rows have y in {1, 3}; comp 2 has y in {2, 4}
  expect_setequal(out$y[out$.component == 1L], c(1, 3))
  expect_setequal(out$y[out$.component == 2L], c(2, 4))
})

test_that("get_data errors when no component has data", {
  gc <- make_two_components()
  expect_error(gc$get_data(), "No data found")
})

test_that("get_PtE combines per-component PtE with a component column", {
  gc <- make_two_components()
  gc$add_observations(
    data = data.frame(component = c(1, 2),
                      edge_number = c(1, 2),
                      distance_on_edge = c(0.2, 0.7),
                      y = c(1, 2)),
    data_coords = "PtE", normalized = TRUE, verbose = 0)
  pte <- gc$get_PtE()
  expect_equal(ncol(pte), 3L)
  expect_setequal(pte[, 1], c(1, 2))
})

test_that("get_groups returns the union of group identifiers", {
  gc <- make_two_components()
  df <- data.frame(
    coord_x = c(0.5, 10.5, 1.5, 11.5),
    coord_y = c(0, 10, 0, 10),
    rep = c("a", "a", "b", "b"),
    y = c(1, 2, 3, 4)
  )
  gc$add_observations(data = df, data_coords = "spatial",
                      group = "rep", verbose = 0)
  g <- gc$get_groups()
  expect_setequal(g, c("a", "b"))
})

test_that("clear_observations wipes data on every component", {
  gc <- make_two_components()
  gc$add_observations(
    data = data.frame(coord_x = c(0.5, 10.5), coord_y = c(0, 10),
                      y = c(1, 2)),
    data_coords = "spatial", verbose = 0)
  gc$clear_observations()
  for (k in seq_len(gc$n)) {
    expect_true(is.null(gc$graphs[[k]]$.__enclos_env__$private$data))
  }
})

## ── coordinates ──────────────────────────────────────────────────────────────

test_that("coordinates(PtE) maps (component, edge, dist) back to XY", {
  gc <- make_two_components()
  PtE_in <- rbind(c(1, 1, 0.5), c(2, 2, 0.5))
  XY_out <- gc$coordinates(PtE = PtE_in, normalized = TRUE)
  expect_equal(XY_out, rbind(c(0.5, 0), c(11.5, 10)))
})

test_that("coordinates(XY) returns (component, edge, dist)", {
  gc <- make_two_components()
  XY <- rbind(c(0.5, 0), c(11.5, 10))
  PtE_out <- gc$coordinates(XY = XY)
  expect_equal(PtE_out[, 1], c(1, 2))
  expect_equal(PtE_out[, 3], c(0.5, 0.5))
})

test_that("coordinates round-trips PtE -> XY -> PtE", {
  gc <- make_two_components()
  PtE_in <- rbind(c(1, 1, 0.3), c(1, 2, 0.7), c(2, 1, 0.4))
  XY <- gc$coordinates(PtE = PtE_in, normalized = TRUE)
  PtE_back <- gc$coordinates(XY = XY)
  expect_equal(PtE_back[, 1], PtE_in[, 1])
  expect_equal(PtE_back[, 2], PtE_in[, 2])
  expect_equal(PtE_back[, 3], PtE_in[, 3], tolerance = 1e-8)
})

test_that("coordinates errors on invalid input", {
  gc <- make_two_components()
  expect_error(gc$coordinates(), "PtE or XY must be provided")
  expect_error(gc$coordinates(PtE = c(1, 1)), "length 3")
  expect_error(gc$coordinates(PtE = matrix(1:4, 2)),
               "three columns")
  expect_error(gc$coordinates(PtE = matrix(c(99, 1, 0.5), 1)),
               "between 1 and")
})

## ── build_mesh / compute_fem ────────────────────────────────────────────────

test_that("build_mesh constructs a mesh on every component", {
  gc <- make_two_components()
  gc$build_mesh(h = 0.25)
  for (k in seq_len(gc$n)) {
    expect_false(is.null(gc$graphs[[k]]$mesh))
    expect_true(nrow(gc$graphs[[k]]$mesh$V) > 0L)
  }
})

test_that("compute_fem populates FEM matrices on every component", {
  gc <- make_two_components()
  gc$build_mesh(h = 0.25)
  gc$compute_fem()
  for (k in seq_len(gc$n)) {
    expect_false(is.null(gc$graphs[[k]]$mesh$C))
    expect_false(is.null(gc$graphs[[k]]$mesh$G))
    expect_equal(nrow(gc$graphs[[k]]$mesh$C),
                 nrow(gc$graphs[[k]]$mesh$V))
  }
})

## ── plot / plot_function ─────────────────────────────────────────────────────

test_that("plot returns a ggplot object for graph-only display", {
  gc <- make_two_components()
  p <- gc$plot()
  expect_s3_class(p, "ggplot")
})

test_that("plot with data argument works and skips empty components", {
  gc <- make_two_components()
  gc$add_observations(
    data = data.frame(coord_x = 0.5, coord_y = 0, y = 5),
    data_coords = "spatial", verbose = 0)
  # Only comp 1 has data; comp 2 must be skipped without error.
  p <- gc$plot(data = "y")
  expect_s3_class(p, "ggplot")
})

test_that("plot_function works using stored data", {
  gc <- make_two_components()
  df <- data.frame(
    coord_x = c(0.2, 0.8, 1.5, 10.2, 10.8, 11.5),
    coord_y = c(0, 0, 0, 10, 10, 10),
    y = c(1, 2, 3, 4, 5, 6)
  )
  gc$add_observations(data = df, data_coords = "spatial", verbose = 0)
  p <- gc$plot_function(data = "y")
  expect_s3_class(p, "ggplot")
})

test_that("plot_function works with newdata containing .component", {
  gc <- make_two_components()
  nd <- data.frame(
    .component = c(1, 1, 2, 2),
    .edge_number = c(1, 2, 1, 2),
    .distance_on_edge = c(0.5, 0.5, 0.5, 0.5),
    y = c(1, 2, 3, 4)
  )
  class(nd) <- c("metric_graph_data", class(nd))
  p <- gc$plot_function(data = "y", newdata = nd)
  expect_s3_class(p, "ggplot")
})

test_that("plot_function errors when newdata lacks .component", {
  gc <- make_two_components()
  nd <- data.frame(
    .edge_number = c(1, 2),
    .distance_on_edge = c(0.5, 0.5),
    y = c(1, 2)
  )
  class(nd) <- c("metric_graph_data", class(nd))
  expect_error(gc$plot_function(data = "y", newdata = nd),
               "\\.component")
})

## ── three-component sanity checks ────────────────────────────────────────────

test_that("three-component routing puts every point in its expected component", {
  gc <- make_three_components()
  expect_equal(gc$n, 3L)
  XY <- rbind(c(2, 0),    # comp 1 (longest)
              c(21, 0),   # comp 2
              c(40.5, 0)) # comp 3
  expect_equal(gc$which_component(XY), c(1L, 2L, 3L))

  df <- data.frame(coord_x = XY[, 1], coord_y = XY[, 2], y = c(1, 2, 3))
  gc$add_observations(data = df, data_coords = "spatial", verbose = 0)
  out <- gc$get_data()
  expect_setequal(out$.component, c(1L, 2L, 3L))
  for (k in 1:3) {
    expect_equal(out$y[out$.component == k], k)
  }
})
