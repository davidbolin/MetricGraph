test_that("point distances preserve duplicate and shared-endpoint rows", {
  graph <- metric_graph$new(edges = list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(1, 1))
  ), verbose = 0)
  points <- rbind(c(1, 0.5), c(2, 0.5), c(1, 0.5),
                  c(1, 1), c(2, 0))
  D <- graph$compute_geodist_PtE(points, include_vertices = FALSE)
  expect_equal(dim(D), c(5L, 5L))
  expect_equal(D, t(D))
  expect_equal(diag(D), rep(0, 5))
  expect_equal(D[1, 3], 0)
  expect_equal(D[4, 5], 0)
  expect_equal(D[1, 2], 1)
  expect_equal(D[4, 2], 0.5)
  expect_equal(dim(graph$compute_geodist_PtE(points)),
               c(graph$nV + 5L, graph$nV + 5L))
})

test_that("absolute positions and disconnected components are handled", {
  graph <- metric_graph$new(edges = list(
    rbind(c(0, 0), c(2, 0)),
    rbind(c(10, 0), c(11, 0))
  ), check_connected = FALSE, verbose = 0)
  points <- rbind(c(1, 0.5), c(1, 1.5), c(2, 0.5))
  D <- graph$compute_geodist_PtE(points, normalized = FALSE,
                                  include_vertices = FALSE)
  expect_equal(D[1, 2], 1)
  expect_true(is.infinite(D[1, 3]))
  expect_equal(diag(D), rep(0, 3))
})

test_that("reversing an edge and its point positions preserves distances", {
  forward <- metric_graph$new(edges = list(rbind(c(0, 0), c(2, 0))),
                              verbose = 0)
  backward <- metric_graph$new(edges = list(rbind(c(2, 0), c(0, 0))),
                               verbose = 0)
  p <- cbind(1, c(0.2, 0.8))
  D1 <- forward$compute_geodist_PtE(p, include_vertices = FALSE)
  D2 <- backward$compute_geodist_PtE(cbind(1, 1 - p[, 2]),
                                     include_vertices = FALSE)
  expect_equal(D1, D2)
})
