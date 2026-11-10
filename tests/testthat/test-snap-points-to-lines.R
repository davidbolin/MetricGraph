# Regression tests for the indexed nearest-edge search that replaced the dense
# point-by-edge distance matrix in `snapPointsToLines()`.

# Reference implementation: the original dense `st_distance()` search, kept here
# so that the indexed version can be checked against it.
snap_reference <- function(points, lines) {
  if (!is.list(lines)) lines <- list(lines)
  points <- as.matrix(points)[, 1:2, drop = FALSE]
  points_sf <- sf::st_as_sf(as.data.frame(points), coords = 1:2)
  lines_sf <- sf::st_sfc(lapply(lines, function(i) sf::st_linestring(i)))
  d <- t(sf::st_distance(points_sf, lines_sf))
  idx <- apply(d, 2, which.min)
  coords <- vapply(seq_len(nrow(points)),
                   function(x) MetricGraph:::nearestPointOnLine(lines[[idx[x]]],
                                                                points[x, ]),
                   FUN.VALUE = c(0, 0))
  list(index = as.integer(idx), coords = coords,
       dist = as.numeric(apply(d, 2, min)))
}

expect_matches_reference <- function(points, lines) {
  ref <- snap_reference(points, lines)
  got <- MetricGraph:::snapPointsToLines(points, lines,
                                         longlat = FALSE, crs = NULL)
  expect_equal(as.integer(got$df$nearest_line_index), ref$index)
  expect_equal(unname(as.matrix(got$coords)), unname(ref$coords))
  expect_equal(as.numeric(got$df$snap_dist), ref$dist, tolerance = 1e-10)
}

make_lines <- function(n = 20) {
  set.seed(7)
  lapply(seq_len(n), function(k) {
    x0 <- runif(1, 0, 10); y0 <- runif(1, 0, 10)
    m <- cbind(x0 + cumsum(rnorm(4, sd = 0.3)), y0 + cumsum(rnorm(4, sd = 0.3)))
    rbind(c(x0, y0), m)
  })
}

test_that("indexed snapping reproduces the dense search on random points", {
  lines <- make_lines(60)
  set.seed(11)
  pts <- cbind(runif(200, -1, 11), runif(200, -1, 11))
  expect_matches_reference(pts, lines)
})

test_that("indexed snapping reproduces the dense search for points on the lines", {
  lines <- make_lines(60)
  all_pts <- do.call(rbind, lines)
  # Vertices of the lines, where several edges are exactly equidistant.
  expect_matches_reference(all_pts, lines)
  # Segment midpoints, i.e. points in the interior of the edges.
  mids <- (all_pts[-nrow(all_pts), , drop = FALSE] +
             all_pts[-1, , drop = FALSE]) / 2
  expect_matches_reference(mids, lines)
})

test_that("indexed snapping handles ties towards the smallest edge index", {
  # Three identical edges: every point is exactly equidistant from all of them.
  e <- rbind(c(0, 0), c(1, 0))
  lines <- list(e, e, e)
  pts <- cbind(c(0.25, 0.5, 2), c(1, -1, 0))
  got <- MetricGraph:::snapPointsToLines(pts, lines, longlat = FALSE, crs = NULL)
  expect_equal(as.integer(got$df$nearest_line_index), rep(1L, 3))
  expect_matches_reference(pts, lines)
})

test_that("indexed snapping handles points far outside the bounding box", {
  lines <- make_lines(60)
  set.seed(13)
  pts <- cbind(runif(50, 500, 600), runif(50, -800, -700))
  expect_matches_reference(pts, lines)
})

test_that("indexed snapping handles a single line and a single point", {
  e <- rbind(c(0, 0), c(1, 1), c(2, 0))
  expect_matches_reference(matrix(c(0.5, 1), 1, 2), e)
  expect_matches_reference(cbind(seq(-1, 3, by = 0.25), 0.4), e)
})

test_that("indexed snapping handles degenerate edges", {
  # A point-like edge (both endpoints identical) must not break the search and
  # must not steal points from the real edge unless it is genuinely closest.
  lines <- list(rbind(c(0, 0), c(0, 0)), rbind(c(0, 1), c(5, 1)))
  pts <- rbind(c(0, 0.4), c(0, 0.6), c(3, 5))
  got <- MetricGraph:::snapPointsToLines(pts, lines, longlat = FALSE, crs = NULL)
  expect_equal(as.integer(got$df$nearest_line_index), c(1L, 2L, 2L))
})

test_that("graph$coordinates() lands on the same edges as the dense search", {
  lines <- make_lines(60)
  g <- metric_graph$new(edges = lines, perform_merges = FALSE,
                        check_connected = FALSE, verbose = 0)
  set.seed(17)
  pts <- do.call(rbind, lines)[sample(240, 120), , drop = FALSE]
  pts <- pts + matrix(rnorm(240, sd = 1e-3), ncol = 2)
  PtE <- g$coordinates(XY = pts, normalized = TRUE)
  ref <- snap_reference(pts, g$edges)
  expect_equal(as.integer(PtE[, 1]), ref$index)
  # Positions must be reproduced by walking back along the edge.
  back <- g$coordinates(PtE = PtE, normalized = TRUE)
  expect_equal(unname(back), unname(t(ref$coords)), tolerance = 1e-8)
})

test_that("indexed snapping matches a brute-force search on random configurations", {
  # Brute force with the same per-segment arithmetic, so that ties resolve the
  # same way and the comparison is exact rather than up to rounding.
  brute <- function(points, lines) {
    idx <- integer(nrow(points))
    for (i in seq_len(nrow(points))) {
      best <- Inf
      for (k in seq_along(lines)) {
        nk <- nrow(lines[[k]])
        if (nk < 2) next
        segs <- vapply(2:nk, function(x)
          MetricGraph:::nearestPointOnSegment(lines[[k]][(x - 1):x, , drop = FALSE],
                                              points[i, ]),
          FUN.VALUE = c(0, 0, 0))
        d <- min(segs[3, ])
        if (d < best) { best <- d; idx[i] <- k }
      }
    }
    idx
  }

  set.seed(2024)
  for (rep in 1:25) {
    nl <- sample(2:12, 1)
    lines <- lapply(seq_len(nl), function(k) {
      type <- sample(1:4, 1)
      if (type == 1) {            # degenerate point-edge
        p <- c(sample(0:3, 1), sample(0:3, 1))
        rbind(p, p)
      } else if (type == 2) {     # axis-aligned, often on the domain border
        y <- sample(0:3, 1)
        rbind(c(0, y), c(3, y))
      } else if (type == 3) {     # vertical
        x <- sample(0:3, 1)
        rbind(c(x, 0), c(x, 3))
      } else {                    # random polyline
        cbind(runif(4, 0, 3), runif(4, 0, 3))
      }
    })
    pts <- rbind(cbind(runif(30, -1, 4), runif(30, -1, 4)),
                 do.call(rbind, lines))
    got <- MetricGraph:::nearest_edge_cpp(lines, pts)
    expect_equal(as.integer(got$index), brute(pts, lines),
                 info = paste("replicate", rep))
  }
})

test_that("vertex-to-edge snapping splits edges at the projected position", {
  # A vertical edge whose lower endpoint sits slightly off a horizontal edge:
  # the constructor must snap it onto the horizontal edge and split it there.
  lines <- list(
    rbind(c(0, 0), c(10, 0)),
    rbind(c(3, 0.02), c(3, 5))
  )
  g <- metric_graph$new(edges = lines, perform_merges = TRUE,
                        tolerance = list(vertex_vertex = 1e-8,
                                         vertex_edge = 0.1,
                                         edge_edge = 0),
                        check_connected = FALSE, verbose = 0)
  # The horizontal edge is split in two at x = 3, the vertical one is extended
  # down to the split point.
  expect_equal(g$nE, 3L)
  expect_true(any(abs(g$V[, 1] - 3) < 1e-12 & abs(g$V[, 2]) < 1e-12))
  expect_equal(sort(as.numeric(g$edge_lengths)), c(3, 5 - 0.02, 7),
               tolerance = 1e-10)
})

test_that("vertex-to-edge snapping is unaffected by the number of edges", {
  # The snapping result for a given pair of edges must not depend on how many
  # unrelated edges are present, which is what the sparse candidate inversion
  # has to guarantee.
  base <- list(rbind(c(0, 0), c(10, 0)), rbind(c(3, 0.02), c(3, 5)))
  filler <- lapply(1:200, function(k) rbind(c(100 + k, 0), c(100 + k, 1)))
  tol <- list(vertex_vertex = 1e-8, vertex_edge = 0.1, edge_edge = 0)
  g1 <- metric_graph$new(edges = base, perform_merges = TRUE,
                         tolerance = tol, check_connected = FALSE, verbose = 0)
  g2 <- metric_graph$new(edges = c(base, filler), perform_merges = TRUE,
                         tolerance = tol, check_connected = FALSE, verbose = 0)
  expect_equal(sort(as.numeric(g1$edge_lengths), decreasing = TRUE),
               sort(as.numeric(g2$edge_lengths), decreasing = TRUE)[1:3])
  expect_equal(g2$nE, 3L + length(filler))
})
