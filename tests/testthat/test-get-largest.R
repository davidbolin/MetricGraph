# `get_largest()` must return exactly the graph `get_components()[[1]]` returns,
# without building the other components.

canon <- function(g) {
  list(nV = g$nV, nE = g$nE, V = g$V,
       E = matrix(as.integer(g$E), ncol = 2),
       el = as.numeric(g$edge_lengths),
       edges = lapply(g$edges, function(x) {
         m <- unclass(x); attributes(m) <- NULL; matrix(m, ncol = 2)
       }))
}

disconnected_edges <- function() {
  # One large piece (a 4x4 grid), one medium (a triangle), one small (an edge).
  big <- make_grid_edges(4)
  mid <- list(rbind(c(50, 50), c(52, 50)),
              rbind(c(52, 50), c(51, 51.5)),
              rbind(c(51, 51.5), c(50, 50)))
  small <- list(rbind(c(100, 100), c(100.3, 100)))
  c(big, mid, small)
}

test_that("get_largest returns the same graph as get_components()[[1]]", {
  g <- metric_graph$new(edges = disconnected_edges(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  expect_true(g$is_disconnected())
  expect_equal(canon(g$get_largest()), canon(g$get_components()[[1]]))
})

test_that("get_largest does not populate the full component cache", {
  g <- metric_graph$new(edges = disconnected_edges(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  lg <- g$get_largest()
  expect_s3_class(lg, "metric_graph")
  # The per-component list must still be unbuilt at this point.
  expect_null(g$.__enclos_env__$private$components_cache)
  # And asking again must be consistent.
  expect_equal(canon(g$get_largest()), canon(lg))
})

test_that("get_largest returns the graph itself when connected", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  expect_false(g$is_disconnected())
  expect_identical(g$get_largest(), g)
})

test_that("get_largest routes observations to the largest component", {
  g <- metric_graph$new(edges = disconnected_edges(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  df <- data.frame(y = c(1, 2, 3),
                   x = c(1.5, 2.5, 50.5),
                   yy = c(1, 1, 50))
  suppressWarnings(g$add_observations(df, coord_x = "x", coord_y = "yy",
                                      data_coords = "spatial", verbose = 0,
                                      suppress_warnings = TRUE))
  lg <- g$get_largest()
  ref <- g$get_components()[[1]]
  expect_equal(canon(lg), canon(ref))
  expect_equal(as.numeric(lg$get_data()[["y"]]),
               as.numeric(ref$get_data()[["y"]]))
})

test_that("get_largest cache is invalidated when observations change", {
  g <- metric_graph$new(edges = disconnected_edges(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  lg1 <- g$get_largest()
  expect_s3_class(lg1, "metric_graph")
  expect_false(is.null(g$.__enclos_env__$private$largest_cache))
  df <- data.frame(y = 1, x = 1.5, yy = 1)
  suppressWarnings(g$add_observations(df, coord_x = "x", coord_y = "yy",
                                      data_coords = "spatial", verbose = 0,
                                      suppress_warnings = TRUE))
  expect_null(g$.__enclos_env__$private$largest_cache)
  lg2 <- g$get_largest()
  expect_equal(nrow(as.data.frame(lg2$get_data())), 1L)
})

test_that("get_largest returns NULL for a graph with no edges", {
  g <- metric_graph$new(edges = make_grid_edges(3), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  g$.__enclos_env__$self$nE <- 0L
  expect_null(g$get_largest())
})
