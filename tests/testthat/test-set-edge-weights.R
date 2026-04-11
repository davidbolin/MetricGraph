test_that("set_edge_weights accepts a numeric vector and attaches to each edge", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  w <- seq_len(g$nE) * 0.1
  g$set_edge_weights(weights = w, verbose = 0)
  # Each edge should have a "weight" attribute equal to its vector entry
  for (i in seq_len(g$nE)) {
    expect_equal(attr(g$edges[[i]], "weight"), w[i],
                 info = sprintf("edge %d weight mismatch", i))
  }
})

test_that("set_edge_weights accepts a data.frame and attaches to each edge", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  df_w <- data.frame(w1 = seq_len(g$nE), w2 = runif(g$nE))
  g$set_edge_weights(weights = df_w, verbose = 0)
  for (i in seq_len(g$nE)) {
    w_attr <- attr(g$edges[[i]], "weight")
    expect_s3_class(w_attr, "data.frame")
    expect_equal(w_attr$w1, as.numeric(df_w$w1[i]))
    expect_equal(w_attr$w2, df_w$w2[i])
  }
})

test_that("set_edge_weights preserves id, length, PtE, longlat attributes", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  # Call set_edge_weights a second time with new weights
  g$set_edge_weights(weights = rep(1, g$nE), verbose = 0)
  for (i in seq_len(g$nE)) {
    e <- g$edges[[i]]
    expect_equal(attr(e, "id"), i)
    expect_equal(attr(e, "length"), g$edge_lengths[i], tolerance = 1e-12)
    pte <- attr(e, "PtE")
    expect_equal(pte[1], 0, tolerance = 1e-12)
    expect_equal(pte[length(pte)], 1, tolerance = 1e-12)
    expect_false(is.null(attr(e, "longlat")))
  }
})

test_that("set_edge_weights re-assigns weights correctly after repeated calls", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  w1 <- rep(0.5, g$nE)
  w2 <- seq_len(g$nE) * 1.0
  g$set_edge_weights(weights = w1, verbose = 0)
  for (i in seq_len(g$nE)) {
    expect_equal(attr(g$edges[[i]], "weight"), w1[i])
  }
  g$set_edge_weights(weights = w2, verbose = 0)
  for (i in seq_len(g$nE)) {
    expect_equal(attr(g$edges[[i]], "weight"), w2[i])
  }
})

test_that("set_edge_weights tags every edge with class 'metric_graph_edge'", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  g$set_edge_weights(weights = rep(1, g$nE), verbose = 0)
  for (i in seq_len(g$nE)) {
    expect_true(inherits(g$edges[[i]], "metric_graph_edge"))
  }
  expect_true(inherits(g$edges, "metric_graph_edges"))
})
