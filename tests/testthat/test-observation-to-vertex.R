test_that("observation_to_vertex assigns every observation to a vertex", {
  g <- metric_graph$new(edges = make_grid_edges(6), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  set.seed(1)
  N <- 200
  obs <- data.frame(
    y = rnorm(N),
    edge_number = sample.int(g$nE, N, replace = TRUE),
    distance_on_edge = runif(N, 0.01, 0.99)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)

  # Every observation should now have a vertex assignment
  expect_equal(length(g$PtV), N)
  expect_true(all(!is.na(g$PtV)))
  expect_true(all(g$PtV >= 1L & g$PtV <= g$nV))
})

test_that("observation_to_vertex preserves total edge length", {
  g <- metric_graph$new(edges = make_grid_edges(6), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  set.seed(2)
  N <- 150
  obs <- data.frame(
    y = rnorm(N),
    edge_number = sample.int(g$nE, N, replace = TRUE),
    distance_on_edge = runif(N, 0.01, 0.99)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  total_before <- sum(g$edge_lengths)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  total_after <- sum(g$edge_lengths)
  expect_equal(total_after, total_before, tolerance = 1e-10)
})

test_that("observation_to_vertex handles observations exactly at edge endpoints", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  nE_before <- g$nE
  nV_before <- g$nV
  obs <- data.frame(
    y = c(1, 2, 3, 4),
    edge_number = c(1L, 1L, 2L, 2L),
    distance_on_edge = c(0, 1, 0, 1)  # all at endpoints
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  # No edges should have been split since all obs were at existing vertices
  expect_equal(g$nE, nE_before)
  expect_equal(g$nV, nV_before)
  expect_true(all(!is.na(g$PtV)))
})

test_that("observation_to_vertex splits one edge per interior observation", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  nE_before <- g$nE
  nV_before <- g$nV
  # Place exactly one observation strictly inside each of 3 distinct edges
  obs <- data.frame(
    y = c(1, 2, 3),
    edge_number = c(1L, 2L, 3L),
    distance_on_edge = c(0.5, 0.5, 0.5)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  # Each interior split adds one vertex and one edge
  expect_equal(g$nV, nV_before + 3L)
  expect_equal(g$nE, nE_before + 3L)
})

test_that("observation_to_vertex handles many observations on the same edge", {
  g <- metric_graph$new(edges = make_grid_edges(3), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  nE_before <- g$nE
  nV_before <- g$nV
  N_on_one <- 20L
  obs <- data.frame(
    y = rnorm(N_on_one),
    edge_number = rep(1L, N_on_one),
    distance_on_edge = sort(runif(N_on_one, 0.01, 0.99))
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  total_before <- sum(g$edge_lengths)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  # One edge split N_on_one times adds N_on_one new vertices and edges
  expect_equal(g$nV, nV_before + N_on_one)
  expect_equal(g$nE, nE_before + N_on_one)
  expect_equal(sum(g$edge_lengths), total_before, tolerance = 1e-10)
  expect_true(all(all_pte_valid(g)))
})

test_that("observation_to_vertex preserves vector edge_weights", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  # Set distinctive vector weights
  set.seed(3)
  w <- runif(g$nE)
  g$set_edge_weights(weights = w, verbose = 0)
  nE_before <- g$nE

  set.seed(4)
  obs <- data.frame(
    y = rnorm(30),
    edge_number = sample.int(nE_before, 30, replace = TRUE),
    distance_on_edge = runif(30, 0.1, 0.9)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)

  new_weights <- g$.__enclos_env__$private$edge_weights
  # Every weight should still be in the original set (new edges inherit parent weight)
  expect_true(all(new_weights %in% w))
  expect_equal(length(new_weights), g$nE)
})

test_that("observation_to_vertex preserves data.frame edge_weights with multiple columns", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  set.seed(5)
  df_w <- data.frame(
    w1 = runif(g$nE),
    w2 = rnorm(g$nE),
    w3 = seq_len(g$nE) * 0.5
  )
  g$set_edge_weights(weights = df_w, verbose = 0)

  # set_edge_weights appends a ".weights" column for the default Kirchhoff
  # weight, so the stored data.frame has one more column than the input.
  stored_cols <- names(g$.__enclos_env__$private$edge_weights)
  expect_true(all(c("w1", "w2", "w3") %in% stored_cols))

  set.seed(6)
  obs <- data.frame(
    y = rnorm(30),
    edge_number = sample.int(g$nE, 30, replace = TRUE),
    distance_on_edge = runif(30, 0.1, 0.9)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)

  new_df <- g$.__enclos_env__$private$edge_weights
  expect_s3_class(new_df, "data.frame")
  expect_equal(nrow(new_df), g$nE)
  # Column structure must be preserved through OTV
  expect_equal(names(new_df), stored_cols)
  expect_true(all(c("w1", "w2", "w3") %in% names(new_df)))
})

test_that("observation_to_vertex output is deterministic across repeated runs", {
  make_and_run <- function() {
    g <- metric_graph$new(edges = make_grid_edges(5), perform_merges = TRUE,
                           check_connected = FALSE, verbose = 0)
    set.seed(42)
    obs <- data.frame(
      y = rnorm(100),
      edge_number = sample.int(g$nE, 100, replace = TRUE),
      distance_on_edge = runif(100, 0.01, 0.99)
    )
    g$add_observations(data = obs, edge_number = "edge_number",
                        distance_on_edge = "distance_on_edge",
                        normalized = TRUE, verbose = 0)
    g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
    list(nV = g$nV, nE = g$nE,
         lens_sorted = sort(g$edge_lengths),
         PtV_sorted = sort(g$PtV))
  }
  r1 <- make_and_run()
  r2 <- make_and_run()
  expect_equal(r1$nV, r2$nV)
  expect_equal(r1$nE, r2$nE)
  expect_equal(r1$lens_sorted, r2$lens_sorted)
  expect_equal(r1$PtV_sorted, r2$PtV_sorted)
})

test_that("observation_to_vertex preserves PtE validity on all new edges", {
  g <- metric_graph$new(edges = make_grid_edges(5), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  set.seed(7)
  obs <- data.frame(
    y = rnorm(80),
    edge_number = sample.int(g$nE, 80, replace = TRUE),
    distance_on_edge = runif(80, 0.01, 0.99)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  expect_true(all(all_pte_valid(g)))
})

test_that("observation_to_vertex on lon/lat conserves edge length exactly", {
  skip_if_not_installed("sf")
  g <- metric_graph$new(edges = make_longlat_grid(5), perform_merges = TRUE,
                         tolerance = list(vertex_vertex = 1e-5,
                                           vertex_edge = 0,
                                           edge_edge = 0),
                         check_connected = FALSE, verbose = 0)
  set.seed(8)
  obs <- data.frame(
    y = rnorm(60),
    edge_number = sample.int(g$nE, 60, replace = TRUE),
    distance_on_edge = runif(60, 0.01, 0.99)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  total_before <- sum(g$edge_lengths)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  total_after <- sum(g$edge_lengths)
  expect_equal(total_after, total_before, tolerance = 1e-10)
  expect_true(all(all_pte_valid(g)))
  expect_true(all(!is.na(g$PtV)))
})

test_that("observation_to_vertex called twice in a row is idempotent on PtV and graph state", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  set.seed(9)
  obs <- data.frame(
    y = rnorm(40),
    edge_number = sample.int(g$nE, 40, replace = TRUE),
    distance_on_edge = runif(40, 0.05, 0.95)
  )
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  nV1 <- g$nV
  nE1 <- g$nE
  lens1 <- g$edge_lengths
  PtV1 <- g$PtV
  # All observations are now at vertices, so a second call should be a no-op
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  expect_equal(g$nV, nV1)
  expect_equal(g$nE, nE1)
  expect_equal(g$edge_lengths, lens1, tolerance = 1e-12)
  expect_equal(sort(g$PtV), sort(PtV1))
})

test_that("observation_to_vertex on data with multiple groups preserves per-group structure", {
  g <- metric_graph$new(edges = make_grid_edges(4), perform_merges = TRUE,
                         check_connected = FALSE, verbose = 0)
  set.seed(10)
  n_per_group <- 20
  # Two groups with the same underlying locations (same edge/distance) but
  # different y values
  base <- data.frame(
    edge_number = sample.int(g$nE, n_per_group, replace = TRUE),
    distance_on_edge = runif(n_per_group, 0.1, 0.9)
  )
  grp1 <- cbind(y = rnorm(n_per_group), base, replicate = "a",
                stringsAsFactors = FALSE)
  grp2 <- cbind(y = rnorm(n_per_group), base, replicate = "b",
                stringsAsFactors = FALSE)
  obs <- rbind(grp1, grp2)
  g$add_observations(data = obs, edge_number = "edge_number",
                      distance_on_edge = "distance_on_edge",
                      group = "replicate",
                      normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE, verbose = 0)
  # Both groups should have their observations attached
  priv_data <- g$.__enclos_env__$private$data
  expect_equal(length(priv_data[[".group"]]), 2 * n_per_group)
  expect_equal(length(unique(priv_data[[".group"]])), 2L)
  # All PtV entries are valid
  expect_true(all(g$PtV >= 1L & g$PtV <= g$nV))
})
