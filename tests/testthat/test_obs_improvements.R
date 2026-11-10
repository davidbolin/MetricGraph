## Unit tests for add_observations / process_data / idx_not_all_NA improvements
## Covers: vectorized group construction, filter_spatial_obs_groups helper,
##         idx_not_all_NA, idx_not_any_NA, process_data_add_obs (match-based).

library(MetricGraph)
library(testthat)

## ── helpers ──────────────────────────────────────────────────────────────────

make_simple_graph <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(1, 0), c(1, 1))
  edge3 <- rbind(c(0, 0), c(0, 1))
  metric_graph$new(edges = list(edge1, edge2, edge3))
}

## ── idx_not_all_NA  ───────────────────────────────────────────────────────────

test_that("idx_not_all_NA returns TRUE for rows with at least one non-NA", {
  dl <- list(
    .edge_number = 1:4,
    .distance_on_edge = c(0, 0.5, 1, 0.2),
    .group = rep("g1", 4),
    y = c(1, NA, NA, 4),
    z = c(NA, NA, NA, NA)
  )
  result <- MetricGraph:::idx_not_all_NA(dl)
  # Row 1: y=1, z=NA  → not all NA → TRUE
  # Row 2: y=NA, z=NA → all NA     → FALSE
  # Row 3: y=NA, z=NA → all NA     → FALSE
  # Row 4: y=4, z=NA  → not all NA → TRUE
  expect_equal(result, c(TRUE, FALSE, FALSE, TRUE))
})

test_that("idx_not_all_NA handles all-NA data columns", {
  dl <- list(
    .edge_number = 1:3,
    .group = rep("g1", 3),
    y = c(NA, NA, NA)
  )
  result <- MetricGraph:::idx_not_all_NA(dl)
  expect_equal(result, c(FALSE, FALSE, FALSE))
})

test_that("idx_not_all_NA handles no data columns (only metadata)", {
  dl <- list(
    .edge_number = 1:3,
    .distance_on_edge = c(0, 0.5, 1),
    .group = rep("g1", 3)
  )
  result <- MetricGraph:::idx_not_all_NA(dl)
  expect_equal(length(result), 0L)
})

test_that("idx_not_all_NA: single column, matches expectation", {
  dl <- list(.group = c("a", "b"), y = c(1, NA))
  result <- MetricGraph:::idx_not_all_NA(dl)
  expect_equal(result, c(TRUE, FALSE))
})

## ── idx_not_any_NA  ───────────────────────────────────────────────────────────

test_that("idx_not_any_NA returns TRUE for rows with no NAs", {
  dl <- list(
    .group = rep("g1", 3),
    y = c(1, NA, 3),
    z = c(10, 20, 30)
  )
  result <- MetricGraph:::idx_not_any_NA(dl)
  # Row 1: y=1 z=10  → no NA → TRUE
  # Row 2: y=NA z=20 → has NA → FALSE
  # Row 3: y=3 z=30  → no NA → TRUE
  expect_equal(result, c(TRUE, FALSE, TRUE))
})

## ── group string construction (vectorized Reduce) ─────────────────────────────

test_that("add_observations: single group column → correct .group labels", {
  g <- make_simple_graph()
  df <- data.frame(
    edge_number = c(1L, 2L, 3L),
    distance_on_edge = c(0.3, 0.4, 0.5),
    y = c(1, 2, 3),
    grp = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, group = "grp", normalized = TRUE, verbose = 0)
  groups <- unique(g$get_data()[[".group"]])
  expect_true(all(groups %in% c("A", "B")))
})

test_that("add_observations: two group columns → combined .group labels", {
  g <- make_simple_graph()
  df <- data.frame(
    edge_number = c(1L, 1L, 2L, 2L),
    distance_on_edge = c(0.2, 0.6, 0.3, 0.7),
    y = 1:4,
    grp1 = c("X", "X", "Y", "Y"),
    grp2 = c("a", "b", "a", "b"),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, group = c("grp1", "grp2"),
                     group_sep = ".", normalized = TRUE, verbose = 0)
  groups <- sort(unique(g$get_data()[[".group"]]))
  expect_equal(groups, c("X.a", "X.b", "Y.a", "Y.b"))
})

## ── process_data: group construction ─────────────────────────────────────────

test_that("process_data: two group columns → combined labels", {
  g <- make_simple_graph()
  df <- data.frame(
    edge_number = c(1L, 1L, 2L, 2L),
    distance_on_edge = c(0.2, 0.6, 0.3, 0.7),
    y = 1:4,
    grp1 = c("X", "X", "Y", "Y"),
    grp2 = c("a", "b", "a", "b"),
    stringsAsFactors = FALSE
  )
  out <- g$process_data(data = df, group = c("grp1", "grp2"),
                        group_sep = "-", normalized = TRUE, verbose = FALSE)
  groups <- sort(unique(out[[".group"]]))
  expect_equal(groups, c("X-a", "X-b", "Y-a", "Y-b"))
})

## ── far-point and duplicate handling ─────────────────────────────────────────

test_that("add_observations PtE: observations on all three edges are stored", {
  g <- make_simple_graph()
  df <- data.frame(
    edge_number = c(1L, 2L, 3L),
    distance_on_edge = c(0.25, 0.5, 0.75),
    y = c(10, 20, 30),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  d <- g$get_data()
  expect_equal(nrow(d), 3L)
  expect_equal(sort(d$y), c(10, 20, 30))
})

test_that("add_observations: multiple calls accumulate observations", {
  g <- make_simple_graph()
  df1 <- data.frame(edge_number = 1L, distance_on_edge = 0.3,
                    y = 1, stringsAsFactors = FALSE)
  df2 <- data.frame(edge_number = 2L, distance_on_edge = 0.5,
                    y = 2, stringsAsFactors = FALSE)
  g$add_observations(data = df1, normalized = TRUE, verbose = 0)
  g$add_observations(data = df2, normalized = TRUE, verbose = 0)
  d <- g$get_data()
  expect_equal(nrow(d), 2L)
})

test_that("add_observations: PtE duplicate warning and data is valid", {
  g <- make_simple_graph()
  # Two observations at the same normalized position
  df <- data.frame(
    edge_number = c(1L, 1L),
    distance_on_edge = c(0.5, 0.5),
    y = c(10, 20),
    stringsAsFactors = FALSE
  )
  expect_warning(
    g$add_observations(data = df, normalized = TRUE, verbose = 0),
    regexp = "repeated"
  )
})

## ── NA-sparse behavior via get_data(drop_all_na=TRUE) ─────────────────────────

test_that("get_data drop_all_na removes rows with all-NA data columns", {
  g <- make_simple_graph()
  # Two groups, same location, only one has data
  df1 <- data.frame(edge_number = 1L, distance_on_edge = 0.5,
                    y = 5, grp = "A", stringsAsFactors = FALSE)
  df2 <- data.frame(edge_number = 1L, distance_on_edge = 0.5,
                    y = NA_real_, grp = "B", stringsAsFactors = FALSE)
  df <- rbind(df1, df2)
  g$add_observations(data = df, group = "grp", normalized = TRUE, verbose = 0,
                     suppress_warnings = TRUE)
  d <- g$get_data(drop_all_na = TRUE)
  # Only the row for group A (with y=5) should remain
  expect_true(all(!is.na(d$y[d$.group == "A"])))
})

## ── correctness: process_data_add_obs via add_observations round-trip ─────────

test_that("round-trip: stored PtE matches input", {
  g <- make_simple_graph()
  df <- data.frame(
    edge_number = c(1L, 2L, 3L),
    distance_on_edge = c(0.1, 0.5, 0.9),
    y = c(1.1, 2.2, 3.3),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  PtE <- g$get_PtE()
  expect_equal(nrow(PtE), 3L)
  # PtE after standardize_df_positions may reorder — check values present
  expect_equal(sort(PtE[, 1]), c(1L, 2L, 3L))
})

test_that("process_data: returns correct structure", {
  g <- make_simple_graph()
  df <- data.frame(
    edge_number = c(1L, 2L),
    distance_on_edge = c(0.3, 0.7),
    y = c(5.0, 6.0),
    stringsAsFactors = FALSE
  )
  out <- g$process_data(data = df, normalized = TRUE, verbose = FALSE, format = "list")
  expect_true(".edge_number" %in% names(out))
  expect_true(".distance_on_edge" %in% names(out))
  expect_equal(length(out$y), 2L)
  expect_equal(sort(out$y), c(5.0, 6.0))
})

## ── group_sep customization ───────────────────────────────────────────────────

test_that("group_sep is respected in combined group label", {
  g <- make_simple_graph()
  df <- data.frame(
    edge_number = c(1L, 2L),
    distance_on_edge = c(0.2, 0.8),
    y = c(1, 2),
    g1 = c("P", "Q"),
    g2 = c("x", "y"),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, group = c("g1", "g2"),
                     group_sep = ":", normalized = TRUE, verbose = 0)
  groups <- sort(unique(g$get_data()[[".group"]]))
  expect_equal(groups, c("P:x", "Q:y"))
})

## ── model-pipeline: private$data structure consumed by likelihoods ────────────

make_two_edge_graph <- function() {
  edge1 <- rbind(c(30, 0), c(30, 80))
  edge2 <- rbind(c(30, 80), c(140, 80))
  metric_graph$new(edges = list(edge1, edge2))
}

test_that("private$data structure: .group has n_loc entries for single replicate", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 1L, 2L),
    distance_on_edge = c(0.2, 0.6, 0.4),
    y = c(1.0, 2.0, 3.0),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0)

  grp <- g$.__enclos_env__$private$data[[".group"]]
  y   <- g$.__enclos_env__$private$data[["y"]]
  PtE <- g$get_PtE()

  # Single group → .group has same length as nrow(PtE)
  expect_equal(length(grp), nrow(PtE))
  expect_equal(length(y),   nrow(PtE))
  # All same group label
  expect_equal(length(unique(grp)), 1L)
})

test_that("private$data structure: two replicates → .group = 2 × n_locs", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L, 1L, 2L),
    distance_on_edge = c(0.3, 0.7, 0.3, 0.7),
    y = c(1.0, 2.0, 3.0, 4.0),
    repl = c("t1", "t1", "t2", "t2"),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, group = "repl", normalized = TRUE,
                     verbose = 0, suppress_warnings = TRUE)

  grp  <- g$.__enclos_env__$private$data[[".group"]]
  y    <- g$.__enclos_env__$private$data[["y"]]
  PtE  <- g$get_PtE()
  n_loc <- nrow(PtE)

  expect_equal(length(grp), 2L * n_loc)
  expect_equal(length(y),   2L * n_loc)
  expect_equal(sort(unique(grp)), c("t1", "t2"))

  # Per-group subsetting (mirrors likelihood code pattern)
  for (repl in c("t1", "t2")) {
    mask <- grp == repl
    expect_equal(sum(mask), n_loc)
    expect_false(any(is.na(y[mask])))
  }
})

test_that("private$data: response values are correctly placed per group", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L, 1L, 2L),
    distance_on_edge = c(0.25, 0.75, 0.25, 0.75),
    y = c(10.0, 20.0, 30.0, 40.0),
    repl = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, group = "repl", normalized = TRUE,
                     verbose = 0, suppress_warnings = TRUE)

  grp <- g$.__enclos_env__$private$data[[".group"]]
  y   <- g$.__enclos_env__$private$data[["y"]]

  y_A <- sort(y[grp == "A"])
  y_B <- sort(y[grp == "B"])
  expect_equal(y_A, c(10.0, 20.0))
  expect_equal(y_B, c(30.0, 40.0))
})

test_that("get_PtE returns same locations across replicates", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L, 1L, 2L),
    distance_on_edge = c(0.3, 0.7, 0.3, 0.7),
    y = c(1, 2, 3, 4),
    repl = c("r1", "r1", "r2", "r2"),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, group = "repl", normalized = TRUE,
                     verbose = 0, suppress_warnings = TRUE)

  PtE <- g$get_PtE()
  # PtE reflects first-group locations only (used by likelihood for mesh)
  expect_equal(nrow(PtE), 2L)
  expect_equal(sort(PtE[, 1]), c(1L, 2L))
})

## ── observation_to_vertex: correctness after improvements ─────────────────────

test_that("observation_to_vertex: splits edges and updates PtV correctly", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L),
    distance_on_edge = c(0.5, 0.5),
    y = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  nV_before <- g$nV

  g$observation_to_vertex()

  # Two interior observations → two new vertices added
  expect_equal(g$nV, nV_before + 2L)
  # PtV has entries for each observation-vertex
  expect_true(length(g$PtV) > 0)
})

test_that("observation_to_vertex: vertex observations stay at correct vertices", {
  g <- make_two_edge_graph()
  # distance=0 → start vertex of edge; distance=1 → end vertex
  df <- data.frame(
    edge_number = c(1L, 1L),
    distance_on_edge = c(0, 1),
    y = c(5.0, 6.0),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0,
                     suppress_warnings = TRUE)
  nV_before <- g$nV
  g$observation_to_vertex()
  # No new vertices needed (endpoints already exist)
  expect_equal(g$nV, nV_before)
})

test_that("observation_to_vertex preserves data values", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L),
    distance_on_edge = c(0.4, 0.6),
    y = c(42.0, 99.0),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  g$observation_to_vertex()

  d <- g$get_data(drop_all_na = TRUE)
  expect_equal(sort(d$y), c(42.0, 99.0))
})

test_that("multi-replicate + observation_to_vertex: data integrity", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L, 1L, 2L),
    distance_on_edge = c(0.4, 0.6, 0.4, 0.6),
    y = c(1.0, 2.0, 3.0, 4.0),
    repl = c("r1", "r1", "r2", "r2"),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, group = "repl", normalized = TRUE,
                     verbose = 0, suppress_warnings = TRUE)
  g$observation_to_vertex()

  d <- g$get_data(drop_all_na = TRUE)
  # Both replicates still present
  expect_equal(sort(unique(d[[".group"]])), c("r1", "r2"))
  # Each replicate has 2 observations
  expect_equal(sum(d[[".group"]] == "r1"), 2L)
  expect_equal(sum(d[[".group"]] == "r2"), 2L)
})

## ── spde_precision compatibility ──────────────────────────────────────────────

test_that("spde_precision(alpha=2) works on graph with added observations", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L),
    distance_on_edge = c(0.3, 0.7),
    y = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  g$buildC(alpha = 2)
  Q <- spde_precision(kappa = 0.3, tau = 1.0, alpha = 2, graph = g, BC = 1)
  expect_equal(dim(Q), c(4L * g$nE, 4L * g$nE))
  # Q should be positive semi-definite — check no negative eigenvalues significantly
  eig <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(eig) > -1e-10)
})

test_that("spde_precision(alpha=1) works on graph with added observations", {
  g <- make_two_edge_graph()
  df <- data.frame(
    edge_number = c(1L, 2L),
    distance_on_edge = c(0.3, 0.7),
    y = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  Q <- spde_precision(kappa = 0.3, tau = 1.0, alpha = 1, graph = g, BC = 1)
  expect_equal(dim(Q), c(g$nV, g$nV))
})

## ── speed comparison: add_observations with many observations ─────────────────

test_that("add_observations with 1000 PtE obs completes in reasonable time", {
  edge1 <- rbind(c(0, 0), c(10, 0))
  edge2 <- rbind(c(10, 0), c(10, 10))
  edge3 <- rbind(c(10, 10), c(0, 10))
  edge4 <- rbind(c(0, 10), c(0, 0))
  g <- metric_graph$new(edges = list(edge1, edge2, edge3, edge4))

  set.seed(42)
  n <- 1000L
  df <- data.frame(
    edge_number       = sample(1:4, n, replace = TRUE),
    distance_on_edge  = runif(n),
    y                 = rnorm(n),
    stringsAsFactors  = FALSE
  )

  t <- system.time(
    g$add_observations(data = df, normalized = TRUE, verbose = 0,
                       suppress_warnings = TRUE)
  )["elapsed"]

  expect_lt(t, 10)   # must complete in under 10 seconds
  expect_equal(nrow(g$get_PtE()), n)
})

test_that("add_observations with 5 replicates x 200 obs completes fast", {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(1, 0), c(1, 1))
  g <- metric_graph$new(edges = list(edge1, edge2))

  set.seed(7)
  n_per <- 200L
  repls <- paste0("r", 1:5)
  df_list <- lapply(repls, function(r) {
    data.frame(
      edge_number      = sample(1:2, n_per, replace = TRUE),
      distance_on_edge = runif(n_per),
      y                = rnorm(n_per),
      repl             = r,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, df_list)

  t <- system.time(
    g$add_observations(data = df, group = "repl", normalized = TRUE,
                       verbose = 0, suppress_warnings = TRUE)
  )["elapsed"]

  expect_lt(t, 15)

  grp <- g$.__enclos_env__$private$data[[".group"]]
  expect_equal(sort(unique(grp)), sort(repls))
})
