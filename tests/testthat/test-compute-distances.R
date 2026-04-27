# --------------------------------------------------------------------------- #
#  Tests for compute_geodist_PtE, compute_geodist, compute_geodist_mesh,
#            compute_resdist_PtE, compute_resdist
# --------------------------------------------------------------------------- #

# ---- helpers -------------------------------------------------------------- #

# Simple L-shaped graph: edge 1 from (0,0)->(1,0), edge 2 from (1,0)->(1,1)
# Both edges have length 1. Vertex 2 = (1,0) is the shared corner.
#
# `make_L_graph()` builds a fresh instance — used by tests that mutate the
# graph (add observations, build a mesh, call `compute_geodist`/`compute_resdist`,
# which write back to fields on `self`).
#
# `shared_L_graph()` returns a single cached instance — safe for tests that
# only read (`g$nE`, `g$nV`, `g$edge_lengths`) or call `compute_geodist_PtE`
# / `compute_resdist_PtE` (those clone internally and don't mutate `self`).
# Sharing skips the per-test construction cost across ~20 read-only tests.
make_L_graph <- function() {
  edges <- list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(1, 1))
  )
  metric_graph$new(edges = edges, perform_merges = TRUE,
                   check_connected = FALSE, verbose = 0)
}

.shared_L_graph <- NULL
shared_L_graph <- function() {
  if (is.null(.shared_L_graph)) {
    .shared_L_graph <<- make_L_graph()
  }
  .shared_L_graph
}

add_obs_to_graph <- function(g, PtE) {
  df <- data.frame(y = rnorm(nrow(PtE)),
                   edge_number = PtE[, 1],
                   distance_on_edge = PtE[, 2])
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
}

# ============================================================================ #
# compute_geodist_PtE
# ============================================================================ #

# ---- basic properties ----------------------------------------------------- #

test_that("geodist_PtE: basic shape and metric properties", {
  # Bundle the cheap sanity checks (dim, diag-zero, symmetry, non-negativity,
  # include_vertices flag) into one block so we pay the
  # compute_geodist_PtE cost (~1s) once instead of five times.
  g <- shared_L_graph()

  PtE_a <- cbind(c(1, 1, 2, 2), c(0.1, 0.9, 0.2, 0.7))
  D_a <- g$compute_geodist_PtE(PtE_a, normalized = TRUE,
                               include_vertices = FALSE, verbose = 0)
  expect_equal(dim(D_a), c(4, 4))
  expect_equal(diag(D_a), rep(0, 4))
  expect_equal(D_a, t(D_a))
  expect_true(all(D_a >= 0))

  # Include vertices: front block is g$nV vertices; back block is the points.
  PtE_b <- cbind(c(1, 2), c(0.5, 0.5))
  D_b <- g$compute_geodist_PtE(PtE_b, normalized = TRUE,
                               include_vertices = TRUE, verbose = 0)
  expect_equal(dim(D_b), c(nrow(PtE_b) + g$nV, nrow(PtE_b) + g$nV))
})

# ---- known distances on a simple graph ------------------------------------ #

test_that("geodist_PtE: correct on a single edge", {
  edges <- list(rbind(c(0, 0), c(2, 0)))
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  # Two points on the single edge (length = 2)
  PtE <- cbind(c(1, 1), c(0.25, 0.75))
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  # Distance = 2 * |0.75 - 0.25| = 1.0
  expect_equal(D[1, 2], 1.0, tolerance = 1e-10)
  expect_equal(D[2, 1], 1.0, tolerance = 1e-10)
})

test_that("geodist_PtE: across edges via shared vertex", {
  g <- shared_L_graph()
  # Point A at midpoint of edge 1, point B at midpoint of edge 2.
  # Geodesic = 0.5 + 0.5 = 1.0
  PtE <- cbind(c(1, 2), c(0.5, 0.5))
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(D[1, 2], 1.0, tolerance = 1e-10)
})

test_that("geodist_PtE: triangle respects shortest path", {
  g <- metric_graph$new(edges = make_triangle(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  # Two points very close to the same vertex but on different edges.
  eps <- 0.01
  PtE <- cbind(c(1, 2), c(1 - eps, eps))
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  e1_len <- g$edge_lengths[1]
  e2_len <- g$edge_lengths[2]
  expected <- eps * e1_len + eps * e2_len
  expect_equal(D[1, 2], expected, tolerance = 1e-8)
})

# ---- triangle inequality ------------------------------------------------- #

test_that("geodist_PtE: triangle inequality holds", {
  g <- metric_graph$new(edges = make_grid_edges(3), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  set.seed(42)
  PtE <- cbind(
    sample(seq_len(g$nE), 8, replace = TRUE),
    runif(8)
  )
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  # Reduce one assertion per triple (used to be n^3 = 512 expects in a
  # nested loop). The worst-case gap over all triples is a single number,
  # so summarise it once.
  n <- nrow(D)
  worst_gap <- 0
  for (k in seq_len(n)) {
    via_k <- outer(D[, k], D[k, ], `+`)
    worst_gap <- max(worst_gap, max(D - via_k))
  }
  expect_true(worst_gap <= 1e-10)
})

# ---- duplicate handling --------------------------------------------------- #

test_that("geodist_PtE: true duplicates produce warning and reduced matrix", {
  g <- shared_L_graph()
  PtE <- cbind(c(1, 1, 2), c(0.3, 0.3, 0.6))
  expect_warning(
    D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  expect_equal(dim(D), c(2, 2))
})

# ---- boundary duplicate regression --------------------------------------- #
# Points that are unique in PtE but map to the same vertex after
# standardize_df_positions().

test_that("geodist_PtE: boundary points at shared vertex detected as duplicates", {
  g <- shared_L_graph()
  # Vertex 2 = (1,0) is the end of edge 1 (pos=1) and start of edge 2 (pos=0).
  PtE <- cbind(c(1, 2, 1), c(1.0, 0.0, 0.5))
  expect_warning(
    D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  expect_equal(dim(D), c(2, 2))
})

test_that("geodist_PtE: unique boundary points at different vertices give full matrix", {
  g <- shared_L_graph()
  # Start of edge 1 (vertex 1) and end of edge 2 (vertex 3): different vertices
  PtE <- cbind(c(1, 2, 1), c(0.0, 1.0, 0.5))
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(dim(D), c(3, 3))
})

test_that("geodist_PtE: star tips are truly unique boundary points", {
  g <- metric_graph$new(edges = make_star(5), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  # Points at the tip of each ray (pos=1, all different vertices)
  PtE <- cbind(seq_len(g$nE), rep(1.0, g$nE))
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(dim(D), c(g$nE, g$nE))
})

test_that("geodist_PtE: star center duplicates are collapsed", {
  g <- metric_graph$new(edges = make_star(5), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  # 5 points all at the center vertex (pos=0 on each edge) + 1 at a tip
  PtE <- cbind(c(seq_len(g$nE), 1), c(rep(0.0, g$nE), 0.5))
  expect_warning(
    D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  expect_equal(dim(D), c(2, 2))
})

test_that("geodist_PtE: nrow equals unique standardized points", {
  g <- metric_graph$new(edges = make_triangle(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  set.seed(123)
  n <- 20
  PtE_interior <- cbind(
    sample(seq_len(g$nE), n, replace = TRUE),
    runif(n, 0.01, 0.99)
  )
  # 3 pairs of boundary points that each map to the same vertex
  PtE_boundary <- cbind(
    c(1, 2, 2, 3, 3, 1),
    c(1, 0, 1, 0, 1, 0)
  )
  PtE <- rbind(PtE_interior, PtE_boundary)
  expect_warning(
    D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  # 20 interior + 6 boundary that collapse to 3 unique vertices = 23
  expect_equal(nrow(D), n + 3)
  expect_equal(ncol(D), n + 3)
})

# ---- single point -------------------------------------------------------- #

test_that("geodist_PtE: single point returns a zero distance", {
  g <- shared_L_graph()
  PtE <- cbind(1, 0.5)
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(as.numeric(D), 0)
})

# ---- include_vertices ordering ------------------------------------------- #

test_that("geodist_PtE: include_vertices adds vertices at front of matrix", {
  g <- shared_L_graph()
  PtE <- cbind(c(1, 2), c(0.5, 0.5))
  D_with <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                                   include_vertices = TRUE, verbose = 0)
  D_without <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                                      include_vertices = FALSE, verbose = 0)
  nv <- g$nV
  np <- nrow(PtE)
  expect_equal(D_with[(nv + 1):(nv + np), (nv + 1):(nv + np)],
               D_without, tolerance = 1e-10)
})

# ============================================================================ #
# compute_geodist
# ============================================================================ #

test_that("compute_geodist (obs=FALSE) returns vertex distance matrix", {
  g <- make_L_graph()
  g$compute_geodist(obs = FALSE, verbose = 0)
  D <- g$geo_dist[[".vertices"]]
  expect_equal(dim(D), c(g$nV, g$nV))
  expect_equal(diag(D), rep(0, g$nV))
  expect_true(isSymmetric(D))
  # V1=(0,0), V2=(1,0), V3=(1,1): d(V1,V2)=1, d(V2,V3)=1, d(V1,V3)=2
  expect_equal(D[1, 2], 1, tolerance = 1e-10)
  expect_equal(D[2, 3], 1, tolerance = 1e-10)
  expect_equal(D[1, 3], 2, tolerance = 1e-10)
})

test_that("compute_geodist (all_groups) uses observations", {
  g <- make_L_graph()
  PtE <- cbind(c(1, 1, 2, 2), c(0.2, 0.7, 0.3, 0.8))
  add_obs_to_graph(g, PtE)
  g$compute_geodist(obs = TRUE, all_groups = TRUE,
                    include_vertices = FALSE, verbose = 0)
  D <- g$geo_dist[[".complete"]]
  expect_equal(dim(D), c(nrow(PtE), nrow(PtE)))
  expect_equal(diag(D), rep(0, nrow(PtE)))
  expect_true(isSymmetric(D))
})

test_that("compute_geodist per-group returns correct dimensions", {
  g <- make_L_graph()
  PtE <- cbind(c(1, 2), c(0.3, 0.6))
  df <- data.frame(y = c(1, 2),
                   edge_number = PtE[, 1],
                   distance_on_edge = PtE[, 2])
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  g$compute_geodist(obs = TRUE, all_groups = FALSE,
                    include_vertices = FALSE, verbose = 0)
  grps <- names(g$geo_dist)
  expect_true(length(grps) >= 1)
  for (nm in grps) {
    D <- g$geo_dist[[nm]]
    expect_true(all(D >= 0))
    expect_true(isSymmetric(D))
  }
})

test_that("compute_geodist include_vertices adds vertex rows/cols", {
  g <- make_L_graph()
  PtE <- cbind(c(1, 2), c(0.5, 0.5))
  add_obs_to_graph(g, PtE)

  g$compute_geodist(obs = TRUE, all_groups = TRUE,
                    include_vertices = TRUE, verbose = 0)
  D_with <- g$geo_dist[[".complete"]]

  g$compute_geodist(obs = TRUE, all_groups = TRUE,
                    include_vertices = FALSE, verbose = 0)
  D_without <- g$geo_dist[[".complete"]]

  nv <- g$nV
  np <- nrow(PtE)
  expect_equal(dim(D_with), c(np + nv, np + nv))
  expect_equal(dim(D_without), c(np, np))
  # Bottom-right block should match
  expect_equal(D_with[(nv + 1):(nv + np), (nv + 1):(nv + np)],
               D_without, tolerance = 1e-10)
})

# ---- compute_geodist: boundary duplicate regression ----------------------- #

test_that("compute_geodist handles boundary duplicates across edges", {
  g <- make_L_graph()
  # Place observations at the shared vertex via different edges
  PtE <- cbind(c(1, 2, 1), c(1.0, 0.0, 0.5))
  add_obs_to_graph(g, PtE)
  # After standardization inside compute_geodist_PtE, (1,1.0) and (2,0.0)
  # both map to vertex 2; only 2 unique locations remain.
  expect_warning(
    g$compute_geodist(obs = TRUE, all_groups = TRUE,
                      include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  D <- g$geo_dist[[".complete"]]
  expect_equal(dim(D), c(2, 2))
})

# ============================================================================ #
# compute_geodist_mesh
# ============================================================================ #

test_that("compute_geodist_mesh: basic shape and metric properties", {
  # Bundled sanity check (dim, diag-zero, symmetry, non-negativity) — one
  # mesh build + distance computation instead of three.
  g <- make_L_graph()
  g$build_mesh(h = 0.5)
  g$compute_geodist_mesh()
  D <- g$mesh$geo_dist
  expect_equal(dim(D), c(nrow(g$mesh$V), nrow(g$mesh$V)))
  expect_equal(diag(D), rep(0, nrow(D)))
  expect_true(isSymmetric(D))
  expect_true(all(D >= 0))
})

test_that("compute_geodist_mesh distances are consistent with edge lengths", {
  edges <- list(rbind(c(0, 0), c(3, 0)))
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  g$build_mesh(h = 1.0)
  g$compute_geodist_mesh()
  D <- g$mesh$geo_dist
  # The first and last mesh vertices should be graph endpoints
  # Total edge length is 3
  expect_equal(max(D), 3.0, tolerance = 1e-10)
})

test_that("compute_geodist_mesh triangle inequality holds", {
  g <- metric_graph$new(edges = make_triangle(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  g$build_mesh(h = 0.3)
  g$compute_geodist_mesh()
  D <- g$mesh$geo_dist
  # Vectorised worst-gap check across all triples — one assertion instead
  # of 20 per-sample expectations.
  n <- nrow(D)
  worst_gap <- 0
  for (k in seq_len(n)) {
    via_k <- outer(D[, k], D[k, ], `+`)
    worst_gap <- max(worst_gap, max(D - via_k))
  }
  expect_true(worst_gap <= 1e-10)
})

# ============================================================================ #
# compute_resdist_PtE
# ============================================================================ #

test_that("compute_resdist_PtE: basic shape and metric properties", {
  # Bundled sanity check (dim, diag-zero, symmetry, non-negativity,
  # include_vertices flag) — one compute call instead of five.
  g <- make_L_graph()

  PtE_a <- cbind(c(1, 1, 2, 2), c(0.1, 0.9, 0.2, 0.7))
  R_a <- g$compute_resdist_PtE(PtE_a, normalized = TRUE,
                               include_vertices = FALSE, verbose = 0)
  expect_equal(dim(R_a), c(4, 4))
  expect_equal(diag(R_a), rep(0, 4), tolerance = 1e-10)
  expect_equal(R_a, t(R_a), tolerance = 1e-10)
  expect_true(all(R_a >= -1e-10))

  PtE_b <- cbind(c(1, 2), c(0.5, 0.5))
  R_b <- g$compute_resdist_PtE(PtE_b, normalized = TRUE,
                               include_vertices = TRUE, verbose = 0)
  expect_equal(dim(R_b), c(nrow(PtE_b) + g$nV, nrow(PtE_b) + g$nV))
})

test_that("compute_resdist_PtE on a line equals geodesic distance", {
  # On a tree (no cycles), resistance distance equals geodesic distance
  edges <- list(rbind(c(0, 0), c(2, 0)))
  g <- metric_graph$new(edges = edges, perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  PtE <- cbind(c(1, 1), c(0.25, 0.75))
  R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(R[1, 2], D[1, 2], tolerance = 1e-8)
})

test_that("compute_resdist_PtE on a tree equals geodesic distance", {
  # Star graph is a tree — resistance == geodesic
  g <- metric_graph$new(edges = make_star(3), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  PtE <- cbind(c(1, 2, 3), c(0.5, 0.5, 0.5))
  R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(as.numeric(R), as.numeric(D), tolerance = 1e-6)
})

test_that("compute_resdist_PtE on cycle is less than geodesic", {
  # On a graph with cycles, resistance distance < geodesic distance
  g <- metric_graph$new(edges = make_triangle(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  PtE <- cbind(c(1, 2), c(0.5, 0.5))
  R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_true(R[1, 2] < D[1, 2])
})

# ---- compute_resdist_PtE: duplicate handling ------------------------------ #

test_that("compute_resdist_PtE warns on true duplicate PtE", {
  g <- shared_L_graph()
  PtE <- cbind(c(1, 1, 2), c(0.3, 0.3, 0.6))
  expect_warning(
    R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  expect_equal(dim(R), c(2, 2))
})

# ---- compute_resdist_PtE: boundary duplicate regression ------------------- #

test_that("compute_resdist_PtE detects boundary duplicates across edges", {
  g <- shared_L_graph()
  # (1, 1.0) and (2, 0.0) both map to vertex 2 after standardization
  PtE <- cbind(c(1, 2, 1), c(1.0, 0.0, 0.5))
  expect_warning(
    R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  expect_equal(dim(R), c(2, 2))
})

test_that("compute_resdist_PtE unique boundary points give full matrix", {
  g <- shared_L_graph()
  # Start of edge 1 (vertex 1) and end of edge 2 (vertex 3): different vertices
  PtE <- cbind(c(1, 2, 1), c(0.0, 1.0, 0.5))
  R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(dim(R), c(3, 3))
})

test_that("compute_resdist_PtE star center duplicates are collapsed", {
  g <- metric_graph$new(edges = make_star(5), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  # 5 points at center (all map to same vertex) + 1 interior point
  PtE <- cbind(c(seq_len(g$nE), 1), c(rep(0.0, g$nE), 0.5))
  expect_warning(
    R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  expect_equal(dim(R), c(2, 2))
})

test_that("compute_resdist_PtE dim matches unique standardized points", {
  g <- metric_graph$new(edges = make_triangle(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  set.seed(123)
  n <- 10
  PtE_interior <- cbind(
    sample(seq_len(g$nE), n, replace = TRUE),
    runif(n, 0.01, 0.99)
  )
  # 3 pairs of boundary points that each map to the same vertex
  PtE_boundary <- cbind(
    c(1, 2, 2, 3, 3, 1),
    c(1, 0, 1, 0, 1, 0)
  )
  PtE <- rbind(PtE_interior, PtE_boundary)
  expect_warning(
    R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                                include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  expect_equal(nrow(R), n + 3)
  expect_equal(ncol(R), n + 3)
})

# ---- compute_resdist_PtE: include_vertices ordering ---------------------- #

test_that("compute_resdist_PtE include_vertices adds front block", {
  g <- shared_L_graph()
  PtE <- cbind(c(1, 2), c(0.5, 0.5))
  R_with <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                                   include_vertices = TRUE, verbose = 0)
  R_without <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                                      include_vertices = FALSE, verbose = 0)
  nv <- g$nV
  np <- nrow(PtE)
  expect_equal(R_with[(nv + 1):(nv + np), (nv + 1):(nv + np)],
               R_without, tolerance = 1e-10)
})

# ============================================================================ #
# compute_resdist (wrapper)
# ============================================================================ #

test_that("compute_resdist (obs=FALSE) returns vertex resistance matrix", {
  g <- shared_L_graph()
  g$compute_resdist(obs = FALSE, verbose = 0)
  R <- g$res_dist[[".vertices"]]
  expect_equal(dim(R), c(g$nV, g$nV))
  expect_true(all(abs(diag(R)) < 1e-10))
  expect_true(all(R >= -1e-10))
})

test_that("compute_resdist (full) uses observations", {
  g <- shared_L_graph()
  PtE <- cbind(c(1, 1, 2), c(0.2, 0.7, 0.4))
  add_obs_to_graph(g, PtE)
  g$compute_resdist(full = TRUE, include_vertices = FALSE, verbose = 0)
  R <- g$res_dist[[".complete"]]
  expect_equal(dim(R), c(nrow(PtE), nrow(PtE)))
  expect_true(all(abs(diag(R)) < 1e-10))
})

test_that("compute_resdist per-group returns correct dimensions", {
  g <- shared_L_graph()
  PtE <- cbind(c(1, 2), c(0.3, 0.6))
  df <- data.frame(y = c(1, 2),
                   edge_number = PtE[, 1],
                   distance_on_edge = PtE[, 2])
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  g$compute_resdist(full = FALSE, verbose = 0)
  grps <- names(g$res_dist)
  expect_true(length(grps) >= 1)
  for (nm in grps) {
    R <- g$res_dist[[nm]]
    expect_true(all(R >= -1e-10))
  }
})

test_that("compute_resdist handles boundary duplicates", {
  g <- make_L_graph()
  PtE <- cbind(c(1, 2, 1), c(1.0, 0.0, 0.5))
  add_obs_to_graph(g, PtE)
  expect_warning(
    g$compute_resdist(full = TRUE, include_vertices = FALSE, verbose = 0),
    "Duplicated locations"
  )
  R <- g$res_dist[[".complete"]]
  expect_equal(dim(R), c(2, 2))
})

# ============================================================================ #
# Cross-function consistency
# ============================================================================ #

test_that("geodist via compute_geodist matches compute_geodist_PtE directly", {
  g <- make_L_graph()
  PtE <- cbind(c(1, 1, 2), c(0.2, 0.8, 0.4))
  add_obs_to_graph(g, PtE)
  g$compute_geodist(obs = TRUE, all_groups = TRUE,
                    include_vertices = FALSE, verbose = 0)
  D_indirect <- g$geo_dist[[".complete"]]
  D_direct <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                                     include_vertices = FALSE, verbose = 0)
  expect_equal(D_indirect, D_direct, tolerance = 1e-10)
})

test_that("resdist via compute_resdist matches compute_resdist_PtE directly", {
  g <- make_L_graph()
  PtE <- cbind(c(1, 1, 2), c(0.2, 0.8, 0.4))
  add_obs_to_graph(g, PtE)
  g$compute_resdist(full = TRUE, include_vertices = FALSE, verbose = 0)
  R_indirect <- g$res_dist[[".complete"]]
  R_direct <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                                     include_vertices = FALSE, verbose = 0)
  expect_equal(R_indirect, R_direct, tolerance = 1e-10)
})

test_that("resistance <= geodesic on graph with cycles", {
  g <- metric_graph$new(edges = make_triangle(), perform_merges = TRUE,
                        check_connected = FALSE, verbose = 0)
  PtE <- cbind(c(1, 2, 3), c(0.3, 0.5, 0.7))
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  # On a graph with cycles, resistance <= geodesic for all pairs
  expect_true(all(R <= D + 1e-10))
})

test_that("resistance == geodesic on tree", {
  g <- make_L_graph()
  PtE <- cbind(c(1, 1, 2, 2), c(0.2, 0.7, 0.3, 0.8))
  D <- g$compute_geodist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  R <- g$compute_resdist_PtE(PtE, normalized = TRUE,
                              include_vertices = FALSE, verbose = 0)
  expect_equal(as.numeric(R), as.numeric(D), tolerance = 1e-6)
})
