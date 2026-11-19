build_all_constraint_test_cases <- function(graph) {
  for (alpha in 1:3) graph$build_constraints(alpha = alpha)
  invisible(NULL)
}

constraint_test_graph <- function(E) {
  # Explicit V/E retains edge orientation, loops and repeated edges.
  nV <- max(E)
  V <- cbind(seq_len(nV), seq_len(nV) %% 2)
  metric_graph$new(V = V, E = E, manual_edge_lengths = rep(1, nrow(E)),
                   auto_remove_point_edges = FALSE, check_connected = FALSE, verbose = 0)
}

constraint_vertex_covariance <- function(setup, precision) {
  Q <- setup$Tc %*% precision %*% t(setup$Tc)
  as.matrix(setup$A %*% solve(Q, t(setup$A)))
}

test_that("integer constraint setup is opt-in and preserves the legacy caches", {
  g <- constraint_test_graph(matrix(c(1L, 2L), ncol = 2))
  expect_null(g$constraints)
  expect_null(g$C)
  expect_null(g$CoB)
  expect_error(g$build_constraints(alpha = 1), "PtV is NULL")
  expect_null(g$constraints)

  g$PtV <- c(2L, 1L, 2L)
  g$buildC(alpha = 2, edge_constraint = FALSE)
  before <- serialize(list(g$C, g$CoB, g$E, g$PtV), NULL)
  expect_invisible(build_all_constraint_test_cases(g))
  expect_identical(before, serialize(list(g$C, g$CoB, g$E, g$PtV), NULL))
  expect_named(g$constraints, paste0("alpha", 1:3))
  expect_error(g$constraints <- list(), "read-only")

  # Rebuilding the legacy cache must not overwrite the new cache either.
  cached <- g$constraints
  g$buildC(alpha = 2, edge_constraint = TRUE)
  expect_identical(g$constraints, cached)
})

test_that("all integer bases impose the correct endpoint conditions", {
  graphs <- list(
    interval = matrix(c(1L, 2L), ncol = 2),
    branch = matrix(c(1L, 2L, 2L, 3L, 2L, 4L), ncol = 2, byrow = TRUE),
    loop = matrix(c(1L, 1L), ncol = 2),
    tadpole = matrix(c(1L, 2L, 2L, 2L), ncol = 2, byrow = TRUE),
    parallel = matrix(c(1L, 2L, 2L, 1L), ncol = 2, byrow = TRUE),
    disconnected = matrix(c(1L, 2L, 3L, 4L), ncol = 2, byrow = TRUE)
  )
  for (E in graphs) {
    g <- constraint_test_graph(E)
    g$PtV <- rev(seq_len(g$nV))
    build_all_constraint_test_cases(g)
    incidence <- as.vector(t(g$E))
    degree <- tabulate(incidence, nbins = g$nV)
    nc <- sum(pmax(degree - 1L, 0L))
    for (alpha in 1:3) {
      setup <- g$constraints[[alpha]]
      N <- 2L * alpha * g$nE
      expected_rows <- c(nc, nc + g$nV, 2L * nc + g$nV)[alpha]
      expect_s4_class(setup$K, "sparseMatrix")
      expect_equal(dim(setup$K), c(expected_rows, N))
      expect_equal(as.matrix(tcrossprod(setup$T)), diag(N), tolerance = 1e-10)
      expect_equal(nrow(setup$Tc), N - expected_rows)
      if (expected_rows > 0L) {
        expect_equal(as.integer(Matrix::rankMatrix(setup$K)), expected_rows)
        expect_equal(as.matrix(setup$K %*% t(setup$Tc)),
                     matrix(0, expected_rows, nrow(setup$Tc)), tolerance = 1e-10)
      }

      # Independently inspect a free endpoint vector using the vertex laws.
      endpoint_values <- as.vector(t(setup$Tc) %*% seq_len(nrow(setup$Tc)))
      for (v in seq_len(g$nV)) {
        ends <- which(incidence == v)
        value_indices <- alpha * (ends - 1L) + 1L
        expect_equal(diff(endpoint_values[value_indices]),
                     rep(0, length(ends) - 1L), tolerance = 1e-10)
        if (alpha >= 2L) {
          signs <- ifelse(ends %% 2L == 1L, 1, -1)
          expect_equal(sum(signs * endpoint_values[value_indices + 1L]),
                       0, tolerance = 1e-10)
        }
        if (alpha == 3L) {
          expect_equal(diff(endpoint_values[value_indices + 2L]),
                       rep(0, length(ends) - 1L), tolerance = 1e-10)
        }
      }
    }
  }
})

test_that("vertex selection prefers starts and preserves repeated PtV order", {
  E <- matrix(c(2L, 1L, 3L, 2L, 2L, 2L), ncol = 2, byrow = TRUE)
  setups <- lapply(1:3, function(alpha) {
    build_constraint_matrices(E, 3L, c(1L, 2L, 3L, 2L), alpha)
  })
  # Vertex 1 only ends an edge; vertex 2 starts the first edge, even though
  # it is also an end elsewhere; vertex 3 starts the second edge.
  ends <- c(2L, 1L, 3L, 1L)
  for (alpha in 1:3) {
    setup <- setups[[alpha]]
    expect_equal(setup$indices, alpha * (ends - 1L) + 1L)
    z <- seq_len(nrow(setup$Tc))
    full <- as.vector(t(setup$Tc) %*% z)
    expect_equal(as.vector(setup$A %*% z), full[setup$indices])
  }
})

test_that("alpha-two covariance agrees with buildC including Neumann conditions", {
  E <- matrix(c(1L, 2L, 2L, 3L, 2L, 2L), ncol = 2, byrow = TRUE)
  g <- constraint_test_graph(E)
  g$PtV <- c(3L, 1L, 2L, 2L)
  build_all_constraint_test_cases(g)
  setup <- g$constraints$alpha2
  g$buildC(alpha = 2, edge_constraint = TRUE)
  T <- g$CoB$T
  keep <- setdiff(seq_len(nrow(T)), seq_along(g$CoB$S))
  Tc <- T[keep, , drop = FALSE]
  reference <- list(Tc = Tc, A = t(Tc)[setup$indices, , drop = FALSE])
  expect_equal(as.integer(Matrix::rankMatrix(rbind(setup$K, g$C))), nrow(g$C))
  set.seed(812)
  N <- ncol(setup$K)
  precision <- crossprod(matrix(rnorm(N * N), N, N)) + diag(N)
  expect_equal(constraint_vertex_covariance(setup, precision),
               constraint_vertex_covariance(reference, precision), tolerance = 1e-9)
})

test_that("constraint cache is reused and invalidates on changes to its inputs", {
  g <- constraint_test_graph(matrix(c(1L, 2L, 2L, 3L), ncol = 2, byrow = TRUE))
  g$PtV <- c(1L, 3L)
  build_all_constraint_test_cases(g)
  cached <- g$constraints
  # A cached call must not invoke the constructor again.
  local_mocked_bindings(
    build_constraint_matrices = function(...) stop("unexpected rebuild"),
    .package = "MetricGraph"
  )
  expect_invisible(build_all_constraint_test_cases(g))
  expect_identical(g$constraints, cached)
  g$set_manual_edge_lengths(g$edge_lengths * 2)
  expect_identical(g$constraints, cached)
  clone <- g$clone(deep = TRUE)
  clone$PtV <- rev(clone$PtV)
  expect_null(clone$constraints)
  expect_identical(g$constraints, cached)
  g$E[1, ] <- rev(g$E[1, ])
  expect_null(g$constraints)
})

test_that("observation changes invalidate constraint maps", {
  g <- constraint_test_graph(matrix(c(1L, 2L), ncol = 2))
  g$add_observations(data = data.frame(edge_number = 1, distance_on_edge = 0.25, y = 1),
                     normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE)
  build_all_constraint_test_cases(g)
  expect_type(g$constraints, "list")
  g$add_observations(data = data.frame(edge_number = 1, distance_on_edge = 0.5, y = 2),
                     normalized = TRUE, verbose = 0)
  g$observation_to_vertex(mesh_warning = FALSE)
  expect_null(g$constraints)
  build_all_constraint_test_cases(g)
  expect_type(g$constraints, "list")
  g$clear_observations()
  expect_null(g$constraints)
})

test_that("empty selections and invalid inputs are handled explicitly", {
  E <- matrix(c(1L, 2L), ncol = 2)
  for (alpha in 1:3) {
    setup <- build_constraint_matrices(E, 2L, integer(0), alpha)
    expect_equal(nrow(setup$A), 0L)
    expect_equal(ncol(setup$A), nrow(setup$Tc))
    empty <- build_constraint_matrices(matrix(integer(0), ncol = 2), 0L,
                                       integer(0), alpha)
    expect_equal(dim(empty$A), c(0L, 0L))
  }
  expect_error(build_constraint_matrices(E, 3L, 3L, 1), "incident to an edge")
  expect_error(build_constraint_matrices(E, 2L, c(NA, 1), 1), "PtV must be")
  expect_error(build_constraint_matrices(E, 2L, 1.5, 1), "PtV must be")
  expect_error(build_kirchhoff_constraint_matrix(E, 2L, 4), "alpha must be")
  expect_error(build_kirchhoff_constraint_matrix(E * 2, 2L, 1), "vertex indices")
})
