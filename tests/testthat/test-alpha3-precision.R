alpha3_precision_test_graph <- function(E = matrix(c(1L, 2L, 2L, 3L, 2L, 4L),
                                                  ncol = 2, byrow = TRUE)) {
  V <- cbind(seq_len(max(E)), seq_len(max(E)) %% 2)
  metric_graph$new(V = V, E = E, manual_edge_lengths = 0.7 + 0.3 * seq_len(nrow(E)),
                   auto_remove_point_edges = FALSE, check_connected = FALSE, verbose = 0)
}

# Independent closed-form Matérn nu=5/2 covariance and its derivatives.
alpha3_stationary_joint <- function(s, t, kappa, tau) {
  h <- s - t
  a <- abs(h)
  derivatives <- c(
    1 + kappa * a + kappa^2 * a^2 / 3,
    -kappa^2 * a / 3 - kappa^3 * a^2 / 3,
    -kappa^2 / 3 - kappa^3 * a / 3 + kappa^4 * a^2 / 3,
    kappa^4 * a - kappa^5 * a^2 / 3,
    kappa^4 - 5 * kappa^5 * a / 3 + kappa^6 * a^2 / 3
  ) * exp(-kappa * a) * 3 / (16 * tau^2 * kappa^5)
  if (h < 0) derivatives[c(2, 4)] <- -derivatives[c(2, 4)]
  outer(0:2, 0:2, function(i, j) (-1)^j * derivatives[i + j + 1L])
}

test_that("alpha-three precision agrees with the closed-form endpoint covariance", {
  g <- alpha3_precision_test_graph(matrix(c(1L, 2L), ncol = 2))
  kappa <- 1.4
  tau <- 0.7
  R00 <- alpha3_stationary_joint(0, 0, kappa, tau)
  R01 <- alpha3_stationary_joint(0, g$edge_lengths[1], kappa, tau)
  stationary <- rbind(cbind(R00, R01), cbind(t(R01), R00))
  expected <- solve(stationary) - as.matrix(Matrix::bdiag(solve(R00), solve(R00))) / 2
  Q <- spde_precision(kappa, tau, alpha = 3, graph = g)
  expect_s4_class(Q, "sparseMatrix")
  expect_equal(dim(Q), c(6L, 6L))
  expect_equal(as.matrix(Q), expected, tolerance = 1e-10)
})

test_that("alpha three is exact and does not request rational coefficients", {
  g <- alpha3_precision_test_graph()
  before <- serialize(list(g$E, g$PtV, g$C, g$CoB, g$constraints), NULL)
  local_mocked_bindings(
    interp_rational_coefficients = function(...) stop("unexpected rational approximation"),
    .package = "rSPDE"
  )
  Q <- spde_precision(1.4, 0.7, 3, g)
  expect_equal(dim(Q), rep(6L * g$nE, 2))
  expect_identical(Q, spde_precision(1.4, 0.7, 3, g, m = NA,
                                    type_rational_approx = "unused", type_interp = "unused"))
  expect_equal(spde_precision(1.4, 1.4, 3, g), 4 * Q, tolerance = 1e-12)
  expect_identical(Q, spde_precision(1.4, 0.7, 3, g, BC = 0))
  expect_identical(before, serialize(list(g$E, g$PtV, g$C, g$CoB, g$constraints), NULL))
})

test_that("alpha-three triplets are built without sparse matrix assembly", {
  g <- alpha3_precision_test_graph()
  Q <- spde_precision(1.4, 0.7, 3, g)
  entries <- spde_precision(1.4, 0.7, 3, g, build = FALSE)
  expect_named(entries, c("i", "j", "x", "dims"))
  expect_equal(as.matrix(Matrix::sparseMatrix(i = entries$i, j = entries$j,
                                              x = entries$x, dims = entries$dims)), as.matrix(Q))
  local_mocked_bindings(
    sparseMatrix = function(...) stop("unexpected sparse construction"),
    bdiag = function(...) stop("unexpected block-matrix construction"),
    .package = "Matrix"
  )
  expect_identical(spde_precision(1.4, 0.7, 3, g, build = FALSE), entries)
})

test_that("alpha-three raw precision can be constrained externally without changing caches", {
  g <- alpha3_precision_test_graph()
  g$PtV <- c(rev(seq_len(g$nV)), 1L)
  g$build_constraints(alpha = 3)
  g$buildC(2, edge_constraint = FALSE)
  before <- serialize(list(g$C, g$CoB, g$constraints), NULL)
  Q <- spde_precision(1.4, 0.7, 3, g)
  setup <- g$constraints$alpha3
  reduced <- setup$Tc %*% Q %*% t(setup$Tc)
  Sigma <- setup$A %*% solve(reduced, t(setup$A))
  expect_true(all(is.finite(Sigma)))
  expect_equal(as.matrix(Sigma), t(as.matrix(Sigma)), tolerance = 1e-10)
  expect_gt(min(eigen(as.matrix(reduced), symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_identical(before, serialize(list(g$C, g$CoB, g$constraints), NULL))
})

test_that("alpha-three edge assembly supports loops and disconnected components", {
  for (E in list(matrix(c(1L, 1L), ncol = 2),
                 matrix(c(1L, 2L, 2L, 2L), ncol = 2, byrow = TRUE),
                 matrix(c(1L, 2L, 2L, 1L), ncol = 2, byrow = TRUE),
                 matrix(c(1L, 2L, 3L, 4L), ncol = 2, byrow = TRUE))) {
    g <- alpha3_precision_test_graph(E)
    Q <- spde_precision(1.4, 0.7, 3, g)
    expect_equal(dim(Q), rep(6L * g$nE, 2))
    expect_equal(as.matrix(Q), t(as.matrix(Q)), tolerance = 1e-12)
  }
  components <- suppressWarnings(graph_components$new(
    edges = list(rbind(c(0, 0), c(1, 0)), rbind(c(3, 0), c(4, 0))), verbose = 0
  ))
  expected <- Matrix::bdiag(lapply(components$graphs, function(g) spde_precision(1.4, 0.7, 3, g)))
  expect_identical(spde_precision(1.4, 0.7, 3, components), expected)
  expect_error(spde_precision(1.4, 0.7, 3, components, build = FALSE),
               "build = FALSE is not supported")
})

test_that("alpha-three raw assembly validates parameters and graph inputs", {
  g <- alpha3_precision_test_graph()
  expect_error(spde_precision(0, 0.7, 3, g), "kappa must")
  expect_error(spde_precision(1.4, 0, 3, g), "tau must")
  expect_error(spde_precision(1.4, 1e300, 3, g), "variance scale")
  g$edge_lengths[1] <- NA_real_
  expect_error(spde_precision(1.4, 0.7, 3, g), "edge_lengths must")
})
