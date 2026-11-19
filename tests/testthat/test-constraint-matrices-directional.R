test_that("build_constraints constructs only the requested case", {
  g <- make_directional_test_graph()
  expect_error(g$build_constraints(), "alpha")
  expect_error(g$build_constraints(alpha = 2, directional = TRUE), "only supported")
  expect_error(g$build_constraints(alpha = 1, directional = NA), "TRUE or FALSE")
  expect_error(g$build_constraints(alpha = 1, directional = 1), "TRUE or FALSE")
  expect_null(g$constraints)
  expect_null(g$PtV)
  expect_invisible(g$build_constraints(alpha = 1, directional = TRUE))
  expect_named(g$constraints, "alpha1_directional")
  expect_named(g$constraints$alpha1_directional, c("K", "T", "S", "keep", "Tc"))
  g$PtV <- rev(seq_len(g$nV))
  directional <- g$constraints$alpha1_directional
  expect_invisible(g$build_constraints(alpha = 2))
  expect_named(g$constraints, c("alpha1_directional", "alpha2"))
  expect_identical(g$constraints$alpha1_directional, directional)
  expect_null(g$constraints$alpha1)
  expect_null(g$constraints$alpha3)
  alpha2 <- g$constraints$alpha2
  expect_invisible(g$build_constraints(alpha = 3))
  expect_identical(g$constraints$alpha2, alpha2)
  expect_invisible(g$build_constraints(alpha = 1))
  expect_setequal(names(g$constraints), c("alpha1", "alpha2", "alpha3", "alpha1_directional"))
})

test_that("non-directional requests construct only their selected order", {
  g <- make_directional_test_graph()
  g$PtV <- rev(seq_len(g$nV))
  original <- build_kirchhoff_constraint_matrix
  local_mocked_bindings(
    build_kirchhoff_constraint_matrix = function(E, nV, alpha) {
      if (alpha != 2) stop("unexpected order")
      original(E, nV, alpha)
    },
    build_directional_constraint_matrices = function(...) stop("unexpected directional setup"),
    .package = "MetricGraph"
  )
  expect_no_error(g$build_constraints(alpha = 2))
  expect_named(g$constraints, "alpha2")
})

test_that("directional cache construction never invokes a legacy builder", {
  g <- make_directional_test_graph()
  # Fail immediately if a directional request attempts non-directional work.
  local_mocked_bindings(
    build_constraint_matrices = function(...) stop("unexpected non-directional setup"),
    build_kirchhoff_constraint_matrix = function(...) stop("unexpected non-directional K"),
    .package = "MetricGraph"
  )
  for (method in c("buildDirectionalConstraints", "buildC")) {
    unlockBinding(method, g)
    g[[method]] <- function(...) stop("unexpected legacy method")
    lockBinding(method, g)
  }
  expect_no_error(g$build_constraints(alpha = 1, directional = TRUE))
  expect_null(g$C)
  expect_null(g$CoB)
})

test_that("directional setups match the legacy basis and endpoint covariance", {
  graphs <- list(make_directional_test_graph(), make_reversed_directional_test_graph(),
                 make_nondendritic_directional_graph(),
                 metric_graph$new(edges = list(rbind(c(0, 0), c(1, 0))), verbose = 0),
                 metric_graph$new(edges = list(rbind(c(0, 0), c(1, 0)),
                                                rbind(c(3, 0), c(4, 0))),
                                  check_connected = FALSE, verbose = 0))
  for (g in graphs) {
    # Custom callbacks exercise the same generalized weighting path as the
    # unchanged method, including source forks and merging incoming edges.
    for (custom in c(FALSE, TRUE)) {
      if (custom) {
        g$DirectionalWeightFunction_in <- function(w) sqrt(w) / sum(sqrt(w))
        g$DirectionalWeightFunction_out <- function(w) -0.75
      }
      g$buildC(alpha = 2, edge_constraint = FALSE)
      snapshot <- serialize(list(g$C, g$CoB, g$E, g$PtV), NULL)
      reference <- g$clone(deep = TRUE)
      reference$buildDirectionalConstraints(alpha = 1)
      g$build_constraints(alpha = 1, directional = TRUE)
      expect_identical(snapshot, serialize(list(g$C, g$CoB, g$E, g$PtV), NULL))
      setup <- g$constraints$alpha1_directional
      expect_equal(setup$K, reference$C, tolerance = 1e-13)
      expect_equal(as.matrix(setup$T), as.matrix(reference$CoB$T), tolerance = 1e-13)
      expect_equal(as.numeric(setup$S), as.numeric(reference$CoB$S), tolerance = 1e-13)
      keep <- setdiff(seq_len(nrow(reference$CoB$T)), seq_along(reference$CoB$S))
      old_Tc <- reference$CoB$T[keep, , drop = FALSE]
      expect_equal(as.matrix(setup$Tc), as.matrix(old_Tc), tolerance = 1e-13)
      if (nrow(setup$K)) {
        expect_lt(max(abs(setup$K %*% t(setup$Tc))), 1e-10)
      }
      Q <- Qalpha1_edges(c(0.8, 1.2), g, w = 0, BC = 1, build = TRUE)
      reduced <- setup$Tc %*% Q %*% t(setup$Tc)
      old_reduced <- old_Tc %*% Q %*% t(old_Tc)
      expect_equal(as.matrix(reduced), as.matrix(old_reduced), tolerance = 1e-12)
      Sigma <- t(setup$Tc) %*% solve(reduced, setup$Tc)
      old_Sigma <- t(old_Tc) %*% solve(old_reduced, old_Tc)
      expect_equal(as.matrix(Sigma), as.matrix(old_Sigma), tolerance = 1e-11)
    }
  }
})

test_that("directional cache is reused without reevaluating weighting callbacks", {
  g <- make_directional_test_graph()
  calls <- 0L
  g$DirectionalWeightFunction_out <- function(w) { calls <<- calls + 1L; -1 }
  g$build_constraints(alpha = 1, directional = TRUE)
  initial_calls <- calls
  expect_gt(initial_calls, 0)
  cached <- g$constraints
  local_mocked_bindings(
    build_directional_constraint_matrices = function(...) stop("unexpected rebuild"),
    .package = "MetricGraph"
  )
  expect_invisible(g$build_constraints(alpha = 1, directional = TRUE))
  expect_identical(g$constraints, cached)
  expect_identical(calls, initial_calls)
})

test_that("cache entries invalidate independently according to their inputs", {
  g <- make_directional_test_graph()
  g$PtV <- seq_len(g$nV)
  g$build_constraints(alpha = 2)
  g$build_constraints(alpha = 1, directional = TRUE)
  directional <- g$constraints$alpha1_directional
  g$PtV <- rev(g$PtV)
  expect_named(g$constraints, "alpha1_directional")
  expect_identical(g$constraints$alpha1_directional, directional)
  g$build_constraints(alpha = 2)
  nondirectional <- g$constraints$alpha2
  weights <- g$get_edge_weights()
  weights$w[1] <- weights$w[1] + 0.5
  g$set_edge_weights(weights, directional_weights = "w", verbose = 0)
  expect_named(g$constraints, "alpha2")
  expect_identical(g$constraints$alpha2, nondirectional)
  g$build_constraints(alpha = 1, directional = TRUE)
  g$setDirectionalWeightFunction(f_in = function(w) rep(1 / length(w), length(w)))
  expect_named(g$constraints, "alpha2")
  expect_identical(g$constraints$alpha2, nondirectional)
  g$build_constraints(alpha = 1, directional = TRUE)
  g$set_manual_edge_lengths(g$edge_lengths * 2)
  expect_setequal(names(g$constraints), c("alpha2", "alpha1_directional"))
  clone <- g$clone(deep = TRUE)
  clone$DirectionalWeightFunction_out <- function(w) -2
  expect_named(clone$constraints, "alpha2")
  expect_setequal(names(g$constraints), c("alpha2", "alpha1_directional"))
  g$E[1, ] <- rev(g$E[1, ])
  expect_null(g$constraints)
})

test_that("directional weight selection and callback failures do not corrupt caches", {
  g <- make_directional_test_graph()
  g$PtV <- seq_len(g$nV)
  g$build_constraints(alpha = 3)
  keep <- g$constraints$alpha3
  weights <- g$get_edge_weights()
  weights$other <- rev(weights$w) + 1
  g$set_edge_weights(weights, directional_weights = "w", verbose = 0)
  g$build_constraints(alpha = 1, directional = TRUE)
  g$set_edge_weights(weights, directional_weights = "other", verbose = 0)
  expect_named(g$constraints, "alpha3")
  expect_identical(g$constraints$alpha3, keep)
  g$buildC(2, edge_constraint = TRUE)
  snapshot <- serialize(list(g$C, g$CoB), NULL)
  g$DirectionalWeightFunction_out <- function(w) stop("callback failed")
  expect_error(g$build_constraints(alpha = 1, directional = TRUE), "callback failed")
  expect_identical(snapshot, serialize(list(g$C, g$CoB), NULL))
  expect_named(g$constraints, "alpha3")
  expect_identical(g$constraints$alpha3, keep)
})
