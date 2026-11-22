fractional_constraint_graph <- function(E = matrix(c(1L, 2L, 2L, 3L, 2L, 4L),
                                                  ncol = 2, byrow = TRUE)) {
  nV <- max(E)
  g <- metric_graph$new(V = cbind(seq_len(nV), seq_len(nV) %% 2), E = E,
                        manual_edge_lengths = rep(1, nrow(E)),
                        auto_remove_point_edges = FALSE,
                        check_connected = FALSE, verbose = 0)
  g$PtV <- c(rev(seq_len(g$nV)), 1L)
  g
}

test_that("fractional requests store only two templates for the selected interval", {
  for (alpha in c(0.7, 1.2, 2.2)) {
    g <- fractional_constraint_graph()
    q <- floor(alpha)
    key <- paste0("alpha", q, q + 1L)
    interval <- if (q == 0) "0.5 < alpha < 1" else paste0(q, " < alpha < ", q + 1L)
    expect_message(g$build_constraints(alpha), paste0(interval, ".*", key))
    expect_named(g$constraints, key)
    templates <- g$constraints[[key]]
    expect_named(templates, c("Tc_joint", "A"))
    expect_named(templates$Tc_joint, c("Tc1", "Tc2"))
    expect_named(templates$A, c("A1", "A2"))
    expect_null(g$constraints$alpha1)
    expect_null(g$constraints$alpha2)
    expect_null(g$constraints$alpha3)
    expect_equal(nrow(templates$A$A1), length(g$PtV))
    expect_equal(nrow(templates$A$A2), length(g$PtV))
    expect_equal(ncol(templates$A$A1), nrow(templates$Tc_joint$Tc1))
    expect_equal(ncol(templates$A$A2), nrow(templates$Tc_joint$Tc2))
  }
})

test_that("fractional requests reuse existing integer bases and avoid joint assembly", {
  g <- fractional_constraint_graph()
  suppressMessages(g$build_constraints(alpha = 1))
  suppressMessages(g$build_constraints(alpha = 2))
  before <- g$constraints
  local_mocked_bindings(
    build_constraint_matrices = function(...) stop("unexpected integer rebuild"),
    .package = "MetricGraph"
  )
  local_mocked_bindings(bdiag = function(...) stop("unexpected joint matrix"),
                         .package = "Matrix")
  expect_message(g$build_constraints(alpha = 1.2), "alpha12")
  expect_identical(g$constraints$alpha1, before$alpha1)
  expect_identical(g$constraints$alpha2, before$alpha2)
  expect_identical(g$constraints$alpha12$Tc_joint$Tc1, before$alpha1$Tc)
  expect_identical(g$constraints$alpha12$Tc_joint$Tc2, before$alpha2$Tc)
  expect_identical(g$constraints$alpha12$A$A1, before$alpha1$A)
  expect_identical(g$constraints$alpha12$A$A2, before$alpha2$A)
})

test_that("only the integer orders needed by a fractional request are computed", {
  g <- fractional_constraint_graph()
  original <- build_constraint_matrices
  seen <- integer(0)
  local_mocked_bindings(
    build_constraint_matrices = function(E, nV, PtV, alpha) {
      seen <<- c(seen, alpha)
      original(E, nV, PtV, alpha)
    },
    .package = "MetricGraph"
  )
  suppressMessages(g$build_constraints(alpha = 2.4))
  expect_setequal(seen, c(2L, 3L))
  expect_length(seen, 2L)
  expect_named(g$constraints, "alpha23")
})

test_that("white-noise templates select vertices directly in PtV order", {
  g <- fractional_constraint_graph()
  suppressMessages(g$build_constraints(alpha = 0.7))
  setup <- g$constraints$alpha01
  expect_equal(as.matrix(setup$Tc_joint$Tc1), diag(g$nV))
  values <- seq_len(g$nV) * 0.3
  expect_equal(as.vector(setup$A$A1 %*% values), values[g$PtV])

  # Compare with the private selector in its normal observation-data setting.
  observed <- metric_graph$new(edges = list(rbind(c(0, 0), c(1, 0)),
                                             rbind(c(1, 0), c(2, 0))), verbose = 0)
  observed$add_observations(data = data.frame(y = c(1, 2), edge_number = c(1, 2),
                                               distance_on_edge = c(0.3, 0.6)),
                             normalized = TRUE, verbose = 0)
  observed$observation_to_vertex(mesh_warning = FALSE)
  observed$PtV <- c(rev(observed$PtV), observed$PtV[1])
  suppressMessages(observed$build_constraints(alpha = 0.7))
  expect_equal(observed$constraints$alpha01$A$A1,
               observed$.__enclos_env__$private$A())
})

test_that("template pairs reconstruct joint covariance for different rational orders", {
  graphs <- list(
    matrix(c(1L, 2L), ncol = 2),
    matrix(c(1L, 2L, 2L, 3L, 2L, 4L), ncol = 2, byrow = TRUE),
    matrix(c(1L, 2L, 2L, 2L), ncol = 2, byrow = TRUE),
    matrix(c(1L, 2L, 2L, 1L), ncol = 2, byrow = TRUE),
    matrix(c(1L, 2L, 3L, 4L), ncol = 2, byrow = TRUE)
  )
  for (E in graphs) {
    g <- fractional_constraint_graph(E)
    for (alpha in c(0.75, 1.5, 2.5)) {
      suppressMessages(g$build_constraints(alpha))
      q <- floor(alpha)
      setup <- g$constraints[[paste0("alpha", q, q + 1L)]]
      original <- setup
      T1 <- setup$Tc_joint$Tc1
      T2 <- setup$Tc_joint$Tc2
      A1 <- setup$A$A1
      A2 <- setup$A$A2
      for (m in 1:3) {
        joint <- Matrix::bdiag(c(list(T1), rep(list(T2), m)))
        A <- cbind(A1, do.call(cbind, rep(list(A2), m)))
        raw <- spde_precision(1.7, 0.8, alpha, g, BC = 0, m = m)
        expect_equal(ncol(joint), sum(vapply(raw, nrow, integer(1))))
        raw_joint <- Matrix::bdiag(raw)
        expect_equal(ncol(A), nrow(joint))
        reduced <- joint %*% raw_joint %*% t(joint)
        Sigma <- A %*% solve(reduced, t(A))
        # Independent components contribute additive observed covariances.
        Q0 <- T1 %*% raw$Q0 %*% t(T1)
        sum_cov <- A1 %*% solve(Q0, t(A1))
        for (i in seq_len(m)) {
          Qi <- T2 %*% raw[[i + 1L]] %*% t(T2)
          sum_cov <- sum_cov + A2 %*% solve(Qi, t(A2))
        }
        expect_equal(as.matrix(Sigma), as.matrix(sum_cov), tolerance = 1e-8)
        expect_true(all(is.finite(Sigma)))
        expect_identical(g$constraints[[paste0("alpha", q, q + 1L)]], original)
      }
    }
  }
})

test_that("fractional entries are reused across alpha values within an interval", {
  g <- fractional_constraint_graph()
  for (alpha in c(0.7, 1.2, 2.2)) suppressMessages(g$build_constraints(alpha))
  cached <- g$constraints
  local_mocked_bindings(
    build_fractional_constraint_templates = function(...) stop("unexpected rebuild"),
    .package = "MetricGraph"
  )
  for (alpha in c(0.9, 1.9, 2.9)) expect_message(g$build_constraints(alpha), "Constraint templates")
  expect_identical(g$constraints, cached)
})

test_that("fractional template invalidation preserves unrelated directional entries", {
  g <- fractional_constraint_graph()
  suppressMessages(g$build_constraints(alpha = 1, directional = TRUE))
  directional <- g$constraints$alpha1_directional
  for (alpha in c(0.7, 1.2, 2.2)) suppressMessages(g$build_constraints(alpha))
  g$buildC(2, edge_constraint = FALSE)
  legacy <- serialize(list(g$C, g$CoB), NULL)
  g$PtV <- rev(g$PtV)
  expect_named(g$constraints, "alpha1_directional")
  expect_identical(g$constraints$alpha1_directional, directional)
  expect_identical(legacy, serialize(list(g$C, g$CoB), NULL))
  suppressMessages(g$build_constraints(alpha = 1.8))
  cached <- g$constraints$alpha12
  g$set_manual_edge_lengths(g$edge_lengths * 2)
  expect_identical(g$constraints$alpha12, cached)
  g$DirectionalWeightFunction_out <- function(w) -2
  expect_named(g$constraints, "alpha12")
  expect_identical(g$constraints$alpha12, cached)
  g$E[1, ] <- rev(g$E[1, ])
  expect_null(g$constraints)
})

test_that("fractional construction preserves legacy caches and integer formats", {
  g <- fractional_constraint_graph()
  g$buildDirectionalConstraints(alpha = 1)
  before <- serialize(list(g$C, g$CoB, g$E, g$PtV), NULL)
  for (alpha in c(0.7, 1.2, 2.2)) suppressMessages(g$build_constraints(alpha))
  expect_identical(before, serialize(list(g$C, g$CoB, g$E, g$PtV), NULL))
  fractional <- g$constraints
  for (alpha in 1:3) {
    suppressMessages(g$build_constraints(alpha))
    expect_named(g$constraints[[paste0("alpha", alpha)]],
                 c("K", "T", "S", "keep", "Tc", "indices", "A"))
  }
  for (key in names(fractional)) expect_identical(g$constraints[[key]], fractional[[key]])
})

test_that("fractional templates handle empty selections and validate the range", {
  g <- fractional_constraint_graph()
  g$PtV <- integer(0)
  for (alpha in c(0.7, 1.2, 2.2)) {
    suppressMessages(g$build_constraints(alpha))
    setup <- g$constraints[[paste0("alpha", floor(alpha), ceiling(alpha))]]
    expect_equal(nrow(setup$A$A1), 0L)
    expect_equal(nrow(setup$A$A2), 0L)
  }
  before <- g$constraints
  for (bad in list(0, 0.5, 3.1, NA_real_, Inf, NaN, "1.2", c(1.2, 1.3))) {
    expect_error(g$build_constraints(bad), "alpha must")
  }
  expect_error(g$build_constraints(1.2, directional = TRUE), "only supported")
  expect_error(g$build_constraints(0.7, m = 3), "unused argument")
  expect_identical(g$constraints, before)
  expect_error(build_kirchhoff_constraint_matrix(g$E, g$nV, 1.2), "integer")
  g$PtV <- NULL
  expect_error(g$build_constraints(0.7), "PtV is NULL")
  expect_null(g$constraints)
})
