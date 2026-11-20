fractional_precision_test_graph <- function() {
  g <- metric_graph$new(
    edges = list(rbind(c(0, 0), c(1, 0)), rbind(c(1, 0), c(2, 1)),
                 rbind(c(1, 0), c(2, -1))), verbose = 0
  )
  g$set_manual_edge_lengths(c(1, 0.7, 1))
  g
}

test_that("integer precision matrices and triplets retain their existing formats", {
  g <- fractional_precision_test_graph()
  for (alpha in 1:2) {
    assembler <- if (alpha == 1) Qalpha1 else Qalpha2
    for (BC in 0:3) {
      for (build in c(TRUE, FALSE)) {
        reference <- assembler(c(0.8, 1.7), g, BC = BC, build = build)
        expect_identical(spde_precision(1.7, 0.8, alpha, g, BC, build), reference)
        # Fractional-only settings must not affect integer calls.
        expect_identical(spde_precision(1.7, 0.8, alpha, g, BC, build,
                                         m = NA, type_rational_approx = "unused",
                                         type_interp = "unused"), reference)
      }
    }
  }
  expect_identical(spde_precision(0, 0.8, 1, g), Qalpha1(c(0.8, 0), g))
})

test_that("fractional precision returns separate scaled matrices or triplets", {
  g <- fractional_precision_test_graph()
  before <- serialize(list(g$E, g$V, g$PtV, g$C, g$CoB, g$constraints), NULL)
  for (alpha in c(0.75, 1.5, 2.5)) {
    for (m in c(1L, 2L, 3L)) {
      q <- floor(alpha)
      n0 <- if (q == 0) g$nV else 2L * q * g$nE
      ni <- 2L * ceiling(alpha) * g$nE
      Q <- spde_precision(1.7, 0.8, alpha, g, m = m)
      entries <- spde_precision(1.7, 0.8, alpha, g, m = m, build = FALSE)
      scaled <- spde_precision(1.7, 1.6, alpha, g, m = m)
      expect_type(Q, "list")
      expect_named(Q, paste0("Q", 0:m))
      expect_named(entries, names(Q))
      for (i in seq_along(Q)) {
        expect_s4_class(Q[[i]], "sparseMatrix")
        expect_equal(dim(Q[[i]]), rep(if (i == 1L) n0 else ni, 2))
        expect_equal(as.matrix(Q[[i]]), t(as.matrix(Q[[i]])), tolerance = 1e-12)
        expect_equal(scaled[[i]], 4 * Q[[i]], tolerance = 1e-12)
        component <- entries[[i]]
        expect_named(component, c("i", "j", "x", "dims"))
        rebuilt <- Matrix::sparseMatrix(i = component$i, j = component$j,
                                        x = component$x, dims = component$dims)
        expect_identical(rebuilt, Q[[i]])
      }
    }
  }
  expect_identical(before, serialize(list(g$E, g$V, g$PtV, g$C, g$CoB, g$constraints), NULL))
})

test_that("zeroth fractional blocks agree with the integer spectral component", {
  g <- fractional_precision_test_graph()
  kappa <- 1.7
  tau <- 0.8
  for (alpha in c(0.75, 1.5, 2.5)) {
    coefficients <- rSPDE:::interp_rational_coefficients(
      order = 2, type_rational_approx = "chebfun", type_interp = "spline", alpha = alpha
    )
    q <- floor(alpha)
    if (q == 0) {
      expected <- Matrix::Diagonal(g$nV, tau^2 * kappa^(2 * alpha) / coefficients$k)
    } else {
      integer_raw <- if (q == 1) {
        Qalpha1_edges(c(1, kappa), g, w = 0.5, stationary_points = "none")
      } else {
        Qalpha2(c(1, kappa), g, BC = 0)
      }
      expected <- integer_raw * tau^2 * kappa^(2 * (alpha - q)) / coefficients$k
    }
    Q <- spde_precision(kappa, tau, alpha, g)
    expect_equal(as.matrix(Q$Q0), as.matrix(expected), tolerance = 1e-10)
  }
})

test_that("fractional computation preserves existing bases and uses raw BC semantics", {
  g <- fractional_precision_test_graph()
  g$PtV <- rev(seq_len(g$nV))
  g$buildC(2, edge_constraint = FALSE)
  for (alpha in 1:3) g$build_constraints(alpha = alpha)
  before <- serialize(list(g$E, g$PtV, g$C, g$CoB, g$constraints), NULL)
  for (alpha in c(0.75, 1.5, 2.5)) {
    Q <- spde_precision(1.7, 0.8, alpha, g)
    expect_identical(Q, spde_precision(1.7, 0.8, alpha, g, BC = 0))
  }
  expect_identical(before, serialize(list(g$E, g$PtV, g$C, g$CoB, g$constraints), NULL))
})

test_that("fractional option combinations retain component matrices", {
  g <- fractional_precision_test_graph()
  for (method in c("chebfun", "brasil", "chebfunLB")) {
    for (interpolation in c("linear", "log", "spline", "logspline")) {
      Q <- spde_precision(1.7, 0.8, 1.5, g,
                          type_rational_approx = method, type_interp = interpolation)
      expect_named(Q, c("Q0", "Q1", "Q2"))
      for (component in Q) expect_s4_class(component, "sparseMatrix")
    }
  }
})

test_that("disconnected fractional precision keeps components and graphs separate", {
  # Unequal graph sizes make accidental concatenation or transposition visible.
  edges <- list(rbind(c(0, 0), c(1, 0)), rbind(c(3, 0), c(4, 0)),
                rbind(c(4, 0), c(5, 1)))
  components <- suppressWarnings(graph_components$new(edges = edges, verbose = 0))
  before <- lapply(components$graphs, function(g) {
    serialize(list(g$E, g$PtV, g$C, g$CoB, g$constraints), NULL)
  })
  for (alpha in c(0.75, 1.5, 2.5)) {
    for (m in c(1L, 3L)) {
      expected <- lapply(components$graphs, function(g) {
        spde_precision(1.7, 0.8, alpha, g, m = m)
      })
      Q <- spde_precision(1.7, 0.8, alpha, components, m = m)
      expect_named(Q, paste0("Q", 0:m))
      for (component in names(Q)) {
        expect_named(Q[[component]], paste0("G", seq_along(expected)))
        for (i in seq_along(expected)) {
          expect_identical(Q[[component]][[i]], expected[[i]][[component]])
        }
      }
    }
  }
  after <- lapply(components$graphs, function(g) {
    serialize(list(g$E, g$PtV, g$C, g$CoB, g$constraints), NULL)
  })
  expect_identical(after, before)
  for (alpha in c(0.75, 1, 1.5, 2, 2.5, 3)) {
    expect_error(spde_precision(1.7, 0.8, alpha, components, build = FALSE),
                 "build = FALSE is not supported")
  }
})

test_that("invalid arguments fail before precision assembly", {
  g <- fractional_precision_test_graph()
  call_precision <- function(...) {
    args <- modifyList(list(kappa = 1.7, tau = 0.8, alpha = 1.5, graph = g), list(...),
                       keep.null = TRUE)
    do.call(spde_precision, args)
  }
  for (bad in list(NULL, NA, Inf, NaN, c(1, 2), "1", 1+1i, 0.5, 3.5, 4)) {
    expect_error(call_precision(alpha = bad), "alpha must")
  }
  for (bad in list(NULL, NA, Inf, c(1, 2), "1", 1+1i, -1, 0)) {
    expect_error(call_precision(kappa = bad), "kappa must")
    expect_error(call_precision(tau = bad), "tau must")
  }
  for (bad in list(NULL, NA, Inf, c(1, 2), "1", 1+1i, 0, -1, 1.5)) {
    expect_error(call_precision(m = bad), "m must")
  }
  for (bad in list(NULL, NA, 0, 1, c(TRUE, FALSE))) {
    expect_error(call_precision(build = bad), "build must")
  }
  for (alpha in 1:2) {
    for (bad in list(NULL, NA, NaN, Inf, c(0, 1), "1", 1+1i)) {
      expect_error(call_precision(alpha = alpha, BC = bad), "BC must")
    }
  }
  expect_error(call_precision(graph = list()), "graph must")
  expect_error(call_precision(type_rational_approx = "wl2"), "type_rational_approx must")
  expect_error(call_precision(type_rational_approx = c("chebfun", "brasil")),
               "type_rational_approx must")
  expect_error(call_precision(type_interp = "other"), "type_interp must")
  expect_error(call_precision(m = 9999), "Could not obtain rational coefficients")

  g_bad <- g$clone(deep = TRUE)
  g_bad$edge_lengths[1] <- 0
  expect_error(call_precision(graph = g_bad), "edge_lengths must")
  g_bad$edge_lengths <- g_bad$edge_lengths[-1]
  expect_error(call_precision(graph = g_bad), "edge_lengths must")
  g_bad <- g$clone(deep = TRUE)
  g_bad$E[1, 1] <- g_bad$nV + 1L
  expect_error(call_precision(graph = g_bad), "valid vertex indices")
})


test_that("degenerate rational coefficients and numerical overflow are rejected", {
  g <- fractional_precision_test_graph()
  expect_error(spde_precision(1.7, 1e300, 1.5, g), "variance scale")
  coefficients <- rSPDE:::interp_rational_coefficients(
    order = 2, type_rational_approx = "chebfun", type_interp = "spline", alpha = 1.5
  )
  supplied <- coefficients
  local_mocked_bindings(
    interp_rational_coefficients = function(...) supplied,
    .package = "rSPDE"
  )
  for (field in c("k", "r", "p")) {
    supplied <- coefficients
    supplied[[field]][1] <- NA_real_
    expect_error(spde_precision(1.7, 0.8, 1.5, g), "Rational coefficients must")
    supplied <- coefficients
    supplied[[field]][1] <- 0
    expect_error(spde_precision(1.7, 0.8, 1.5, g), "Rational coefficients must")
  }
})


test_that("fractional components and build FALSE remain triplets throughout", {
  g <- fractional_precision_test_graph()
  # Any sparse construction here would defeat the purpose of this path.
  local_mocked_bindings(
    sparseMatrix = function(...) stop("unexpected sparse construction"),
    bdiag = function(...) stop("unexpected block-matrix construction"),
    .package = "Matrix"
  )
  for (alpha in c(0.75, 1.5, 2.5)) {
    components <- fractional_spde_precision_components(1.7, 0.8, alpha, g, m = 2)
    expect_named(components, c("Q0", "Q1", "Q2"))
    for (component in components) {
      expect_type(component, "list")
      expect_named(component, c("i", "j", "x", "dims"))
      expect_length(component$i, length(component$x))
      expect_length(component$j, length(component$x))
      expect_true(all(is.finite(component$x)))
    }
    entries <- spde_precision(1.7, 0.8, alpha, g, build = FALSE)
    expect_identical(entries, components)
  }
})

test_that("fractional matrix output does not assemble a joint block matrix", {
  g <- fractional_precision_test_graph()
  local_mocked_bindings(
    bdiag = function(...) stop("unexpected block-matrix construction"),
    .package = "Matrix"
  )
  for (alpha in c(0.75, 1.5, 2.5)) {
    Q <- spde_precision(1.7, 0.8, alpha, g)
    expect_named(Q, c("Q0", "Q1", "Q2"))
    for (component in Q) expect_s4_class(component, "sparseMatrix")
  }
})

test_that("triplet assembly retains zero blocks and their coordinate offsets", {
  blocks <- list(
    list(i = integer(0), j = integer(0), x = numeric(0), dims = c(2L, 2L)),
    list(i = 1L, j = 1L, x = 4, dims = c(3L, 3L)),
    list(i = integer(0), j = integer(0), x = numeric(0), dims = c(1L, 1L))
  )
  entries <- combine_spde_precision_triplets(blocks)
  expect_identical(entries, list(i = 3L, j = 3L, x = 4, dims = c(6L, 6L)))
  expect_identical(combine_spde_precision_triplets(list()),
                   list(i = integer(0), j = integer(0), x = numeric(0), dims = c(0L, 0L)))
})


test_that("BC validation rejects codes outside the package boundary options", {
  g <- fractional_precision_test_graph()
  for (alpha in 1:2) {
    for (BC in c(-1, -0.5, 0.5, 1.5, 2.5, 4, 3000)) {
      expect_error(spde_precision(1.7, 0.8, alpha, g, BC = BC),
                   "BC must be one of 0, 1, 2, or 3")
    }
  }
})

test_that("BC is not validated or used by raw alpha-three and fractional backends", {
  g <- fractional_precision_test_graph()
  for (alpha in c(0.75, 1.5, 2.5, 3)) {
    expected <- spde_precision(1.7, 0.8, alpha, g)
    for (ignored in list(NULL, NA, Inf, c(0, 1), "unused")) {
      expect_identical(spde_precision(1.7, 0.8, alpha, g, BC = ignored), expected)
    }
  }
})


test_that("fractional components retain their values across derivative units", {
  g <- metric_graph$new(edges = list(rbind(c(0, 0), c(1, 0))), verbose = 0)
  for (alpha in c(0.75, 1.4, 2.6)) {
    g$set_manual_edge_lengths(1)
    reference <- spde_precision(kappa = 1, tau = 0.8, alpha = alpha, graph = g)
    for (kappa in c(1e-4, 1e4)) {
      g$set_manual_edge_lengths(1 / kappa)
      Q <- spde_precision(kappa, tau = 0.8, alpha = alpha, graph = g)
      expect_named(Q, names(reference))
      for (i in seq_along(Q)) {
        expect_true(all(is.finite(Q[[i]])))
        q <- floor(alpha)
        if (i == 1L && q == 0L) {
          dimensionless <- as.matrix(Q[[i]]) / kappa^(2 * alpha)
        } else {
          n_state <- if (i == 1L) q else q + 1L
          scales <- rep(kappa^(0:(n_state - 1L)), 2)
          dimensionless <- as.matrix(Q[[i]]) * outer(scales, scales) /
            kappa^(2 * alpha - 1)
        }
        expect_equal(dimensionless, as.matrix(reference[[i]]), tolerance = 1e-9)
      }
    }
  }
})
