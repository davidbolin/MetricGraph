# Independent references for nonzero rational components.
# Closed-form physical Green functions are independent of the production
# endpoint helper, matern.p.* functions, and sigma2 normalization.
# Moderate edge lengths keep this test separate from short-edge conditioning.

fractional_ref_triplet_matrix <- function(entries) {
  as.matrix(Matrix::sparseMatrix(i = entries$i, j = entries$j,
                                x = entries$x, dims = entries$dims))
}
fractional_ref_relative_error <- function(actual, expected) {
  max(abs(actual - expected)) / max(abs(expected))
}
fractional_ref_expect_lt <- function(object, expected, info = NULL) {
  expect_true(object < expected, info = info)
}
# These tests verify coefficient use, not the approximation tables themselves.
fractional_ref_coefficients <- function(alpha, m) {
  rSPDE:::interp_rational_coefficients(
    order = m, type_rational_approx = "chebfun",
    type_interp = "spline", alpha = alpha
  )
}
fractional_ref_graph <- function() {
  g <- metric_graph$new(
    edges = list(rbind(c(0, 0), c(0.8, 0)),
                 rbind(c(0.8, 0), c(2, 0))), verbose = 0
  )
  # Exercise observation ordering as well as the degree-two vertex.
  g$PtV <- c(3L, 1L, 2L)
  g
}
fractional_ref_cases <- expand.grid(alpha = c(0.75, 1.4, 2.6), m = c(2L, 3L))
fractional_ref_cases$kappa <- ifelse(fractional_ref_cases$m == 2L, 1.3, 0.9)
fractional_ref_cases$tau <- ifelse(fractional_ref_cases$m == 2L, 0.7, 1.15)

# Derivatives of exp(-rate*h) times a polynomial, for h >= 0.
fractional_ref_exp_polynomial_derivative <- function(polynomial, rate, order, h) {
  for (j in seq_len(order)) {
    derivative <- if (length(polynomial) == 1L) 0 else {
      c(seq_len(length(polynomial) - 1L) * polynomial[-1L], 0)
    }
    polynomial <- derivative - rate * polynomial
  }
  sum(polynomial * h^(seq_along(polynomial) - 1L)) * exp(-rate * h)
}

# Real-line Green function of
# beta / ((kappa^2 + w^2)^q * (a^2 + w^2)), a=kappa*sqrt(1-p).
# These physical coefficients come directly from the covariance operator:
# beta = r_i * tau^-2 * kappa^(2*(q+1-alpha)).
# In particular this reference never computes sigma2 or calls matern.p.*.
fractional_ref_endpoint_block <- function(length, kappa, tau, alpha, pole, residue) {
  q <- floor(alpha)
  n <- q + 1L
  a <- kappa * sqrt(1 - pole)
  delta <- a^2 - kappa^2
  beta <- residue * kappa^(2 * (q + 1 - alpha)) / tau^2
  terms <- switch(as.character(q),
    `0` = list(list(rate = a, polynomial = beta / (2 * a))),
    `1` = list(
      list(rate = kappa, polynomial = beta / (2 * kappa * delta)),
      list(rate = a, polynomial = -beta / (2 * a * delta))
    ),
    `2` = list(
      list(rate = kappa, polynomial = beta * c(
        -1 / (2 * kappa * delta^2) + 1 / (4 * kappa^3 * delta),
        1 / (4 * kappa^2 * delta)
      )),
      list(rate = a, polynomial = beta / (2 * a * delta^2))
    )
  )
  fractional_ref_precision_from_terms(length, n, terms)
}

fractional_ref_precision_from_terms <- function(length, n, terms) {
  joint <- function(distance) {
    outer(0:(n - 1L), 0:(n - 1L), Vectorize(function(i, j) {
      order <- i + j
      if (distance == 0 && order %% 2L == 1L) return(0)
      value <- sum(vapply(terms, function(term) {
        fractional_ref_exp_polynomial_derivative(term$polynomial, term$rate,
                                      order, abs(distance))
      }, numeric(1)))
      (-1)^j * if (distance < 0 && order %% 2L == 1L) -value else value
    }))
  }
  C0 <- joint(0)
  C0l <- joint(-length)
  C <- rbind(cbind(C0, C0l), cbind(t(C0l), C0))
  correction <- as.matrix(Matrix::bdiag(solve(C0) / 2, solve(C0) / 2))
  solve(C) - correction
}
fractional_ref_raw_component <- function(g, kappa, tau, alpha, pole, residue) {
  as.matrix(Matrix::bdiag(lapply(g$edge_lengths, function(length) {
    fractional_ref_endpoint_block(length, kappa, tau, alpha, pole, residue)
  })))
}

# Exact Neumann resolvent kernel on [0,L], using stable image terms.
# power=2 is -(1/(2*rate)) times the derivative of the power=1 kernel.
fractional_ref_neumann_kernel <- function(points, L, rate, power = 1L) {
  denominator <- -expm1(-2 * rate * L)
  outer(points, points, Vectorize(function(s, t) {
    distances <- c(abs(s - t), s + t, 2 * L - s - t,
                   2 * L - abs(s - t))
    terms <- exp(-rate * distances)
    if (power == 1L) return(sum(terms) / (2 * rate * denominator))
    sum(terms * (distances + 1 / rate +
                  2 * L * exp(-2 * rate * L) / denominator)) /
      (4 * rate^2 * denominator)
  }))
}
fractional_ref_interval_component <- function(points, L, kappa, tau, alpha, pole, residue) {
  q <- floor(alpha)
  a <- kappa * sqrt(1 - pole)
  delta <- a^2 - kappa^2
  beta <- residue * kappa^(2 * (q + 1 - alpha)) / tau^2
  Ga <- fractional_ref_neumann_kernel(points, L, a)
  if (q == 0L) return(beta * Ga)
  Gk <- fractional_ref_neumann_kernel(points, L, kappa)
  if (q == 1L) return(beta * (Gk - Ga) / delta)
  Hk <- fractional_ref_neumann_kernel(points, L, kappa, power = 2L)
  beta * ((Ga - Gk) / delta^2 + Hk / delta)
}
fractional_ref_zero_interval_component <- function(g, kappa, tau, alpha, constant) {
  q <- floor(alpha)
  # For alpha<1 the package Q0 is a discrete vertex white-noise block.
  # This does not interpret continuous white noise as pointwise values.
  if (q == 0L) return(diag(g$nV) * constant / (tau^2 * kappa^(2 * alpha)))
  points <- g$V[g$PtV, 1L]
  constant * kappa^(2 * (q - alpha)) / tau^2 *
    fractional_ref_neumann_kernel(points, sum(g$edge_lengths), kappa, power = q)
}
fractional_ref_constrained_covariance <- function(Q, setup) {
  Q <- as.matrix(Q)
  Tc <- as.matrix(setup$Tc)
  A <- as.matrix(setup$A)
  A %*% solve(Tc %*% Q %*% t(Tc), t(A))
}


test_that("nonzero fractional components match physical endpoint Green functions", {
  g <- fractional_ref_graph()
  for (case in seq_len(nrow(fractional_ref_cases))) {
    args <- as.list(fractional_ref_cases[case, ])
    coeff <- fractional_ref_coefficients(args$alpha, args$m)
    Q <- do.call(spde_precision, c(args, list(graph = g)))
    entries <- do.call(spde_precision, c(args, list(graph = g, build = FALSE)))
    expect_named(Q, paste0("Q", 0:args$m))
    expect_named(entries, names(Q))
    for (i in seq_len(args$m)) {
      component <- paste0("Q", i)
      expected <- fractional_ref_raw_component(
        g, args$kappa, args$tau, args$alpha, coeff$p[i], coeff$r[i]
      )
      info <- sprintf("alpha=%g, m=%d, %s", args$alpha, args$m, component)
      fractional_ref_expect_lt(
        fractional_ref_relative_error(as.matrix(Q[[component]]), expected),
        2e-8, info = paste("matrix", info)
      )
      fractional_ref_expect_lt(
        fractional_ref_relative_error(fractional_ref_triplet_matrix(entries[[component]]),
                                      expected),
        2e-8, info = paste("triplets", info)
      )
    }
  }
})

test_that("fractional constrained covariance matches independent interval resolvents", {
  g <- fractional_ref_graph()
  points <- g$V[g$PtV, 1L]
  L <- sum(g$edge_lengths)
  for (case in seq_len(nrow(fractional_ref_cases))) {
    args <- as.list(fractional_ref_cases[case, ])
    q <- floor(args$alpha)
    coeff <- fractional_ref_coefficients(args$alpha, args$m)
    components <- do.call(spde_precision, c(args, list(graph = g)))
    suppressMessages(g$build_constraints(alpha = args$alpha))
    templates <- g$constraints[[paste0("alpha", q, q + 1L)]]
    setup0 <- list(Tc = templates$Tc_joint$Tc1, A = templates$A$A1)
    setupi <- list(Tc = templates$Tc_joint$Tc2, A = templates$A$A2)
    expected <- fractional_ref_zero_interval_component(
      g, args$kappa, args$tau, args$alpha, coeff$k
    )
    fractional_ref_expect_lt(fractional_ref_relative_error(
      fractional_ref_constrained_covariance(components$Q0, setup0), expected), 2e-8)
    for (i in seq_len(args$m)) {
      reference <- fractional_ref_interval_component(
        points, L, args$kappa, args$tau, args$alpha, coeff$p[i], coeff$r[i]
      )
      actual <- fractional_ref_constrained_covariance(components[[i + 1L]], setupi)
      fractional_ref_expect_lt(fractional_ref_relative_error(actual, reference), 2e-8,
        info = sprintf("alpha=%g, m=%d, constrained Q%d", args$alpha, args$m, i))
      expected <- expected + reference
    }
    # Assemble in the same rational-component order as the public API and
    # external covariance wrappers. The interval reference uses no basis.
    setups <- c(list(setup0), rep(list(setupi), args$m))
    Tc_joint <- as.matrix(Matrix::bdiag(lapply(setups, function(s) s$Tc)))
    A_joint <- do.call(cbind, lapply(setups, function(s) as.matrix(s$A)))
    Qraw <- as.matrix(Matrix::bdiag(components))
    Q <- Tc_joint %*% Qraw %*% t(Tc_joint)
    actual <- A_joint %*% solve(Q, t(A_joint))
    fractional_ref_expect_lt(fractional_ref_relative_error(actual, expected), 2e-8,
      info = sprintf("alpha=%g, m=%d, joint covariance", args$alpha, args$m))
  }
})


test_that("disconnected metric_graph keeps original coordinates for every fractional component", {
  # Two intervals, with interleaved edges A1, B1, A2, B2. The last edge is
  # reversed, all lengths differ, and PtV is reordered with a repeated vertex.
  # References below use these explicit fixture coordinates, never
  # get_components(), production endpoint kernels, or production output order.
  # Vertex IDs follow first endpoint occurrence, matching the constructor's
  # canonical numbering without grouping vertices by connected piece.
  vertices <- rbind(c(0, 0), c(0.8, 0), c(5, 0),
                    c(6.1, 0), c(2.4, 0), c(6.75, 0))
  edges <- matrix(c(1L, 2L, 3L, 4L, 2L, 5L, 6L, 4L), ncol = 2, byrow = TRUE)
  lengths <- c(0.8, 1.1, 1.6, 0.65)
  observed <- c(6L, 2L, 4L, 1L, 5L, 3L, 2L)
  g <- metric_graph$new(V = vertices, E = edges, manual_edge_lengths = lengths,
                        auto_remove_point_edges = FALSE,
                        check_connected = FALSE, verbose = 0)
  g$PtV <- observed
  expect_true(g$is_disconnected())
  expect_equal(g$E, edges)
  expect_equal(unname(g$V), vertices)
  expect_equal(g$edge_lengths, lengths)

  interval_vertices <- list(c(1L, 2L, 5L), c(3L, 4L, 6L))
  interval_points <- list(c(0, 0.8, 2.4), c(0, 1.1, 1.75))
  for (case in seq_len(nrow(fractional_ref_cases))) {
    args <- as.list(fractional_ref_cases[case, ])
    q <- floor(args$alpha)
    coeff <- fractional_ref_coefficients(args$alpha, args$m)
    beta0 <- coeff$k * args$kappa^(2 * (q - args$alpha)) / args$tau^2
    matrices <- do.call(spde_precision, c(args, list(graph = g)))
    entries <- do.call(spde_precision, c(args, list(graph = g, build = FALSE)))
    rebuilt <- lapply(entries, fractional_ref_triplet_matrix)
    expect_named(matrices, paste0("Q", 0:args$m))
    expect_named(entries, names(matrices))
    suppressMessages(g$build_constraints(alpha = args$alpha))
    templates <- g$constraints[[paste0("alpha", q, q + 1L)]]
    setup0 <- list(Tc = templates$Tc_joint$Tc1, A = templates$A$A1)
    setupi <- list(Tc = templates$Tc_joint$Tc2, A = templates$A$A2)
    expected_joint <- matrix(0, length(observed), length(observed))

    for (i in 0:args$m) {
      component <- paste0("Q", i)
      info <- sprintf("disconnected alpha=%g, m=%d, %s", args$alpha, args$m, component)
      if (i == 0L && q == 0L) {
        # Discrete vertex white noise is ordered by the original vertices.
        variance <- coeff$k / (args$tau^2 * args$kappa^(2 * args$alpha))
        expected_raw <- diag(nrow(vertices)) / variance
        expected_cov <- outer(observed, observed, "==") * variance
      } else {
        n_state <- if (i == 0L) q else q + 1L
        expected_raw <- matrix(0, 2L * n_state * nrow(edges),
                               2L * n_state * nrow(edges))
        for (edge in seq_along(lengths)) {
          if (i == 0L) {
            # Integer Green kernel for beta0 * (kappa^2 - d^2/ds^2)^(-q).
            polynomial <- if (q == 1L) beta0 / (2 * args$kappa) else {
              beta0 * c(1, args$kappa) / (4 * args$kappa^3)
            }
            block <- fractional_ref_precision_from_terms(lengths[edge], n_state,
              list(list(rate = args$kappa, polynomial = polynomial)))
          } else {
            block <- fractional_ref_endpoint_block(lengths[edge], args$kappa,
              args$tau, args$alpha, coeff$p[i], coeff$r[i])
          }
          ids <- 2L * n_state * (edge - 1L) + seq_len(2L * n_state)
          expected_raw[ids, ids] <- block
        }
        # Neumann kernels on the two intervals, placed directly in the
        # original full graph's PtV order. Cross-piece covariance is zero.
        expected_cov <- matrix(0, length(observed), length(observed))
        for (piece in seq_along(interval_vertices)) {
          ids <- which(observed %in% interval_vertices[[piece]])
          points <- interval_points[[piece]][match(observed[ids], interval_vertices[[piece]])]
          L <- tail(interval_points[[piece]], 1L)
          expected_cov[ids, ids] <- if (i == 0L) {
            beta0 * fractional_ref_neumann_kernel(points, L, args$kappa, power = q)
          } else {
            fractional_ref_interval_component(points, L, args$kappa, args$tau,
              args$alpha, coeff$p[i], coeff$r[i])
          }
        }
      }
      setup <- if (i == 0L) setup0 else setupi
      for (output in list(matrices, rebuilt)) {
        fractional_ref_expect_lt(fractional_ref_relative_error(
          as.matrix(output[[component]]), expected_raw), 2e-8, info = paste("raw", info))
        actual_cov <- fractional_ref_constrained_covariance(output[[component]], setup)
        fractional_ref_expect_lt(fractional_ref_relative_error(actual_cov, expected_cov),
                                 2e-8, info = paste("constrained", info))
      }
      expected_joint <- expected_joint + expected_cov
    }
    setups <- c(list(setup0), rep(list(setupi), args$m))
    Tc <- as.matrix(Matrix::bdiag(lapply(setups, function(s) s$Tc)))
    A <- do.call(cbind, lapply(setups, function(s) as.matrix(s$A)))
    for (output in list(matrices, rebuilt)) {
      raw <- as.matrix(Matrix::bdiag(output))
      actual_joint <- A %*% solve(Tc %*% raw %*% t(Tc), t(A))
      fractional_ref_expect_lt(fractional_ref_relative_error(actual_joint, expected_joint),
                               2e-8, info = sprintf("joint disconnected alpha=%g, m=%d",
                                                    args$alpha, args$m))
    }
    expect_equal(g$E, edges)
    expect_identical(g$PtV, observed)
    expect_equal(g$edge_lengths, lengths)
  }
})
