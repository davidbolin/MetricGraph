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
