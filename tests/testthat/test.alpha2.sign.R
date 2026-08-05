# Sign-convention regression tests for the alpha = 2 (Whittle-Matern, nu = 3/2)
# code paths.
#
# Background: the observation-/bridge-side covariance blocks used to treat the
# endpoint-derivative states as -u' while the prior precision (Q00 / Qalpha2)
# and the vertex constraints treat them as +u'. Because both enter the same
# quadratic forms, the alpha = 2 covariances, likelihoods, posteriors and exact
# samples deviated from the true model. The blocks were made consistent with the
# plain +u' convention (see the sign-convention fix). These tests pin the
# behaviour to the exact Whittle-Matern covariance on a circle, where a closed
# form / Fourier reference is available and no boundary-condition ambiguity
# enters (both vertices have degree 2).

# --- helpers ---------------------------------------------------------------

# Two parallel edges e1 = (v1, v2) of length a and e2 = (v2, v1) of length b
# form a circle of circumference L = a + b. e2 is drawn as a bent polyline so
# that its geometric length is b (!= a). Both vertices have degree 2.
make_circle_graph <- function(a, b) {
  stopifnot(b > a)
  yb <- sqrt((b / 2)^2 - (a / 2)^2)
  e1 <- rbind(c(0, 0), c(a, 0))
  e2 <- rbind(c(a, 0), c(a / 2, yb), c(0, 0))
  g <- metric_graph$new(edges = list(e1, e2), verbose = 0)
  stopifnot(g$nV == 2, g$nE == 2)
  stopifnot(max(abs(sort(g$edge_lengths) - sort(c(a, b)))) < 1e-8)
  g
}

# Circle coordinate of a point (edge, t_norm): e1 is parameterised v1 -> v2 so
# the coordinate is t_norm * a; e2 is parameterised v2 -> v1 so it is a + t_norm * b.
circle_coord <- function(PtE, a, b) ifelse(PtE[, 1] == 1, PtE[, 2] * a, a + PtE[, 2] * b)

# Exact stationary Whittle-Matern covariance on a circle of circumference L,
# spectral density proportional to (kappa^2 + (2 pi k / L)^2)^(-alpha).
# Terms decay like k^(-2 alpha); K = 1e5 is far into the tail for alpha = 2.
C_wm_circle <- function(d, kappa, tau, L, alpha, K = 1e5) {
  k <- -K:K
  vapply(d, function(dd)
    sum((kappa^2 + (2 * pi * k / L)^2)^(-alpha) * cos(2 * pi * k * dd / L)) /
      (tau^2 * L),
    numeric(1))
}

# Exact OU-on-circle covariance (alpha = 1) in closed form (no truncation).
C_ou_circle <- function(d, kappa, tau, L) {
  dd <- d %% L
  dd <- pmin(dd, L - dd)
  cosh(kappa * (L / 2 - dd)) / (2 * kappa * tau^2 * sinh(kappa * L / 2))
}


# --- 1. spde_covariance (alpha = 2) vs exact Fourier series ----------------

test_that("alpha=2 spde_covariance matches the exact circle covariance", {
  a <- 0.9; b <- 1.4; L <- a + b
  g <- make_circle_graph(a, b)
  g$build_mesh(h = 0.25)
  PtE <- g$mesh$PtE
  # spde_covariance returns [v1, v2, interior mesh points in PtE order]
  circ_all <- c(0, a, circle_coord(PtE, a, b))

  for (par in list(c(kappa = 1.1, tau = 0.8), c(kappa = 0.7, tau = 1.3))) {
    kappa <- par[["kappa"]]; tau <- par[["tau"]]
    for (P in list(c(1, 0.4), c(2, 0.65))) {
      Pcirc <- circle_coord(matrix(P, 1), a, b)
      Cpkg <- MetricGraph:::spde_covariance(P = P, kappa = kappa, tau = tau,
                                            alpha = 2, graph = g)
      Cex <- C_wm_circle(Pcirc - circ_all, kappa, tau, L, alpha = 2)
      expect_equal(as.numeric(Cpkg), Cex, tolerance = 1e-7)
    }
  }
})


# --- 2. likelihood_alpha2 (+ precomputed variant) vs dense Gaussian --------

test_that("alpha=2 likelihood matches a dense Gaussian from the exact covariance", {
  a <- 0.9; b <- 1.4; L <- a + b
  # asymmetric observation locations on both edges so a sign error cannot
  # cancel by symmetry
  t1 <- c(0.13, 0.37, 0.58, 0.81)
  t2 <- c(0.17, 0.44, 0.72, 0.93)
  set.seed(42)
  yv <- rnorm(length(t1) + length(t2))
  df <- data.frame(y = yv,
                   edge_number = c(rep(1, length(t1)), rep(2, length(t2))),
                   distance_on_edge = c(t1, t2))

  for (par in list(c(kappa = 1.1, tau = 0.8, sigma_e = 0.35),
                   c(kappa = 0.7, tau = 1.3, sigma_e = 0.20))) {
    kappa <- par[["kappa"]]; tau <- par[["tau"]]; sigma_e <- par[["sigma_e"]]

    g <- make_circle_graph(a, b)
    g$add_observations(data = df, normalized = TRUE, verbose = 0)
    g$buildC(2, FALSE)

    # spde parameterisation: theta = (log sigma_e, log(1/tau), log kappa)
    theta <- c(log(sigma_e), log(1 / tau), log(kappa))

    ll_pkg <- MetricGraph:::likelihood_alpha2(theta = theta, graph = g,
                 data_name = "y", X_cov = NULL, repl = NULL, BC = 1,
                 parameterization = "spde")

    pc <- MetricGraph:::precompute_alpha2(g, data_name = "y")
    ll_pc <- MetricGraph:::likelihood_alpha2_precompute(theta = theta,
                 precomputed_data = pc, BC = 1, parameterization = "spde")

    # dense reference: y ~ N(0, Sigma + sigma_e^2 I) at the observation coords,
    # ordered to match graph$get_PtE() / graph$get_data().
    circ_g <- circle_coord(g$get_PtE(), a, b)
    y_g <- g$get_data()[["y"]]
    Sig <- outer(circ_g, circ_g,
                 function(s, t) C_wm_circle(s - t, kappa, tau, L, alpha = 2))
    Sy <- Sig; diag(Sy) <- diag(Sy) + sigma_e^2
    n <- length(y_g)
    ll_dense <- as.numeric(-0.5 * (n * log(2 * pi) +
                  determinant(Sy)$modulus + t(y_g) %*% solve(Sy, y_g)))

    expect_equal(as.numeric(ll_pkg), ll_dense, tolerance = 1e-7)
    expect_equal(as.numeric(ll_pc),  ll_dense, tolerance = 1e-7)
  }
})


# --- 3. alpha = 1 regression: unaffected by the fix ------------------------

test_that("alpha=1 spde_covariance is unchanged (matches exact OU on circle)", {
  a <- 0.9; b <- 1.4; L <- a + b
  g <- make_circle_graph(a, b)
  g$build_mesh(h = 0.25)
  circ_all <- c(0, a, circle_coord(g$mesh$PtE, a, b))

  for (par in list(c(kappa = 1.1, tau = 0.8), c(kappa = 0.7, tau = 1.3))) {
    kappa <- par[["kappa"]]; tau <- par[["tau"]]
    Pcirc <- circle_coord(matrix(c(1, 0.4), 1), a, b)
    Cpkg <- MetricGraph:::spde_covariance(P = c(1, 0.4), kappa = kappa,
                                          tau = tau, alpha = 1, graph = g)
    Cex <- C_ou_circle(Pcirc - circ_all, kappa, tau, L)
    expect_equal(as.numeric(Cpkg), Cex, tolerance = 1e-8)
  }
})
