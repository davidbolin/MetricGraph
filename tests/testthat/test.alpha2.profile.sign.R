# Regression tests for the alpha = 2 (Whittle-Matern, nu = 3/2) *profiled*
# likelihood code path in R/graph_likelihoods_v2.R.
#
# Background: profile_lik_core_alpha2() (the workhorse behind
# likelihood_alpha2_profile / likelihood_alpha2_profile_precompute) builds
# its own copy of the endpoint-derivative cross-covariance block
#   S[d.index, -d.index] <- r_2(D[1:2, ], ..., deriv = 1)
# independently of likelihood_alpha2_precompute in graph_likelihoods.R. A
# spurious leading minus sign on that line was fixed there (matching the
# sign-convention fix David Bolin made to graph_likelihoods.R / other files
# earlier), but graph_likelihoods_v2.R was added afterwards and silently
# carried the same bug for the profiled path until it was independently
# fixed. The existing "profile equals joint at beta_hat" tests in
# test_profile_likelihood.R only caught this because they cross-check v2
# against the (already-fixed) v1 likelihood_alpha2. These tests instead pin
# likelihood_alpha2_profile directly against an *exact* closed-form
# reference, so they do not depend on likelihood_alpha2's own correctness.
#
# Same circle-graph device as test.alpha2.sign.R: two parallel edges between
# the same two degree-2 vertices, so there is no boundary-condition
# ambiguity, plus a closed-form (truncated Fourier) Whittle-Matern
# covariance to check against. See that file for more background.

# --- helpers (adapted from test.alpha2.sign.R) ------------------------------

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

# Builds a circle graph with asymmetric observations on both edges and
# returns it together with the exact dense reference log-likelihood for the
# given (kappa, tau, sigma_e).
circle_setup <- function(a = 0.9, b = 1.4, seed = 42) {
  L <- a + b
  # asymmetric observation locations on both edges so a sign error cannot
  # cancel by symmetry
  t1 <- c(0.13, 0.37, 0.58, 0.81)
  t2 <- c(0.17, 0.44, 0.72, 0.93)
  set.seed(seed)
  yv <- rnorm(length(t1) + length(t2))
  df <- data.frame(y = yv,
                   edge_number = c(rep(1, length(t1)), rep(2, length(t2))),
                   distance_on_edge = c(t1, t2))

  g <- make_circle_graph(a, b)
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  g$buildC(2, FALSE)

  list(graph = g, a = a, b = b, L = L)
}

dense_reference_loglik <- function(setup, kappa, tau, sigma_e) {
  circ_g <- circle_coord(setup$graph$get_PtE(), setup$a, setup$b)
  y_g <- setup$graph$get_data()[["y"]]
  Sig <- outer(circ_g, circ_g,
              function(s, t) C_wm_circle(s - t, kappa, tau, setup$L, alpha = 2))
  Sy <- Sig; diag(Sy) <- diag(Sy) + sigma_e^2
  n <- length(y_g)
  as.numeric(-0.5 * (n * log(2 * pi) +
              determinant(Sy)$modulus + t(y_g) %*% solve(Sy, y_g)))
}

param_grid <- list(c(kappa = 1.1, tau = 0.8, sigma_e = 0.35),
                   c(kappa = 0.7, tau = 1.3, sigma_e = 0.20))

# --- 1. likelihood_alpha2_profile (no covariates) vs dense Gaussian --------

test_that("alpha=2 profile likelihood with no covariates matches the exact circle covariance", {
  setup <- circle_setup()

  for (par in param_grid) {
    kappa <- par[["kappa"]]; tau <- par[["tau"]]; sigma_e <- par[["sigma_e"]]

    # spde parameterisation: theta = (log sigma_e, log(1/tau), log kappa)
    theta <- c(log(sigma_e), log(1 / tau), log(kappa))

    ll_prof <- MetricGraph:::likelihood_alpha2_profile(theta = theta,
                 graph = setup$graph, data_name = "y", X_cov = NULL,
                 repl = NULL, BC = 1, parameterization = "spde")

    ll_dense <- dense_reference_loglik(setup, kappa, tau, sigma_e)

    expect_equal(as.numeric(ll_prof), ll_dense, tolerance = 1e-7)
  }
})

# --- 2. precompute path agrees with the non-precompute path ----------------

test_that("alpha=2 profile likelihood agrees between precompute and non-precompute paths", {
  setup <- circle_setup()

  for (par in param_grid) {
    kappa <- par[["kappa"]]; tau <- par[["tau"]]; sigma_e <- par[["sigma_e"]]
    theta <- c(log(sigma_e), log(1 / tau), log(kappa))

    ll_prof <- MetricGraph:::likelihood_alpha2_profile(theta = theta,
                 graph = setup$graph, data_name = "y", X_cov = NULL,
                 repl = NULL, BC = 1, parameterization = "spde")

    pc <- MetricGraph:::precompute_alpha2(setup$graph, data_name = "y")
    ll_prof_pc <- MetricGraph:::likelihood_alpha2_profile_precompute(
                 theta = theta, precomputed_data = pc, BC = 1,
                 parameterization = "spde")

    expect_equal(as.numeric(ll_prof_pc), as.numeric(ll_prof), tolerance = 1e-10)

    ll_dense <- dense_reference_loglik(setup, kappa, tau, sigma_e)
    expect_equal(as.numeric(ll_prof_pc), ll_dense, tolerance = 1e-7)
  }
})

# --- 3. cross-check with covariates: profile matches joint at beta_hat -----

test_that("alpha=2 profile likelihood with covariates matches the joint likelihood at beta_hat on the circle graph", {
  setup <- circle_setup()
  n <- nrow(setup$graph$get_PtE())
  set.seed(99)
  X_cov <- cbind(1, rnorm(n))

  for (par in param_grid) {
    kappa <- par[["kappa"]]; tau <- par[["tau"]]; sigma_e <- par[["sigma_e"]]
    theta <- c(log(sigma_e), log(1 / tau), log(kappa))

    bhat <- MetricGraph:::profile_beta_estimate(theta, model = "alpha2",
                 graph = setup$graph, data_name = "y", X_cov = X_cov,
                 repl = NULL, BC = 1, parameterization = "spde")$beta

    ll_prof <- MetricGraph:::likelihood_alpha2_profile(theta = theta,
                 graph = setup$graph, data_name = "y", X_cov = X_cov,
                 repl = NULL, BC = 1, parameterization = "spde")

    ll_joint <- MetricGraph:::likelihood_alpha2(theta = c(theta, bhat),
                 graph = setup$graph, data_name = "y", X_cov = X_cov,
                 repl = NULL, BC = 1, parameterization = "spde")

    expect_equal(as.numeric(ll_prof), as.numeric(ll_joint), tolerance = 1e-8)
  }
})
