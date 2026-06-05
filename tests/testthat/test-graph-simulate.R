## Tests for simulate.metric_graph, simulate_parallel, and C++ equivalence

# ---------------------------------------------------------------------------
# Helper: a small planar square graph (4 vertices, 4 edges, edge length = 1)
# ---------------------------------------------------------------------------
make_square_graph <- function() {
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  suppressMessages(metric_graph$new(V = V, E = E))
}

# 100 points spread over the square graph
make_PtE <- function(g, n = 5) {
  t_norm <- seq(0.1, 0.9, length.out = n)
  do.call(rbind, lapply(seq_len(g$nE), function(e) cbind(e, t_norm)))
}

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("bad alpha raises error", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  expect_error(simulate(g, alpha = 3, method = "direct",
                        kappa = 1, tau = 1, PtE = PtE))
})

test_that("bad method raises error", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  expect_error(simulate(g, alpha = 1, method = "foobar",
                        kappa = 1, tau = 1, PtE = PtE))
})

test_that("missing kappa/tau and range/sigma raises error", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  expect_error(simulate(g, alpha = 1, method = "direct", PtE = PtE))
})

test_that("bad type raises error", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  expect_error(simulate(g, alpha = 1, method = "direct",
                        kappa = 1, tau = 1, PtE = PtE, type = "bogus"))
})

test_that("manual type without PtE raises error", {
  g <- make_square_graph()
  expect_error(simulate(g, alpha = 1, method = "direct",
                        kappa = 1, tau = 1, type = "manual"))
})

# ---------------------------------------------------------------------------
# 2. Output shape
# ---------------------------------------------------------------------------

test_that("nsim = 1 returns a vector", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  u   <- simulate(g, alpha = 1, method = "direct",
                  kappa = 1, tau = 1, PtE = PtE, seed = 1)
  expect_true(is.numeric(u) && is.null(dim(u)))
  expect_equal(length(u), nrow(PtE))
})

test_that("nsim > 1 returns a matrix with nsim columns", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  U   <- simulate(g, nsim = 4, alpha = 1, method = "direct",
                  kappa = 1, tau = 1, PtE = PtE, seed = 1)
  expect_true(is.matrix(U))
  expect_equal(dim(U), c(nrow(PtE), 4L))
})

test_that("range/sigma parameterisation gives same result as equivalent kappa/tau", {
  g    <- make_square_graph()
  PtE  <- make_PtE(g)
  nu   <- 0.5
  rng  <- 1.5
  sig  <- 2.0
  kap  <- sqrt(8 * nu) / rng
  tau_ <- sqrt(gamma(nu) / (sig^2 * kap^(2 * nu) *
                              (4 * pi)^0.5 * gamma(nu + 0.5)))
  set.seed(7); u1 <- simulate(g, alpha = 1, method = "direct",
                               kappa = kap, tau = tau_, PtE = PtE)
  set.seed(7); u2 <- simulate(g, alpha = 1, method = "direct",
                               range = rng, sigma = sig, PtE = PtE)
  expect_equal(u1, u2, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. Method equivalence — same distribution (empirical variance)
# ---------------------------------------------------------------------------

test_that("direct and kriging have equal marginal variance: alpha = 1", {
  skip_on_cran()
  g   <- make_square_graph()
  PtE <- make_PtE(g, n = 3)
  n_s <- 2000L
  kap <- 1; tau_ <- 1
  r0  <- 1 / (2 * kap * tau_^2)

  set.seed(10)
  UA <- do.call(cbind, lapply(seq_len(n_s), function(s)
    simulate(g, alpha = 1, method = "direct",  kappa = kap, tau = tau_, PtE = PtE)))
  set.seed(10)
  UB <- do.call(cbind, lapply(seq_len(n_s), function(s)
    simulate(g, alpha = 1, method = "kriging", kappa = kap, tau = tau_, PtE = PtE)))

  varA <- rowMeans(UA^2)
  varB <- rowMeans(UB^2)
  # Theoretical marginal variance ≈ r_1(0, kappa, tau) for interior points
  expect_equal(varA, rep(r0, length(varA)), tolerance = 0.1)
  expect_equal(varB, rep(r0, length(varB)), tolerance = 0.1)
  # Also check that both methods agree empirically
  expect_equal(mean(varA), mean(varB), tolerance = 0.05)
})

test_that("direct and kriging have equal marginal variance: alpha = 2", {
  skip_on_cran()
  g   <- make_square_graph()
  PtE <- make_PtE(g, n = 3)
  n_s <- 2000L
  kap <- 1; tau_ <- 1
  # r_2(0, kap, tau_) = 1/(4*kap^3*tau_^2)
  r0  <- 1 / (4 * kap^3 * tau_^2)

  set.seed(20)
  UA <- do.call(cbind, lapply(seq_len(n_s), function(s)
    simulate(g, alpha = 2, method = "direct",  kappa = kap, tau = tau_, PtE = PtE)))
  set.seed(20)
  UB <- do.call(cbind, lapply(seq_len(n_s), function(s)
    simulate(g, alpha = 2, method = "kriging", kappa = kap, tau = tau_, PtE = PtE)))

  varA <- rowMeans(UA^2)
  varB <- rowMeans(UB^2)
  expect_equal(varA, rep(r0, length(varA)), tolerance = 0.1)
  expect_equal(varB, rep(r0, length(varB)), tolerance = 0.1)
  expect_equal(mean(varA), mean(varB), tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# 4. Seed reproducibility
# ---------------------------------------------------------------------------

test_that("same seed gives identical output (serial, alpha = 1)", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  u1  <- simulate(g, alpha = 1, method = "direct",  kappa = 1, tau = 1, PtE = PtE, seed = 99)
  u2  <- simulate(g, alpha = 1, method = "direct",  kappa = 1, tau = 1, PtE = PtE, seed = 99)
  expect_equal(u1, u2)
})

test_that("same seed gives identical output (serial, alpha = 2)", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  u1  <- simulate(g, alpha = 2, method = "kriging", kappa = 1, tau = 1, PtE = PtE, seed = 7)
  u2  <- simulate(g, alpha = 2, method = "kriging", kappa = 1, tau = 1, PtE = PtE, seed = 7)
  expect_equal(u1, u2)
})

# ---------------------------------------------------------------------------
# 5. Serial vs parallel reproducibility
# ---------------------------------------------------------------------------

test_that("serial and parallel give identical output for the same seed", {
  skip_on_cran()
  skip_on_os("windows")   # fork-based parallel not available on Windows CRAN
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  u_ser <- simulate(g, alpha = 1, method = "kriging",
                    kappa = 1, tau = 1, PtE = PtE,
                    seed = 42, parallel = FALSE)
  u_par <- simulate_parallel(g, alpha = 1, method = "kriging",
                              kappa = 1, tau = 1, PtE = PtE,
                              seed = 42, n_cores = 2)
  expect_equal(u_ser, u_par)
})

test_that("simulate_parallel is a wrapper around simulate.metric_graph", {
  skip_on_cran()
  skip_on_os("windows")
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  u1  <- simulate_parallel(g, alpha = 2, method = "direct",
                            kappa = 1, tau = 1, PtE = PtE, seed = 5, n_cores = 2)
  u2  <- simulate(g, alpha = 2, method = "direct",
                  kappa = 1, tau = 1, PtE = PtE, seed = 5, parallel = TRUE, n_cores = 2)
  expect_equal(u1, u2)
})

# ---------------------------------------------------------------------------
# 6. PtE ordering contract — output aligned to input row order
# ---------------------------------------------------------------------------

test_that("output rows match PtE input ordering", {
  g    <- make_square_graph()
  PtE  <- make_PtE(g)
  # Shuffle PtE deliberately
  idx  <- sample(nrow(PtE))
  PtE_s <- PtE[idx, ]
  u_s   <- simulate(g, alpha = 1, method = "direct",
                    kappa = 1, tau = 1, PtE = PtE_s, seed = 123)
  u_r   <- simulate(g, alpha = 1, method = "direct",
                    kappa = 1, tau = 1, PtE = PtE,   seed = 123)
  # u_s[i] should equal u_r[idx[i]]
  expect_equal(u_s, u_r[idx], tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 7. R vs C++ numerical equivalence (same seed → same draws)
# ---------------------------------------------------------------------------

test_that("R and C++ draw_edge_direct match: alpha = 1", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  set.seed(5)
  u_R <- simulate(g, alpha = 1, method = "direct", impl = "R",
                  kappa = 1, tau = 1, PtE = PtE)
  set.seed(5)
  u_C <- simulate(g, alpha = 1, method = "direct", impl = "cpp",
                  kappa = 1, tau = 1, PtE = PtE)
  expect_equal(u_R, u_C, tolerance = 1e-10)
})

test_that("R and C++ draw_edge_direct match: alpha = 2", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  set.seed(6)
  u_R <- simulate(g, alpha = 2, method = "direct", impl = "R",
                  kappa = 1, tau = 1, PtE = PtE)
  set.seed(6)
  u_C <- simulate(g, alpha = 2, method = "direct", impl = "cpp",
                  kappa = 1, tau = 1, PtE = PtE)
  expect_equal(u_R, u_C, tolerance = 1e-10)
})

test_that("R and C++ draw_edge_kriging match: alpha = 1", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  set.seed(7)
  u_R <- simulate(g, alpha = 1, method = "kriging", impl = "R",
                  kappa = 1, tau = 1, PtE = PtE)
  set.seed(7)
  u_C <- simulate(g, alpha = 1, method = "kriging", impl = "cpp",
                  kappa = 1, tau = 1, PtE = PtE)
  expect_equal(u_R, u_C, tolerance = 1e-10)
})

test_that("R and C++ draw_edge_kriging match: alpha = 2", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  set.seed(8)
  u_R <- simulate(g, alpha = 2, method = "kriging", impl = "R",
                  kappa = 1, tau = 1, PtE = PtE)
  set.seed(8)
  u_C <- simulate(g, alpha = 2, method = "kriging", impl = "cpp",
                  kappa = 1, tau = 1, PtE = PtE)
  expect_equal(u_R, u_C, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 8. method = "extended" — basic shape and distributional checks
# ---------------------------------------------------------------------------

test_that("extended method returns correct shape: alpha = 1", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  u   <- simulate(g, alpha = 1, method = "extended",
                  kappa = 1, tau = 1, PtE = PtE, seed = 11)
  expect_true(is.numeric(u) && is.null(dim(u)))
  expect_equal(length(u), nrow(PtE))
})

test_that("extended method returns correct shape: alpha = 2", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  u   <- simulate(g, alpha = 2, method = "extended",
                  kappa = 1, tau = 1, PtE = PtE, seed = 12)
  expect_true(is.numeric(u) && is.null(dim(u)))
  expect_equal(length(u), nrow(PtE))
})

test_that("extended and direct have equal marginal variance: alpha = 1", {
  skip_on_cran()
  g   <- make_square_graph()
  PtE <- make_PtE(g, n = 3)
  n_s <- 2000L
  kap <- 1; tau_ <- 1
  r0  <- 1 / (2 * kap * tau_^2)

  set.seed(30)
  UA <- simulate(g, nsim = n_s, alpha = 1, method = "direct",
                 kappa = kap, tau = tau_, PtE = PtE)
  set.seed(30)
  UE <- simulate(g, nsim = n_s, alpha = 1, method = "extended",
                 kappa = kap, tau = tau_, PtE = PtE)

  expect_equal(mean(rowMeans(UA^2)), r0, tolerance = 0.1)
  expect_equal(mean(rowMeans(UE^2)), r0, tolerance = 0.1)
  expect_equal(mean(rowMeans(UA^2)), mean(rowMeans(UE^2)), tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# 9. impl parameter validation
# ---------------------------------------------------------------------------

test_that("bad impl raises error", {
  g   <- make_square_graph()
  PtE <- make_PtE(g)
  expect_error(simulate(g, alpha = 1, method = "direct",
                        kappa = 1, tau = 1, PtE = PtE, impl = "fortran"))
})
