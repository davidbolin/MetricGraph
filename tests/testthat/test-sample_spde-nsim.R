## Regression tests for sample_spde() with nsim > 1 (method = "conditional").
## The recursive nsim > 1 branch used to forward range/sigma even when only
## kappa/tau were supplied, failing with 'argument "range" is missing'.

make_line_graph <- function() {
  metric_graph$new(edges = list(rbind(c(0, 0), c(1, 0))), verbose = 0)
}

make_square_graph <- function() {
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  suppressMessages(metric_graph$new(V = V, E = E, verbose = 0))
}

kappa_tau_from_range_sigma <- function(range, sigma, alpha) {
  nu <- alpha - 0.5
  kappa <- sqrt(8 * nu) / range
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
                             (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
  list(kappa = kappa, tau = tau)
}

for (alpha in 1:2) {

  test_that(paste0("sample_spde nsim > 1 works with kappa/tau: alpha = ", alpha), {
    g <- make_line_graph()
    PtE <- cbind(1, c(0.2, 0.5))
    u <- sample_spde(kappa = 1, tau = 1, alpha = alpha, graph = g,
                     PtE = PtE, type = "manual", nsim = 10)
    expect_true(is.matrix(u))
    expect_equal(dim(u), c(nrow(PtE), 10L))
    expect_true(all(is.finite(u)))
  })

  test_that(paste0("sample_spde nsim > 1 works with range/sigma: alpha = ", alpha), {
    g <- make_line_graph()
    PtE <- cbind(1, c(0.2, 0.5))
    u <- sample_spde(range = 1, sigma = 1, alpha = alpha, graph = g,
                     PtE = PtE, type = "manual", nsim = 10)
    expect_true(is.matrix(u))
    expect_equal(dim(u), c(nrow(PtE), 10L))
    expect_true(all(is.finite(u)))
  })

  test_that(paste0("sample_spde nsim > 1 matches repeated nsim = 1 draws: alpha = ", alpha), {
    g <- make_square_graph()
    PtE <- cbind(1:4, c(0.2, 0.5, 0.7, 0.3))
    nsim <- 3
    # BC = 0 is non-default, so this also checks that BC is forwarded.
    set.seed(1)
    U <- sample_spde(kappa = 2, tau = 0.5, alpha = alpha, graph = g,
                     PtE = PtE, type = "manual", nsim = nsim, BC = 0)
    set.seed(1)
    U_ref <- sapply(seq_len(nsim), function(i) {
      sample_spde(kappa = 2, tau = 0.5, alpha = alpha, graph = g,
                  PtE = PtE, type = "manual", nsim = 1, BC = 0)
    })
    expect_equal(U, U_ref)
  })

  test_that(paste0("sample_spde nsim > 1: range/sigma equals equivalent kappa/tau: alpha = ", alpha), {
    g <- make_square_graph()
    PtE <- cbind(1:4, c(0.2, 0.5, 0.7, 0.3))
    par <- kappa_tau_from_range_sigma(range = 1.5, sigma = 2, alpha = alpha)
    set.seed(2)
    U_rs <- sample_spde(range = 1.5, sigma = 2, alpha = alpha, graph = g,
                        PtE = PtE, type = "manual", nsim = 4)
    set.seed(2)
    U_kt <- sample_spde(kappa = par$kappa, tau = par$tau, alpha = alpha,
                        graph = g, PtE = PtE, type = "manual", nsim = 4)
    expect_equal(U_rs, U_kt, tolerance = 1e-12)
  })
}
