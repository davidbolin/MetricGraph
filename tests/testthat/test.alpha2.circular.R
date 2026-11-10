# alpha = 2 on graphs with circular edges (edges that start and end in the same
# vertex). A circular edge keeps separate states for its two ends, which are
# identified by the Kirchhoff constraints u(0) = u(l), u'(0) = u'(l) at the
# shared vertex. The references are the exact Whittle-Matern covariance on a
# circle, and the equivalent graph where the circular edge is split in two by a
# degree-2 vertex.

# --- helpers ---------------------------------------------------------------

# Exact stationary Whittle-Matern covariance on a circle of circumference L.
C_wm_circle <- function(d, kappa, tau, L, alpha, K = 1e5) {
  k <- -K:K
  vapply(d, function(dd)
    sum((kappa^2 + (2 * pi * k / L)^2)^(-alpha) * cos(2 * pi * k * dd / L)) /
      (tau^2 * L),
    numeric(1))
}

# Unit circle as a single circular edge (one vertex of degree 2).
make_loop_graph <- function(n = 401) {
  th <- seq(0, 2 * pi, length.out = n)
  e <- cbind(cos(th), sin(th))
  e[n, ] <- e[1, ]
  g <- metric_graph$new(edges = list(e), verbose = 0)
  stopifnot(g$nV == 1, g$nE == 1, g$E[1, 1] == g$E[1, 2])
  g
}

# "Lollipop": a line (-1, 0) -> (0, 0) with a unit circle attached at (0, 0).
# With split = FALSE the circle is one circular edge; with split = TRUE it is
# two arcs joined at (2, 0). Both polylines use the same nodes, so the two
# graphs are the same metric graph.
make_lollipop_graph <- function(split) {
  line <- rbind(c(-1, 0), c(0, 0))
  if (!split) {
    th <- seq(pi, 3 * pi, length.out = 401)
    loop <- cbind(1 + cos(th), sin(th))
    loop[c(1, 401), ] <- 0
    edges <- list(line, loop)
  } else {
    th1 <- seq(pi, 2 * pi, length.out = 201)
    th2 <- seq(2 * pi, 3 * pi, length.out = 201)
    a1 <- cbind(1 + cos(th1), sin(th1))
    a2 <- cbind(1 + cos(th2), sin(th2))
    a1[1, ] <- 0
    a2[201, ] <- 0
    a2[1, ] <- a1[201, ]
    edges <- list(line, a1, a2)
  }
  metric_graph$new(edges = edges, verbose = 0)
}

lollipop_data <- function(seed = 2) {
  set.seed(seed)
  ang <- pi + c(0.4, 1.3, 2.2, 3.5, 4.9, 5.8)
  xy <- rbind(cbind(c(-0.8, -0.35), 0), cbind(1 + cos(ang), sin(ang)))
  data.frame(y = rnorm(nrow(xy)), coord_x = xy[, 1], coord_y = xy[, 2])
}


# --- 1. single circular edge vs exact circle covariance --------------------

test_that("alpha=2 precision and constraints on a circular edge give the exact circle covariance", {
  g <- make_loop_graph()
  L <- g$edge_lengths
  kappa <- 1.1; tau <- 0.8

  expect_no_warning(
    Q <- MetricGraph:::Qalpha2(c(tau, kappa), g, BC = 1)
  )
  g$buildC(2, FALSE)
  Tc <- g$CoB$T[-seq_along(g$CoB$S), , drop = FALSE]
  Sig <- as.matrix(t(Tc) %*% solve(Tc %*% Q %*% t(Tc)) %*% Tc)

  # both ends of the edge are the same point
  expect_equal(Sig[1, 1], Sig[3, 3])
  expect_equal(Sig[1, 1], Sig[1, 3])
  expect_equal(Sig[2, 2], Sig[2, 4])
  expect_equal(Sig[1, 1], C_wm_circle(0, kappa, tau, L, alpha = 2),
               tolerance = 1e-8)

  g$build_mesh(h = 0.3)
  circ_all <- c(0, g$mesh$PtE[, 2] * L)
  for (P in list(c(1, 0.37), c(1, 0.81))) {
    Cpkg <- MetricGraph:::spde_covariance(P = P, kappa = kappa, tau = tau,
                                          alpha = 2, graph = g)
    Cex <- C_wm_circle(P[2] * L - circ_all, kappa, tau, L, alpha = 2)
    expect_equal(as.numeric(Cpkg), Cex, tolerance = 1e-7)
  }
})

test_that("alpha=2 likelihoods on a circular edge match a dense Gaussian from the exact covariance", {
  kappa <- 1.1; tau <- 0.8; sigma_e <- 0.3
  theta <- c(log(sigma_e), log(1 / tau), log(kappa))
  t_obs <- c(0.05, 0.21, 0.33, 0.58, 0.71, 0.94)
  set.seed(1)
  df <- data.frame(y = rnorm(length(t_obs)), edge_number = 1,
                   distance_on_edge = t_obs)

  g <- make_loop_graph()
  L <- g$edge_lengths
  g$add_observations(data = df, normalized = TRUE, verbose = 0)
  g$buildC(2, FALSE)

  expect_no_warning(
    ll_pkg <- MetricGraph:::likelihood_alpha2(theta = theta, graph = g,
                 data_name = "y", X_cov = NULL, repl = NULL, BC = 1,
                 parameterization = "spde")
  )
  pc <- MetricGraph:::precompute_alpha2(g, data_name = "y")
  expect_no_warning(
    ll_pc <- MetricGraph:::likelihood_alpha2_precompute(theta = theta,
                 precomputed_data = pc, BC = 1, parameterization = "spde")
  )

  s <- g$get_PtE()[, 2] * L
  y <- g$get_data()[["y"]]
  Sy <- outer(s, s, function(a, b) C_wm_circle(a - b, kappa, tau, L, alpha = 2))
  diag(Sy) <- diag(Sy) + sigma_e^2
  ll_dense <- as.numeric(-0.5 * (length(y) * log(2 * pi) +
                determinant(Sy)$modulus + t(y) %*% solve(Sy, y)))

  expect_equal(as.numeric(ll_pkg), ll_dense, tolerance = 1e-7)
  expect_equal(as.numeric(ll_pc), ll_dense, tolerance = 1e-7)
})


# --- 2. circular edge attached to a line vs split circle -------------------

test_that("alpha=2 likelihood and posterior mean agree for a circular edge and the split circle", {
  kappa <- 1.1; tau <- 0.8; sigma_e <- 0.3
  theta <- c(log(sigma_e), log(1 / tau), log(kappa))
  df <- lollipop_data()

  graphs <- lapply(c(FALSE, TRUE), function(split) {
    g <- make_lollipop_graph(split)
    g$add_observations(data = df, data_coords = "spatial", verbose = 0,
                       tolerance = 0.05)
    g$buildC(2, FALSE)
    g
  })
  g_loop <- graphs[[1]]
  g_split <- graphs[[2]]
  expect_equal(nrow(g_loop$E), 2)
  expect_equal(g_loop$E[2, 1], g_loop$E[2, 2])
  expect_equal(nrow(g_split$E), 3)

  ll <- sapply(graphs, function(g) {
    MetricGraph:::likelihood_alpha2(theta = theta, graph = g, data_name = "y",
                                    X_cov = NULL, repl = NULL, BC = 1,
                                    parameterization = "spde")
  })
  ll_pc <- sapply(graphs, function(g) {
    pc <- MetricGraph:::precompute_alpha2(g, data_name = "y")
    MetricGraph:::likelihood_alpha2_precompute(theta = theta,
                                               precomputed_data = pc, BC = 1,
                                               parameterization = "spde")
  })
  expect_equal(ll[1], ll[2], tolerance = 1e-10)
  expect_equal(ll_pc[1], ll_pc[2], tolerance = 1e-10)
  expect_equal(ll[1], ll_pc[1], tolerance = 1e-10)

  # Edge 1 (the line) and the start of edge 2 (the circle leaving (0, 0) in
  # the same direction) have the same states in both graphs.
  post <- lapply(graphs, function(g) {
    as.numeric(MetricGraph:::posterior_mean_alpha2(
      theta = c(sigma_e, tau, kappa), graph = g,
      resp = g$get_data()[["y"]], PtE_resp = g$get_PtE()))
  })
  expect_equal(post[[1]][1:6], post[[2]][1:6], tolerance = 1e-8)
})


# --- 3. prediction on a circular edge vs exact kriging ---------------------

test_that("alpha=2 prediction on a circular edge matches exact kriging", {
  g <- make_loop_graph()
  L <- g$edge_lengths
  kappa <- 1.1; tau <- 0.8; sigma_e <- 0.3
  t_obs <- seq(0.03, 0.97, length.out = 15)
  u <- simulate(g, seed = 5, alpha = 2, method = "extended", kappa = kappa,
                tau = tau, PtE = cbind(1, t_obs))
  set.seed(6)
  g$add_observations(data = data.frame(y = u + sigma_e * rnorm(length(u)),
                                       edge_number = 1,
                                       distance_on_edge = t_obs),
                     normalized = TRUE, verbose = 0)
  g$buildC(2, FALSE)
  PtE_obs <- g$get_PtE()
  y <- g$get_data()[["y"]]
  t_pred <- c(0, 0.12, 0.4, 0.66, 0.99)

  so <- PtE_obs[, 2] * L
  sp <- t_pred * L
  krige <- function(kappa, tau, sigma_e) {
    C_fun <- function(s1, s2) {
      outer(s1, s2, function(a, b) C_wm_circle(a - b, kappa, tau, L, alpha = 2))
    }
    Soo <- C_fun(so, so)
    diag(Soo) <- diag(Soo) + sigma_e^2
    K <- C_fun(sp, so)
    list(mean = as.numeric(K %*% solve(Soo, y)),
         var = diag(C_fun(sp, sp) - K %*% solve(Soo, t(K))))
  }

  mu <- MetricGraph:::posterior_mean_obs_alpha2(
    c(sigma_e, tau, kappa), graph = g, resp = y, PtE_resp = PtE_obs,
    PtE_pred = cbind(1, t_pred))
  expect_equal(as.numeric(mu), krige(kappa, tau, sigma_e)$mean,
               tolerance = 1e-8)

  # graph_lme + predict, compared with exact kriging at the fitted parameters
  fit <- suppressWarnings(
    graph_lme(y ~ -1, graph = g, model = "WM2", parallel = FALSE)
  )
  ex <- krige(fit$coeff$random_effects[["kappa"]],
              fit$coeff$random_effects[["tau"]],
              fit$coeff$measurement_error[[1]])
  nd <- data.frame(edge_number = 1, distance_on_edge = t_pred)
  p_mean <- predict(fit, newdata = nd, normalized = TRUE)
  p_var <- predict(fit, newdata = nd, normalized = TRUE,
                   compute_variances = TRUE)
  expect_equal(as.numeric(p_mean$mean), ex$mean, tolerance = 1e-8)
  expect_equal(as.numeric(p_var$mean), ex$mean, tolerance = 1e-8)
  expect_equal(as.numeric(p_var$variance), ex$var, tolerance = 1e-8)
})


# --- 4. sampling on a circular edge ----------------------------------------

test_that("alpha=2 samplers on a circular edge reproduce the exact circle covariance", {
  g <- make_loop_graph()
  L <- g$edge_lengths
  kappa <- 1.1; tau <- 0.8

  # the sampled edge state satisfies u(0) = u(l), u'(0) = u'(l)
  set.seed(3)
  b_e <- MetricGraph:::.draw_vertex_state_wm(g, kappa, tau, 2L, BC = 1)
  expect_equal(b_e[1], b_e[3], tolerance = 1e-10)
  expect_equal(b_e[2], b_e[4], tolerance = 1e-10)

  PtE <- cbind(1, c(0.1, 0.3, 0.55, 0.8))
  s <- PtE[, 2] * L
  Sig <- outer(s, s, function(a, b) C_wm_circle(a - b, kappa, tau, L, alpha = 2))
  n <- 2000
  # standard errors of the Gaussian sample covariance entries
  se <- sqrt((outer(diag(Sig), diag(Sig)) + Sig^2) / n)
  max_z <- function(U) max(abs((cov(t(U)) - Sig) / se))

  set.seed(7)
  U <- replicate(n, sample_spde(kappa = kappa, tau = tau, alpha = 2,
                                graph = g, PtE = PtE, type = "manual"))
  expect_lt(max_z(U), 5)

  for (m in c("direct", "kriging", "extended")) {
    U <- simulate(g, nsim = n, seed = 10, alpha = 2, method = m,
                  kappa = kappa, tau = tau, PtE = PtE)
    expect_lt(max_z(U), 5)
  }
})
