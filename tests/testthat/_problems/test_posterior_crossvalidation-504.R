# Extracted from test_posterior_crossvalidation.R:504

# prequel ----------------------------------------------------------------------
.loo_reference <- function(fit) {
  graph <- fit$graph$clone()
  response_name <- as.character(fit$response_var)
  data <- graph$.__enclos_env__$private$data
  y_all   <- data[[response_name]]
  grp_all <- data[[".group"]]

  sigma_e <- fit$coeff$measurement_error
  tau     <- fit$coeff$random_effects[1]
  kappa   <- fit$coeff$random_effects[2]

  # Build the full covariance matrix at the (within-replicate) observation
  # locations. The graph is the same for every replicate, so we build it
  # once and re-use it per replicate.
  alpha <- fit$latent_model$alpha
  g_v <- graph$clone()
  g_v$observation_to_vertex()
  if (alpha == 1) {
    Q <- MetricGraph:::Qalpha1(c(tau, kappa), g_v, BC = fit$BC)
    Sigma <- as.matrix(solve(Q))[g_v$PtV, g_v$PtV]
  } else if (alpha == 2) {
    g_v$buildC(2, fit$BC == 0)
    n.c <- 1:length(g_v$CoB$S)
    Q <- MetricGraph:::Qalpha2(c(tau, kappa), g_v, BC = fit$BC)
    Qtilde <- (g_v$CoB$T) %*% Q %*% t(g_v$CoB$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.full <- t(g_v$CoB$T[-n.c, ]) %*% solve(Qtilde) %*% (g_v$CoB$T[-n.c, ])
    PtE <- g_v$get_PtE()
    index.obs <- 4 * (PtE[, 1] - 1) +
      (1 * (PtE[, 2] == 0)) + (3 * (PtE[, 2] != 0))
    Sigma <- as.matrix(Sigma.full[index.obs, index.obs])
  } else {
    stop("alpha must be 1 or 2 for the reference")
  }
  Sigma_o <- Sigma
  diag(Sigma_o) <- diag(Sigma_o) + sigma_e^2

  # Per-replicate LOO using the same Sigma (graph is shared) but only the
  # rows/columns of observations in that replicate.
  n_total <- length(y_all)
  mu  <- rep(NA_real_, n_total)
  var <- rep(NA_real_, n_total)
  for (r in unique(grp_all)) {
    pos_r <- which(grp_all == r)
    # Within-replicate observations sit at the first length(pos_r) entries
    # of Sigma — every replicate uses the same set of observation
    # locations on the (shared) graph.
    nr <- length(pos_r)
    S_r   <- Sigma[seq_len(nr), seq_len(nr), drop = FALSE]
    So_r  <- Sigma_o[seq_len(nr), seq_len(nr), drop = FALSE]
    y_r   <- y_all[pos_r]
    for (i in seq_len(nr)) {
      mu[pos_r[i]]  <- S_r[i, -i, drop = FALSE] %*%
                         solve(So_r[-i, -i, drop = FALSE], y_r[-i])
      var[pos_r[i]] <- So_r[i, i] - S_r[i, -i, drop = FALSE] %*%
                         solve(So_r[-i, -i, drop = FALSE],
                               t(S_r[i, -i, drop = FALSE]))
    }
  }
  list(mu = mu, var = var)
}
.fit_small <- function(alpha = 1, n_repl = 1, seed = 1L,
                       obs_per_edge = 4L) {
  set.seed(seed)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)

  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, obs_per_edge), runif(obs_per_edge)))
  }

  kappa <- 5; tau <- 1; sigma_e <- 0.3
  df_all <- NULL
  for (r in seq_len(n_repl)) {
    u <- sample_spde(kappa = kappa, tau = tau, alpha = alpha,
                     graph = graph, PtE = PtE)
    y <- u + sigma_e * rnorm(length(u))
    df_all <- rbind(df_all,
      data.frame(y = y, edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
  group_var <- if (n_repl > 1) "repl" else NULL
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = group_var)

  # No intercept so that the reference (which treats y as zero-mean)
  # matches posterior_crossvalidation byte-for-byte.
  graph_lme(y ~ -1, graph = graph,
            model = list(type = "WhittleMatern", alpha = alpha))
}

# test -------------------------------------------------------------------------
skip_if_not_installed("rSPDE")
set.seed(2)
V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
graph <- metric_graph$new(V = V, E = E, verbose = 0)
graph$build_mesh(h = 0.1)
PtE <- NULL
for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, 4), runif(4)))
  }
df_all <- NULL
for (r in 1:2) {
    u <- sample_spde(kappa = 5, tau = 1, alpha = 1, graph = graph, PtE = PtE)
    y <- u + 0.3 * rnorm(length(u))
    df_all <- rbind(df_all,
      data.frame(y = y, edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2], repl = r))
  }
graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = "repl")
fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE))
expect_true(inherits(fit, "rspde_lme"))
i_na <- 5
cv_model <- MetricGraph:::update_graph_lme_with_na(fit, i_na)
gd <- fit$graph$.__enclos_env__$private$data
nd <- list(
    .edge_number      = gd[[".edge_number"]][i_na],
    .distance_on_edge = gd[[".distance_on_edge"]][i_na],
    .group            = "1",
    y                 = NA_real_
  )
pred <- predict(cv_model,
                  newdata = nd,
                  edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge",
                  normalized = TRUE, compute_variances = TRUE)
