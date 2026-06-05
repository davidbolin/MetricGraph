# Extracted from test_fem_predict_manual.R:1021

# prequel ----------------------------------------------------------------------
.fem_manual_predict <- function(fit, loc_pred, which_repl = NULL,
                                na_idx = NULL) {
  sigma_e <- as.numeric(fit$coeff$measurement_error)

  coeff_random <- fit$coeff$random_effects
  update_params <- list()
  for (param_name in names(coeff_random)) {
    update_params[[gsub(" \\(fixed\\)$", "", param_name)]] <-
      coeff_random[[param_name]]
  }
  update_params$check_stationarity <- FALSE
  new_obj <- do.call(update, c(list(fit$latent_model), update_params))
  Q <- new_obj$Q

  alpha_param <- grep("^alpha", names(coeff_random), value = TRUE)
  nu_param    <- grep("^nu",    names(coeff_random), value = TRUE)
  if (length(alpha_param) > 0) {
    alpha_val <- coeff_random[[alpha_param]]
  } else {
    alpha_val <- coeff_random[[nu_param]] + fit$latent_model$d / 2
  }
  is_integer_alpha <- (alpha_val %% 1 == 0)

  A_pred <- rSPDE::make_A(new_obj, loc_pred)
  if (!is_integer_alpha) {
    A_pred <- kronecker(matrix(1, 1, fit$rspde_order + 1), A_pred)
  }
  mu_corr <- if (fit$mean_correction) new_obj$mean_correction(full = TRUE) else 0

  mm <- as.matrix(fit$model_matrix)
  Y_all <- mm[, 1]
  if (ncol(mm) > 1) {
    X_fit <- mm[, 2:ncol(mm), drop = FALSE]
    Y_all <- Y_all - as.numeric(X_fit %*% fit$coeff$fixed_effects)
  }
  Y_all <- as.numeric(Y_all)
  if (!is.null(na_idx)) Y_all[na_idx] <- NA

  repl_vec <- fit$repl
  u_repl <- if (is.null(which_repl)) unique(repl_vec) else unique(which_repl)

  out_mean <- list(); out_var <- list()
  for (r in u_repl) {
    idx_r    <- repl_vec == r
    y_r      <- Y_all[idx_r]
    obs_mask <- !is.na(y_r)
    y_r      <- y_r[obs_mask]
    A_r      <- fit$A_list[[r]]
    if (nrow(A_r) == length(obs_mask) && any(!obs_mask)) {
      A_r <- A_r[obs_mask, , drop = FALSE]
    }
    if (fit$mean_correction) y_r <- y_r - as.numeric(A_r %*% mu_corr)

    Q_xy   <- Matrix::t(A_r) %*% A_r / sigma_e^2 + Q
    mu_lat <- solve(Q_xy, as.vector(Matrix::t(A_r) %*% y_r / sigma_e^2))
    mu_pp  <- as.numeric(A_pred %*% mu_lat)

    fe_pred <- 0
    if (ncol(mm) > 1) {
      # Fixed-effects intercept is in fit$coeff$fixed_effects; we built
      # loc_pred without covariates so the intercept (if any) is the
      # only fe contribution. Match predict.rspde_lme's behaviour of
      # adding X_pred %*% beta with X_pred = matrix(1, nrow=n_pred, 1)
      # when no covariate values were provided alongside the loc.
      fe_pred <- as.numeric(fit$coeff$fixed_effects[1])
    }
    if (fit$mean_correction) mu_pp <- mu_pp + as.numeric(A_pred %*% mu_corr)

    out_mean[[as.character(r)]] <- mu_pp + fe_pred
    out_var[[as.character(r)]]  <-
      pmax(diag(as.matrix(A_pred %*% solve(Q_xy, Matrix::t(A_pred)))), 0)
  }
  list(mean = out_mean, variance = out_var)
}
.build_fem_multirep_graph <- function(n_repl, alpha = 1, obs_per_edge = 5,
                                      seed = 11, mesh_h = 0.05,
                                      shared_locations = TRUE) {
  set.seed(seed)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  graph$build_mesh(h = mesh_h)
  df_all <- NULL
  if (shared_locations) {
    PtE <- NULL
    for (i in seq_len(graph$nE)) {
      PtE <- rbind(PtE, cbind(rep(i, obs_per_edge), runif(obs_per_edge)))
    }
  }
  for (r in seq_len(n_repl)) {
    if (!shared_locations) {
      PtE <- NULL
      for (i in seq_len(graph$nE)) {
        PtE <- rbind(PtE, cbind(rep(i, obs_per_edge), runif(obs_per_edge)))
      }
    }
    u <- sample_spde(kappa = 5, tau = 1, alpha = alpha, graph = graph, PtE = PtE)
    df_all <- rbind(df_all,
      data.frame(y = u + 0.3 * rnorm(length(u)),
                 edge_number = PtE[, 1],
                 distance_on_edge = PtE[, 2],
                 repl = r))
  }
  group_var <- if (n_repl > 1L) "repl" else NULL
  graph$add_observations(data = df_all, normalized = TRUE, verbose = 0,
                         group = group_var)
  graph
}

# test -------------------------------------------------------------------------
skip_if_not_installed("rSPDE")
graph <- .build_fem_multirep_graph(n_repl = 3L, alpha = 1, seed = 11)
n_mesh <- nrow(graph$mesh$VtE)
spatial_cov <- graph$mesh$VtE[, 2] - 0.5
B.range <- cbind(rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh), spatial_cov)
B.sigma <- cbind(rep(0, n_mesh), rep(0, n_mesh), rep(1, n_mesh), rep(0, n_mesh))
fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1, fem = TRUE,
                                B.range = B.range, B.sigma = B.sigma))
expect_true(inherits(fit, "rspde_lme"))
expect_true(any(grepl("^theta", names(fit$coeff$random_effects))))
res_pre <- posterior_crossvalidation(fit, mode = "loo", true_CV = FALSE,
                                       use_precomputed = TRUE,
                                       scores = c("logscore", "rmse"))
