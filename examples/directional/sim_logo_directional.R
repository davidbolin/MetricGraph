rm(list = ls())
library(MetricGraph)
library(ggplot2)

set.seed(2024)

nrepl         <- 30      # Monte-Carlo replicates per scenario
obs_per_edge  <- 5       # training observations per edge
test_per_edge <- 5       # held-out test locations per edge
sigma_true    <- 1
range_true    <- 1.0     # practical correlation range
sigma_e_true  <- 0.1
alpha         <- 1       # only alpha = 1 is supported in the directional model

graph <- metric_graph$new(edges           = logo_lines(),
                          perform_merges  = TRUE,
                          remove_circles  = TRUE,
                          verbose         = 0)

## Training and test locations are fixed across replicates; only the field
## (and noise) change. all_loc stacks train rows on top of test rows so we
## can draw a single joint sample of u from sample_spde() and split it.
obs_loc  <- do.call(rbind, lapply(seq_len(graph$nE), function(e)
                                  cbind(e, sort(runif(obs_per_edge)))))
test_loc <- do.call(rbind, lapply(seq_len(graph$nE), function(e)
                                  cbind(e, sort(runif(test_per_edge)))))
all_loc  <- rbind(obs_loc, test_loc)
n_obs    <- nrow(obs_loc)
n_test   <- nrow(test_loc)
cat(sprintf("Using %d training and %d test locations.\n", n_obs, n_test))

## Helpers

## Gaussian CRPS in closed form. Vectorised over (mu, sd, x).
crps_gauss <- function(mu, sd, x) {
  z <- (x - mu) / sd
  sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
}

## Simulate one replicate from a (directional/non-directional) WM field at
## all_loc, returning the (train, test) split of the latent field u and the
## noisy training observations y.
sim_one <- function(directional) {
  u_all <- sample_spde(range       = range_true,
                       sigma       = sigma_true,
                       alpha       = alpha,
                       directional = directional,
                       graph       = graph,
                       PtE         = all_loc)
  u_train <- u_all[seq_len(n_obs)]
  u_test  <- u_all[n_obs + seq_len(n_test)]
  y_train <- u_train + rnorm(n_obs, sd = sigma_e_true)
  y_test  <- u_test  + rnorm(n_test, sd = sigma_e_true)
  list(y_train = y_train, u_test = u_test, y_test = y_test)
}

## Fit a model on training data and predict y at the held-out test
## locations. Returns parameter estimates plus prediction scores. RMSE
## and MAE are computed against the noise-free latent field u_test (so
## they isolate the kriging signal); CRPS is computed against y_test
fit_and_score <- function(sim, model_label) {
  df_train <- data.frame(y                = as.vector(sim$y_train),
                         edge_number      = obs_loc[, 1],
                         distance_on_edge = obs_loc[, 2])
  graph$add_observations(data = df_train, normalized = TRUE,
                         clear_obs = TRUE, suppress_warnings = TRUE,
                         verbose = 0)
  fit <- tryCatch(
    graph_lme(y ~ 1, graph = graph, model = model_label),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(c(sigma = NA_real_, range = NA_real_, sigma_e = NA_real_,
             rmse = NA_real_, mae = NA_real_, crps = NA_real_, ok = 0))
  }

  df_test <- data.frame(edge_number      = test_loc[, 1],
                        distance_on_edge = test_loc[, 2])
  pred <- tryCatch(
    predict(fit, newdata = df_test, normalized = TRUE,
            compute_variances = TRUE),
    error = function(e) e)
  if (inherits(pred, "error")) {
    rmse <- mae <- crps <- NA_real_
  } else {
    mu  <- as.vector(pred$mean)
    rmse <- sqrt(mean((sim$u_test - mu)^2))
    mae  <- mean(abs(sim$u_test - mu))
    if (!is.null(pred$variance)) {
      sdv  <- sqrt(pmax(as.vector(pred$variance), .Machine$double.eps))
      crps <- mean(crps_gauss(mu, sdv, sim$y_test))
    } else {
      crps <- NA_real_
    }
  }

  m <- fit$matern_coeff$random_effects
  c(sigma   = unname(m["sigma"]),
    range   = unname(m["range"]),
    sigma_e = unname(fit$coeff$measurement_error),
    rmse    = rmse,
    mae     = mae,
    crps    = crps,
    ok      = 1)
}

## Run all (true_model x fitted_model) combinations for one replicate.
one_replicate <- function(rep_id) {
  out <- list()
  for (true_dir in c(FALSE, TRUE)) {
    sim <- sim_one(directional = true_dir)
    for (fit_label in c("WM1", "WMD1")) {
      est <- fit_and_score(sim, fit_label)
      out[[length(out) + 1L]] <- data.frame(
        rep        = rep_id,
        true_model = if (true_dir) "WMD1 (true)" else "WM1 (true)",
        fit_model  = fit_label,
        sigma      = est["sigma"],
        range      = est["range"],
        sigma_e    = est["sigma_e"],
        rmse       = est["rmse"],
        mae        = est["mae"],
        crps       = est["crps"],
        ok         = est["ok"])
    }
  }
  do.call(rbind, out)
}

## Main loop

cat(sprintf("Running %d replicates x 2 true models x 2 fits = %d fits.\n",
            nrepl, nrepl * 4))
t0 <- Sys.time()
results <- do.call(rbind, lapply(seq_len(nrepl), function(i) {
  if (i %% 5 == 0 || i == 1)
    cat(sprintf("  replicate %d / %d (%.1fs elapsed)\n",
                i, nrepl, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  one_replicate(i)
}))
cat(sprintf("Done in %.1fs.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Summarise

truth <- c(sigma = sigma_true, range = range_true, sigma_e = sigma_e_true)

summary_tbl <- do.call(rbind, lapply(split(results,
                                           list(results$true_model,
                                                results$fit_model)),
  function(d) {
    if (nrow(d) == 0) return(NULL)
    data.frame(
      true_model   = d$true_model[1],
      fit_model    = d$fit_model[1],
      n_ok         = sum(d$ok == 1),
      bias_sigma   = mean(d$sigma   - truth["sigma"],   na.rm = TRUE),
      bias_range   = mean(d$range   - truth["range"],   na.rm = TRUE),
      bias_sigma_e = mean(d$sigma_e - truth["sigma_e"], na.rm = TRUE),
      mean_rmse    = mean(d$rmse, na.rm = TRUE),
      mean_mae     = mean(d$mae,  na.rm = TRUE),
      mean_crps    = mean(d$crps, na.rm = TRUE))
  }))
rownames(summary_tbl) <- NULL

cat("\n--- Truth ---\n"); print(truth)
cat("\n--- Bias of parameters and mean prediction scores ---\n")
print(summary_tbl, row.names = FALSE, digits = 3)

results$specification <- ifelse(
  (results$true_model == "WMD1 (true)" & results$fit_model == "WMD1") |
  (results$true_model == "WM1 (true)"  & results$fit_model == "WM1"),
  "correctly specified", "misspecified")

# Plots

param_long <- rbind(
  data.frame(results, parameter = "sigma",   estimate = results$sigma,
             truth = truth["sigma"]),
  data.frame(results, parameter = "range",   estimate = results$range,
             truth = truth["range"]),
  data.frame(results, parameter = "sigma_e", estimate = results$sigma_e,
             truth = truth["sigma_e"]))

p_params <- ggplot(param_long, aes(x = fit_model, y = estimate,
                                   fill = specification)) +
  geom_boxplot(outlier.size = 0.7) +
  geom_hline(aes(yintercept = truth), linetype = 2) +
  facet_grid(parameter ~ true_model, scales = "free_y") +
  labs(x = "Fitted model", y = "Estimate",
       title = "Parameter estimates by generative and fitted model",
       subtitle = sprintf("Logo graph, %d edges, %d obs/edge, %d replicates",
                          graph$nE, obs_per_edge, nrepl),
       caption = "Dashed line: true value. Misspecified columns reveal bias.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
print(p_params)

## Prediction scores by (true, fit) model. Lower is better for all three.
score_long <- rbind(
  data.frame(results, score = "RMSE", value = results$rmse),
  data.frame(results, score = "MAE",  value = results$mae),
  data.frame(results, score = "CRPS", value = results$crps))
score_long$score <- factor(score_long$score, levels = c("RMSE", "MAE", "CRPS"))

p_scores <- ggplot(score_long, aes(x = fit_model, y = value,
                                   fill = specification)) +
  geom_boxplot(outlier.size = 0.7) +
  facet_grid(score ~ true_model, scales = "free_y") +
  labs(x = "Fitted model", y = "Score (lower is better)",
       title = "Prediction accuracy at held-out test locations",
       subtitle = sprintf("Logo graph, %d test pts/edge (%d total), %d replicates",
                          test_per_edge, n_test, nrepl),
       caption = "Scores compare predicted latent u to true u.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
print(p_scores)

invisible(list(results  = results,
               summary  = summary_tbl,
               p_params = p_params,
               p_scores = p_scores))
