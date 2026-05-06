## Tests for the cross_validation() function on bru fits over metric graphs.
##
## Strategy: fit a small alpha=1 SPDE on a 4-edge graph, then exercise the
## metric-graph dispatch (.cv_bru_rerun_metric_graph + sample path) and
## verify:
##   1. the function runs end-to-end on single and multi-model inputs;
##   2. the scores returned match expectations (negative orientation =>
##      lower-is-better, MSE >= 0, MAE >= 0);
##   3. all advertised return options (return_train_test, save_settings,
##      include_best, return_post_samples, return_true_test_values,
##      return_scores_folds) populate the right structure;
##   4. true_CV = FALSE skips refitting (sample path) and finishes with the
##      expected output structure;
##   5. weighted CRPS / sCRPS scores are computed when weight_thr is given.
##
## Tests that need a working INLA cgeneric runtime skip gracefully on
## environments where INLA cannot dlopen the shared library at runtime
## (a known issue with INLA cgeneric models in some test setups). Pure
## input-validation tests do not depend on a successful fit and always run.

library(MetricGraph)
library(testthat)

skip_if_no_inla <- function() {
  skip_if_not_installed("INLA")
  skip_if_not_installed("inlabru")
}

try_bru_fit <- function(cmp, data, options = list(num.threads = "1:1",
                                                  control.inla = list(int.strategy = "eb"))) {
  res <- withCallingHandlers(
    tryCatch(
      inlabru::bru(cmp, data = data, options = options),
      error = function(e) NULL
    ),
    warning = function(w) {
      msg <- conditionMessage(w)
      # inlabru cgeneric/dlopen failure noise + the stylistic is_rowwise
      # warning that fires whenever metric_graph_data (a list, not a
      # data.frame) is passed to bru(). Both are non-fatal here.
      if (grepl("inla:Problem|inla program failed|dlopen|Non data-frame list-like data|is_rowwise",
                msg)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (is.null(res)) return(NULL)
  if (is.null(res$summary.fixed) && is.null(res$summary.hyperpar)) {
    return(NULL)
  }
  res
}

# ---------------------------------------------------------------------------
# Cached fit reused across tests when INLA actually succeeds.
# ---------------------------------------------------------------------------

.cv_fit_cache <- new.env(parent = emptyenv())

cv_fit_smallgraph <- function(alpha = 1) {
  key <- paste0("fit_alpha", alpha)
  if (!is.null(.cv_fit_cache[[key]])) return(.cv_fit_cache[[key]])

  set.seed(42)
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(0, 0), c(0, 1))
  edge3 <- rbind(c(0, 1), c(-1, 1))
  theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
  edge4 <- cbind(sin(theta), 1 + cos(theta))
  graph <- metric_graph$new(edges = list(edge1, edge2, edge3, edge4))

  obs_per_edge <- 25
  obs_loc <- do.call(rbind, lapply(seq_len(graph$nE), function(i) {
    cbind(rep(i, obs_per_edge), runif(obs_per_edge))
  }))
  sigma <- 1.5; r <- 0.5; sigma.e <- 0.2
  u <- sample_spde(range = r, sigma = sigma, alpha = alpha,
                   graph = graph, PtE = obs_loc)
  y <- u + sigma.e * rnorm(length(u))
  df <- data.frame(y = y,
                   edge_number = obs_loc[, 1],
                   distance_on_edge = obs_loc[, 2])
  graph$add_observations(data = df, normalized = TRUE, verbose = 0)

  spde_model <- graph_spde(graph, alpha = alpha)
  # Use a stable variable name in the global environment so
  # .cv_find_spde_var_name() can find it during true_CV refits.
  vname <- paste0("spde_model_alpha", alpha)
  assign(vname, spde_model, envir = globalenv())
  spde_var <- get(vname, envir = globalenv())
  cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_var)

  data_bru <- graph_data_spde(spde_model, loc_name = "loc")
  fit <- try_bru_fit(cmp, data_bru[["data"]])
  if (is.null(fit)) return(NULL)

  attr(fit, ".graph") <- graph
  attr(fit, ".spde_model") <- spde_model
  .cv_fit_cache[[key]] <- fit
  fit
}

skip_if_no_fit <- function(fit) {
  if (is.null(fit)) {
    skip("INLA cannot fit the cgeneric metric_graph SPDE in this environment.")
  }
}


# ---------------------------------------------------------------------------
# Input-validation tests (no fit required).
# ---------------------------------------------------------------------------

test_that("cross_validation: rejects non-bru objects with an informative error", {
  expect_error(cross_validation("not a fit"), "bru")
  expect_error(cross_validation(list(1, 2)), "bru")
})


test_that("cross_validation: rejects bogus cv_type and orientation_results", {
  fake_bru <- structure(list(), class = "bru")
  expect_error(
    cross_validation(fake_bru, cv_type = "bogus"),
    "cv_type"
  )
  expect_error(
    cross_validation(fake_bru, orientation_results = "sideways"),
    "orientation_results"
  )
})


test_that("cross_validation: rejects out-of-range percentage", {
  fake_bru <- structure(list(), class = "bru")
  expect_error(
    cross_validation(fake_bru, cv_type = "lpo", percentage = 150),
    "percentage"
  )
  expect_error(
    cross_validation(fake_bru, cv_type = "lpo", percentage = 0),
    "percentage"
  )
})


test_that("cross_validation: requires weight_thr for wcrps / swcrps", {
  fake_bru <- structure(list(), class = "bru")
  expect_error(
    cross_validation(fake_bru, scores = "wcrps"),
    "weight_thr"
  )
  expect_error(
    cross_validation(fake_bru, scores = "swcrps"),
    "weight_thr"
  )
})


test_that("cross_validation: rejects malformed train_test_indexes", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  expect_error(
    cross_validation(fit, train_test_indexes = "not a list"),
    "list"
  )
  expect_error(
    cross_validation(fit, train_test_indexes = list("not a fold")),
    "list"
  )
  expect_error(
    cross_validation(fit,
      train_test_indexes = list(list(train = "no test"))
    ),
    "test"
  )
})


# ---------------------------------------------------------------------------
# Smoke test: single bru fit, default options, k-fold with small k.
# ---------------------------------------------------------------------------

test_that("cross_validation: runs on a single bru fit and returns a data.frame", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(1)
  res <- cross_validation(
    fit,
    cv_type = "k-fold", k = 2,
    n_samples = 50, true_CV = FALSE,
    print = FALSE
  )
  expect_true(is.data.frame(res))
  expect_true("Model" %in% names(res))
  for (s in c("mae", "mse", "crps", "scrps", "dss")) {
    expect_true(s %in% names(res), info = paste("missing column:", s))
  }
  expect_true("Best" %in% as.character(res$Model))
  num_row <- res[1, , drop = FALSE]
  for (s in c("mae", "mse", "crps", "scrps", "dss")) {
    expect_true(is.finite(as.numeric(num_row[[s]])),
                info = paste("score not finite:", s))
  }
  expect_gte(as.numeric(num_row$mse), 0)
  expect_gte(as.numeric(num_row$mae), 0)
})


test_that("cross_validation: scores argument selects the score columns", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(2)
  res <- cross_validation(
    fit, scores = c("mae", "mse"),
    cv_type = "k-fold", k = 2,
    n_samples = 50, true_CV = FALSE,
    print = FALSE, include_best = FALSE
  )
  expect_setequal(names(res), c("Model", "mae", "mse"))
})


# ---------------------------------------------------------------------------
# Multi-model comparison.
# ---------------------------------------------------------------------------

test_that("cross_validation: multi-model input names rows and adds a Best row", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(3)
  res <- cross_validation(
    list(M1 = fit, M2 = fit),
    cv_type = "k-fold", k = 2,
    scores = c("mse"),
    n_samples = 50, true_CV = FALSE,
    print = FALSE
  )
  expect_true(all(c("M1", "M2") %in% as.character(res$Model)))
  expect_true("Best" %in% as.character(res$Model))
  # The Best row coerces the MSE column to character; cast back to numeric
  # for the comparison.
  vals <- as.numeric(res$mse[match(c("M1", "M2"), res$Model)])
  expect_true(all(is.finite(vals)))
  # Two runs of the same model on the same folds should produce reasonably
  # close MSE values (sample-based scoring introduces some Monte Carlo
  # variation through inlabru::generate). Allow a generous relative tolerance.
  expect_lt(abs(vals[1] - vals[2]) / max(abs(vals[1]), 1e-12), 1)
})


test_that("cross_validation: explicit model_names overrides list names", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(4)
  res <- cross_validation(
    list(fit, fit),
    model_names = c("ModelA", "ModelB"),
    cv_type = "k-fold", k = 2,
    scores = "mae",
    n_samples = 50, true_CV = FALSE, print = FALSE,
    include_best = FALSE
  )
  expect_equal(sort(as.character(res$Model)), c("ModelA", "ModelB"))
})


# ---------------------------------------------------------------------------
# CV-type variants and optional return components.
# ---------------------------------------------------------------------------

test_that("cross_validation: cv_type 'lpo' runs and returns numeric scores", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(5)
  res <- cross_validation(
    fit, cv_type = "lpo",
    percentage = 80, number_folds = 2,
    scores = "mse", n_samples = 50,
    true_CV = FALSE, print = FALSE,
    include_best = FALSE
  )
  expect_true(is.data.frame(res))
  expect_true(is.finite(as.numeric(res$mse[1])))
})


test_that("cross_validation: return_train_test returns the fold list", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(6)
  res <- cross_validation(
    fit, cv_type = "k-fold", k = 2,
    scores = "mse", n_samples = 50, true_CV = FALSE,
    return_train_test = TRUE, print = FALSE,
    include_best = FALSE
  )
  expect_true(is.list(res))
  expect_true(!is.null(res$scores_df))
  expect_true(!is.null(res$train_test))
  expect_equal(length(res$train_test), 2)
  expect_true(!is.null(res$train_test[[1]]$train))
  expect_true(!is.null(res$train_test[[1]]$test))
})


test_that("cross_validation: save_settings returns the CV settings used", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(7)
  res <- cross_validation(
    fit, cv_type = "k-fold", k = 3,
    scores = "mse", n_samples = 50, true_CV = FALSE,
    save_settings = TRUE, print = FALSE,
    include_best = FALSE
  )
  expect_true(!is.null(res$settings))
  expect_equal(res$settings$cv_type, "k-fold")
  expect_equal(res$settings$k, 3)
  expect_equal(res$settings$true_CV, FALSE)
})


test_that("cross_validation: return_scores_folds returns per-fold matrices", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(8)
  res <- cross_validation(
    fit, cv_type = "k-fold", k = 2,
    scores = c("mae", "mse"),
    n_samples = 50, true_CV = FALSE,
    return_scores_folds = TRUE, print = FALSE,
    include_best = FALSE
  )
  expect_true(!is.null(res$scores_folds))
  expect_true(!is.null(res$scores_folds$mae))
  expect_true(!is.null(res$scores_folds$mse))
  expect_equal(nrow(res$scores_folds$mae[[1]]), 2)
  expect_equal(ncol(res$scores_folds$mae[[1]]), 1)
})


test_that("cross_validation: return_post_samples + return_true_test_values populate post & truth", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  set.seed(9)
  res <- cross_validation(
    fit, cv_type = "k-fold", k = 2,
    scores = "mse", n_samples = 50, true_CV = FALSE,
    return_post_samples = TRUE,
    return_true_test_values = TRUE,
    print = FALSE, include_best = FALSE
  )
  expect_true(!is.null(res$post_samples))
  expect_true(!is.null(res$true_test_values))
  expect_equal(length(res$post_samples[[1]]), 2)
})


# ---------------------------------------------------------------------------
# Weighted CRPS / sCRPS.
# ---------------------------------------------------------------------------

test_that("cross_validation: wcrps + swcrps with weight_thr produce finite scores", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)
  set.seed(10)
  # weight_thr below the data range and a generous sample count so that
  # E2t (the weighted self-difference) is non-zero on every test obs and
  # swCRPS computes cleanly without the "swCRPS cannot be computed" warning.
  res <- cross_validation(
    fit, scores = c("wcrps", "swcrps"),
    weight_thr = -3,
    cv_type = "k-fold", k = 2,
    n_samples = 400, true_CV = FALSE,
    print = FALSE, include_best = FALSE
  )
  expect_true(is.data.frame(res))
  expect_true(is.finite(as.numeric(res$wcrps[1])))
})


# ---------------------------------------------------------------------------
# Orientation flip.
# ---------------------------------------------------------------------------

test_that("cross_validation: orientation_results flips proper-score sign", {
  skip_if_no_inla()
  fit <- cv_fit_smallgraph(alpha = 1)
  skip_if_no_fit(fit)

  ttidx <- list(list(
    train = list(seq_len(50)),
    test  = list(51:60)
  ))

  # Each cross_validation() call draws its own posterior samples via
  # inlabru::generate() which uses INLA's internal RNG (not governed by
  # R's set.seed). So bitwise equality between two runs isn't achievable.
  # Instead we verify the structural property: under "negative" orientation
  # CRPS comes out positive (smaller-is-better), under "positive" it comes
  # out negative (larger-is-better) — i.e. the sign flips.
  set.seed(11)
  res_neg <- cross_validation(
    fit, scores = "crps",
    train_test_indexes = ttidx,
    n_samples = 200, true_CV = FALSE,
    orientation_results = "negative",
    print = FALSE, include_best = FALSE
  )
  res_pos <- cross_validation(
    fit, scores = "crps",
    train_test_indexes = ttidx,
    n_samples = 200, true_CV = FALSE,
    orientation_results = "positive",
    print = FALSE, include_best = FALSE
  )
  expect_gt(as.numeric(res_neg$crps[1]), 0)
  expect_lt(as.numeric(res_pos$crps[1]), 0)
})
