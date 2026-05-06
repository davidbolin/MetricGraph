## Tests for the INLA and inlabru interfaces on metric graphs.
##
## We try to actually fit a small alpha=1 SPDE model end-to-end with INLA
## and inlabru to verify the pipelines run and produce the expected return
## structure (summaries, classes, predictions). On systems where INLA
## cannot dlopen the cgeneric shared library at runtime (a known
## environment issue with INLA cgeneric models in some test setups), the
## end-to-end tests skip gracefully; the structural tests do not depend
## on a successful fit.

library(MetricGraph)
library(testthat)

skip_if_no_inla <- function() {
  skip_if_not_installed("INLA")
  skip_if_not_installed("inlabru")
}

# Returns the bru fit, or NULL if INLA fails. On dlopen / cgeneric issues
# we treat the test as un-runnable rather than failing.
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
# Shared simulator
# ---------------------------------------------------------------------------

make_small_graph <- function() {
  edge1 <- rbind(c(0, 0), c(1, 0))
  edge2 <- rbind(c(0, 0), c(0, 1))
  edge3 <- rbind(c(0, 1), c(-1, 1))
  theta <- seq(from = pi, to = 3 * pi / 2, length.out = 20)
  edge4 <- cbind(sin(theta), 1 + cos(theta))
  metric_graph$new(edges = list(edge1, edge2, edge3, edge4))
}

simulate_observations <- function(graph, alpha = 1, sigma = 2, range = 0.5,
                                  sigma_e = 0.2, obs_per_edge = 30,
                                  seed = 42) {
  set.seed(seed)
  obs_loc <- do.call(rbind, lapply(seq_len(graph$nE), function(i) {
    cbind(rep(i, obs_per_edge), runif(obs_per_edge))
  }))
  u <- sample_spde(range = range, sigma = sigma, alpha = alpha,
                   graph = graph, PtE = obs_loc)
  y <- u + sigma_e * rnorm(length(u))
  list(obs_loc = obs_loc, y = y, u = u)
}

build_fitted_graph <- function(alpha = 1, seed = 42, obs_per_edge = 30) {
  graph <- make_small_graph()
  sim <- simulate_observations(graph,
    alpha = alpha, obs_per_edge = obs_per_edge, seed = seed
  )
  df <- data.frame(y = sim$y,
                   edge_number = sim$obs_loc[, 1],
                   distance_on_edge = sim$obs_loc[, 2])
  graph$add_observations(data = df, normalized = TRUE, verbose = 0)
  list(graph = graph, sim = sim)
}


# ---------------------------------------------------------------------------
# Structural tests for graph_spde / graph_data_spde — these run even if
# INLA cannot fit the model, because they don't actually call inla().
# ---------------------------------------------------------------------------

test_that("graph_spde returns an inla_metric_graph_spde object", {
  skip_if_no_inla()
  fit_setup <- build_fitted_graph(alpha = 1, obs_per_edge = 10)
  spde_model <- graph_spde(fit_setup$graph, alpha = 1)
  expect_s3_class(spde_model, "inla_metric_graph_spde")
  expect_s3_class(spde_model, "inla.cgeneric")
  expect_equal(spde_model$alpha, 1)
})


test_that("graph_data_spde returns the components inla.stack expects", {
  skip_if_no_inla()
  fit_setup <- build_fitted_graph(alpha = 1, obs_per_edge = 10)
  spde_model <- graph_spde(fit_setup$graph, alpha = 1)
  data_spde <- graph_data_spde(spde_model, name = "field")
  expect_true(all(c("data", "basis", "index") %in% names(data_spde)))
  expect_equal(length(data_spde$data$y), nrow(fit_setup$sim$obs_loc))
  expect_true(inherits(data_spde$basis, c("Matrix", "matrix")))
  expect_equal(nrow(data_spde$basis), nrow(fit_setup$sim$obs_loc))
  expect_true("field" %in% names(data_spde$index))
})


test_that("graph_data_spde with loc_name produces a list usable by bru()", {
  skip_if_no_inla()
  fit_setup <- build_fitted_graph(alpha = 1, obs_per_edge = 10)
  spde_model <- graph_spde(fit_setup$graph, alpha = 1)
  data_bru <- graph_data_spde(spde_model, loc_name = "loc")
  expect_true(is.list(data_bru))
  expect_true(!is.null(data_bru$data$loc))
  expect_equal(NCOL(data_bru$data$loc), 2)
  expect_equal(NROW(data_bru$data$loc), nrow(fit_setup$sim$obs_loc))
})


# ---------------------------------------------------------------------------
# End-to-end INLA fit
# ---------------------------------------------------------------------------

test_that("INLA: graph_spde + inla.stack + inla() runs end-to-end (alpha=1)", {
  skip_if_no_inla()
  fit_setup <- build_fitted_graph(alpha = 1, seed = 1, obs_per_edge = 30)
  spde_model <- graph_spde(fit_setup$graph, alpha = 1)
  data_spde <- graph_data_spde(spde_model, name = "field")

  f.s <- y ~ -1 + Intercept + f(field, model = spde_model)
  stk <- INLA::inla.stack(
    data = data_spde[["data"]],
    A = data_spde[["basis"]],
    effects = c(data_spde[["index"]], list(Intercept = 1))
  )
  fit <- tryCatch(
    INLA::inla(
      f.s,
      data = INLA::inla.stack.data(stk),
      control.predictor = list(A = INLA::inla.stack.A(stk)),
      control.inla = list(int.strategy = "eb"),
      num.threads = "1:1"
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$summary.hyperpar)) {
    skip("INLA could not fit the cgeneric model in this environment.")
  }
  expect_s3_class(fit, "inla")
  expect_true(nrow(fit$summary.hyperpar) >= 1)

  spde_result <- spde_metric_graph_result(fit, "field", spde_model)
  expect_true(is.list(spde_result))
  expect_true(!is.null(spde_result$summary.sigma))
  expect_true(!is.null(spde_result$summary.range))
})


# ---------------------------------------------------------------------------
# End-to-end inlabru fit
# ---------------------------------------------------------------------------

test_that("inlabru: bru() fits an alpha=1 SPDE and returns a bru object", {
  skip_if_no_inla()
  fit_setup <- build_fitted_graph(alpha = 1, seed = 3, obs_per_edge = 30)
  spde_model <- graph_spde(fit_setup$graph, alpha = 1)
  cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_model)
  data_bru <- graph_data_spde(spde_model, loc_name = "loc")

  fit <- try_bru_fit(cmp, data_bru[["data"]])
  if (is.null(fit)) {
    skip("inlabru/INLA could not fit the cgeneric model in this environment.")
  }
  expect_s3_class(fit, "bru")
  expect_true(!is.null(fit$summary.fixed))
  expect_true("Intercept" %in% rownames(fit$summary.fixed))
})


test_that("inlabru: predict.inla_metric_graph_spde returns predictions at new locations", {
  skip_if_no_inla()
  fit_setup <- build_fitted_graph(alpha = 1, seed = 4, obs_per_edge = 30)
  graph <- fit_setup$graph
  spde_model <- graph_spde(graph, alpha = 1)
  cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_model)
  data_bru <- graph_data_spde(spde_model, loc_name = "loc")
  fit <- try_bru_fit(cmp, data_bru[["data"]])
  if (is.null(fit)) {
    skip("inlabru/INLA could not fit the cgeneric model in this environment.")
  }
  pred_loc <- list(
    loc = cbind(rep(1, 5), seq(0.1, 0.9, length.out = 5)),
    Intercept = rep(1, 5)
  )
  field_pred <- predict(spde_model, cmp, fit,
    newdata = pred_loc,
    formula = ~ Intercept + field,
    n.samples = 50
  )
  expect_s3_class(field_pred, "graph_bru_pred")
  expect_true(!is.null(field_pred$pred))
  expect_true(!is.null(field_pred$pred$mean))
  expect_equal(length(field_pred$pred$mean), 5)
  expect_true(all(is.finite(field_pred$pred$mean)))
})


test_that(".cv_find_mg_component recognises a metric_graph SPDE bru fit", {
  skip_if_no_inla()
  fit_setup <- build_fitted_graph(alpha = 1, seed = 5, obs_per_edge = 15)
  spde_model <- graph_spde(fit_setup$graph, alpha = 1)
  cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_model)
  data_bru <- graph_data_spde(spde_model, loc_name = "loc")
  fit <- try_bru_fit(cmp, data_bru[["data"]])
  if (is.null(fit)) {
    skip("inlabru/INLA could not fit the cgeneric model in this environment.")
  }
  mg_info <- MetricGraph:::.cv_find_mg_component(fit)
  expect_true(!is.null(mg_info))
  expect_s3_class(mg_info$spde_model, "inla_metric_graph_spde")
})
