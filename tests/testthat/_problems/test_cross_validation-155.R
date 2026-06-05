# Extracted from test_cross_validation.R:155

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
fake_bru <- structure(list(), class = "bru")
expect_error(
    cross_validation(fake_bru, scores = "wcrps"),
    "weight_thr"
  )
expect_error(
    cross_validation(fake_bru, scores = "swcrps"),
    "weight_thr"
  )
