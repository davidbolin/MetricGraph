# The graph-Laplacian likelihood used to build its observation matrix once with
# `group = ".all"`, which is the block-diagonal operator over every replicate,
# and then use it inside the per-replicate loop against a single-replicate
# precision matrix. That conforms only when there is one replicate, so the
# likelihood errored for replicated data. These tests pin both the single- and
# multi-replicate cases against the independent covariance formulation.

build_graph <- function(n_repl) {
  set.seed(1)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1),
             c(-1, 1), c(-1, 0), c(0, -1))
  E <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5),
             c(5, 6), c(6, 1), c(4, 1), c(1, 7))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)

  n_per_edge <- 5
  PtE <- NULL
  for (i in seq_len(graph$nE)) {
    PtE <- rbind(PtE, cbind(rep(i, n_per_edge), runif(n_per_edge)))
  }
  n_row <- nrow(PtE)
  df <- data.frame(y = rnorm(n_row * n_repl),
                   edge_number      = rep(PtE[, 1], times = n_repl),
                   distance_on_edge = rep(PtE[, 2], times = n_repl),
                   repl             = rep(seq_len(n_repl), each = n_row))
  suppressWarnings(graph$add_observations(
    data = df, normalized = TRUE,
    group = if (n_repl > 1) "repl" else NULL,
    verbose = 0, suppress_warnings = TRUE))
  graph$observation_to_vertex(verbose = 0)
  graph$compute_laplacian(full = TRUE)
  graph
}

test_that("graph-Laplacian likelihood works with replicates", {
  for (n_repl in c(1L, 2L, 4L)) {
    graph <- build_graph(n_repl)
    y <- graph$get_data()$y
    for (alpha in c(1, 2)) {
      lik <- MetricGraph:::likelihood_graph_laplacian(
        graph, alpha = alpha, y_graph = y, repl = NULL, X_cov = NULL,
        parameterization = "spde")
      val <- lik(log(c(0.1, 20, 10)))
      expect_true(is.finite(as.numeric(val)),
                  info = paste("n_repl", n_repl, "alpha", alpha))
    }
  }
})

test_that("Laplacian and covariance formulations agree with replicates", {
  theta <- c(0.1, 20, 10)
  for (n_repl in c(1L, 2L, 4L)) {
    graph <- build_graph(n_repl)
    y <- graph$get_data()$y
    for (alpha in c(1, 2)) {
      lik_cov <- MetricGraph:::likelihood_graph_covariance(
        graph, model = paste0("GL", alpha), log_scale = FALSE,
        y_graph = y, repl = NULL, X_cov = NULL)(theta)
      lik_lap <- MetricGraph:::likelihood_graph_laplacian(
        graph, alpha = alpha, y_graph = y, repl = NULL, X_cov = NULL,
        parameterization = "spde")(log(theta))
      expect_equal(as.numeric(lik_lap), as.numeric(lik_cov), tolerance = 1e-8,
                   info = paste("n_repl", n_repl, "alpha", alpha))
    }
  }
})

test_that("Laplacian likelihood handles NA observations with replicates", {
  graph <- build_graph(3L)
  y <- graph$get_data()$y
  y[seq(1, length(y), by = 7)] <- NA
  for (alpha in c(1, 2)) {
    val <- MetricGraph:::likelihood_graph_laplacian(
      graph, alpha = alpha, y_graph = y, repl = NULL, X_cov = NULL,
      parameterization = "spde")(log(c(0.1, 20, 10)))
    expect_true(is.finite(as.numeric(val)), info = paste("alpha", alpha))
  }
})
