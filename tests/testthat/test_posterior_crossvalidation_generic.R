## posterior_crossvalidation() is the S3 generic from rSPDE. MetricGraph
## re-exports it and provides the graph_lme method, so attaching rSPDE after
## MetricGraph does not mask it.

test_that("MetricGraph re-exports rSPDE's posterior_crossvalidation generic", {
  expect_identical(MetricGraph::posterior_crossvalidation,
                   rSPDE::posterior_crossvalidation)
})

test_that("a list of graph_lme fits is dispatched to the graph_lme method", {
  set.seed(1)
  V <- rbind(c(0, 0), c(1, 0), c(1, 1))
  E <- rbind(c(1, 2), c(2, 3))
  graph <- metric_graph$new(V = V, E = E, verbose = 0)
  obs_per_edge <- 10
  PtE <- cbind(rep(seq_len(graph$nE), each = obs_per_edge),
               runif(graph$nE * obs_per_edge))
  u <- sample_spde(kappa = 5, tau = 1, alpha = 1, graph = graph, PtE = PtE)
  graph$add_observations(data = data.frame(y = u + 0.3 * rnorm(length(u)),
                                           edge_number = PtE[, 1],
                                           distance_on_edge = PtE[, 2]),
                         normalized = TRUE, verbose = 0)
  fit <- graph_lme(y ~ -1, graph = graph,
                   model = list(type = "WhittleMatern", alpha = 1))

  single <- posterior_crossvalidation(fit, mode = "loo")
  both <- posterior_crossvalidation(list(A = fit, B = fit), mode = "loo")

  expect_equal(both$scores$Model, c("A", "B"))
  expect_equal(both$mu$A, single$mu)
  expect_equal(both$scores$mae, rep(single$scores$mae, 2))
})
