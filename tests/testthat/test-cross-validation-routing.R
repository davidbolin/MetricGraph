test_that("cross_validation routes identifiable rSPDE fits", {
  fake_fit <- function(model_class) {
    structure(list(bru_info = list(model = list(effects = list(
      field = list(main = list(model = structure(list(), class = model_class)))
    )))), class = "bru")
  }

  rspde_fit <- fake_fit("inla_rspde")
  graph_fit <- fake_fit("inla_metric_graph_spde")
  expect_identical(MetricGraph:::.cv_model_origin(rspde_fit), "rSPDE")
  expect_identical(MetricGraph:::.cv_model_origin(list(rspde_fit, graph_fit)),
                   "MetricGraph")
  # Validation wording identifies the delegated engine without requiring INLA fits.
  expect_error(MetricGraph::cross_validation(rspde_fit, cv_type = "bogus"),
               "possible options for cv_type")
  expect_error(MetricGraph::cross_validation(list(rspde_fit, graph_fit),
                                             cv_type = "bogus"),
               "cv_type must be one of")
})
