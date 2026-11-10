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
  # rSPDE fits are passed on to rSPDE, so this check does not need INLA.
  expect_error(MetricGraph::cross_validation(rspde_fit, cv_type = "bogus"),
               "possible options for cv_type")
})

test_that("cross_validation validates mixed fits in MetricGraph", {
  # MetricGraph::cross_validation() checks for INLA before validating cv_type.
  skip_if_no_safe_inla()
  fake_fit <- function(model_class) {
    structure(list(bru_info = list(model = list(effects = list(
      field = list(main = list(model = structure(list(), class = model_class)))
    )))), class = "bru")
  }

  rspde_fit <- fake_fit("inla_rspde")
  graph_fit <- fake_fit("inla_metric_graph_spde")
  expect_error(MetricGraph::cross_validation(list(rspde_fit, graph_fit),
                                             cv_type = "bogus"),
               "cv_type must be one of")
})
