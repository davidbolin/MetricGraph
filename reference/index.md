# Package index

## MetricGraph package

- [`MetricGraph`](https://davidbolin.github.io/MetricGraph/reference/MetricGraph-package.md)
  [`MetricGraph-package`](https://davidbolin.github.io/MetricGraph/reference/MetricGraph-package.md)
  : Gaussian processes on metric graphs
- [`logo_lines()`](https://davidbolin.github.io/MetricGraph/reference/logo_lines.md)
  : Create lines for package name

## Metric graph constructor

- [`metric_graph`](https://davidbolin.github.io/MetricGraph/reference/metric_graph.md)
  : Metric graph
- [`graph_components`](https://davidbolin.github.io/MetricGraph/reference/graph_components.md)
  **\[deprecated\]** : Connected components of metric graph
- [`update_graph()`](https://davidbolin.github.io/MetricGraph/reference/update_graph.md)
  : Update an older version metric graph to the current package version

## INLA and SPDE approach on metric graphs

- [`graph_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md)
  : 'INLA' implementation of Whittle-Matérn fields for metric graphs

- [`spde_metric_graph_result()`](https://davidbolin.github.io/MetricGraph/reference/spde_metric_graph_result.md)
  : Metric graph SPDE result extraction from 'INLA' estimation results

- [`summary(`*`<metric_graph_spde_result>`*`)`](https://davidbolin.github.io/MetricGraph/reference/summary.metric_graph_spde_result.md)
  :

  Summary for posteriors of field parameters for an `inla_rspde` model
  from a `rspde.result` object

- [`gg_df(`*`<metric_graph_spde_result>`*`)`](https://davidbolin.github.io/MetricGraph/reference/gg_df.metric_graph_spde_result.md)
  : Data frame for metric_graph_spde_result objects to be used in
  'ggplot2'

- [`graph_data_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_data_spde.md)
  : Data extraction for 'spde' models

- [`graph_spde_basis()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde_basis.md)
  : Deprecated - Observation/prediction matrices for 'SPDE' models

- [`graph_spde_make_A()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde_make_A.md)
  : Deprecated - Observation/prediction matrices for 'SPDE' models

## inlabru and SPDE approach on metric graphs

- [`bru_get_mapper(`*`<inla_metric_graph_spde>`*`)`](https://davidbolin.github.io/MetricGraph/reference/bru_mapper.inla_metric_graph_spde.md)
  [`ibm_n(`*`<bru_mapper_inla_metric_graph_spde>`*`)`](https://davidbolin.github.io/MetricGraph/reference/bru_mapper.inla_metric_graph_spde.md)
  [`ibm_values(`*`<bru_mapper_inla_metric_graph_spde>`*`)`](https://davidbolin.github.io/MetricGraph/reference/bru_mapper.inla_metric_graph_spde.md)
  [`ibm_jacobian(`*`<bru_mapper_inla_metric_graph_spde>`*`)`](https://davidbolin.github.io/MetricGraph/reference/bru_mapper.inla_metric_graph_spde.md)
  : Metric graph 'inlabru' mapper

- [`predict(`*`<inla_metric_graph_spde>`*`)`](https://davidbolin.github.io/MetricGraph/reference/predict.inla_metric_graph_spde.md)
  : Predict method for 'inlabru' fits on Metric Graphs

- [`predict(`*`<rspde_metric_graph>`*`)`](https://davidbolin.github.io/MetricGraph/reference/predict.rspde_metric_graph.md)
  : Predict method for 'inlabru' fits on Metric Graphs for 'rSPDE'
  models

- [`plot(`*`<graph_bru_pred>`*`)`](https://davidbolin.github.io/MetricGraph/reference/plot.graph_bru_pred.md)
  : Plot of predicted values with 'inlabru'

- [`plot(`*`<graph_bru_proc_pred>`*`)`](https://davidbolin.github.io/MetricGraph/reference/plot.graph_bru_proc_pred.md)
  : Plot of processed predicted values with 'inlabru'

- [`graph_bru_process_data()`](https://davidbolin.github.io/MetricGraph/reference/graph_bru_process_data.md)
  : Prepare data frames or data lists to be used with 'inlabru' in
  metric graphs

- [`process_rspde_predictions()`](https://davidbolin.github.io/MetricGraph/reference/process_rspde_predictions.md)
  :

  Process predictions of `rspde_metric_graph` objects obtained by using
  `inlabru`

- [`cross_validation()`](https://davidbolin.github.io/MetricGraph/reference/cross_validation.md)
  : Perform cross-validation on a list of fitted inlabru models on
  metric graphs.

## Linear mixed-effects models

- [`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
  : Metric graph linear mixed effects models

- [`predict(`*`<graph_lme>`*`)`](https://davidbolin.github.io/MetricGraph/reference/predict.graph_lme.md)
  : Prediction for a mixed effects regression model on a metric graph

- [`simulate(`*`<graph_lme>`*`)`](https://davidbolin.github.io/MetricGraph/reference/simulate.graph_lme.md)
  : Simulation of models on metric graphs

- [`summary(`*`<graph_lme>`*`)`](https://davidbolin.github.io/MetricGraph/reference/summary.graph_lme.md)
  :

  Summary Method for `graph_lme` Objects

- [`posterior_crossvalidation(`*`<graph_lme>`*`)`](https://davidbolin.github.io/MetricGraph/reference/posterior_crossvalidation.graph_lme.md)
  :

  Cross-validation for `graph_lme` models assuming observations at the
  vertices of metric graphs

- [`posterior_crossvalidation_loo()`](https://davidbolin.github.io/MetricGraph/reference/posterior_crossvalidation_loo.md)
  :

  Leave-one-out pseudo-crossvalidation for `graph_lme` models assuming
  observations at the vertices of metric graphs

- [`graph_starting_values()`](https://davidbolin.github.io/MetricGraph/reference/graph_starting_values.md)
  : Starting values for random field models on metric graphs

- [`glance(`*`<graph_lme>`*`)`](https://davidbolin.github.io/MetricGraph/reference/glance.graph_lme.md)
  :

  Glance at a `graph_lme` object

- [`augment(`*`<graph_lme>`*`)`](https://davidbolin.github.io/MetricGraph/reference/augment.graph_lme.md)
  :

  Augment data with information from a `graph_lme` object

## Log-Cox Gaussian processes

- [`graph_lgcp_sim()`](https://davidbolin.github.io/MetricGraph/reference/graph_lgcp_sim.md)
  : Simulate log-Gaussian Cox processes on metric graphs
- [`lgcp_graph()`](https://davidbolin.github.io/MetricGraph/reference/lgcp_graph.md)
  : Fit log-Gaussian Cox process models on metric graphs
- [`precompute_lgcp_graph()`](https://davidbolin.github.io/MetricGraph/reference/precompute_lgcp_graph.md)
  : Precompute expensive quantities for efficient LGCP model fitting

## Space-time models

- [`simulate_spacetime()`](https://davidbolin.github.io/MetricGraph/reference/simulate_spacetime.md)
  : space-time simulation based on implicit Euler discretization in time
- [`make_Q_euler()`](https://davidbolin.github.io/MetricGraph/reference/make_Q_euler.md)
  : Space-time precision operator Euler discretization
- [`make_Q_spacetime()`](https://davidbolin.github.io/MetricGraph/reference/make_Q_spacetime.md)
  : Space-time precision operator discretization

## Sampling SPDE on metric graphs

- [`sample_spde()`](https://davidbolin.github.io/MetricGraph/reference/sample_spde.md)
  : Samples a Whittle-Matérn field on a metric graph
- [`simulate(`*`<metric_graph>`*`)`](https://davidbolin.github.io/MetricGraph/reference/simulate.metric_graph.md)
  : Simulate a Whittle-Matérn field on a metric graph
- [`simulate_parallel()`](https://davidbolin.github.io/MetricGraph/reference/simulate_parallel.md)
  : Parallel simulation of a Whittle-Matérn field on a metric graph

## Data manipulation on metric graphs

- [`select(`*`<metric_graph_data>`*`)`](https://davidbolin.github.io/MetricGraph/reference/select.metric_graph_data.md)
  :

  A version of
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
  function for datasets on metric graphs

- [`filter(`*`<metric_graph_data>`*`)`](https://davidbolin.github.io/MetricGraph/reference/filter.metric_graph_data.md)
  :

  A version of
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  function for datasets on metric graphs

- [`mutate(`*`<metric_graph_data>`*`)`](https://davidbolin.github.io/MetricGraph/reference/mutate.metric_graph_data.md)
  :

  A version of
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
  function for datasets on metric graphs

- [`summarise(`*`<metric_graph_data>`*`)`](https://davidbolin.github.io/MetricGraph/reference/summarise.metric_graph_data.md)
  :

  A version of
  [`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
  function for datasets on metric graphs

- [`drop_na(`*`<metric_graph_data>`*`)`](https://davidbolin.github.io/MetricGraph/reference/drop_na.metric_graph_data.md)
  :

  A version of
  [`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html)
  function for datasets on metric graphs

- [`match_mesh_data()`](https://davidbolin.github.io/MetricGraph/reference/match_mesh_data.md)
  : Match Data Frame Rows to Graph Mesh Order

## Methods for metric graphs

- [`summary(`*`<metric_graph>`*`)`](https://davidbolin.github.io/MetricGraph/reference/summary.metric_graph.md)
  :

  Summary Method for `metric_graph` Objects

## Precision matrices

- [`spde_precision()`](https://davidbolin.github.io/MetricGraph/reference/spde_precision.md)
  : Precision matrix for Whittle-Matérn fields

## Covariance functions

- [`exp_covariance()`](https://davidbolin.github.io/MetricGraph/reference/exp_covariance.md)
  : Exponential covariance function
- [`spde_covariance()`](https://davidbolin.github.io/MetricGraph/reference/spde_covariance.md)
  : Covariance function for Whittle-Matérn fields
- [`spde_variance()`](https://davidbolin.github.io/MetricGraph/reference/spde_variance.md)
  : Variancefor Whittle-Matérn fields
- [`directional_ou_covariance()`](https://davidbolin.github.io/MetricGraph/reference/directional_ou_covariance.md)
  : Closed-form covariance of the directional OU model
- [`directional_ou_variance()`](https://davidbolin.github.io/MetricGraph/reference/directional_ou_variance.md)
  : Closed-form marginal variance of the directional OU model

## Auxiliary constructors

- [`linnet.to.graph()`](https://davidbolin.github.io/MetricGraph/reference/linnet.to.graph.md)
  :

  Convert a `linnet` object to a metric graph object

- [`psp.to.graph()`](https://davidbolin.github.io/MetricGraph/reference/psp.to.graph.md)
  :

  Convert a `psp` object to a metric graph object

- [`stlpp.to.graph()`](https://davidbolin.github.io/MetricGraph/reference/stlpp.to.graph.md)
  :

  Convert an `stlpp` object to a metric graph object

- [`fetch_osm()`](https://davidbolin.github.io/MetricGraph/reference/fetch_osm.md)
  : Fetch OpenStreetMap data via the Overpass API

- [`metric_graph_from_osm()`](https://davidbolin.github.io/MetricGraph/reference/metric_graph_from_osm.md)
  : Build a metric_graph directly from an OpenStreetMap query

## Misc functions

- [`selected_inv()`](https://davidbolin.github.io/MetricGraph/reference/selected_inv.md)
  : Selected Inverse Calculation

## Datasets

- [`pems`](https://davidbolin.github.io/MetricGraph/reference/pems.md) :
  Traffic speed data from San Jose, California
- [`pems_repl`](https://davidbolin.github.io/MetricGraph/reference/pems_repl.md)
  : Traffic speed data with replicates from San Jose, California
- [`columbia_main_component`](https://davidbolin.github.io/MetricGraph/reference/columbia_main_component.md)
  : Largest connected component of the Mid-Columbia stream network
