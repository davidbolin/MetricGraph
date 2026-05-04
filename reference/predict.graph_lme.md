# Prediction for a mixed effects regression model on a metric graph

Prediction for a mixed effects regression model on a metric graph

## Usage

``` r
# S3 method for class 'graph_lme'
predict(
  object,
  newdata = NULL,
  mesh = FALSE,
  mesh_h = 0.01,
  which_repl = NULL,
  compute_variances = FALSE,
  compute_pred_variances = FALSE,
  posterior_samples = FALSE,
  pred_samples = FALSE,
  n_samples = 100,
  edge_number = "edge_number",
  distance_on_edge = "distance_on_edge",
  normalized = FALSE,
  no_nugget = FALSE,
  return_as_list = FALSE,
  return_original_order = TRUE,
  check_euclidean = TRUE,
  advanced_options = list(),
  ...,
  data = deprecated()
)
```

## Arguments

- object:

  The fitted object with the
  [`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
  function.

- newdata:

  A `data.frame` or a `list` containing the covariates, the edge number
  and the distance on edge for the locations to obtain the prediction.
  Observe that you should not provide the locations for each replicate.
  Only a single set of locations and covariates, and the predictions for
  the different replicates will be obtained for this same set of
  locations.

- mesh:

  Obtain predictions for mesh nodes? The graph must have a mesh and
  should not have covariates.

- mesh_h:

  If the graph does not have a mesh, one will be created with this value
  of 'h'.

- which_repl:

  Which replicates to obtain the prediction. If `NULL` predictions will
  be obtained for all replicates. Default is `NULL`.

- compute_variances:

  Set to TRUE to compute the kriging variances.

- compute_pred_variances:

  Set to TRUE to compute the prediction variances. Will only be computed
  if newdata is `NULL`.

- posterior_samples:

  If `TRUE`, posterior samples for the random effect will be returned.

- pred_samples:

  If `TRUE`, prediction samples for the response variable will be
  returned. Will only be computed if newdata is `NULL`.

- n_samples:

  Number of samples to be returned. Will only be used if `sampling` is
  `TRUE`.

- edge_number:

  Name of the variable that contains the edge number, the default is
  `edge_number`.

- distance_on_edge:

  Name of the variable that contains the distance on edge, the default
  is `distance_on_edge`.

- normalized:

  Are the distances on edges normalized?

- no_nugget:

  Should the prediction be carried out without the nugget?

- return_as_list:

  Should the means of the predictions and the posterior samples be
  returned as a list, with each replicate being an element?

- return_original_order:

  Should the results be return in the original (input) order or in the
  order inside the graph?

- check_euclidean:

  Check if the graph used to compute the resistance distance has
  Euclidean edges? The graph used to compute the resistance distance has
  the observation locations as vertices.

- advanced_options:

  Advanced options for internal use only. This parameter is intended to
  be used by the cross-validation function and should not be used
  otherwise.

- ...:

  Not used.

- data:

  **\[deprecated\]** Use `newdata` instead.

## Value

A list with elements `mean`, which contains the means of the
predictions, `fe_mean`, which is the prediction for the fixed effects,
`re_mean`, which is the prediction for the random effects, `variance`
(if `compute_variance` is `TRUE`), which contains the posterior
variances of the random effects, `samples` (if `posterior_samples` is
`TRUE`), which contains the posterior samples.
