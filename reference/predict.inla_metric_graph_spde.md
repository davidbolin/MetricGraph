# Predict method for 'inlabru' fits on Metric Graphs

Auxiliar function to obtain predictions of the field using 'inlabru'.

## Usage

``` r
# S3 method for class 'inla_metric_graph_spde'
predict(
  object,
  cmp,
  bru_fit,
  newdata = NULL,
  formula = NULL,
  data_coords = c("PtE", "euclidean"),
  normalized = TRUE,
  repl = NULL,
  repl_col = NULL,
  group = NULL,
  group_col = NULL,
  n.samples = 100,
  seed = 0L,
  probs = c(0.025, 0.5, 0.975),
  return_original_order = TRUE,
  num.threads = NULL,
  include = NULL,
  exclude = NULL,
  drop = FALSE,
  tolerance_merge = 1e-05,
  ...,
  data = deprecated()
)
```

## Arguments

- object:

  An `inla_metric_graph_spde` object built with the
  [`graph_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md)
  function.

- cmp:

  The 'inlabru' component used to fit the model.

- bru_fit:

  A fitted model using 'inlabru' or 'INLA'.

- newdata:

  A data.frame of covariates needed for the prediction. The locations
  must be normalized PtE.

- formula:

  A formula where the right hand side defines an R expression to
  evaluate for each generated sample. If NULL, the latent and
  hyperparameter states are returned as named list elements. See Details
  for more information.

- data_coords:

  It decides which coordinate system to use. If `PtE`, the user must
  provide the locations as a data frame with the first column being the
  edge number and the second column as the distance on edge, otherwise
  if `euclidean`, the user must provide a data frame with the first
  column being the `x` Euclidean coordinates and the second column being
  the `y` Euclidean coordinates.

- normalized:

  if `TRUE`, then the distances in distance on edge are assumed to be
  normalized to (0,1). Default TRUE. Will not be used if `data_coords`
  is `euclidean`.

- repl:

  Which replicates? If there is no replicates, one can set `repl` to
  `NULL`. If one wants all replicates, then one sets to `repl` to
  `.all`.

- repl_col:

  Column containing the replicates. If the replicate is the internal
  group variable, set the replicates to ".group". If not replicates, set
  to `NULL`.

- group:

  Which groups? If there is no groups, one can set `group` to `NULL`. If
  one wants all groups, then one sets to `group` to `.all`.

- group_col:

  Which "column" of the data contains the group variable?

- n.samples:

  Integer setting the number of samples to draw in order to calculate
  the posterior statistics. The default is rather low but provides a
  quick approximate result.

- seed:

  Random number generator seed passed on to `inla.posterior.sample()`

- probs:

  A numeric vector of probabilities with values in the standard unit
  interval to be passed to stats::quantile

- return_original_order:

  Should the predictions be returned in the original order?

- num.threads:

  Specification of desired number of threads for parallel computations.
  Default NULL, leaves it up to 'INLA'. When seed != 0, overridden to
  "16:1"

- include:

  Character vector of component labels that are needed by the predictor
  expression; Default: NULL (include all components that are not
  explicitly excluded)

- exclude:

  Character vector of component labels that are not used by the
  predictor expression. The exclusion list is applied to the list as
  determined by the include parameter; Default: NULL (do not remove any
  components from the inclusion list)

- drop:

  logical; If keep=FALSE, data is a SpatialDataFrame, and the prediciton
  summary has the same number of rows as data, then the output is a
  SpatialDataFrame object. Default FALSE.

- tolerance_merge:

  Tolerance for merging prediction points into original points to
  increase stability.

- ...:

  Additional arguments passed on to `inla.posterior.sample()`.

- data:

  **\[deprecated\]** Use `newdata` instead.

## Value

A list with predictions.
