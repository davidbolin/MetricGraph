# Predict method for 'inlabru' fits on Metric Graphs for 'rSPDE' models

Auxiliar function to obtain predictions of the field using 'inlabru' and
'rSPDE'.

## Usage

``` r
# S3 method for class 'rspde_metric_graph'
predict(
  object,
  cmp,
  bru_fit,
  newdata = NULL,
  formula = NULL,
  data_coords = c("PtE", "euclidean"),
  normalized = TRUE,
  n.samples = 100,
  seed = 0L,
  probs = c(0.025, 0.5, 0.975),
  num.threads = NULL,
  include = NULL,
  exclude = NULL,
  drop = FALSE,
  ...,
  data = deprecated()
)
```

## Arguments

- object:

  An `rspde_metric_graph` object built with the `rspde.metric_graph()`
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

- n.samples:

  Integer setting the number of samples to draw in order to calculate
  the posterior statistics. The default is rather low but provides a
  quick approximate result.

- seed:

  Random number generator seed passed on to inla.posterior.sample

- probs:

  A numeric vector of probabilities with values in the standard unit
  interval to be passed to stats::quantile.

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

- ...:

  Additional arguments passed on to inla.posterior.sample.

- data:

  **\[deprecated\]** Use `newdata` instead.

## Value

A list with predictions.
