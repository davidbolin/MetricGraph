# Augment data with information from a `graph_lme` object

Augment accepts a model object and a dataset and adds information about
each observation in the dataset. It includes predicted values in the
`.fitted` column, residuals in the `.resid` column, and standard errors
for the fitted values in a `.se.fit` column. It also contains the New
columns always begin with a . prefix to avoid overwriting columns in the
original dataset.

## Usage

``` r
# S3 method for class 'graph_lme'
augment(
  x,
  newdata = NULL,
  which_repl = NULL,
  sd_post_re = FALSE,
  se_fit = FALSE,
  conf_int = FALSE,
  pred_int = FALSE,
  level = 0.95,
  edge_number = "edge_number",
  distance_on_edge = "distance_on_edge",
  coord_x = "coord_x",
  coord_y = "coord_y",
  data_coords = c("PtE", "spatial"),
  normalized = FALSE,
  no_nugget = FALSE,
  check_euclidean = FALSE,
  ...
)
```

## Arguments

- x:

  A `graph_lme` object.

- newdata:

  A `data.frame` or a `list` containing the covariates, the edge number
  and the distance on edge for the locations to obtain the prediction.
  If `NULL`, the fitted values will be given for the original locations
  where the model was fitted.

- which_repl:

  Which replicates to obtain the prediction. If `NULL` predictions will
  be obtained for all replicates. Default is `NULL`.

- sd_post_re:

  Logical indicating whether or not a .sd_post_re column should be added
  to the augmented output containing the posterior standard deviations
  of the random effects.

- se_fit:

  Logical indicating whether or not a .se_fit column should be added to
  the augmented output containing the standard errors of the fitted
  values. If `TRUE`, the posterior standard deviations of the random
  effects will also be returned.

- conf_int:

  Logical indicating whether or not confidence intervals for the
  posterior mean of the random effects should be built.

- pred_int:

  Logical indicating whether or not prediction intervals for the fitted
  values should be built. If `TRUE`, the confidence intervals for the
  posterior random effects will also be built.

- level:

  Level of confidence and prediction intervals if they are constructed.

- edge_number:

  Name of the variable that contains the edge number, the default is
  `edge_number`.

- distance_on_edge:

  Name of the variable that contains the distance on edge, the default
  is `distance_on_edge`.

- coord_x:

  Column (or entry on the list) of the `data` that contains the x
  coordinate. If not supplied, the column with name "coord_x" will be
  chosen. Will not be used if `Spoints` is not `NULL` or if
  `data_coords` is `PtE`.

- coord_y:

  Column (or entry on the list) of the `data` that contains the y
  coordinate. If not supplied, the column with name "coord_x" will be
  chosen. Will not be used if `Spoints` is not `NULL` or if
  `data_coords` is `PtE`.

- data_coords:

  To be used only if `Spoints` is `NULL`. It decides which coordinate
  system to use. If `PtE`, the user must provide `edge_number` and
  `distance_on_edge`, otherwise if `spatial`, the user must provide
  `coord_x` and `coord_y`.

- normalized:

  Are the distances on edges normalized?

- no_nugget:

  Should the prediction be done without nugget?

- check_euclidean:

  Check if the graph used to compute the resistance distance has
  Euclidean edges? The graph used to compute the resistance distance has
  the observation locations as vertices.

- ...:

  Additional arguments.

## Value

A
[`tidyr::tibble()`](https://tidyr.tidyverse.org/reference/reexports.html)
with columns:

- `.fitted` Fitted or predicted value.

- `.relwrconf` Lower bound of the confidence interval of the random
  effects, if conf_int = TRUE

- `.reuprconf` Upper bound of the confidence interval of the random
  effects, if conf_int = TRUE

- `.fittedlwrpred` Lower bound of the prediction interval, if conf_int =
  TRUE

- `.fitteduprpred` Upper bound of the prediction interval, if conf_int =
  TRUE

- `.fixed` Prediction of the fixed effects.

- `.random` Prediction of the random effects.

- `.resid` The ordinary residuals, that is, the difference between
  observed and fitted values.

- `.std_resid` The standardized residuals, that is, the ordinary
  residuals divided by the standard error of the fitted values (by the
  prediction standard error), if se_fit = TRUE or pred_int = TRUE.

- `.se_fit` Standard errors of fitted values, if se_fit = TRUE.

- `.sd_post_re` Standard deviation of the posterior mean of the random
  effects, if se_fit = TRUE.

## See also

[glance.graph_lme](https://davidbolin.github.io/MetricGraph/reference/glance.graph_lme.md)
