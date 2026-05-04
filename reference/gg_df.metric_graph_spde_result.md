# Data frame for metric_graph_spde_result objects to be used in 'ggplot2'

Returns a 'ggplot2'-friendly data-frame with the marginal posterior
densities.

## Usage

``` r
# S3 method for class 'metric_graph_spde_result'
gg_df(
  result,
  parameter = result$params,
  transform = TRUE,
  restrict_x_axis = parameter,
  restrict_quantiles = list(sigma = c(0, 1), range = c(0, 1), kappa = c(0, 1), sigma =
    c(0, 1)),
  ...
)
```

## Arguments

- result:

  A metric_graph_spde_result object.

- parameter:

  Vector. Which parameters to get the posterior density in the
  data.frame? The options are `sigma`, `range` or `kappa`.

- transform:

  Should the posterior density be given in the original scale?

- restrict_x_axis:

  Variables to restrict the range of x axis based on quantiles.

- restrict_quantiles:

  List of quantiles to restrict x axis.

- ...:

  Not being used.

## Value

A `data.frame` containing the posterior densities.
