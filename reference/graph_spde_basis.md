# Deprecated - Observation/prediction matrices for 'SPDE' models

Constructs observation/prediction weight matrices for metric graph
models.

## Usage

``` r
graph_spde_basis(graph_spde, repl = NULL, drop_na = FALSE, drop_all_na = TRUE)
```

## Arguments

- graph_spde:

  An `inla_metric_graph_spde` object built with the
  [`graph_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md)
  function.

- repl:

  Which replicates? If there is no replicates, or to use all replicates,
  one can set to `NULL`.

- drop_na:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- drop_all_na:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

## Value

The observation matrix.
