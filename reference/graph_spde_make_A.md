# Deprecated - Observation/prediction matrices for 'SPDE' models

Constructs observation/prediction weight matrices for metric graph
models.

## Usage

``` r
graph_spde_make_A(graph_spde, repl = NULL)
```

## Arguments

- graph_spde:

  An `inla_metric_graph_spde` object built with the
  [`graph_spde()`](https://davidbolin.github.io/MetricGraph/reference/graph_spde.md)
  function.

- repl:

  Which replicates? If there is no replicates, or to use all replicates,
  one can set to `NULL`.

## Value

The observation matrix.
