# Process predictions of `rspde_metric_graph` objects obtained by using `inlabru`

Auxiliar function to transform the predictions of the field into a plot
friendly object.

## Usage

``` r
process_rspde_predictions(pred, graph, PtE = NULL)
```

## Arguments

- pred:

  The predictions of the field obtained by using `inlabru`

- graph:

  The original `metric_graph` object in which the predictions were
  obtained.

- PtE:

  Normalized locations of the points on the edge.

## Value

A list with predictions.
