# Summary Method for `metric_graph` Objects

Function providing a summary of several informations/characteristics of
a metric graph object.

## Usage

``` r
# S3 method for class 'metric_graph'
summary(
  object,
  messages = FALSE,
  compute_characteristics = NULL,
  check_euclidean = NULL,
  check_distance_consistency = NULL,
  ...
)
```

## Arguments

- object:

  an object of class `metric_graph`.

- messages:

  Should message explaining how to build the results be given for
  missing quantities?

- compute_characteristics:

  Should the characteristics of the graph be computed? If `NULL` it will
  be determined based on the size of the graph.

- check_euclidean:

  Check if the graph has Euclidean edges? If `NULL` it will be
  determined based on the size of the graph.

- check_distance_consistency:

  Check the distance consistency assumption?#' If `NULL` it will be
  determined based on the size of the graph.

- ...:

  not used.

## Value

An object of class `summary_graph_lme` containing information about a
*metric_graph* object.
