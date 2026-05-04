# Prepare data frames or data lists to be used with 'inlabru' in metric graphs

Prepare data frames or data lists to be used with 'inlabru' in metric
graphs

## Usage

``` r
graph_bru_process_data(
  data,
  edge_number = "edge_number",
  distance_on_edge = "distance_on_edge",
  loc = "loc"
)
```

## Arguments

- data:

  A `data.frame` or a `list` containing the covariates, the edge number
  and the distance on edge for the locations to obtain the prediction.

- edge_number:

  Name of the variable that contains the edge number, the default is
  `edge_number`.

- distance_on_edge:

  Name of the variable that contains the distance on edge, the default
  is `distance_on_edge`.

- loc:

  character. Name of the locations to be used in 'inlabru' component.

## Value

A list containing the processed data to be used in a user-friendly
manner by 'inlabru'.
