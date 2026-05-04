# Plot of processed predicted values with 'inlabru'

Auxiliary function to obtain plots of the processed predictions of the
field using 'inlabru'.

## Usage

``` r
# S3 method for class 'graph_bru_proc_pred'
plot(x, y = NULL, vertex_size = 0, ...)
```

## Arguments

- x:

  A processed predicted object obtained with the
  `process_rspde_predictions` function.

- y:

  Not used.

- vertex_size:

  Size of the vertices.

- ...:

  Additional parameters to be passed to plot_function.

## Value

A 'ggplot2' object.
