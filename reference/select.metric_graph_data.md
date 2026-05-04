# A version of `dplyr::select()` function for datasets on metric graphs

Selects columns on metric graphs, while keeps the spatial positions.

## Usage

``` r
# S3 method for class 'metric_graph_data'
select(.data, ...)
```

## Arguments

- .data:

  The data list or
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  obtained from a metric graph object.

- ...:

  Additional parameters to be passed to
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html).

## Value

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with the resulting selected columns.
