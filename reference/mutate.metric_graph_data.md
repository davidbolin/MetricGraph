# A version of `dplyr::mutate()` function for datasets on metric graphs

Applies
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
function for datasets obtained from a metric graph object.

## Usage

``` r
# S3 method for class 'metric_graph_data'
mutate(.data, ...)
```

## Arguments

- .data:

  The data list or
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  obtained from a metric graph object.

- ...:

  Additional parameters to be passed to
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).

## Value

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with the resulting selected columns.
