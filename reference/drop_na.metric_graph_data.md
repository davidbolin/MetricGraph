# A version of `tidyr::drop_na()` function for datasets on metric graphs

Applies
[`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html)
function for datasets obtained from a metric graph object.

## Usage

``` r
# S3 method for class 'metric_graph_data'
drop_na(data, ...)
```

## Arguments

- data:

  The data list or
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  obtained from a metric graph object.

- ...:

  Additional parameters to be passed to
  [`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html).

## Value

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with the resulting selected columns.
