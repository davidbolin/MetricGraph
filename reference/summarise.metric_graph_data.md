# A version of `dplyr::summarise()` function for datasets on metric graphs

Creates summaries, while keeps the spatial positions.

## Usage

``` r
# S3 method for class 'metric_graph_data'
summarise(.data, ..., .include_graph_groups = FALSE, .groups = NULL)
```

## Arguments

- .data:

  The data list or
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  obtained from a metric graph object.

- ...:

  Additional parameters to be passed to
  [`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).

- .include_graph_groups:

  Should the internal graph groups be included in the grouping
  variables? The default is `FALSE`. This means that, when summarising,
  the data will be grouped by the internal group variable together with
  the spatial locations.

- .groups:

  A vector of strings containing the names of the columns to be
  additionally grouped, when computing the summaries. The default is
  `NULL`.

## Value

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with the resulting selected columns.
