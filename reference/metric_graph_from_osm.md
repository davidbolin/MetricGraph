# Build a metric_graph directly from an OpenStreetMap query

Convenience wrapper that fetches OSM data with
[`fetch_osm()`](https://davidbolin.github.io/MetricGraph/reference/fetch_osm.md)
and builds a
[metric_graph](https://davidbolin.github.io/MetricGraph/reference/metric_graph.md)
from the returned linestrings in one call. `longlat = TRUE` is set
automatically because OSM coordinates are always in WGS84 latitude /
longitude. Other graph-construction arguments (`tolerance`,
`perform_merges`, ...) are forwarded to
[metric_graph](https://davidbolin.github.io/MetricGraph/reference/metric_graph.md)'s
`$new()` method without overriding its defaults.

## Usage

``` r
metric_graph_from_osm(
  bbox,
  key = "highway",
  value = NULL,
  endpoint = "https://overpass-api.de/api/interpreter",
  timeout = 180,
  user_agent = "MetricGraph-R-package",
  cache_path = NULL,
  retries = 2,
  ...
)
```

## Arguments

- bbox:

  A bounding box, passed to
  [`osmdata::opq()`](https://docs.ropensci.org/osmdata/reference/opq.html).
  May be a length-4 numeric vector `c(xmin, ymin, xmax, ymax)`, an
  `sf::bbox` object, or any other form accepted by
  [`osmdata::opq()`](https://docs.ropensci.org/osmdata/reference/opq.html).

- key:

  OSM key to filter on, e.g. `"highway"`. Passed to
  [`osmdata::add_osm_feature()`](https://docs.ropensci.org/osmdata/reference/add_osm_feature.html).
  Use `NULL` to skip the feature filter and fetch everything in the
  bbox.

- value:

  Optional character vector of values for `key`, e.g.
  `c("motorway", "motorway_link")`.

- endpoint:

  Overpass endpoint URL, or a character vector of URLs that are tried in
  order. Defaults to the main `overpass-api.de` endpoint. Mirror URLs
  (e.g. `https://overpass.kumi.systems/api/interpreter`) work too.

- timeout:

  Server-side query timeout in seconds. Forwarded to both
  [`osmdata::opq()`](https://docs.ropensci.org/osmdata/reference/opq.html)
  and the HTTP request.

- user_agent:

  String sent as the `User-Agent` HTTP header.

- cache_path:

  Optional path. If provided and the file already exists, the download
  is skipped and the cached response is parsed instead. If provided and
  the file does not exist, the response is saved there for re-use.

- retries:

  Number of times a request to each endpoint is retried after a
  transient failure (HTTP 429 or 5xx, or a connection error), with
  exponential backoff between attempts.

- ...:

  Forwarded to
  [metric_graph](https://davidbolin.github.io/MetricGraph/reference/metric_graph.md)'s
  `$new()`. Use this to pass `tolerance`, `perform_merges`,
  `check_connected`, `which_longlat`, etc.

## Value

A
[metric_graph](https://davidbolin.github.io/MetricGraph/reference/metric_graph.md)
object built from the OSM ways.

## See also

[`fetch_osm()`](https://davidbolin.github.io/MetricGraph/reference/fetch_osm.md)
for the lower-level fetch + parse step.

## Examples

``` r
if (FALSE) { # \dontrun{
## Build a freeway graph for a small bbox around San Jose, splitting
## at intersections via the edge_edge tolerance.
g <- metric_graph_from_osm(
  bbox           = c(-122.10, 37.22, -121.80, 37.45),
  value          = c("motorway", "motorway_link"),
  perform_merges = TRUE,
  tolerance      = list(edge_edge = 1e-5)
)
g$plot()
} # }
```
