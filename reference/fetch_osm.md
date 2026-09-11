# Fetch OpenStreetMap data via the Overpass API

A drop-in replacement for
[`osmdata::osmdata_sf()`](https://docs.ropensci.org/osmdata/reference/osmdata_sf.html)
that sidesteps the Overpass server-status check `osmdata` performs
before each download.`fetch_osm()` builds the Overpass query with
`osmdata`'s helpers, POSTs it directly with the `curl` package, and then
hands the raw response back to `osmdata` for parsing.

## Usage

``` r
fetch_osm(
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
  [`osmdata::add_osm_feature()`](https://docs.ropensci.org/osmdata/reference/add_osm_feature.html).

## Value

An `osmdata` list with elements `$osm_points`, `$osm_lines`,
`$osm_polygons`, `$osm_multilines`, `$osm_multipolygons` (same structure
as
[`osmdata::osmdata_sf()`](https://docs.ropensci.org/osmdata/reference/osmdata_sf.html)).
Tags appear as columns of each `sf` data frame.

## See also

[`metric_graph_from_osm()`](https://davidbolin.github.io/MetricGraph/reference/metric_graph_from_osm.md)
which wraps this and returns a ready-to-use `metric_graph`.

## Examples

``` r
if (FALSE) { # \dontrun{
## All highway ways in a small bbox around San Jose.
osm <- fetch_osm(
  bbox  = c(-122.10, 37.22, -121.80, 37.45),
  key   = "highway",
  value = c("motorway", "motorway_link", "primary", "secondary")
)
head(osm$osm_lines[, c("highway", "oneway", "geometry")])
} # }
```
