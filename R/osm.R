## ============================================================================
## Helpers for fetching OpenStreetMap data and turning it into a
## metric_graph. The Overpass API is queried directly so we can
## sidestep the rate-limited server-status check that osmdata performs
## before each download.
## ============================================================================

#' Fetch OpenStreetMap data via the Overpass API
#'
#' @description
#' A drop-in replacement for [osmdata::osmdata_sf()] that sidesteps
#' the Overpass server-status check `osmdata` performs before each
#' download.`fetch_osm()` builds the Overpass query
#' with `osmdata`'s helpers, POSTs it directly with the `curl`
#' package, and then hands the raw response back to `osmdata` for
#' parsing.
#'
#' @param bbox A bounding box, passed to [osmdata::opq()]. May be a
#'   length-4 numeric vector `c(xmin, ymin, xmax, ymax)`, an `sf::bbox`
#'   object, or any other form accepted by [osmdata::opq()].
#' @param key OSM key to filter on, e.g. `"highway"`. Passed to
#'   [osmdata::add_osm_feature()]. Use `NULL` to skip the feature
#'   filter and fetch everything in the bbox.
#' @param value Optional character vector of values for `key`, e.g.
#'   `c("motorway", "motorway_link")`.
#' @param endpoint Overpass endpoint URL, or a character vector of URLs
#'   that are tried in order. Defaults to the main `overpass-api.de`
#'   endpoint. Mirror URLs (e.g.
#'   `https://overpass.kumi.systems/api/interpreter`) work too.
#' @param timeout Server-side query timeout in seconds. Forwarded to
#'   both [osmdata::opq()] and the HTTP request.
#' @param user_agent String sent as the `User-Agent` HTTP header.
#' @param cache_path Optional path. If provided and the file already
#'   exists, the download is skipped and the cached response is parsed
#'   instead. If provided and the file does not exist, the response is
#'   saved there for re-use.
#' @param retries Number of times a request to each endpoint is retried
#'   after a transient failure (HTTP 429 or 5xx, or a connection error),
#'   with exponential backoff between attempts.
#' @param ... Forwarded to [osmdata::add_osm_feature()].
#'
#' @return An `osmdata` list with elements `$osm_points`,
#'   `$osm_lines`, `$osm_polygons`, `$osm_multilines`,
#'   `$osm_multipolygons` (same structure as [osmdata::osmdata_sf()]).
#'   Tags appear as columns of each `sf` data frame.
#'
#' @seealso [metric_graph_from_osm()] which wraps this and returns a
#'   ready-to-use `metric_graph`.
#'
#' @examples
#' \dontrun{
#' ## All highway ways in a small bbox around San Jose.
#' osm <- fetch_osm(
#'   bbox  = c(-122.10, 37.22, -121.80, 37.45),
#'   key   = "highway",
#'   value = c("motorway", "motorway_link", "primary", "secondary")
#' )
#' head(osm$osm_lines[, c("highway", "oneway", "geometry")])
#' }
#'
#' @export
fetch_osm <- function(bbox,
                      key        = "highway",
                      value      = NULL,
                      endpoint   = "https://overpass-api.de/api/interpreter",
                      timeout    = 180,
                      user_agent = "MetricGraph-R-package",
                      cache_path = NULL,
                      retries    = 2,
                      ...) {
  if (!requireNamespace("osmdata", quietly = TRUE)) {
    stop("Package 'osmdata' is required for fetch_osm(). ",
         "Install with install.packages(\"osmdata\").",
         call. = FALSE)
  }

  query <- osmdata::opq(bbox = bbox, timeout = timeout)
  if (!is.null(key)) {
    query <- osmdata::add_osm_feature(query, key = key, value = value, ...)
  }

  out_path <- if (is.null(cache_path)) tempfile(fileext = ".osm") else cache_path
  need_download <- is.null(cache_path) || !file.exists(out_path)

  if (need_download) {
    if (!is.null(cache_path)) {
      dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    }
    body <- osmdata::opq_string(query)
    ## Overpass signals overload with 429 / 5xx; these are worth retrying.
    transient_status <- c(429L, 500L, 502L, 503L, 504L)
    failures <- character(0)
    success <- FALSE
    for (url in endpoint) {
      for (attempt in seq_len(retries + 1L)) {
        h <- curl::new_handle()
        curl::handle_setopt(h,
                            customrequest = "POST",
                            postfields    = body,
                            useragent     = user_agent,
                            timeout       = timeout)
        resp <- tryCatch(curl::curl_fetch_disk(url, out_path, handle = h),
                         error = function(e) e)
        if (!inherits(resp, "error") && resp$status_code == 200L) {
          success <- TRUE
          break
        }
        if (inherits(resp, "error")) {
          msg <- conditionMessage(resp)
          transient <- TRUE
        } else {
          head_lines <- tryCatch(
            paste(readLines(out_path, n = 5, warn = FALSE), collapse = "\n"),
            error = function(e) "")
          msg <- sprintf("HTTP %d. Response head:\n%s",
                         resp$status_code, substr(head_lines, 1, 500))
          transient <- resp$status_code %in% transient_status
        }
        failures <- c(failures, sprintf("%s (attempt %d): %s", url, attempt, msg))
        if (!transient) break
        if (attempt <= retries) Sys.sleep(min(60, 5 * 2^(attempt - 1)))
      }
      if (success) break
    }
    if (!success) {
      ## Do not leave an error page behind where a cached response is expected.
      unlink(out_path)
      stop("Overpass request failed:\n", paste(failures, collapse = "\n"),
           call. = FALSE)
    }
  }

  osmdata::osmdata_sf(query, doc = out_path)
}


#' Build a metric_graph directly from an OpenStreetMap query
#'
#' @description
#' Convenience wrapper that fetches OSM data with [fetch_osm()] and
#' builds a [metric_graph] from the returned linestrings in one call.
#' `longlat = TRUE` is set automatically because OSM coordinates are
#' always in WGS84 latitude / longitude. Other graph-construction
#' arguments (`tolerance`, `perform_merges`, ...) are forwarded to
#' [metric_graph]'s `$new()` method without overriding its defaults.
#'
#' @inheritParams fetch_osm
#' @param ... Forwarded to [metric_graph]'s `$new()`. Use this to
#'   pass `tolerance`, `perform_merges`, `check_connected`,
#'   `which_longlat`, etc.
#'
#' @return A [metric_graph] object built from the OSM ways.
#'
#' @seealso [fetch_osm()] for the lower-level fetch + parse step.
#'
#' @examples
#' \dontrun{
#' ## Build a freeway graph for a small bbox around San Jose, splitting
#' ## at intersections via the edge_edge tolerance.
#' g <- metric_graph_from_osm(
#'   bbox           = c(-122.10, 37.22, -121.80, 37.45),
#'   value          = c("motorway", "motorway_link"),
#'   perform_merges = TRUE,
#'   tolerance      = list(edge_edge = 1e-5)
#' )
#' g$plot()
#' }
#'
#' @export
metric_graph_from_osm <- function(bbox,
                                  key        = "highway",
                                  value      = NULL,
                                  endpoint   = "https://overpass-api.de/api/interpreter",
                                  timeout    = 180,
                                  user_agent = "MetricGraph-R-package",
                                  cache_path = NULL,
                                  retries    = 2,
                                  ...) {
  osm <- fetch_osm(bbox       = bbox,
                   key        = key,
                   value      = value,
                   endpoint   = endpoint,
                   timeout    = timeout,
                   user_agent = user_agent,
                   cache_path = cache_path,
                   retries    = retries)

  if (is.null(osm$osm_lines) || nrow(osm$osm_lines) == 0L) {
    stop("OSM query returned no linestrings.", call. = FALSE)
  }

  metric_graph$new(edges   = osm$osm_lines,
                   longlat = TRUE,
                   ...)
}
