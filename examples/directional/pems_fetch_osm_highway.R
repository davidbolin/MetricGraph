## --------------------------------------------------------------------------
## Fetch OSM `highway=*` tags for the pems road network and save a small
## per-edge classification next to this script (pems_osm_highway.rds).
##
## The bundled `pems` dataset only carries geometry and unit edge weights,
## so to get real road-class information we have to re-download the
## graph from OpenStreetMap. We use the Overpass API directly via curl
## (one HTTP POST, no rate-limited status checks like `osmdata`) and
## match OSM ways to pems edges by nearest-feature distance (every pems
## edge in the current data set matches an OSM way within ~3 m).
##
## Run this script once; it produces pems_osm_highway.rds, which the
## main LOO comparison loads if present. Re-run if the pems dataset or
## OSM data change.
##
## Requires: an internet connection, jsonlite, sf.
## --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(jsonlite)
  library(sf)
})
library(MetricGraph)
data(pems, package = "MetricGraph")

## Build a padded bbox around the pems graph.
bb <- sf::st_bbox(pems$edges)
bb["xmin"] <- bb["xmin"] - 0.005
bb["xmax"] <- bb["xmax"] + 0.005
bb["ymin"] <- bb["ymin"] - 0.005
bb["ymax"] <- bb["ymax"] + 0.005
cat("Padded bbox:\n"); print(bb)

## Overpass QL: fetch ways with `highway` in {motorway, trunk, primary,
## secondary} plus their _link variants. `out geom;` returns each way's
## node geometry inline so we don't have to issue a second node query.
q <- sprintf(
  paste(
    "[out:json][timeout:90];",
    "(way(%s,%s,%s,%s)",
    "  [highway~\"^(motorway|motorway_link|trunk|trunk_link|primary|primary_link|secondary|secondary_link)$\"];);",
    "out geom;", sep = "\n"),
  bb["ymin"], bb["xmin"], bb["ymax"], bb["xmax"])
tmp_q <- tempfile(fileext = ".txt")
writeLines(q, tmp_q)
tmp_j <- tempfile(fileext = ".json")

ua <- "MetricGraph-research/1.0 (research)"
url <- "https://overpass-api.de/api/interpreter"
cat(sprintf("POSTing query to %s ...\n", url))
status <- system2("curl",
                  c("-s", "-o", shQuote(tmp_j),
                    "-w", "%{http_code}",
                    "-A", shQuote(ua),
                    "--max-time", "180",
                    "--data-binary", paste0("@", tmp_q),
                    shQuote(url)),
                  stdout = TRUE)
cat(sprintf("HTTP %s, %.1f MB downloaded.\n",
            status, file.info(tmp_j)$size / 1024^2))

cat("Parsing JSON...\n")
osm <- jsonlite::fromJSON(tmp_j, simplifyVector = FALSE)
ways <- osm$elements
ways <- ways[vapply(ways, function(w) identical(w$type, "way"), logical(1))]

highway_tag <- vapply(ways, function(w) {
  h <- w$tags$highway
  if (is.null(h)) NA_character_ else as.character(h)
}, character(1))

geoms <- lapply(ways, function(w) {
  g <- w$geometry
  if (is.null(g) || length(g) < 2L) return(NULL)
  m <- do.call(rbind, lapply(g, function(p) c(p$lon, p$lat)))
  sf::st_linestring(m)
})
ok <- !vapply(geoms, is.null, logical(1))
ways_sf <- sf::st_sf(highway  = highway_tag[ok],
                     geometry = sf::st_sfc(geoms[ok], crs = 4326))
cat(sprintf("Built %d OSM linestrings.\n", nrow(ways_sf)))

## Build pems-edges sf and project to UTM zone 10N for metric distances.
pg <- metric_graph$new(edges = pems$edges, verbose = 0)
edges_sf <- sf::st_sf(
  edge_id  = seq_along(pg$edges),
  geometry = sf::st_sfc(
    lapply(pg$edges, function(m) {
      m <- unclass(m); class(m) <- NULL
      sf::st_linestring(unname(m))
    }),
    crs = 4326))
crs_m   <- 32610
ways_p  <- sf::st_transform(ways_sf,  crs_m)
edges_p <- sf::st_transform(edges_sf, crs_m)

## Nearest OSM way per pems edge, with distance.
nidx <- sf::st_nearest_feature(edges_p, ways_p)
ndist <- as.numeric(sf::st_distance(edges_p, ways_p[nidx, ],
                                    by_element = TRUE))
match_tag <- ways_sf$highway[nidx]
match_tag[ndist > 50] <- "unmatched"

cat("\nMatch distance summary (m):\n");      print(summary(ndist))
cat("\nMatched highway class (threshold 50 m):\n"); print(table(match_tag))

out <- data.frame(edge_id = seq_along(match_tag),
                  highway = match_tag,
                  dist_m  = ndist)
out_path <- file.path(dirname(sys.frame(1)$ofile %||% "."),
                      "pems_osm_highway.rds")
if (is.null(out_path) || out_path == "") {
  out_path <- "pems_osm_highway.rds"
}
saveRDS(out, out_path)
cat(sprintf("\nSaved per-edge tags to %s (%d rows).\n",
            out_path, nrow(out)))
