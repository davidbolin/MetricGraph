## ============================================================================
## Build an OSM-direction-corrected pems graph.
##
## Outputs (saved next to this script):
##   * pems_osm.rds          - per-edge OSM metadata (highway class,
##                             oneway tag, match distance, reversed,
##                             auth_flip).
##   * pems_osmdir_edges.rds - list of edge coordinate matrices that
##                             defines `pems_graph_osmdir`; load with
##                             metric_graph$new(edges = ., longlat = TRUE).
##
## ============================================================================

library(sf)
library(MetricGraph)
data(pems, package = "MetricGraph")

HIGHWAY_CLASSES <- c("motorway", "motorway_link",
                     "trunk",    "trunk_link",
                     "primary",  "primary_link",
                     "secondary","secondary_link")
## Implicit one-way road classes (in addition to anything explicitly
## tagged `oneway = yes`).
IMPLICIT_ONEWAY <- c("motorway", "motorway_link", "trunk_link",
                     "primary_link", "secondary_link")
## How close an OSM way must be to count as a match (metres). Every
## pems edge in the current data matches within ~3 m.
MATCH_TOL_M <- 50
## UTM zone 10N -- used to compute distances in metres.
CRS_METRIC  <- 32610


# Download OSM data
bbox <- sf::st_bbox(pems$edges) + c(-0.005, -0.005, 0.005, 0.005)
osm_ways <- fetch_osm(bbox = bbox, key = "highway",
                      value = HIGHWAY_CLASSES)$osm_lines

#sf views of pems and OSM, both projected to metres
pems_graph <- metric_graph$new(edges = pems$edges, verbose = 0)
pems_edges_sf <- sf::st_sf(
  edge_id  = seq_along(pems_graph$edges),
  geometry = sf::st_sfc(
    lapply(pems_graph$edges, function(m) {
      m <- unclass(m); class(m) <- NULL
      sf::st_linestring(unname(m))
    }), crs = 4326))

pems_edges_m <- sf::st_transform(pems_edges_sf, CRS_METRIC)
osm_ways_m   <- sf::st_transform(osm_ways,      CRS_METRIC)


# Match each pems edge to its nearest OSM way
match_idx  <- sf::st_nearest_feature(pems_edges_m, osm_ways_m)
match_dist <- as.numeric(sf::st_distance(
  pems_edges_m, osm_ways_m[match_idx, ], by_element = TRUE))

## Compare orientations: each linestring's "direction" is the vector
## from its first coordinate to its last. If pems and OSM point in
## opposite directions, their dot product is negative.
direction_vector <- function(linestring) {
  coords <- sf::st_coordinates(linestring)[, 1:2, drop = FALSE]
  coords[nrow(coords), ] - coords[1, ]
}
pems_dirs <- vapply(seq_len(nrow(pems_edges_sf)),
                    function(i) direction_vector(pems_edges_sf$geometry[[i]]),
                    numeric(2))
osm_dirs  <- vapply(match_idx,
                    function(j) direction_vector(osm_ways$geometry[[j]]),
                    numeric(2))
reversed <- colSums(pems_dirs * osm_dirs) < 0
reversed[match_dist > MATCH_TOL_M] <- FALSE

## A reversal is authoritative when the road is one-way in OSM:
## either explicitly tagged or one of the implicit one-way classes.
highway_class <- ifelse(match_dist > MATCH_TOL_M, "unmatched",
                        as.character(osm_ways$highway[match_idx]))
oneway_tag    <- osm_ways$oneway[match_idx]
is_oneway     <- highway_class %in% IMPLICIT_ONEWAY |
                 oneway_tag %in% c("yes", "true", "1")
auth_flip     <- reversed & is_oneway

osm_meta <- data.frame(
  edge_id   = seq_len(pems_graph$nE),
  highway   = highway_class,
  oneway    = oneway_tag,
  dist_m    = match_dist,
  reversed  = reversed,
  auth_flip = auth_flip)

saveRDS(osm_meta, "pems_osm.rds")

# Helpers for orientation correction

## Reverse the coordinate order of a subset of edges in `edge_list`.
reverse_edges <- function(edge_list, idx) {
  for (i in idx) {
    edge_list[[i]] <- edge_list[[i]][nrow(edge_list[[i]]):1, , drop = FALSE]
  }
  edge_list
}

## Walk every degree-2 chain in `g` and return a logical vector of
## edges whose orientation must be flipped to make every pass-through
## vertex consistent (exactly one inbound and one outbound edge).
propagate_orientation <- function(g) {
  adjacent <- vector("list", g$nV)
  for (e in seq_len(g$nE)) {
    adjacent[[g$E[e, 1]]] <- c(adjacent[[g$E[e, 1]]], e)
    adjacent[[g$E[e, 2]]] <- c(adjacent[[g$E[e, 2]]], e)
  }
  degree  <- lengths(adjacent)
  flipped <- logical(g$nE)
  visited <- logical(g$nE)
  endpoints <- function(e) {
    if (flipped[e]) c(g$E[e, 2], g$E[e, 1]) else g$E[e, ]
  }
  for (seed in seq_len(g$nE)) {
    if (visited[seed]) next
    queue <- seed
    while (length(queue) > 0L) {
      cur <- queue[1L]; queue <- queue[-1L]
      if (visited[cur]) next
      visited[cur] <- TRUE
      cur_ep <- endpoints(cur)
      for (v in cur_ep) {
        if (degree[v] != 2L) next
        other <- setdiff(adjacent[[v]], cur)
        if (length(other) != 1L || visited[other]) next
        ## At v, exactly one edge should be inbound and one outbound.
        ## If both `cur` and `other` come into v (or both leave v),
        ## flip `other`.
        if ((v == cur_ep[2L]) == (v == endpoints(other)[2L])) {
          flipped[other] <- !flipped[other]
        }
        queue <- c(queue, other)
      }
    }
  }
  flipped
}


# Apply OSM flips and BFS cascade
edge_list <- lapply(pems_graph$edges, function(m) {
  m <- unclass(m); class(m) <- NULL
  unname(m)
})
edge_list <- reverse_edges(edge_list, which(auth_flip))

## Diagnostic: count pass-through (deg-2) vertices that became
## inconsistent because of OSM flips on adjacent edges. These are the
## inversions the BFS cascade is about to resolve.
g_intermediate <- metric_graph$new(edges = edge_list, longlat = TRUE,
                                   verbose = 0)
cascade <- propagate_orientation(g_intermediate)
edge_list <- reverse_edges(edge_list, which(cascade))

# Save final graph
osm_graph <- metric_graph$new(edges = edge_list, longlat = TRUE, verbose = 0)

saveRDS(edge_list, "pems_osmdir_edges.rds")

