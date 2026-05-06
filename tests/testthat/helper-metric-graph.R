# Build a planar n x n square grid as a list of edges (2x2 matrices).
make_grid_edges <- function(n) {
  xs <- ys <- seq_len(n)
  out <- vector("list", 2L * n * (n - 1L))
  k <- 1L
  for (j in seq_len(n)) {
    for (i in seq_len(n - 1L)) {
      out[[k]] <- rbind(c(xs[i], ys[j]), c(xs[i + 1L], ys[j]))
      k <- k + 1L
    }
  }
  for (i in seq_len(n)) {
    for (j in seq_len(n - 1L)) {
      out[[k]] <- rbind(c(xs[i], ys[j]), c(xs[i], ys[j + 1L]))
      k <- k + 1L
    }
  }
  out
}

# Build a small lon/lat grid as an sf object with CRS 4326.
make_longlat_grid <- function(n, lon0 = 2.34, lon1 = 2.36,
                              lat0 = 48.85, lat1 = 48.87) {
  lons <- seq(lon0, lon1, length.out = n)
  lats <- seq(lat0, lat1, length.out = n)
  edges <- list()
  k <- 1L
  for (j in seq_len(n)) {
    for (i in seq_len(n - 1L)) {
      edges[[k]] <- rbind(c(lons[i], lats[j]), c(lons[i + 1L], lats[j]))
      k <- k + 1L
    }
  }
  for (i in seq_len(n)) {
    for (j in seq_len(n - 1L)) {
      edges[[k]] <- rbind(c(lons[i], lats[j]), c(lons[i], lats[j + 1L]))
      k <- k + 1L
    }
  }
  ls_list <- lapply(edges, function(e) sf::st_linestring(unname(e)))
  sf::st_sf(geometry = sf::st_sfc(ls_list, crs = 4326))
}

# Build an open chain graph: a single run of n connected edges v1-v2-v3-...-v(n+1)
# All interior vertices are degree 2 — exactly the case prune_vertices simplifies.
make_open_chain <- function(n) {
  lapply(seq_len(n), function(i) rbind(c(i - 1, 0), c(i, 0)))
}

# Build a triangle (3 edges, all deg-2 — closed loop, tests fallback path).
make_triangle <- function() {
  list(
    rbind(c(0, 0), c(1, 0)),
    rbind(c(1, 0), c(0.5, 1)),
    rbind(c(0.5, 1), c(0, 0))
  )
}

# Build a 5-point star (one center vertex, 5 leaves, 5 edges — no deg-2 to prune).
make_star <- function(n_rays = 5) {
  lapply(seq_len(n_rays), function(k) {
    ang <- 2 * pi * (k - 1) / n_rays
    rbind(c(0, 0), c(cos(ang), sin(ang)))
  })
}

# Helper: check that every edge in the graph has a valid PtE attribute
# (starts at 0, ends at 1, monotone non-decreasing, length = nrow(edge)).
all_pte_valid <- function(g, tol = 1e-12) {
  vapply(g$edges, function(e) {
    pte <- attr(e, "PtE")
    !is.null(pte) &&
      length(pte) == nrow(e) &&
      abs(pte[1]) < tol &&
      abs(pte[length(pte)] - 1) < tol &&
      all(diff(pte) >= -tol)
  }, logical(1))
}

# Helper: compute Euclidean length of a polyline.
polyline_length <- function(coords) {
  if (is.null(coords) || nrow(coords) < 2) return(0)
  d <- coords[-1, , drop = FALSE] - coords[-nrow(coords), , drop = FALSE]
  sum(sqrt(rowSums(d * d)))
}
