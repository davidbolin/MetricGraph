
#' Convert an `stlpp` object to a metric graph object
#'
#' This function converts an `stlpp` object (from the `stlnpp` package) into a metric graph object.
#'
#' @param stlpp.obj An `stlpp` object to be converted.
#' @param ... Additional arguments to be passed to the `metric_graph` constructor.

#' @return A metric graph object
#' @export
stlpp.to.graph <- function(stlpp.obj,...) {
  return(linnet.to.graph(stlpp.obj$domain,...))
}

#' Convert a `linnet` object to a metric graph object
#'
#' This function converts a `linnet` object (from the `spatstat` package) into a metric graph object.
#'
#' @param linnet.object A `linnet` object to be converted.
#' @param crs The coordinate reference system of the graph.
#' @param ... Additional arguments to be passed to the `metric_graph` constructor.
#'
#' @return A metric graph object with edges defined by the network.
#' @export
#'
linnet.to.graph <- function(linnet.object, crs, ...) {
  # Extract number of edges and vertices
  n <- length(linnet.object$from)
  vertices <- as.data.frame(linnet.object$vertices)
  lines <- vector("list", n)

  # Create LINESTRING geometries for each edge
  for (i in 1:n) {
    coords <- rbind(
      vertices[linnet.object$from[i], 1:2],
      vertices[linnet.object$to[i], 1:2]
    )
    lines[[i]] <- sf::st_linestring(as.matrix(coords))
  }

  # Convert to an sf object
  edges_sf <- sf::st_sf(
    geometry = sf::st_sfc(lines, crs = crs)
  )

  return(metric_graph$new(edges = edges_sf, ...))
}

#' Convert a `psp` object to a metric graph object
#'
#' This function converts a `psp` object (from the `spatstat` package) into a metric graph object.
#'
#' @param psp.object A `psp` object to be converted.
#'
#' @return A metric graph object with edges defined by the segments in the `psp` object.
#' @export
psp.to.graph <- function(psp.object) {
  n <- nrow(psp.object$ends)
  lines <- vector("list", n)

  for (i in 1:n) {
    lines[[i]] <- rbind(
      as.matrix(psp.object$ends[i, 1:2]),
      as.matrix(psp.object$ends[i, 3:4])
    )
  }

  graph <- metric_graph$new(edges = lines)
  return(graph)
}
