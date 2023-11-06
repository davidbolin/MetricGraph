#' Traffic speed data from San Jose, California
#'
#' Data set of traffic speed observations on highways in the city of San Jose,
#' California.
#'
#' @format ## `pems`
#' A list with three elements:
#' \describe{
#'   \item{edges}{A `list` object containing the coordinates of the road segments.}
#'   \item{data}{Locations of the observations on the road segments as a
#'   `data.frame` with 325 rows and 2 columns. The first column indicates the edge
#'   number and the second column indicates the distance on edge of the position.}
#'   \item{Y}{Observations of traffic speed. Consists of a vector with
#'   325 observations.}
#' }
#' @source https://www.openstreetmap.org
#' @source https://github.com/spbu-math-cs/Graph-Gaussian-Processes/blob/main/examples/data/PEMS.zip
#' @references Chen, C., K. Petty, A. Skabardonis, P. Varaiya, and Z. Jia (2001). Freeway performance measurement system: mining loop detector data. Transportation Research Record 1748(1), 96–102.
#' @references OpenStreetMap contributors (2017). Planet dump retrieved from https://planet.osm.org. https://www.openstreetmap.org.
"pems"
