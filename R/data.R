#' Traffic speed data from San Jose, California
#'
#' Data set of traffic speed observations on highways in the city of San Jose,
#' California.
#'
#' @format ## `pems`
#' A list with two elements:
#' \describe{
#'   \item{edges}{A `list` object containing the coordinates of the road segments.}
#'   \item{data}{Locations of the observations on the road segments as a
#'   `data.frame` with 325 rows and 3 columns. The first column indicates the edge
#'   number, the second column indicates the distance on edge of the position,
#'  and the third column indicates the average speed observed.}
#' }
#' @source https://www.openstreetmap.org
#' @source https://github.com/spbu-math-cs/Graph-Gaussian-Processes/blob/main/examples/data/PEMS.zip
#' @references Chen, C., K. Petty, A. Skabardonis, P. Varaiya, and Z. Jia (2001). Freeway performance measurement system: mining loop detector data. Transportation Research Record 1748(1), 96-102.
#' @references OpenStreetMap contributors (2017). Planet dump retrieved from https://planet.osm.org. https://www.openstreetmap.org.
"pems"

#' Traffic speed data with replicates from San Jose, California
#'
#' Data set of traffic speed observations on highways in the city of San Jose,
#' California.
#'
#' @format ## `pems_repl`
#' A list with two elements:
#' \describe{
#'   \item{edges}{A `list` object containing the coordinates of the road segments.}
#'   \item{data}{Locations of the observations on the road segments as a
#'   `data.frame` with 325 rows and 4 columns. The first column indicates the observed speed,
#' the second column indicates the edge
#'   number, the third column indicates the distance on edge of the position,
#'  and the fourth column indicates the replicate number.}
#' }
#' @source https://www.openstreetmap.org
#' @source https://github.com/spbu-math-cs/Graph-Gaussian-Processes/blob/main/examples/data/PEMS.zip
#' @references Chen, C., K. Petty, A. Skabardonis, P. Varaiya, and Z. Jia (2001). Freeway performance measurement system: mining loop detector data. Transportation Research Record 1748(1), 96-102.
#' @references OpenStreetMap contributors (2017). Planet dump retrieved from https://planet.osm.org. https://www.openstreetmap.org.
"pems_repl"


#' Largest connected component of the Mid-Columbia stream network
#'
#' A prepared full-graph river component for the directional and symmetric
#' Whittle--Matérn examples. The original 9,521 source rows were imported as
#' 2,758 unique ungrouped graph locations; the largest component retains
#' 18,668 edges and 2,080 observations. The transformation preserves stable
#' source-row identifiers in `columbia_obs_id` and records all selection counts
#' and a content fingerprint in `summary`.
#'
#' @format ## `columbia_main_component`
#' A `columbia_main_component` list with five elements:
#' \describe{
#'   \item{edges}{An `sf` object containing the retained stream edges,
#'   `h2oAreaKm2` directional weights, and edge geometry.}
#'   \item{observations}{A `data.frame` containing normalized edge positions,
#'   `STREAM_AUG`, `ELEV`, `SLOPE`, `PRECIP`, and stable
#'   `columbia_obs_id` values.}
#'   \item{summary}{Component-selection counts and the audited fingerprint.}
#'   \item{weights_name}{The edge attribute used for directional weights.}
#'   \item{fingerprint}{A SHA-256 fingerprint of the prepared component.}
#' }
#' @source Jay Ver Hoef (2023), *MidColumbia.zip - a large spatial stream
#'   network data set*, Figshare,
#'   \doi{10.6084/m9.figshare.24132840.v1}. The source data are licensed under
#'   the Creative Commons Attribution 4.0 International license.
#' @references Isaak, D.J. et al. (2017). The NorWeST Database and Modeled
#'   Summer Temperature Scenarios. *Water Resources Research*, 53(11),
#'   9181--9205. \doi{10.1002/2017WR020969}.
#' @references Ver Hoef, J.M., Dumelle, M., Higham, M., Peterson, E.E., and
#'   Isaak, D.J. (2023). Indexing and Partitioning the Spatial Linear Model for
#'   Large Data Sets. *PLOS ONE*, 18(11), e0291906.
"columbia_main_component"
