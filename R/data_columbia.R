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
