# Metric graph

Class representing a general metric graph.

## Value

Object of [`R6Class`](https://r6.r-lib.org/reference/R6Class.html) for
creating metric graphs.

## Details

A graph object created from vertex and edge matrices, or from an
[`sp::SpatialLines`](https://edzer.github.io/sp/reference/SpatialLines.html)
object where each line is representing and edge. For more details, see
the vignette:
[`vignette("metric_graph", package = "MetricGraph")`](https://davidbolin.github.io/MetricGraph/articles/metric_graph.md)

## Public fields

- `V`:

  Matrix with positions in Euclidean space of the vertices of the graph.

- `nV`:

  The number of vertices.

- `E`:

  Matrix with the edges of the graph, where each row represents an edge,
  `E[i,1]` is the vertex at the start of the ith edge and `E[i,2]` is
  the vertex at the end of the edge.

- `nE`:

  The number of edges.

- `edge_lengths`:

  Vector with the lengths of the edges in the graph.

- `C`:

  Constraint matrix used to set Kirchhoff constraints.

- `CoB`:

  Change-of-basis object used for Kirchhoff constraints.

- `PtV`:

  Vector with the indices of the vertices which are observation
  locations.

- `mesh`:

  Mesh object used for plotting.

- `edges`:

  The coordinates of the edges in the graph.

- `DirectionalWeightFunction_in`:

  Function for inwards weights in directional models

- `DirectionalWeightFunction_out`:

  Function for outwards weights in directional models

- `vertices`:

  The coordinates of the vertices in the graph, along with several
  attributes.

- `geo_dist`:

  Geodesic distances between the vertices in the graph.

- `res_dist`:

  Resistance distances between the observation locations.

- `Laplacian`:

  The weighted graph Laplacian of the vertices in the graph. The weights
  are given by the edge lengths.

- `characteristics`:

  List with various characteristics of the graph.

## Methods

### Public methods

- [`metric_graph$new()`](#method-metric_graph-initialize)

- [`metric_graph$remove_small_circles()`](#method-metric_graph-remove_small_circles)

- [`metric_graph$get_edges()`](#method-metric_graph-get_edges)

- [`metric_graph$get_bounding_box()`](#method-metric_graph-get_bounding_box)

- [`metric_graph$get_vertices()`](#method-metric_graph-get_vertices)

- [`metric_graph$export()`](#method-metric_graph-export)

- [`metric_graph$leaflet()`](#method-metric_graph-leaflet)

- [`metric_graph$mapview()`](#method-metric_graph-mapview)

- [`metric_graph$set_edge_weights()`](#method-metric_graph-set_edge_weights)

- [`metric_graph$get_edge_weights()`](#method-metric_graph-get_edge_weights)

- [`metric_graph$get_vertices_incomp_dir()`](#method-metric_graph-get_vertices_incomp_dir)

- [`metric_graph$summary()`](#method-metric_graph-summary)

- [`metric_graph$print()`](#method-metric_graph-print)

- [`metric_graph$compute_characteristics()`](#method-metric_graph-compute_characteristics)

- [`metric_graph$check_euclidean()`](#method-metric_graph-check_euclidean)

- [`metric_graph$check_distance_consistency()`](#method-metric_graph-check_distance_consistency)

- [`metric_graph$compute_geodist()`](#method-metric_graph-compute_geodist)

- [`metric_graph$compute_geodist_PtE()`](#method-metric_graph-compute_geodist_PtE)

- [`metric_graph$compute_geodist_mesh()`](#method-metric_graph-compute_geodist_mesh)

- [`metric_graph$compute_resdist()`](#method-metric_graph-compute_resdist)

- [`metric_graph$compute_resdist_PtE()`](#method-metric_graph-compute_resdist_PtE)

- [`metric_graph$get_degrees()`](#method-metric_graph-get_degrees)

- [`metric_graph$compute_PtE_edges()`](#method-metric_graph-compute_PtE_edges)

- [`metric_graph$compute_resdist_mesh()`](#method-metric_graph-compute_resdist_mesh)

- [`metric_graph$compute_laplacian()`](#method-metric_graph-compute_laplacian)

- [`metric_graph$prune_vertices()`](#method-metric_graph-prune_vertices)

- [`metric_graph$set_manual_edge_lengths()`](#method-metric_graph-set_manual_edge_lengths)

- [`metric_graph$get_groups()`](#method-metric_graph-get_groups)

- [`metric_graph$get_PtE()`](#method-metric_graph-get_PtE)

- [`metric_graph$get_edge_lengths()`](#method-metric_graph-get_edge_lengths)

- [`metric_graph$get_locations()`](#method-metric_graph-get_locations)

- [`metric_graph$observation_to_vertex()`](#method-metric_graph-observation_to_vertex)

- [`metric_graph$edgeweight_to_data()`](#method-metric_graph-edgeweight_to_data)

- [`metric_graph$get_mesh_locations()`](#method-metric_graph-get_mesh_locations)

- [`metric_graph$clear_observations()`](#method-metric_graph-clear_observations)

- [`metric_graph$process_data()`](#method-metric_graph-process_data)

- [`metric_graph$add_observations()`](#method-metric_graph-add_observations)

- [`metric_graph$mutate_weights()`](#method-metric_graph-mutate_weights)

- [`metric_graph$select_weights()`](#method-metric_graph-select_weights)

- [`metric_graph$filter_weights()`](#method-metric_graph-filter_weights)

- [`metric_graph$summarise_weights()`](#method-metric_graph-summarise_weights)

- [`metric_graph$drop_na_weights()`](#method-metric_graph-drop_na_weights)

- [`metric_graph$mutate()`](#method-metric_graph-mutate)

- [`metric_graph$drop_na()`](#method-metric_graph-drop_na)

- [`metric_graph$select()`](#method-metric_graph-select)

- [`metric_graph$filter()`](#method-metric_graph-filter)

- [`metric_graph$summarise()`](#method-metric_graph-summarise)

- [`metric_graph$get_data()`](#method-metric_graph-get_data)

- [`metric_graph$setDirectionalWeightFunction()`](#method-metric_graph-setDirectionalWeightFunction)

- [`metric_graph$buildDirectionalConstraints()`](#method-metric_graph-buildDirectionalConstraints)

- [`metric_graph$buildC()`](#method-metric_graph-buildC)

- [`metric_graph$build_mesh()`](#method-metric_graph-build_mesh)

- [`metric_graph$get_version()`](#method-metric_graph-get_version)

- [`metric_graph$is_disconnected()`](#method-metric_graph-is_disconnected)

- [`metric_graph$get_components()`](#method-metric_graph-get_components)

- [`metric_graph$which_component()`](#method-metric_graph-which_component)

- [`metric_graph$compute_fem()`](#method-metric_graph-compute_fem)

- [`metric_graph$compute_mesh_weights()`](#method-metric_graph-compute_mesh_weights)

- [`metric_graph$mesh_A()`](#method-metric_graph-mesh_A)

- [`metric_graph$fem_basis()`](#method-metric_graph-fem_basis)

- [`metric_graph$VtEfirst()`](#method-metric_graph-VtEfirst)

- [`metric_graph$plot()`](#method-metric_graph-plot)

- [`metric_graph$plot_connections()`](#method-metric_graph-plot_connections)

- [`metric_graph$is_tree()`](#method-metric_graph-is_tree)

- [`metric_graph$plot_function()`](#method-metric_graph-plot_function)

- [`metric_graph$plot_movie()`](#method-metric_graph-plot_movie)

- [`metric_graph$add_mesh_observations()`](#method-metric_graph-add_mesh_observations)

- [`metric_graph$get_initial_graph()`](#method-metric_graph-get_initial_graph)

- [`metric_graph$update_graph()`](#method-metric_graph-update_graph)

- [`metric_graph$coordinates()`](#method-metric_graph-coordinates)

- [`metric_graph$clone()`](#method-metric_graph-clone)

------------------------------------------------------------------------

### `metric_graph$new()`

Create a new `metric_graph` object.

#### Usage

    metric_graph$new(
      edges = NULL,
      V = NULL,
      E = NULL,
      vertex_unit = NULL,
      length_unit = NULL,
      edge_weights = NULL,
      kirchhoff_weights = NULL,
      directional_weights = NULL,
      longlat = NULL,
      crs = NULL,
      proj4string = NULL,
      which_longlat = "sp",
      include_obs = NULL,
      include_edge_weights = NULL,
      project = FALSE,
      project_data = FALSE,
      which_projection = "Winkel tripel",
      manual_edge_lengths = NULL,
      perform_merges = NULL,
      approx_edge_PtE = TRUE,
      tolerance = list(vertex_vertex = 0.001, vertex_edge = 0.001, edge_edge = 0),
      check_connected = TRUE,
      remove_deg2 = FALSE,
      merge_close_vertices = NULL,
      factor_merge_close_vertices = 1,
      remove_circles = FALSE,
      auto_remove_point_edges = TRUE,
      verbose = 1,
      add_obs_options = list(return_removed = FALSE, verbose = verbose),
      lines = deprecated(),
      .assemble = NULL
    )

#### Arguments

- `edges`:

  A list containing coordinates as `m x 2` matrices (that is, of
  `matrix` type) or m x 2 data frames (`data.frame` type) of sequence of
  points connected by straightlines. Alternatively, you can also prove
  an object of type `SSN`, `osmdata_sp`, `osmdata_sf`,
  `SpatialLinesDataFrame` or `SpatialLines` (from `sp` package) or
  `MULTILINESTRING` (from `sf` package).

- `V`:

  n x 2 matrix with Euclidean coordinates of the n vertices. If
  non-NULL, no merges will be performed.

- `E`:

  m x 2 matrix where each row represents one of the m edges. If
  non-NULL, no merges will be performed.

- `vertex_unit`:

  The unit in which the vertices are specified. The options are 'degree'
  (the great circle distance in km), 'km', 'm' and 'miles'. The default
  is `NULL`, which means no unit. However, if you set `length_unit`, you
  need to set `vertex_unit`.

- `length_unit`:

  The unit in which the lengths will be computed. The options are 'km',
  'm' and 'miles'. The default, when longlat is `TRUE`, or an `sf` or
  `sp` objects are provided, is 'km'.

- `edge_weights`:

  Either a number, a numerical vector with length given by the number of
  edges, providing the edge weights, or a `data.frame` with the number
  of rows being equal to the number of edges, where each row gives a
  vector of weights to its corresponding edge. Can be changed by using
  the `set_edge_weights()` method.

- `kirchhoff_weights`:

  If non-null, the name (or number) of the column of `edge_weights` that
  contain the Kirchhoff weights. Must be equal to 1 (or `TRUE`) in case
  `edge_weights` is a single number and those are the Kirchhoff weights.

- `directional_weights`:

  If non-null, the name (or number) of the column of `edge_weights` that
  contain the directional weights. The default is the first column of
  the edge weights.

- `longlat`:

  There are three options: `NULL`, `TRUE` or `FALSE`. If `NULL` (the
  default option), the `edges` argument will be checked to see if there
  is a CRS or proj4string available, if so, `longlat` will be set to
  `TRUE`, otherwise, it will be set to `FALSE`. If `TRUE`, then it is
  assumed that the coordinates are given. in Longitude/Latitude and that
  distances should be computed in meters. If `TRUE` it takes precedence
  over `vertex_unit` and `length_unit`, and is equivalent to
  `vertex_unit = 'degree'` and `length_unit = 'm'`.

- `crs`:

  Coordinate reference system to be used in case `longlat` is set to
  `TRUE` and `which_longlat` is `sf`. Object of class crs. The default
  choice, if the `edges` object does not have CRS nor proj4string, is
  `sf::st_crs(4326)`.

- `proj4string`:

  Projection string of class CRS-class to be used in case `longlat` is
  set to `TRUE` and `which_longlat` is `sp`. The default choice, if the
  `edges` object does not have CRS nor proj4string, is
  `sp::CRS("+proj=longlat +datum=WGS84")`.

- `which_longlat`:

  Compute the distance using which package? The options are `sp` and
  `sf`. The default is `sp`.

- `include_obs`:

  If the object is of class `SSN`, should the observations be added? If
  `NULL` and the edges are of class `SSN`, the data will be
  automatically added. If `FALSE`, the data will not be added.
  Alternatively, one can set this argument to the numbers or names of
  the columns of the observations to be added as observations.

- `include_edge_weights`:

  If the object is of class `SSN`, `osmdata_sp`, `osmdata_sf`,
  `SpatialLinesDataFrame`, `MULTILINESTRING`, `LINESTRING`,
  `sfc_LINESTRING`, `sfc_MULTILINESTRING`, should the edge data (if any)
  be added as edge weights? If `NULL`, the edge data will be added as
  edge weights, if `FALSE` they will not be added. Alternatively, one
  can set this argument to the numbers or names of the columns of the
  edge data to be added as edge weights.

- `project`:

  If `longlat` is `TRUE` should a projection be used to compute the
  distances to be used for the tolerances (see `tolerance` below)? The
  default is `FALSE`. When `TRUE`, the construction of the graph is
  faster.

- `project_data`:

  If `longlat` is `TRUE` should the vertices be project to planar
  coordinates? The default is `FALSE`. When `TRUE`, the construction of
  the graph is faster.

- `which_projection`:

  Which projection should be used in case `project` is `TRUE`? The
  options are `Robinson`, `Winkel tripel` or a proj4string. The default
  is `Winkel tripel`.

- `manual_edge_lengths`:

  If non-NULL, a vector containing the edges lengths, and all the
  quantities related to edge lengths will be computed in terms of these.
  If merges are performed, it is likely that the merges will override
  the manual edge lengths. In such a case, to provide manual edge
  lengths, one should either set the `perform_merges` argument to
  `FALSE` or use the `set_manual_edge_lengths()` method.

- `perform_merges`:

  There are three options, `NULL`, `TRUE` or `FALSE`. The default option
  is `NULL`. If `NULL`, it will be set to `FALSE` unless 'edges', 'V'
  and 'E' are `NULL`, in which case it will be set to `TRUE`. If FALSE,
  this will take priority over the other arguments, and no merges
  (except the optional `merge_close_vertices` below) will be performed.
  Note that the merge on the additional `merge_close_vertices` might
  still be performed, if it is set to `TRUE`.

- `approx_edge_PtE`:

  Should the relative positions on the edges be approximated? The
  default is `TRUE`. If `FALSE`, the speed can be considerably slower,
  especially for large metric graphs.

- `tolerance`:

  List that provides tolerances during the construction of the graph:

  - `vertex_vertex` Vertices that are closer than this number are merged
    (default = 1e-7).

  - `vertex_edge` If a vertex at the end of one edge is closer than this
    number to another edge, this vertex is connected to that edge
    (default = 1e-7). Previously `vertex_line`, which is now deprecated.

  - `edge_edge` If two edges at some point are closer than this number,
    a new vertex is added at that point and the two edges are connected
    (default = 0).

  - `vertex_line`, Deprecated. Use `vertex_edge` instead.

  - `line_line`, Deprecated. Use `edge_edge` instead.

  In case `longlat = TRUE`, the tolerances are given in `length_unit`.

- `check_connected`:

  If `TRUE`, it is checked whether the graph is connected and a warning
  is given if this is not the case.

- `remove_deg2`:

  Set to `TRUE` to remove all vertices of degree 2 in the
  initialization. Default is `FALSE`.

- `merge_close_vertices`:

  Should an additional step to merge close vertices be done? The options
  are `NULL` (the default), `TRUE` or `FALSE`. If `NULL`, it will be
  determined automatically. If `TRUE` this step will be performed even
  if `perfom_merges` is set to `FALSE`.

- `factor_merge_close_vertices`:

  Which factor to be multiplied by tolerance `vertex_vertex` when
  merging close vertices at the additional step?

- `remove_circles`:

  All circlular edges with a length smaller than this number are
  removed. If `TRUE`, the `vertex_vertex` tolerance will be used. If
  `FALSE`, no circles will be removed.

- `auto_remove_point_edges`:

  Should edges of length zero, that is, edges that are actually points,
  be automatically removed?

- `verbose`:

  Print progress of graph creation. There are 3 levels of verbose, level
  0, 1 and 2. In level 0, no messages are printed. In level 1, only
  messages regarding important steps are printed. Finally, in level 2,
  messages detailing all the steps are printed. The default is 1.

- `add_obs_options`:

  List containing additional options to be passed to the
  `add_observations()` method when adding observations from `SSN` data?

- `lines`:

  **\[deprecated\]** Use `edges` instead.

- `.assemble`:

  Expert option for assembling the graph by skipping the entire
  constructor pipeline and instead populate the R6 object directly from
  a pre-computed list. Intended for internal use.

#### Details

A graph object can be initialized in two ways. The first method is to
specify V and E. In this case, all edges are assumed to be straight
lines. The second option is to specify the graph via the `lines` input.
In this case, the vertices are set by the end points of the lines. Thus,
if two lines are intersecting somewhere else, this will not be viewed as
a vertex.

#### Returns

A `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$remove_small_circles()`

Sets the edge weights

#### Usage

    metric_graph$remove_small_circles(tolerance, verbose = 1)

#### Arguments

- `tolerance`:

  Tolerance at which circles with length less than this will be removed.

- `verbose`:

  Print progress of graph creation. There are 3 levels of verbose, level
  0, 1 and 2. In level 0, no messages are printed. In level 1, only
  messages regarding important steps are printed. Finally, in level 2,
  messages detailing all the steps are printed. The default is 1.

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$get_edges()`

Exports the edges of the MetricGraph object as an `sf` or `sp`.

#### Usage

    metric_graph$get_edges(format = c("sf", "sp", "list"))

#### Arguments

- `format`:

  The format for the exported object. The options are `sf` (default),
  `sp` and `list`.

#### Returns

For `format == "sf"`, the function returns an `sf` object of
`LINESTRING` geometries, where the associated data frame includes edge
weights.

For `format == "sp"`, the function returns a `SpatialLinesDataFrame`
where the data frame includes edge weights.

------------------------------------------------------------------------

### `metric_graph$get_bounding_box()`

Bounding box of the metric graph

#### Usage

    metric_graph$get_bounding_box(format = "sf")

#### Arguments

- `format`:

  If the metric graph has a coordinate reference system, the format for
  the exported object. The options are `sf` (default), `sp` and
  `matrix`.

#### Returns

A bounding box of the metric graph

------------------------------------------------------------------------

### `metric_graph$get_vertices()`

Exports the vertices of the MetricGraph object as an `sf`, `sp` or as a
matrix.

#### Usage

    metric_graph$get_vertices(format = c("sf", "sp", "list"))

#### Arguments

- `format`:

  The format for the exported object. The options are `sf` (default),
  `sp` and `matrix`.

#### Returns

For `which_format == "sf"`, the function returns an `sf` object of
`POINT` geometries.

For `which_format == "sp"`, the function returns a
`SpatialPointsDataFrame` object.

------------------------------------------------------------------------

### `metric_graph$export()`

Exports the MetricGraph object as an `sf` or `sp` object.

#### Usage

    metric_graph$export(format = "sf")

#### Arguments

- `format`:

  The format for the exported object. The options are `sf` (default) and
  `sp`.

#### Returns

Returns a list with three elements: `edges`, `vertices`, and `data`.

For `format == "sf"`, `edges` is an `sf` object of `LINESTRING`
geometries with edge weights, and `vertices` and `data` are `sf` objects
with `POINT` geometries.

For `format == "sp"`, `edges` is a `SpatialLinesDataFrame` with edge
weights, and `vertices` and `data` are `SpatialPointsDataFrame`.

------------------------------------------------------------------------

### `metric_graph$leaflet()`

Return the metric graph as a
[`leaflet::leaflet()`](https://rstudio.github.io/leaflet/reference/leaflet.html)
object to be built upon.

#### Usage

    metric_graph$leaflet(
      width = NULL,
      height = NULL,
      padding = 0,
      options = leafletOptions(),
      elementId = NULL,
      sizingPolicy = leafletSizingPolicy(padding = padding)
    )

#### Arguments

- `width`:

  the width of the map

- `height`:

  the height of the map

- `padding`:

  the padding of the map

- `options`:

  the map options

- `elementId`:

  Use an explicit element ID for the widget (rather than an
  automatically generated one).

- `sizingPolicy`:

  htmlwidgets sizing policy object. Defaults to `leafletSizingPolicy()`.

------------------------------------------------------------------------

### `metric_graph$mapview()`

Returns a
[`mapview::mapview()`](https://r-spatial.github.io/mapview/reference/mapView.html)
object of the metric graph

#### Usage

    metric_graph$mapview(...)

#### Arguments

- `...`:

  Additional arguments to be passed to
  [`mapview::mapview()`](https://r-spatial.github.io/mapview/reference/mapView.html).
  The `x` argument of mapview, containing the metric graph is already
  passed internally.

------------------------------------------------------------------------

### `metric_graph$set_edge_weights()`

Sets the edge weights

#### Usage

    metric_graph$set_edge_weights(
      weights = NULL,
      kirchhoff_weights = NULL,
      directional_weights = NULL,
      verbose = 0
    )

#### Arguments

- `weights`:

  Either a number, a numerical vector with length given by the number of
  edges, providing the edge weights, or a `data.frame` with the number
  of rows being equal to the number of edges, where each row gives a
  vector of weights to its corresponding edge.

- `kirchhoff_weights`:

  If non-null, the name (or number) of the column of `weights` that
  contain the Kirchhoff weights. Must be equal to 1 (or `TRUE`) in case
  `weights` is a single number and those are the Kirchhoff weights.

- `directional_weights`:

  If non-null, the name (or number) of the column of `weights` that
  contain the directional weights.

- `verbose`:

  There are 3 levels of verbose, level 0, 1 and 2. In level 0, no
  messages are printed. In level 1, only messages regarding important
  steps are printed. Finally, in level 2, messages detailing all the
  steps are printed. The default is 1.

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$get_edge_weights()`

Gets the edge weights

#### Usage

    metric_graph$get_edge_weights(
      data.frame = FALSE,
      format = c("tibble", "sf", "sp", "list"),
      tibble = deprecated()
    )

#### Arguments

- `data.frame`:

  If the edge weights are given as vectors, should the result be
  returned as a data.frame?

- `format`:

  Which format should the data be returned? The options are `tibble` for
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
  `sf` for `POINT`, `sp` for `SpatialPointsDataFrame` and `list` for the
  internal list format.

- `tibble`:

  **\[deprecated\]** Use `format` instead.

#### Returns

A vector or `data.frame` containing the edge weights.

------------------------------------------------------------------------

### `metric_graph$get_vertices_incomp_dir()`

Gets vertices with incompatible directions

#### Usage

    metric_graph$get_vertices_incomp_dir()

#### Returns

A vector containing the vertices with incompatible directions.

------------------------------------------------------------------------

### `metric_graph$summary()`

Prints a summary of various informations of the graph

#### Usage

    metric_graph$summary(
      messages = FALSE,
      compute_characteristics = NULL,
      check_euclidean = NULL,
      check_distance_consistency = NULL
    )

#### Arguments

- `messages`:

  Should message explaining how to build the results be given for
  missing quantities?

- `compute_characteristics`:

  Should the characteristics of the graph be computed? If `NULL` it will
  be determined based on the size of the graph.

- `check_euclidean`:

  Check if the graph has Euclidean edges? If `NULL` it will be
  determined based on the size of the graph.

- `check_distance_consistency`:

  Check the distance consistency assumption? If `NULL` it will be
  determined based on the size of the graph.

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$print()`

Prints various characteristics of the graph

#### Usage

    metric_graph$print()

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$compute_characteristics()`

Computes various characteristics of the graph

#### Usage

    metric_graph$compute_characteristics(check_euclidean = FALSE)

#### Arguments

- `check_euclidean`:

  Also check if the graph has Euclidean edges? This essentially means
  that the distance consistency check will also be perfomed. If the
  graph does not have Euclidean edges due to another reason rather than
  the distance consistency, then it will already be indicated that the
  graph does not have Euclidean edges.

#### Returns

No return value. Called for its side effects. The computed
characteristics are stored in the `characteristics` element of the
`metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$check_euclidean()`

Check if the graph has Euclidean edges.

#### Usage

    metric_graph$check_euclidean()

#### Returns

Returns `TRUE` if the graph has Euclidean edges, or `FALSE` otherwise.
The result is stored in the `characteristics` element of the
`metric_graph` object. The result is displayed when the graph is
printed.

------------------------------------------------------------------------

### `metric_graph$check_distance_consistency()`

Checks distance consistency of the graph.

#### Usage

    metric_graph$check_distance_consistency()

#### Returns

No return value. The result is stored in the `characteristics` element
of the `metric_graph` object. The result is displayed when the graph is
printed.

------------------------------------------------------------------------

### `metric_graph$compute_geodist()`

Computes shortest path distances between the vertices in the graph

#### Usage

    metric_graph$compute_geodist(
      obs = TRUE,
      include_vertices = FALSE,
      all_groups = FALSE,
      group = NULL,
      verbose = 0,
      full = lifecycle::deprecated()
    )

#### Arguments

- `obs`:

  Should the geodesic distances be computed at the observation locations
  or only at vertices?

- `include_vertices`:

  If `obs` is `TRUE`, should the vertex locations be included in the
  resulting distance matrix?

- `all_groups`:

  Should the geodesic distances be computed for all the available
  locations accross all groups? If `FALSE`, it will be computed
  separately for the locations of each group.

- `group`:

  Vector or list containing which groups to compute the distance for. If
  `NULL`, it will be computed for all groups.

- `verbose`:

  Print progress of the computation of the geodesic distances. There are
  3 levels of verbose, level 0, 1 and 2. In level 0, no messages are
  printed. In level 1, only messages regarding important steps are
  printed. Finally, in level 2, messages detailing all the steps are
  printed. The default is 1.

- `full`:

  **\[deprecated\]** Use `all_groups` instead.

#### Returns

No return value. Called for its side effects. The computed geodesic
distances are stored in the `geo_dist` element of the `metric_graph`
object.

------------------------------------------------------------------------

### `metric_graph$compute_geodist_PtE()`

Computes shortest path distances between the vertices in the graph.

#### Usage

    metric_graph$compute_geodist_PtE(
      PtE,
      normalized = TRUE,
      include_vertices = TRUE,
      verbose = 0
    )

#### Arguments

- `PtE`:

  Points to compute the metric for.

- `normalized`:

  are the locations in PtE in normalized distance?

- `include_vertices`:

  Should the original vertices be included in the distance matrix?

- `verbose`:

  Print progress of the computation of the geodesic distances. There are
  3 levels of verbose, level 0, 1 and 2. In level 0, no messages are
  printed. In level 1, only messages regarding important steps are
  printed. Finally, in level 2, messages detailing all the steps are
  printed. The default is 1.

#### Returns

A matrix containing the geodesic distances.

------------------------------------------------------------------------

### `metric_graph$compute_geodist_mesh()`

Computes shortest path distances between the vertices in the mesh.

#### Usage

    metric_graph$compute_geodist_mesh()

#### Returns

No return value. Called for its side effects. The geodesic distances on
the mesh are stored in `mesh$geo_dist` in the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$compute_resdist()`

Computes the resistance distance between the observation locations.

#### Usage

    metric_graph$compute_resdist(
      full = FALSE,
      obs = TRUE,
      group = NULL,
      check_euclidean = FALSE,
      include_vertices = FALSE,
      verbose = 0
    )

#### Arguments

- `full`:

  Should the resistance distances be computed for all the available
  locations. If `FALSE`, it will be computed separately for the
  locations of each group.

- `obs`:

  Should the resistance distances be computed at the observation
  locations?

- `group`:

  Vector or list containing which groups to compute the distance for. If
  `NULL`, it will be computed for all groups.

- `check_euclidean`:

  Check if the graph used to compute the resistance distance has
  Euclidean edges? The graph used to compute the resistance distance has
  the observation locations as vertices.

- `include_vertices`:

  Should the vertices of the graph be also included in the resulting
  matrix when using `FULL=TRUE`?

- `verbose`:

  Print progress of the computation of the resistance distances. There
  are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are
  printed. In level 1, only messages regarding important steps are
  printed. Finally, in level 2, messages detailing all the steps are
  printed. The default is 1.

#### Returns

No return value. Called for its side effects. The geodesic distances are
stored in the `res_dist` element of the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$compute_resdist_PtE()`

Computes the resistance distance between the observation locations.

#### Usage

    metric_graph$compute_resdist_PtE(
      PtE,
      normalized = TRUE,
      include_vertices = FALSE,
      check_euclidean = FALSE,
      verbose = 0
    )

#### Arguments

- `PtE`:

  Points to compute the metric for.

- `normalized`:

  Are the locations in PtE in normalized distance?

- `include_vertices`:

  Should the original vertices be included in the Laplacian matrix?

- `check_euclidean`:

  Check if the graph used to compute the resistance distance has
  Euclidean edges? The graph used to compute the resistance distance has
  the observation locations as vertices.

- `verbose`:

  Print progress of the computation of the resistance distances. There
  are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are
  printed. In level 1, only messages regarding important steps are
  printed. Finally, in level 2, messages detailing all the steps are
  printed. The default is 1.

#### Returns

A matrix containing the resistance distances.

------------------------------------------------------------------------

### `metric_graph$get_degrees()`

Returns the degrees of the vertices in the metric graph.

#### Usage

    metric_graph$get_degrees(which = "degree")

#### Arguments

- `which`:

  If "degree", returns the degree of the vertex. If "indegree", returns
  the indegree, and if "outdegree", it returns the outdegree.

#### Returns

A vector containing the degrees of the vertices.

------------------------------------------------------------------------

### `metric_graph$compute_PtE_edges()`

Computes the relative positions of the coordinates of the edges and save
it as an attribute to each edge. This improves the quality of plots
obtained by the `plot_function()` method, however it might be costly to
compute.

#### Usage

    metric_graph$compute_PtE_edges(approx = TRUE, verbose = 0)

#### Arguments

- `approx`:

  Should the computation of the relative positions be approximate?
  Default is `TRUE`. If `FALSE`, the speed can be considerably slower,
  especially for large metric graphs.

- `verbose`:

  Level of verbosity, 0, 1 or 2. The default is 0.

#### Returns

No return value, called for its side effects.

------------------------------------------------------------------------

### `metric_graph$compute_resdist_mesh()`

Computes the resistance metric between the vertices in the mesh.

#### Usage

    metric_graph$compute_resdist_mesh()

#### Returns

No return value. Called for its side effects. The geodesic distances on
the mesh are stored in the `mesh$res_dist` element in the `metric_graph`
object.

------------------------------------------------------------------------

### `metric_graph$compute_laplacian()`

Computes the weigthed graph Laplacian for the graph.

#### Usage

    metric_graph$compute_laplacian(
      full = FALSE,
      obs = TRUE,
      group = NULL,
      verbose = 0
    )

#### Arguments

- `full`:

  Should the resistance distances be computed for all the available
  locations. If `FALSE`, it will be computed separately for the
  locations of each group.

- `obs`:

  Should the resistance distances be computed at the observation
  locations? It will only compute for locations in which there is at
  least one observations that is not NA.

- `group`:

  Vector or list containing which groups to compute the Laplacian for.
  If `NULL`, it will be computed for all groups.

- `verbose`:

  Print progress of the computation of the Laplacian. There are 3 levels
  of verbose, level 0, 1 and 2. In level 0, no messages are printed. In
  level 1, only messages regarding important steps are printed. Finally,
  in level 2, messages detailing all the steps are printed. The default
  is 1.

#### Returns

No reutrn value. Called for its side effects. The Laplacian is stored in
the `Laplacian` element in the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$prune_vertices()`

Removes vertices of degree 2 from the metric graph.

#### Usage

    metric_graph$prune_vertices(
      check_weights = TRUE,
      check_circles = TRUE,
      verbose = FALSE
    )

#### Arguments

- `check_weights`:

  If `TRUE` will only prune edges with different weights.

- `check_circles`:

  If `TRUE` will not prune a vertex such that the resulting edge is a
  circle.

- `verbose`:

  Print progress of pruning. There are 3 levels of verbose, level 0, 1
  and 2. In level 0, no messages are printed. In level 1, only messages
  regarding important steps are printed. Finally, in level 2, messages
  detailing all the steps are printed. The default is 1.

#### Details

Vertices of degree 2 are removed as long as the corresponding edges that
would be merged are compatible in terms of direction.

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$set_manual_edge_lengths()`

Gets the groups from the data.

#### Usage

    metric_graph$set_manual_edge_lengths(edge_lengths, unit = NULL)

#### Arguments

- `edge_lengths`:

  edge lengths to be set to the metric graph edges.

- `unit`:

  set or override the edge lengths unit.

#### Returns

does not return anything. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$get_groups()`

Gets the groups from the data.

#### Usage

    metric_graph$get_groups(get_cols = FALSE)

#### Arguments

- `get_cols`:

  Should the names of the columns that created the group variable be
  returned?

#### Returns

A vector containing the available groups in the internal data.

------------------------------------------------------------------------

### `metric_graph$get_PtE()`

Gets PtE from the data.

#### Usage

    metric_graph$get_PtE()

#### Returns

A matrix with two columns, where the first column contains the edge
number and the second column contains the distance on edge of the
observation locations.

------------------------------------------------------------------------

### `metric_graph$get_edge_lengths()`

Gets the edge lengths with the corresponding unit.

#### Usage

    metric_graph$get_edge_lengths(unit = NULL)

#### Arguments

- `unit`:

  If non-NULL, changes from `length_unit` from the graph construction to
  `unit`.

#### Returns

a vector with the length unit (if the graph was constructed with a
length unit).

------------------------------------------------------------------------

### `metric_graph$get_locations()`

Gets the spatial locations from the data.

#### Usage

    metric_graph$get_locations()

#### Returns

A `data.frame` object with observation locations. If `longlat = TRUE`,
the column names are lon and lat, otherwise the column names are x and
y.

------------------------------------------------------------------------

### `metric_graph$observation_to_vertex()`

Adds observation locations as vertices in the graph.

#### Usage

    metric_graph$observation_to_vertex(
      mesh_warning = TRUE,
      verbose = 0,
      tolerance = deprecated()
    )

#### Arguments

- `mesh_warning`:

  Display a warning if the graph structure change and the metric graph
  has a mesh object.

- `verbose`:

  Print progress of the steps when adding observations. There are 3
  levels of verbose, level 0, 1 and 2. In level 0, no messages are
  printed. In level 1, only messages regarding important steps are
  printed. Finally, in level 2, messages detailing all the steps are
  printed. The default is 1.

- `tolerance`:

  **\[deprecated\]**. Not used anymore

- `share_weights`:

  Should the same weight be shared among the split edges? If `FALSE`,
  the weights will be removed, and a common weight given by 1 will be
  given.

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$edgeweight_to_data()`

Turns edge weights into data on the metric graph

#### Usage

    metric_graph$edgeweight_to_data(
      loc = NULL,
      mesh = FALSE,
      data_loc = FALSE,
      weight_col = NULL,
      add = TRUE,
      data_coords = c("PtE", "spatial"),
      normalized = FALSE,
      tibble = FALSE,
      format = c("tibble", "sf", "sp", "list"),
      verbose = 1,
      suppress_warnings = FALSE,
      return = FALSE
    )

#### Arguments

- `loc`:

  A `matrix` or `data.frame` with two columns containing the locations
  to generate the data from the edge weights. If `data_coords` is
  'spatial', the first column must be the x-coordinate of the data, and
  the second column must be the y-coordinate. If `data_coords` is 'PtE',
  the first column must be the edge number and the second column must be
  the distance on edge.

- `mesh`:

  Should the data be generated to the mesh locations? In this case, the
  `loc` argument will be ignored. Observe that the metric graph must
  have a mesh built for one to use this option. CAUTION: To add
  edgeweight to data to both the data locations and mesh locations,
  please, add at the data locations first, then to mesh locations.

- `data_loc`:

  Should the data be generated to the data locations? In this case, the
  `loc` argument will be ignored. Observe that the metric graph must
  have data for one to use this option. CAUTION: To add edgeweight to
  data to both the data locations and mesh locations, please, add at the
  data locations first, then to mesh locations.

- `weight_col`:

  Which columns of the edge weights should be turned into data? If
  `NULL`, all columns will be turned into data.

- `add`:

  Should the data generated be added to the metric graph internal data?

- `data_coords`:

  To be used only if `mesh` is `FALSE`. It decides which coordinate
  system to use. If `PtE`, the user must provide `edge_number` and
  `distance_on_edge`, otherwise if `spatial`, the user must provide
  `coord_x` and `coord_y`.

- `normalized`:

  if TRUE, then the distances in `distance_on_edge` are assumed to be
  normalized to (0,1). Default FALSE.

- `tibble`:

  Should the data be returned in a `tibble` format?

- `format`:

  If `return` is `TRUE`, the format of the output: "tibble", "sf", or
  "sp". Default is "tibble".

- `verbose`:

  Print progress of the steps when adding observations. There are 3
  levels of verbose, level 0, 1 and 2. In level 0, no messages are
  printed. In level 1, only messages regarding important steps are
  printed. Finally, in level 2, messages detailing all the steps are
  printed. The default is 1.

- `suppress_warnings`:

  Suppress warnings related to duplicated observations?

- `return`:

  Should the data be returned? If `return_removed` is `TRUE`, only the
  removed locations will be return (if there is any).

------------------------------------------------------------------------

### `metric_graph$get_mesh_locations()`

Returns a list or a matrix with the mesh locations.

#### Usage

    metric_graph$get_mesh_locations(
      bru = FALSE,
      loc = c(".edge_number", ".distance_on_edge"),
      loc_name = NULL,
      normalized = TRUE
    )

#### Arguments

- `bru`:

  Should an 'inlabru'-friendly list be returned?

- `loc`:

  If `bru` is set to `TRUE`, the column names of the location variables.
  The default name is `c('.edge_number', '.distance_on_edge')`.

- `loc_name`:

  The name of the location variables. Not needed for `rSPDE` models.

- `normalized`:

  If TRUE, then the distances in `distance_on_edge` are assumed to be
  normalized to (0,1). Default TRUE.

#### Returns

A list or a matrix containing the mesh locations.

------------------------------------------------------------------------

### `metric_graph$clear_observations()`

Clear all observations from the `metric_graph` object.

#### Usage

    metric_graph$clear_observations()

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$process_data()`

Process data to the metric graph data format.

#### Usage

    metric_graph$process_data(
      data = NULL,
      edge_number = "edge_number",
      distance_on_edge = "distance_on_edge",
      coord_x = "coord_x",
      coord_y = "coord_y",
      data_coords = c("PtE", "spatial"),
      group = NULL,
      group_sep = ".",
      normalized = FALSE,
      format = c("tibble", "sf", "sp", "list"),
      duplicated_strategy = "closest",
      include_distance_to_graph = TRUE,
      only_return_removed = FALSE,
      tolerance = max(self$edge_lengths)/2,
      verbose = FALSE,
      suppress_warnings = FALSE,
      Spoints = lifecycle::deprecated(),
      tibble = lifecycle::deprecated()
    )

#### Arguments

- `data`:

  A `data.frame` or named list containing the observations. In case of
  groups, the data.frames for the groups should be stacked vertically,
  with a column indicating the index of the group. If `data` is not
  `NULL`, it takes priority over any eventual data in `Spoints`.

- `edge_number`:

  Column (or entry on the list) of the `data` that contains the edge
  numbers. If not supplied, the column with name "edge_number" will be
  chosen. Will not be used if `Spoints` is not `NULL`.

- `distance_on_edge`:

  Column (or entry on the list) of the `data` that contains the edge
  numbers. If not supplied, the column with name "distance_on_edge" will
  be chosen. Will not be used if `Spoints` is not `NULL`.

- `coord_x`:

  Column (or entry on the list) of the `data` that contains the x
  coordinate. If not supplied, the column with name "coord_x" will be
  chosen. Will not be used if `Spoints` is not `NULL` or if
  `data_coords` is `PtE`.

- `coord_y`:

  Column (or entry on the list) of the `data` that contains the y
  coordinate. If not supplied, the column with name "coord_x" will be
  chosen. Will not be used if `Spoints` is not `NULL` or if
  `data_coords` is `PtE`.

- `data_coords`:

  It decides which coordinate system to use. If `PtE`, the user must
  provide `edge_number` and `distance_on_edge`, otherwise if `spatial`,
  the user must provide `coord_x` and `coord_y`. The option `euclidean`
  is **\[deprecated\]**. Use `spatial` instead.

- `group`:

  Vector. If the data is grouped (for example measured at different time
  points), this argument specifies the columns (or entries on the list)
  in which the group variables are stored. It will be stored as a single
  column `.group` with the combined entries.

- `group_sep`:

  separator character for creating the new group variable when grouping
  two or more variables.

- `normalized`:

  if TRUE, then the distances in `distance_on_edge` are assumed to be
  normalized to (0,1). Default FALSE.

- `format`:

  Which format should the data be returned? The options are `tibble` for
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
  `sf` for `POINT`, `sp` for `SpatialPointsDataFrame` and `list` for the
  internal list format.

- `duplicated_strategy`:

  Which strategy to handle observations on the same location on the
  metric graph (that is, if there are two or more observations projected
  at the same location). The options are 'closest' and 'jitter'. If
  'closest', only the closest observation will be used. If 'jitter', a
  small perturbation will be performed on the projected observation
  location. The default is 'closest'.

- `include_distance_to_graph`:

  When `data_coord` is 'spatial', should the distance of the
  observations to the graph be included as a column?

- `only_return_removed`:

  Should the removed data (if it exists) when using 'closest'
  `duplicated_strategy` be returned instead of the processed data?

- `tolerance`:

  Parameter to control a warning when adding observations. If the
  distance of some location and the closest point on the graph is
  greater than the tolerance, the function will display a warning. This
  helps detecting mistakes on the input locations when adding new data.

- `verbose`:

  If `TRUE`, report steps and times.

- `suppress_warnings`:

  Suppress warnings related to duplicated observations?

- `Spoints`:

  **\[deprecated\]** Use `data` instead.

- `tibble`:

  **\[deprecated\]** Use `format` instead.

#### Returns

No return value. Called for its side effects. The observations are
stored in the `data` element of the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$add_observations()`

Add observations to the metric graph.

#### Usage

    metric_graph$add_observations(
      data = NULL,
      edge_number = "edge_number",
      distance_on_edge = "distance_on_edge",
      coord_x = "coord_x",
      coord_y = "coord_y",
      data_coords = c("PtE", "spatial"),
      group = NULL,
      group_sep = ".",
      normalized = FALSE,
      clear_obs = FALSE,
      tibble = FALSE,
      tolerance = max(self$edge_lengths)/2,
      duplicated_strategy = "closest",
      include_distance_to_graph = TRUE,
      return_removed = TRUE,
      tolerance_merge = 0,
      merge_strategy = "merge",
      verbose = 1,
      suppress_warnings = FALSE,
      Spoints = lifecycle::deprecated()
    )

#### Arguments

- `data`:

  A `data.frame` or named list containing the observations. In case of
  groups, the data.frames for the groups should be stacked vertically,
  with a column indicating the index of the group. `data` can also be an
  `sf` object, a `SpatialPointsDataFrame` object or an `SSN` object. in
  which case `data_coords` will automatically be spatial, and there is
  no need to specify the `coord_x` or `coord_y` arguments.

- `edge_number`:

  Column (or entry on the list) of the `data` that contains the edge
  numbers. If not supplied, the column with name "edge_number" will be
  chosen. Will not be used if `Spoints` is not `NULL`.

- `distance_on_edge`:

  Column (or entry on the list) of the `data` that contains the edge
  numbers. If not supplied, the column with name "distance_on_edge" will
  be chosen. Will not be used if `Spoints` is not `NULL`.

- `coord_x`:

  Column (or entry on the list) of the `data` that contains the x
  coordinate. If not supplied, the column with name "coord_x" will be
  chosen. Will not be used if `Spoints` is not `NULL` or if
  `data_coords` is `PtE`.

- `coord_y`:

  Column (or entry on the list) of the `data` that contains the y
  coordinate. If not supplied, the column with name "coord_x" will be
  chosen. Will not be used if `Spoints` is not `NULL` or if
  `data_coords` is `PtE`.

- `data_coords`:

  It decides which coordinate system to use. If `PtE`, the user must
  provide `edge_number` and `distance_on_edge`, otherwise if `spatial`,
  the user must provide `coord_x` and `coord_y`. The option `euclidean`
  is **\[deprecated\]**. Use `spatial` instead.

- `group`:

  Vector. If the data is grouped (for example measured at different time
  points), this argument specifies the columns (or entries on the list)
  in which the group variables are stored. It will be stored as a single
  column `.group` with the combined entries.

- `group_sep`:

  separator character for creating the new group variable when grouping
  two or more variables.

- `normalized`:

  if TRUE, then the distances in `distance_on_edge` are assumed to be
  normalized to (0,1). Default FALSE.

- `clear_obs`:

  Should the existing observations be removed before adding the data?

- `tibble`:

  Should the data be returned as a
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)?

- `tolerance`:

  Parameter to control a warning when adding observations. If the
  distance of some location and the closest point on the graph is
  greater than the tolerance, the function will display a warning. This
  helps detecting mistakes on the input locations when adding new data.

- `duplicated_strategy`:

  Which strategy to handle observations on the same location on the
  metric graph (that is, if there are two or more observations projected
  at the same location). The options are 'closest' and 'jitter'. If
  'closest', only the closest observation will be used. If 'jitter', a
  small perturbation will be performed on the projected observation
  location. The default is 'closest'.

- `include_distance_to_graph`:

  When `data_coord` is 'spatial', should the distance of the
  observations to the graph be included as a column?

- `return_removed`:

  Should the removed data (if it exists) due to being projected to the
  same place when using 'closest' `duplicated_strategy`, or due to some
  merge strategy, be returned?

- `tolerance_merge`:

  tolerance (in edge_length units) for merging points that are very
  close and are on a common edge. By default, this tolerance is zero,
  meaning no merges will be performed.

- `merge_strategy`:

  The strategies to handle observations that are within the tolerance.
  The options are `remove`, `merge`, `average`. The default is `merge`,
  in which one of the observations will be chosen, and the remaining
  will be used to try to fill all columns with non-NA values. The second
  strategy is `remove`, meaning that if two observations are within the
  tolerance one of them will be removed. Finally, `average` will take
  the average over the close observations for numerical variables, and
  will choose one non-NA for non-numerical variables.

- `verbose`:

  Print progress of the steps when adding observations. There are 3
  levels of verbose, level 0, 1 and 2. In level 0, no messages are
  printed. In level 1, only messages regarding important steps are
  printed. Finally, in level 2, messages detailing all the steps are
  printed. The default is 1.

- `suppress_warnings`:

  Suppress warnings related to duplicated observations?

- `Spoints`:

  **\[deprecated\]** Use `data` instead.

#### Returns

No return value. Called for its side effects. The observations are
stored in the `data` element of the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$mutate_weights()`

Use [`dplyr::mutate`](https://dplyr.tidyverse.org/reference/mutate.html)
function on the internal edge weights object.

#### Usage

    metric_graph$mutate_weights(
      ...,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
on the internal edge weights object and return the result in the
requested format.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
`sf` or `sp` object containing the resulting data list after the mutate.

------------------------------------------------------------------------

### `metric_graph$select_weights()`

Use [`dplyr::select`](https://dplyr.tidyverse.org/reference/select.html)
function on the internal edge weights object.

#### Usage

    metric_graph$select_weights(
      ...,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html).

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
on the internal edge weights object and return the result in the
requested format.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
`sf` or `sp` object containing the resulting data list after the select.

------------------------------------------------------------------------

### `metric_graph$filter_weights()`

Use [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)
function on the internal edge weights object.

#### Usage

    metric_graph$filter_weights(
      ...,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
on the internal edge weights object and return the result in the
requested format.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
`sf` or `sp` object containing the resulting data list after the filter.

------------------------------------------------------------------------

### `metric_graph$summarise_weights()`

Use
[`dplyr::summarise`](https://dplyr.tidyverse.org/reference/summarise.html)
function on the internal edge weights object grouped by the edge
numbers.

#### Usage

    metric_graph$summarise_weights(
      ...,
      .groups = NULL,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).

- `.groups`:

  A vector of strings containing the names of the columns to be grouped,
  when computing the summaries. The default is `NULL`.

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
on the internal edge weights object and return the result in the
requested format.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
`sf` or `sp` object containing the resulting data list after the
summarise.

------------------------------------------------------------------------

### `metric_graph$drop_na_weights()`

Use
[`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html)
function on the internal edge weights object.

#### Usage

    metric_graph$drop_na_weights(..., format = "tibble")

#### Arguments

- `...`:

  Arguments to be passed to
  [`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html).

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html)
within the internal edge weights object.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
`sf`, or `sp` object containing the resulting data list after the
drop_na.

------------------------------------------------------------------------

### `metric_graph$mutate()`

Use [`dplyr::mutate`](https://dplyr.tidyverse.org/reference/mutate.html)
function on the internal metric graph data object.

#### Usage

    metric_graph$mutate(
      ...,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
within the internal metric graph data object and return the result in
the requested format.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
`sf`, or `sp` object containing the resulting data list after the
mutate.

------------------------------------------------------------------------

### `metric_graph$drop_na()`

Use
[`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html)
function on the internal metric graph data object.

#### Usage

    metric_graph$drop_na(..., format = "tibble")

#### Arguments

- `...`:

  Arguments to be passed to
  [`tidyr::drop_na()`](https://tidyr.tidyverse.org/reference/drop_na.html).

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use `dplyr::drop_na()` within the internal metric graph
data object.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
object containing the resulting data list after the drop_na.

------------------------------------------------------------------------

### `metric_graph$select()`

Use [`dplyr::select`](https://dplyr.tidyverse.org/reference/select.html)
function on the internal metric graph data object.

#### Usage

    metric_graph$select(
      ...,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html).

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
within the internal metric graph data object. Observe that it is a bit
different from directly using
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
since it does not allow to remove the internal positions that are needed
for the metric_graph methods to work.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
object containing the resulting data list after the selection.

------------------------------------------------------------------------

### `metric_graph$filter()`

Use [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)
function on the internal metric graph data object.

#### Usage

    metric_graph$filter(
      ...,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
within the internal metric graph data object.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
object containing the resulting data list after the filter.

------------------------------------------------------------------------

### `metric_graph$summarise()`

Use
[`dplyr::summarise`](https://dplyr.tidyverse.org/reference/summarise.html)
function on the internal metric graph data object grouped by the spatial
locations and the internal group variable.

#### Usage

    metric_graph$summarise(
      ...,
      .include_graph_groups = FALSE,
      .groups = NULL,
      .drop_na = FALSE,
      .drop_all_na = TRUE,
      format = "tibble"
    )

#### Arguments

- `...`:

  Arguments to be passed to
  [`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).

- `.include_graph_groups`:

  Should the internal graph groups be included in the grouping
  variables? The default is `FALSE`. This means that, when summarising,
  the data will be grouped by the internal group variable together with
  the spatial locations.

- `.groups`:

  A vector of strings containing the names of the columns to be
  additionally grouped, when computing the summaries. The default is
  `NULL`.

- `.drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `.drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `format`:

  The format of the output: "tibble", "sf", or "sp". Default is
  "tibble".

#### Details

A wrapper to use
[`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
within the internal metric graph data object grouped by manually
inserted groups (optional), the internal group variable (optional) and
the spatial locations. Observe that if the integral group variable was
not used as a grouping variable for the summarise, a new column, called
`.group`, will be added, with the same value 1 for all rows.

#### Returns

A [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
object containing the resulting data list after the summarise.

------------------------------------------------------------------------

### `metric_graph$get_data()`

Return the internal data with the option to filter by groups.

#### Usage

    metric_graph$get_data(
      group = NULL,
      format = c("tibble", "sf", "sp", "list"),
      drop_na = FALSE,
      drop_all_na = TRUE,
      tibble = deprecated()
    )

#### Arguments

- `group`:

  A vector contaning which groups should be returned? The default is
  `NULL`, which gives the result for the all groups.

- `format`:

  Which format should the data be returned? The options are `tibble` for
  [`tidyr::tibble`](https://tibble.tidyverse.org/reference/tibble.html),
  `sf` for `POINT`, `sp` for `SpatialPointsDataFrame` and `list` for the
  internal list format.

- `drop_na`:

  Should the rows with at least one NA for one of the columns be
  removed? DEFAULT is `FALSE`.

- `drop_all_na`:

  Should the rows with all variables being NA be removed? DEFAULT is
  `TRUE`.

- `tibble`:

  **\[deprecated\]** Use `format` instead.

------------------------------------------------------------------------

### `metric_graph$setDirectionalWeightFunction()`

Define the columns to be used for creating the directional vertex
weights. Also possible to supply user defined functions for input and
output to create ones own weights.

#### Usage

    metric_graph$setDirectionalWeightFunction(f_in = NULL, f_out = NULL)

#### Arguments

- `f_in`:

  functions for the input vertex (default `w/sum(w)`) uses the columns
  of name_column

- `f_out`:

  functions for the output vertex (deafult `rep(-1,length(w))`) uses the
  columns of name_column

#### Details

For more details see paper (that does not exists yet).

#### Returns

No return value.

------------------------------------------------------------------------

### `metric_graph$buildDirectionalConstraints()`

Build directional ODE constraint matrix from edges.

#### Usage

    metric_graph$buildDirectionalConstraints(alpha = 1)

#### Arguments

- `alpha`:

  how many derivatives the processes has

- `weight`:

  weighting for each vertex used in the constraint (E x 2)

#### Details

Currently not implemented for circles (edges that start and end in the
same vertex)

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$buildC()`

Build Kirchoff constraint matrix from edges.

#### Usage

    metric_graph$buildC(alpha = 2, edge_constraint = FALSE)

#### Arguments

- `alpha`:

  the type of constraint (currently only supports 2)

- `edge_constraint`:

  if TRUE, add constraints on vertices of degree 1

#### Details

Currently not implemented for circles (edges that start and end in the
same vertex)

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$build_mesh()`

Builds mesh object for graph.

#### Usage

    metric_graph$build_mesh(
      h = NULL,
      n = NULL,
      continuous = TRUE,
      continuous.outs = FALSE,
      continuous.deg2 = FALSE
    )

#### Arguments

- `h`:

  Maximum distance between mesh nodes (should be provided if n is not
  provided).

- `n`:

  Maximum number of nodes per edge (should be provided if h is not
  provided).

- `continuous`:

  If `TRUE` (default), the mesh contains only one node per vertex. If
  `FALSE`, each vertex v is split into deg(v) disconnected nodes to
  allow for the creation of discontinuities at the vertices.

- `continuous.outs`:

  If `continuous = FALSE` and `continuous.outs = TRUE`, continuity is
  assumed for the outgoing edges from each vertex.

- `continuous.deg2`:

  If `TRUE`, continuity is assumed at degree 2 vertices.

#### Details

The mesh is a list with the objects:

- `PtE` The mesh locations excluding the original vertices;

- `V` The verties of the mesh;

- `E` The edges of the mesh;

- `n_e` The number of vertices in the mesh per original edge in the
  graph;

- `h_e` The mesh width per edge in the graph;

- `ind` The indices of the vertices in the mesh;

- `VtE` All mesh locations including the original vertices.

#### Returns

No return value. Called for its side effects. The mesh is stored in the
`mesh` element of the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$get_version()`

Get the version of MetricGraph package used to build the graph

#### Usage

    metric_graph$get_version()

#### Returns

A character string with the version number

------------------------------------------------------------------------

### `metric_graph$is_disconnected()`

Does this graph have more than one connected component? The
decomposition is computed and cached at construction time, so this is an
O(1) lookup.

#### Usage

    metric_graph$is_disconnected()

#### Returns

`TRUE` if the graph has two or more connected components, `FALSE`
otherwise.

------------------------------------------------------------------------

### `metric_graph$get_components()`

Return the connected components of the graph as a list of `metric_graph`
objects. For a connected graph this is simply `list(self)`. For a
disconnected graph, the components are returned in order of decreasing
total edge length and any observations stored on `self` are routed to
the appropriate component. The result is cached internally; subsequent
calls are O(1) until observations change (which invalidates the cache).

#### Usage

    metric_graph$get_components(verbose = 0)

#### Arguments

- `verbose`:

  Verbosity level passed to the per-component constructors (default
  `0`).

#### Returns

A list of `metric_graph` objects.

------------------------------------------------------------------------

### `metric_graph$which_component()`

For each spatial point, determine which connected component of the graph
it belongs to. The component is the one whose nearest edge is closest in
Euclidean distance to the point.

#### Usage

    metric_graph$which_component(XY)

#### Arguments

- `XY`:

  An `n x 2` matrix of spatial coordinates (or a length-2 numeric vector
  for a single point).

#### Returns

An integer vector of length `n` with the component index for each point.
Indices match those of `get_components()` (i.e. components are sorted by
total edge length, descending).

------------------------------------------------------------------------

### `metric_graph$compute_fem()`

Build mass and stiffness matrices for given mesh object.

#### Usage

    metric_graph$compute_fem(petrov = FALSE)

#### Arguments

- `petrov`:

  Compute Petrov-Galerkin matrices? (default `FALSE`). These are defined
  as \\Cpet\_{ij} = \<\phi_i, \psi_j\>\\ and \\Gpet\_{ij} = \<d\phi_i,
  \psi_j\>\\, where \\\psi\_{i}\\ are piecewise constant basis functions
  on the edges of the mesh.

#### Details

The function builds: The matrix `C` which is the mass matrix with
elements \\C\_{ij} = \<\phi_i, \phi_j\>\\, the matrix `G` which is the
stiffness matrix with elements \\G\_{ij} = \<d\phi_i, d\phi_j\>\\, the
matrix `B` with elements \\B\_{ij} = \<d\phi_i, \phi_j\>\\, the matrix
`D` with elements \\D\_{ij} = \sum\_{v\in V}\phi_i(v)\phi_j(v)\\, and
the vector with weights \\\<\phi_i, 1\>\\.

#### Returns

No return value. Called for its side effects. The finite element
matrices `C`, `G` and `B` are stored in the `mesh` element in the
`metric_graph` object. If `petrov=TRUE`, the corresponding
Petrov-Galerkin matrices are stored in `Cpet` and `Gpet`.

------------------------------------------------------------------------

### `metric_graph$compute_mesh_weights()`

Compute the weights of the mesh nodes.

#### Usage

    metric_graph$compute_mesh_weights()

#### Details

Compute the weights of the mesh nodes.

#### Returns

No return value. Called for its side effects. The weights are stored in
the `mesh` element in the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$mesh_A()`

Deprecated - Computes observation matrix for mesh.

**\[deprecated\]** in favour of `metric_graph$fem_basis()`.

#### Usage

    metric_graph$mesh_A(PtE)

#### Arguments

- `PtE`:

  Locations given as (edge number in graph, normalized location on edge)

#### Details

For n locations and a mesh with m nodes, `A` is an n x m matrix with
elements \\A\_{ij} = \phi_j(s_i)\\.

#### Returns

The observation matrix.

------------------------------------------------------------------------

### `metric_graph$fem_basis()`

Computes observation matrix for mesh.

#### Usage

    metric_graph$fem_basis(PtE)

#### Arguments

- `PtE`:

  Locations given as (edge number in graph, normalized location on edge)

#### Details

For n locations and a mesh with m nodes, `A` is an n x m matrix with
elements \\A\_{ij} = \phi_j(s_i)\\.

#### Returns

The observation matrix.

------------------------------------------------------------------------

### `metric_graph$VtEfirst()`

Find one edge corresponding to each vertex.

#### Usage

    metric_graph$VtEfirst()

#### Returns

A nV x 2 matrix the first element of the `i`th row is the edge number
corresponding to the `i`th vertex and the second value is 0 if the
vertex is at the start of the edge and 1 if the vertex is at the end of
the edge.

------------------------------------------------------------------------

### `metric_graph$plot()`

Plots the metric graph.

#### Usage

    metric_graph$plot(
      data = NULL,
      newdata = NULL,
      group = 1,
      type = c("ggplot", "plotly", "mapview"),
      interactive = FALSE,
      vertex_size = 3,
      vertex_color = "black",
      edge_width = 0.3,
      edge_color = "black",
      data_size = 1,
      support_width = 0.5,
      support_color = "gray",
      mesh = FALSE,
      X = NULL,
      X_loc = NULL,
      p = NULL,
      degree = FALSE,
      direction = FALSE,
      arrow_size = ggplot2::unit(0.25, "inches"),
      edge_weight = NULL,
      edge_width_weight = NULL,
      scale_color_main = ggplot2::scale_color_viridis_c(option = "D"),
      scale_color_weights = ggplot2::scale_color_viridis_c(option = "C"),
      scale_color_degree = ggplot2::scale_color_viridis_d(option = "D"),
      scale_color_weights_discrete = ggplot2::scale_color_viridis_d(option = "C"),
      scale_color_main_discrete = ggplot2::scale_color_viridis_d(option = "C"),
      add_new_scale_weights = TRUE,
      scale_color_mapview = viridis::viridis(100, option = "D"),
      scale_color_weights_mapview = viridis::viridis(100, option = "C"),
      scale_color_weights_discrete_mapview = NULL,
      scale_color_degree_mapview = NULL,
      plotly = deprecated(),
      components = FALSE,
      ...
    )

#### Arguments

- `data`:

  Which column of the data to plot? If `NULL`, no data will be plotted.

- `newdata`:

  A dataset of class `metric_graph_data`, obtained by any `get_data()`,
  [`mutate()`](https://davidbolin.github.io/MetricGraph/reference/mutate.metric_graph_data.md),
  [`filter()`](https://davidbolin.github.io/MetricGraph/reference/filter.metric_graph_data.md),
  [`summarise()`](https://davidbolin.github.io/MetricGraph/reference/summarise.metric_graph_data.md),
  [`drop_na()`](https://davidbolin.github.io/MetricGraph/reference/drop_na.metric_graph_data.md)
  methods of metric graphs, see the vignette on data manipulation for
  more details.

- `group`:

  If there are groups, which group to plot? If `group` is a number and
  `newdata` is `NULL`, it will be the index of the group as stored
  internally and if `newdata` is provided, it will be the index of the
  group stored in `newdata`. If `group` is a character, then the group
  will be chosen by its name.

- `type`:

  The type of plot to be returned. The options are `ggplot` (the
  default), that uses `ggplot2`; `plotly` that uses `plot_ly` for 3D
  plots, which requires the `plotly` package, and `mapview` that uses
  the `mapview` function, to build interactive plots, which requires the
  `mapview` package.

- `interactive`:

  Only works for 2d plots. If `TRUE`, an interactive plot will be
  displayed. Unfortunately, `interactive` is not compatible with
  `edge_weight` if `add_new_scale_weights` is TRUE.

- `vertex_size`:

  Size of the vertices.

- `vertex_color`:

  Color of vertices.

- `edge_width`:

  Line width for edges. If `edge_width_weight` is not `NULL`, this
  determines the maximum edge width.

- `edge_color`:

  Color of edges.

- `data_size`:

  Size of markers for data.

- `support_width`:

  For 3D plot, width of support lines.

- `support_color`:

  For 3D plot, color of support lines.

- `mesh`:

  Plot the mesh locations?

- `X`:

  Additional values to plot.

- `X_loc`:

  Locations of the additional values in the format (edge, normalized
  distance on edge).

- `p`:

  Existing objects obtained from 'ggplot2' or 'plotly' to add the graph
  to

- `degree`:

  Show the degrees of the vertices?

- `direction`:

  Show the direction of the edges? For `type == "mapview"` the arrows
  are not shown, only the color of the vertices indicating whether they
  are problematic or not.

- `arrow_size`:

  The size of the arrows if direction is TRUE.

- `edge_weight`:

  Which column from edge weights to determine the colors of the edges?
  If `NULL` edge weights are not plotted. To plot the edge weights when
  the metric graph `edge_weights` is a vector instead of a `data.frame`,
  simply set to 1. `edge_weight` is only available for 2d plots. For 3d
  plots with edge weights, please use the `plot_function()` method.

- `edge_width_weight`:

  Which column from edge weights to determine the edges widths? If
  `NULL` edge width will be determined from `edge_width`. Currently it
  is not supported for `type = "mapview"`.

- `scale_color_main`:

  Color scale for the data to be plotted.

- `scale_color_weights`:

  Color scale for the edge weights. Will only be used if
  `add_new_scale_weights` is TRUE.

- `scale_color_degree`:

  Color scale for the degrees.

- `scale_color_weights_discrete`:

  Color scale for discrete edge weights. Will only be used if
  `add_new_scale_weights` is TRUE.

- `scale_color_main_discrete`:

  Color scale for the data to be plotted, for discrete data.

- `add_new_scale_weights`:

  Should a new color scale for the edge weights be created?

- `scale_color_mapview`:

  Color scale to be applied for data when `type = "mapview"`.

- `scale_color_weights_mapview`:

  Color scale to be applied for edge weights when `type = "mapview"`.

- `scale_color_weights_discrete_mapview`:

  Color scale to be applied for degrees when `type = "mapview"`. If
  `NULL` `RColorBrewer::brewer.pal(n = n_weights, "Set1")` will be used
  where `n_weights` is the number of different degrees.

- `scale_color_degree_mapview`:

  Color scale to be applied for degrees when `type = "mapview"`. If
  `NULL` `RColorBrewer::brewer.pal(n = n_degrees, "Set1")` will be used
  where `n_degrees` is the number of different degrees.

- `plotly`:

  **\[deprecated\]** Use `type` instead.

- `components`:

  Color the connected components separately. Either `FALSE` (the
  default; the graph is plotted as-is), `TRUE` (each component is drawn
  in a randomly-chosen color), or an `n x 3` numeric matrix of RGB
  values in `[0, 1]` (one row per component, indices matching
  `get_components()`).

- `...`:

  Additional arguments to pass to `ggplot()` or `plot_ly()`

#### Returns

A `plot_ly` (if `type = "plotly"`) or `ggplot` object.

------------------------------------------------------------------------

### `metric_graph$plot_connections()`

Plots the connections in the graph

#### Usage

    metric_graph$plot_connections()

#### Returns

No return value. Called for its side effects.

------------------------------------------------------------------------

### `metric_graph$is_tree()`

Checks if the graph is a tree (without considering directions)

#### Usage

    metric_graph$is_tree()

#### Returns

TRUE if the graph is a tree and FALSE otherwise.

------------------------------------------------------------------------

### `metric_graph$plot_function()`

Plots continuous function on the graph.

#### Usage

    metric_graph$plot_function(
      data = NULL,
      newdata = NULL,
      group = 1,
      type = c("ggplot", "plotly", "mapview"),
      continuous = TRUE,
      interpolate_plot = TRUE,
      edge_weight = NULL,
      vertex_size = 5,
      vertex_color = "black",
      edge_width = 1,
      edge_color = "black",
      line_width = NULL,
      line_color = "rgb(0,0,200)",
      scale_color = ggplot2::scale_color_viridis_c(option = "D"),
      scale_color_mapview = viridis::viridis(100, option = "D"),
      support_width = 0.5,
      support_color = "gray",
      mapview_caption = "Function",
      p = NULL,
      plotly = deprecated(),
      improve_plot = deprecated(),
      X = deprecated(),
      ...
    )

#### Arguments

- `data`:

  Which column of the data to plot? If `NULL`, no data will be plotted.

- `newdata`:

  A dataset of class `metric_graph_data`, obtained by any `get_data()`,
  [`mutate()`](https://davidbolin.github.io/MetricGraph/reference/mutate.metric_graph_data.md),
  [`filter()`](https://davidbolin.github.io/MetricGraph/reference/filter.metric_graph_data.md),
  [`summarise()`](https://davidbolin.github.io/MetricGraph/reference/summarise.metric_graph_data.md),
  [`drop_na()`](https://davidbolin.github.io/MetricGraph/reference/drop_na.metric_graph_data.md)
  methods of metric graphs, see the vignette on data manipulation for
  more details.

- `group`:

  If there are groups, which group to plot? If `group` is a number, it
  will be the index of the group as stored internally. If `group` is a
  character, then the group will be chosen by its name.

- `type`:

  The type of plot to be returned. The options are `ggplot` (the
  default), that uses `ggplot2`; `plotly` that uses `plot_ly` for 3D
  plots, which requires the `plotly` package, and `mapview` that uses
  the `mapview` function, to build interactive plots, which requires the
  `mapview` package.

- `continuous`:

  Should continuity be assumed when the plot uses `newdata`?

- `interpolate_plot`:

  Should the values to be plotted be interpolated?

- `edge_weight`:

  Which column from edge weights to plot? If `NULL` edge weights are not
  plotted. To plot the edge weights when the metric graph `edge_weights`
  is a vector instead of a `data.frame`, simply set to 1.

- `vertex_size`:

  Size of the vertices.

- `vertex_color`:

  Color of vertices.

- `edge_width`:

  Width for edges.

- `edge_color`:

  For 3D plot, color of edges.

- `line_width`:

  For 3D plot, line width of the function curve.

- `line_color`:

  Color of the function curve.

- `scale_color`:

  Color scale to be used for data and weights.

- `scale_color_mapview`:

  Color scale to be applied for data when `type = "mapview"`.

- `support_width`:

  For 3D plot, width of support lines.

- `support_color`:

  For 3D plot, color of support lines.

- `mapview_caption`:

  Caption for the function if `type = "mapview"`.

- `p`:

  Previous plot to which the new plot should be added.

- `plotly`:

  **\[deprecated\]** Use `type` instead.

- `improve_plot`:

  **\[deprecated\]** Use `interpolate` instead. There is no need to use
  it to improve the edges.

- `X`:

  **\[deprecated\]** Use `newdata` instead.

- `...`:

  Additional arguments for `ggplot()` or `plot_ly()`

#### Returns

Either a `ggplot` (if `plotly = FALSE`) or a `plot_ly` object.

------------------------------------------------------------------------

### `metric_graph$plot_movie()`

Plots a movie of a continuous function evolving on the graph.

#### Usage

    metric_graph$plot_movie(
      X,
      type = "plotly",
      vertex_size = 5,
      vertex_color = "black",
      edge_width = 1,
      edge_color = "black",
      line_width = NULL,
      line_color = "rgb(0,0,200)",
      ...
    )

#### Arguments

- `X`:

  A m x T matrix where the ith column represents the function at the ith
  time, evaluated at the mesh locations.

- `type`:

  Type of plot. Either `"plotly"` or `"ggplot"`.

- `vertex_size`:

  Size of the vertices.

- `vertex_color`:

  Color of vertices.

- `edge_width`:

  Width for edges.

- `edge_color`:

  For 3D plot, color of edges.

- `line_width`:

  For 3D plot, line width of the function curve.

- `line_color`:

  Color of the function curve.

- `...`:

  Additional arguments for ggplot or plot_ly.

#### Returns

Either a `ggplot` (if `plotly=FALSE`) or a `plot_ly` object.

------------------------------------------------------------------------

### `metric_graph$add_mesh_observations()`

Add observations on mesh to the object.

#### Usage

    metric_graph$add_mesh_observations(data = NULL, group = NULL)

#### Arguments

- `data`:

  A `data.frame` or named list containing the observations. In case of
  groups, the data.frames for the groups should be stacked vertically,
  with a column indicating the index of the group. If `data_frame` is
  not `NULL`, it takes priority over any eventual data in `Spoints`.

- `group`:

  If the data_frame contains groups, one must provide the column in
  which the group indices are stored.

#### Returns

No return value. Called for its side effects. The observations are
stored in the `data` element in the `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$get_initial_graph()`

Returns a copy of the initial metric graph.

#### Usage

    metric_graph$get_initial_graph()

#### Returns

A `metric_graph` object.

------------------------------------------------------------------------

### `metric_graph$update_graph()`

Update an older version metric graph to the current package version.

#### Usage

    metric_graph$update_graph(verbose = TRUE)

#### Arguments

- `verbose`:

  Print progress messages. Default is TRUE.

#### Details

This function creates a new metric_graph object from an existing one,
preserving all data, edge weights, mesh, distances, and other components
while ensuring compatibility with the current package version. The new
graph is created with minimal changes by disabling merges and other
transformations.

#### Returns

A new `metric_graph` object compatible with the current package version.

------------------------------------------------------------------------

### `metric_graph$coordinates()`

Convert between locations on the graph and Euclidean coordinates.

#### Usage

    metric_graph$coordinates(PtE = NULL, XY = NULL, normalized = TRUE)

#### Arguments

- `PtE`:

  Matrix with locations on the graph (edge number and normalized
  position on the edge).

- `XY`:

  Matrix with locations in Euclidean space

- `normalized`:

  If `TRUE`, it is assumed that the positions in `PtE` are normalized to
  (0,1), and the object returned if `XY` is specified contains
  normalized locations.

#### Returns

If `PtE` is specified, then a matrix with Euclidean coordinates of the
locations is returned. If `XY` is provided, then a matrix with the
closest locations on the graph is returned. Gets the edge weights
data.frame If the edge weights are given as vectors, should the result
be returned as a data.frame? A vector or `data.frame` containing the
edge weights. data List containing data on the metric graph.

------------------------------------------------------------------------

### `metric_graph$clone()`

The objects of this class are cloneable with this method.

#### Usage

    metric_graph$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
edge1 <- rbind(c(0, 0), c(2, 0))
edge2 <- rbind(c(2, 0), c(1, 1))
edge3 <- rbind(c(1, 1), c(0, 0))
edges <- list(edge1, edge2, edge3)
graph <- metric_graph$new(edges)
#> Starting graph creation...
#> LongLat is set to FALSE
#> Creating edges...
#> Setting edge weights...
#> Computing bounding box...
#> Setting up edges
#> Merging close vertices
#> Total construction time: 0.30 secs
#> Creating and updating vertices...
#> Storing the initial graph...
#> Computing the relative positions of the edges...
graph$plot()

```
