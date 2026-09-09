# MetricGraph (development version)

* Graph construction no longer materializes dense point-by-edge matrices.
  `metric_graph$new()` used to densify the `st_is_within_distance()` result
  into an `nEdges x nVertices` logical matrix before snapping vertices to
  nearby edges, and `snapPointsToLines()` (used by `add_observations()` with
  `data_coords = "spatial"`, `coordinates()` and `which_component()`) built a
  full `nPoints x nEdges` distance matrix. Both are now replaced by sparse
  candidate sets and C++ kernels (`snap_points_to_edges_cpp()` and a
  grid-indexed `nearest_edge_cpp()`), which return the same graphs and the same
  snapped locations. 
* The edge-edge intersection search (`tolerance$edge_edge > 0`) also keeps its
  neighbour list sparse instead of building an `nEdges x nEdges` logical matrix.
* Edge weights are no longer copied once per edge split. `add_vertices()` now
  records which weight row each new edge inherits and materializes the weight
  table once, instead of `rbind`-ing the whole table inside every
  `split_edge()` call. `set_edge_weights()` builds the per-edge weight rows 
  without a `[.data.frame` call per edge. The values are unchanged but row names 
  of weight rows duplicated by a split are now e.g. `"990509915.1"` rather than 
  `"9905099151"`.
* `prune_vertices()` no longer pays a `[.data.frame` call per edge when
  rebuilding the edge attributes, and the serial fallback for closed-loop
  chains compacts the weight table once instead of once per removed vertex.
* Added `metric_graph$get_largest()`, which returns only the largest connected
  component. It is equivalent to `get_components()[[1]]` but constructs just
  that component.
* Added a mirrored closed-form directional OU covariance fast path for
  out-trees, including direction-reversed K1 continuity graphs. Reversed
  river networks now share the K1/K2 C++ pairwise kernel, using a cached
  Euler/RMQ LCA index and source-normalized transfers instead of allocating
  R objects or walking parent chains for every covariance pair.
* Added C++-backed directional alpha-1 edge precision and closed-form
  directional OU covariance calculations. The C++ paths are now the default,
  with pure-R reference implementations retained for validation and fallback.
* Added a packaged-data Mid-Columbia speed example covering K1, K2,
  reversed-continuity and covariance likelihood evaluation.
* Fixed a sign-convention inconsistency in the exact α = 2 codes. The 
  observation-/bridge-side covariance blocks treated the
  endpoint-derivative states as −u′, while the prior precision (`Q00` /
  `Qalpha2`) and the vertex constraints used +u′. Because both enter the same
  quadratic forms, α = 2 covariances, likelihoods, posteriors and exact samples
  were biased on graphs containing cycles. 
* Fixed the same α = 2 sign-convention bug in the profiled-likelihood (v2)
  code path (`profile_lik_core_alpha2` in `R/graph_likelihoods_v2.R`), which
  had reimplemented the affected S-matrix construction independently and was
  not covered by the fix above.
* Added `simulate.metric_graph()` S3 method and `simulate_parallel()` for
  unconditional prior simulation of Whittle-Matérn fields: Method A (`method = "direct"`, 
  O(m³) per edge), Method B (`method = "kriging"`, O(m) per edge), and the extended method
  (`method = "extended"`, single sparse Cholesky on the graph with simulation
  locations promoted to vertices). Both α = 1 and α = 2 are supported.
* C++ (Rcpp/Eigen) implementations of the per-edge bridge draws,
  `draw_edge_direct_cpp` and `draw_edge_kriging_cpp`, are now the default
  (`impl = "cpp"`); the pure-R reference implementations remain available via
  `impl = "R"`.
* Updated `examples/fast_simulation/run_study.R` with warm-up iterations,
  `n_rep = 5` medians, extended `n_pts` sweep {8,...,2048}, and all three
  simulation methods.
* Added `examples/fast_simulation/benchmark.R` for per-edge R vs C++ speedup
  and whole-field method comparison.

# MetricGraph 1.6.0

* `metric_graph` now fully supports disconnected graphs, and `graph_components` 
has therefore been deprecated.
* The internal data storage system has been improved to reduce memory requirements. 
* Several improvements to internal methods to speed up the construction and data 
handling for of very large graphs.
* Various stability improvements and bug fixes. 
* Added a `cross_validation()` function for `inlabru` fits on metric
  graphs, mirroring `rSPDE::cross_validation()`.
* Bumped minimum `rSPDE` version to 2.5.0 and added `rlang` to `Imports`.
* Added various unit tests covering input validation, smoke tests, multi-model 
comparison, all CV types, `inla` and `inlabru` models. 

# MetricGraph 1.5.1
* Added `precompute_lgcp_graph()` to precompute expensive quantities for LGCP models.
* Added `update_graph()` method to update graph objects to newer versions of the package.
* Added lower and upper bounds for hyperparameters in SPDE and LGCP exact models.
* Improved speed and robustness of the `observations_to_vertex` method.
* Converted `fem_basis` computation to C for improved performance.
* Updated SPDE result object for exact models.
* Fixed a prediction bug in `graph_lme`.
* Fixed handling of units in edge weights.
* Fixed documentation macros in `graph_lgcp_sim()` and `lgcp_graph()` man pages.

# MetricGraph 1.5.0
* Several optimization improvements for the models in `graph_lme`.
* Added a `parallel` option to `posterior_crossvalidation`.
* Created an `INLA` interface for LGCP processes on metric graphs.

# MetricGraph 1.4.1
* Adding selected inverse function, for computing the inverse matrix elements only on nonzero entries of the original sparse matrix.
* Adjusts on `add_observations()` and `plot_function()` methods.
* Adding helper functions to use `stlnpp` objects.

# MetricGraph 1.4.0
* Added an INLA implementation for `alpha=2`.
* Added a vignette for handling multiple likelihoods in `R-INLA` and `inlabru`, and updated the `graph_spde_data()` function for such cases.
* Added an INLA implementation for directional models.
* Added support for directional edge weights.
* Added support for creating metric graphs from `SSN`, `osmdata_sp` and `osmdata_sf` objects. In such cases, if available, edge weights and data will be automatically added to the graph.
* Now if one creates the metric graph from `SpatialLinesDataFrame`, `LINESTRING`, `MULTILINESTRING`, etc., if the object contain edge data, they will be automatically added as edge weights.
* Added an option to not perform merges (which is now the default), that makes the graph creation faster and using less memory.
* Added vignettes with a river example and with an example of directional models.
* Updated the vignettes to account for the additions.
* `prune_vertices` now has an option to avoid creating circles when pruning.
* Several updates and quality of life improvements for building very large graphs faster.
* Added the `export()` method, that allows one to export a MetricGraph object as an `sf`, `sp` or `SSN2` object.
* Added wrappers for `leaflet` and `mapview` as methods.
* Added `get_edges()` and `get_vertices()` methods were created, that return the edges and vertices, respectively, in either 'sf', 'sp' or the internal formats.
* Added an option `format` to the `get_data()` method that allows one to also return the data in `sf` or `sp` formats.
* The `plot` method now has an argument `type`, that now also allows one to set `type` to `mapview`, thus it can return `ggplot2`, `plotly` and `mapview` objects.
* Adding methods to do data manipulation on weights, `mutate_weight`, `select_weights`, `filter_weights`, `summarise_weights` and `drop_na_weights`. They have a `format` argument that allows one to also return `sf` or `sp` objects.
* Updated the methods, `mutate`, `filter`, `select`, `drop_na` and `summarise` to have a format argument to also return `sf` or `sp` objects.
* Massive improvement on the `observation_to_vertex` method.
* Massive improvement on the metric graph creation speed.
* Deprecated `improve_plot` option, as now all plots from `plot_function()` method are improved.
* Added `merge_strategy` option for `add_observations()` method, for handling observations that are very close.
* Update the metric graph data vignette for illustrating how to use the new tools for data manipulation. 
* Massive improvement for building constraint matrices for `alpha=2`, and for building directional constraints. 
* Updated starting values to use bounding boxes to be more efficient.

# MetricGraph 1.3.0
* Handlers were added in `add_observations()` for situations where observations are projected at the same location, specifically for the `duplicated_strategy` argument.
* A `simulate` method was added for `graph_lme` objects.
* The possibility of fixing parameters during estimation was added.
* The `Spoints` argument in `add_observations()` has been deprecated. Now, `SpatialPointsDataFrame` can be added directly in the `data` argument.
* `sf` objects containing data can also be directly added using the `add_observations()` method in the `data` argument.
* The option of using a `graph_lme` object to provide starting values when fitting a model using `graph_lme()` was added.
* The option of fitting a directional Whittle-Matérn model with `alpha=1` when using `graph_lme()` was added.
* The `kirchhoff_weights` argument was added to obtain weights for Kirchhoff vertex conditions from `edge_weights`.
* Handling of edge weights was improved. For example, if pruning changes any edge weight, a warning will be given.
* The `edgeweight_to_data()` method was added to turn edge weights into data in the internal metric graph format.
* `edge_weight` and `edge_width_weight` were added to the `plot()` method so that plots on metric graphs can be produced with the weights providing colors to the edges, and also with (possibly different) weights providing the thickness of the edges.
* `edge_weights` were added to `graph_components` so that the connected components will have the correct corresponding edge weights.
* `edge_weight` and `edge_width_weight` were added to the `plot_function()` method, where they work in a similar manner to their counterparts for the `plot()` method. The difference is that the weights are plotted as piecewise constant functions.
* The `prune_vertices` now has an option to not prune vertices whose edges have incompatible edge weights.
* The `plot` method has an `interactive` argument that returns the 2D plot as a plotly object, which is interactive when using `ggplotly`.
* The dependency on the `viridis` package has been removed.

# MetricGraph 1.2.0
* Changed argument `data` to `newdata` in `predict` methods. The argument `data` was deprecated.
* Bugfixes on sample_spde and when adding observations based on Euclidean positions.
* Added options `vertex_unit` and `length_unit` on graph creations. Units are given in edge lengths with the `get_edge_lengths()` method.
* Added a method to check if the graph is a tree.
* The graph construction was thoroughly refactored. The resulting construction is faster and cleaner.
* The graph constructions now accepts list of coordinates (where the coordinates are given as either matrices or data frames), `SpatialLines`, `SpatialLinesDataFrames` or `MULTILINESTRING`.
* Adding two options (`sf` package or `sp` package) for handling `longlat` by using the `which_longlat` option.
* Adding `crs` (if using `sf`) and `proj4string` (if using `sp`) for handling general coordinate reference systems.
* Moving `data` to the `private` environment.
* Several data manipulation helper tools and methods were introduced, together with a vignette with a brief tutorial on these tools.
* The method `mesh_A()` has been deprecated, use `fem_basis()` instead.
* Several quality of life improvements.
* Improved the `plot()` method with the option `plotly=TRUE`.
* Improved the `plot_function()` method to accept `data` and `newdata`.
* Included a `process_data()` method for metric graphs
* Renamed the data internal structure from "__group", "__edge_number", "__distance_on_edge", "__coord_x", "__coord_y" to ".group", ".edge_number", ".distance_on_edge", ".coord_x" and ".coord_y".
* Added an "advanced grouping" option, in which the group variable can be a combination of several columns.
* Improved graph_lme() behavior to avoid having NaN as std.errors.
* Added check for distance consistency and, more generally, check to see if the graph has euclidean edges.
* Added method `get_groups()` to get the unique groups, and also to retrieve the columns that were used to create the group variable.
* Added the `get_data()` method to get the data in a user-friendly manner.
* Added `glance()` and `augment()` methods for `graph_lme()` objects.
* Added `get_vertices_incomp_dir()` method to return vertices with incompatible directions.
* Added `print()`, `summary()`, `compute_characteristics()`, `check_euclidean()`, `check_distance_consistency()` methods.
* Added support for edge weights.
* Created `vertices` element in the metric graph object, containing information such as degrees, indegrees, outdegrees.
* Created `print` methods for `edges`, `vertices`, and for their entries.
* Added the `improve_plot` option on `plot_function`.
* Added support for discontinuous meshes (at the vertices).
* Added support for discontinuous functions (at the vertices) for `plot_function()`.

# MetricGraph 1.1.2
* Adjusts to ensure compatibility with future releases of the Matrix package.

# MetricGraph 1.1.1
* Adjusts on documentation for CRAN.

# MetricGraph 1.1.0
* Improved the documentation.
* Reorganized some functions.

# MetricGraph 1.0.0
* First version of the package.
