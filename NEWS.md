# MetricGraph (development version)
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
