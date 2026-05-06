#' @title Connected components of metric graph
#' @description `r lifecycle::badge("deprecated")` `graph_components` is
#' deprecated. `metric_graph` now handles disconnected graphs directly:
#' construct it with `check_connected = FALSE`, then use
#' `mg$get_components()` to obtain a list of per-component
#' `metric_graph` objects, `mg$which_component(XY)` to route spatial
#' points, and `mg$plot(components = TRUE)` to colour the components.
#' Distance methods (`compute_geodist`, `compute_resdist`,
#' `compute_laplacian`) and SPDE machinery (`spde_precision`,
#' `sample_spde`, `graph_lme`, ...) all work on disconnected
#' `metric_graph` objects natively — precision matrices come out
#' block-diagonal by construction.
#' @details A list of `metric_graph` objects (representing the different
#' connected components in the full graph) created from vertex and edge matrices,
#' or from an sp::SpatialLines object where each line is representing and edge.
#' For more details, see the vignette:
#' \code{vignette("metric_graph", package = "MetricGraph")}
#' @return Object of \code{\link[R6]{R6Class}} for creating metric graph components.
#' @keywords internal
#' @examples
#' library(sp)
#' edge1 <- rbind(c(0, 0), c(1, 0))
#' edge2 <- rbind(c(1, 0), c(2, 0))
#' edge3 <- rbind(c(1, 1), c(2, 1))
#' edges <- list(edge1, edge2, edge3)
#'
#' suppressWarnings(graphs <- graph_components$new(edges))
#' graphs$plot()
#' @export
graph_components <-  R6::R6Class("graph_components",
   public = list(
     #' @field graphs List of the graphs representing the connected components.
     graphs = NULL,

     #' @field n The number of graphs.
     n = 0,

     #' @field sizes Number of vertices for each of the graphs.
     sizes = NULL,

     #' @field lengths Total edge lengths for each of the graphs.
     lengths = NULL,

     #' @field nE_total Total number of edges across all components.
     nE_total = 0L,

     #' @field edge_offsets Integer vector of length `n` giving the
     #' cumulative number of edges in components before each one. The
     #' global edge number of local edge `e` in component `k` is
     #' `edge_offsets[k] + e`.
     edge_offsets = NULL,

     #' Create metric graphs for connected components
     #'
     #' @param edges A list containing coordinates as `m x 2` matrices (that is,
     #' of `matrix` type) or m x 2 data frames (`data.frame` type) of sequence of
     #' points connected by straightlines. Alternatively, you can also prove an
     #' object of type `SpatialLinesDataFrame` or `SpatialLines` (from `sp` package)
     #' or `MULTILINESTRING` (from `sf` package).
     #' @param V n x 2 matrix with Euclidean coordinates of the n vertices.
     #' @param E m x 2 matrix where each row represents an edge.
     #' @param vertex_unit The unit in which the vertices are specified.
     #' The options are 'degree' (the great circle distance in km), 'km', 'm' and
     #' 'miles'. The default is `NULL`, which means no unit. However, if you set
     #' `length_unit`, you need to set `vertex_unit`.
     #' @param length_unit The unit in which the lengths will be computed.
     #' The options are 'km', 'm' and 'miles'. The default is `vertex_unit`.
     #' Observe that if `vertex_unit` is `NULL`, `length_unit` can only be `NULL`.
     #' If `vertex_unit` is 'degree', then the default value for `length_unit` is 'km'.
     #' @param longlat If TRUE, then it is assumed that the coordinates are given.
     #' in Longitude/Latitude and that distances should be computed in meters.
     #' It takes precedence over `vertex_unit` and `length_unit`, and is equivalent
     #' to `vertex_unit = 'degree'` and `length_unit = 'm'`.
     #' @param tolerance Vertices that are closer than this number are merged when
     #' constructing the graph (default = 1e-10). If `longlat = TRUE`, the
     #' tolerance is given in km.
     #' @param by_length Sort the components by total edge length? If `FALSE`,
     #' the components are sorted by the number of vertices.
     #' @param edge_weights Either a number, a numerical vector with length given
     #' by the number of edges, providing the edge weights, or a `data.frame` with
     #' the number of rows being equal to the number of edges, where
     #' @param ... Additional arguments used when specifying the graphs
     #' @param lines `r lifecycle::badge("deprecated")` Use `edges` instead.
     #' @return A `graph_components` object.
     initialize = function(edges = NULL,
                           V = NULL,
                           E = NULL,
                           by_length = TRUE,
                           edge_weights = NULL,
                           ...,
                           lines = deprecated()) {

       lifecycle::deprecate_warn(
         "1.5.2",
         "graph_components$new()",
         details = c(
           "metric_graph now handles disconnected graphs directly.",
           i = "Construct with `metric_graph$new(..., check_connected = FALSE)`.",
           i = "Use `mg$get_components()` for the per-component list, `mg$which_component(XY)` for spatial routing, and `mg$plot(components = TRUE)` for per-component colouring.",
           i = "See `vignette(\"graph_components\", package = \"MetricGraph\")` for the full migration guide."
         )
       )

       if (lifecycle::is_present(lines)) {
         if (is.null(edges)) {
           lifecycle::deprecate_warn("1.2.0",
                                     "graph_components$new(lines)", "graph_components$new(edges)",
                                     details = c("`lines` was provided but not `edges`. Setting `edges <- lines`.")
           )
           edges <- lines
         } else {
           lifecycle::deprecate_warn("1.2.0",
                                     "graph_components$new(lines)", "graph_components$new(edges)",
                                     details = c("Both `edges` and `lines` were provided. Only `edges` will be considered.")
           )
         }
         lines <- NULL
       }

       dots_args <- list(...)
       dots_list <- as.list(dots_args)

       if (is.null(dots_args$verbose)) {
         verbose <- 1
       } else {
         verbose <- dots_args$verbose
       }

       if(!is.null(dots_list[["project_data"]])){
         warning("The argument project_data is not compatible with graph_components. Setting project_data to FALSE.")
         dots_list[["project_data"]] <- FALSE
         dots_list[["edges"]] <- edges
         dots_list[["V"]] <- V
         dots_list[["E"]] <- E
         dots_list[["check_connected"]] <- FALSE
         dots_list[["edge_weights"]] <- edge_weights
         graph <- do.call(metric_graph$new, dots_list)
       } else{
         graph <- metric_graph$new(edges = edges, V = V, E = E,
                                   check_connected = FALSE,
                                   edge_weights = edge_weights,...)
       }

       # Making a combinatorial graph to extract the components

       if(verbose > 0){
         message("Extracting components...")
       }

       g <- make_graph(edges = c(t(graph$E)), directed = FALSE)

       if(!is.null(edge_weights)){
         edge_weights <- graph$get_edge_weights(data.frame=TRUE)
       }

       igraph::E(g)$weight <- graph$edge_lengths
       #  components <- igraph::clusters(g, mode="weak")
       components <- igraph::components(g, mode="weak")

       self$n <- components$no

       if(verbose > 0){
         message(sprintf("Number of components: %d", self$n))
       }

       dots_list[["longlat"]] <- graph$.__enclos_env__$private$longlat
       dots_list[["crs"]] <- graph$.__enclos_env__$private$crs
       dots_list[["proj4string"]] <- graph$.__enclos_env__$private$proj4string
       dots_list[["which_longlat"]] <- graph$.__enclos_env__$private$which_longlat
       dots_list[["check_connected"]] <- FALSE

       if(is.null(edge_weights)){
         edge_weights <- graph$.__enclos_env__$private$edge_weights
       }

       data_tmp <- graph$.__enclos_env__$private$data

       if(verbose > 0){
         message("Constructing graphs...")
       }

       if(self$n > 1) {
         self$graphs <- vector(mode = "list", length = self$n)
         for(k in 1:self$n) {
           if(verbose > 0){
             message(paste("Processing component", k))
           }
           vert_ids <- igraph::V(g)[components$membership == k]
           edge_rem <- NULL
           if(verbose == 2){
             message("Detecting the edges of the component...")
           }
           # Vectorized operation to identify edges to remove
           edge_rem <- which(!(graph$E[, 1] %in% vert_ids) & !(graph$E[, 2] %in% vert_ids))
           if(verbose == 2){
             message("Processing the edges to keep...")
           }
           edge_keep <- setdiff(1:graph$nE, edge_rem)
           ind_keep <- rep(0,graph$nE)
           ind_keep[edge_keep] <- 1
           if(is.null(edge_weights)){
             ew_tmp <- NULL
           } else{
             if(verbose == 2){
               message("Processing the edge weights...")
             }
             if(is.vector(edge_weights)){
               ew_tmp <- edge_weights[which(ind_keep!=0)]
             } else{
               ew_tmp <- edge_weights[which(ind_keep!=0), , drop= FALSE]
             }
           }
           if(!is.null(data_tmp)){
             if(verbose == 2){
               message("Processing the data...")
             }
             add_obs_opts <- dots_list[["add_obs_options"]]
             if(is.null(add_obs_opts)){
               add_obs_opts <- list()
             }
             idx_obs_add <- (data_tmp[[".edge_number"]]%in%edge_keep)
             data_tmp_graph <- lapply(data_tmp, function(dat){dat[idx_obs_add]})
             data_tmp_graph[[".edge_number"]] <- match(data_tmp_graph[[".edge_number"]],
                                                       edge_keep)
             class(data_tmp_graph) <- "metric_graph_data"
             add_obs_opts[["data"]] <- data_tmp_graph
           }
           if(length(graph$edges[which(ind_keep!=0)]) > 0){
             if(verbose > 0){
               message("Starting graph construction...")
             }
             dots_list[["edges"]] <- graph$edges[which(ind_keep!=0)]
             dots_list[["edge_weights"]] <- ew_tmp
             self$graphs[[k]] = do.call(metric_graph$new, dots_list)
             if(!is.null(data_tmp)){
               do.call(self$graphs[[k]]$add_observations, add_obs_opts)
             }
           }
         }
         for(i in self$n:1){
           if(is.null(self$graphs[[i]])){
             self$graphs[[i]] <- NULL
             self$n <- self$n - 1
           }
         }
         self$sizes <- components$csize
         self$lengths <- unlist(lapply(1:self$n,
                                       function(x) sum(self$graphs[[x]]$edge_lengths)))
         if(inherits(self$graphs[[1]]$get_edge_lengths(), "units")){
           units(self$lengths) <- units(self$graphs[[1]]$get_edge_lengths())
         }

         if(by_length) {
           reo <- order(self$lengths, decreasing = TRUE)
         } else {
           reo <- sort(self$sizes, decreasing = TRUE)
         }
         self$graphs <- self$graphs[reo]
         self$lengths <- self$lengths[reo]
         self$sizes <- self$sizes[reo]
       } else {
         self$graphs <- list(graph)
         self$lengths <- sum(graph$edge_lengths)
         self$sizes <- graph$nV
       }

       nE_per <- vapply(self$graphs, function(g) g$nE, integer(1))
       self$nE_total <- as.integer(sum(nE_per))
       if (self$n > 0L) {
         self$edge_offsets <- as.integer(c(0L, cumsum(nE_per)[-self$n]))
       } else {
         self$edge_offsets <- integer(0)
       }
     },

     #' @description Map global edge numbers to component indices.
     #' Each edge in the disconnected graph has a unique global edge
     #' number; this method returns the component each edge belongs to.
     #' @param edge_number Integer vector of global edge numbers in
     #' `1:self$nE_total`.
     #' @return Integer vector of component indices, the same length as
     #' `edge_number`.
     edge_to_component = function(edge_number) {
       edge_number <- as.integer(edge_number)
       if (any(is.na(edge_number)) ||
           any(edge_number < 1L) ||
           any(edge_number > self$nE_total)) {
         stop(sprintf(
           "edge_number must be integers between 1 and %d.",
           self$nE_total))
       }
       findInterval(edge_number, self$edge_offsets + 1L)
     },

     #' @description Returns the largest component in the graph.
     #' @return A `metric_graph` object.
     get_largest = function() {
       return(self$graphs[[1]])
     },

     #' @description Combine all components into a single (disconnected)
     #' `metric_graph` object via fast in-memory stacking.
     #' Edges keep their relative order within each component, and the
     #' `.edge_number` of every observation is shifted by the cumulative
     #' edge count of preceding components. The resulting graph carries
     #' a `disconnected = TRUE` flag (queryable via
     #' `metric_graph$is_disconnected()`); calling
     #' `add_observations()` on the returned graph is disallowed because
     #' the user-facing edge numbering is no longer the per-component
     #' one — add observations to the original `graph_components` and
     #' call `as_metric_graph()` again.
     #'
     #' Call this method explicitly before passing to repeated
     #' downstream operations (e.g., `graph_lme`, `graph_spde`,
     #' `sample_spde`) to amortise the assembly cost across calls.
     #' @return A `metric_graph` object representing the disjoint union
     #' of the components, with mesh and FEM matrices block-diagonally
     #' stacked when every component has them.
     as_metric_graph = function() {
       n_comp <- self$n
       graphs <- self$graphs
       ref <- graphs[[1]]
       priv_ref <- ref$.__enclos_env__$private

       # Vertex / edge / mesh-vertex offsets per component
       nV_per     <- vapply(graphs, function(g) g$nV, integer(1))
       nE_per     <- vapply(graphs, function(g) g$nE, integer(1))
       vert_off   <- c(0L, cumsum(nV_per)[-n_comp])
       edge_off   <- c(0L, cumsum(nE_per)[-n_comp])

       # Public structural fields
       V_combined <- do.call(rbind, lapply(graphs, function(g) g$V))
       E_combined <- do.call(rbind, lapply(seq_len(n_comp), function(k) {
         graphs[[k]]$E + vert_off[k]
       }))
       edges_combined <- unlist(lapply(graphs, function(g) g$edges),
                                recursive = FALSE)
       edge_lengths_combined <- unlist(
         lapply(graphs, function(g) as.numeric(g$edge_lengths)),
         use.names = FALSE
       )
       # Preserve the units attribute of edge_lengths if present
       first_el <- graphs[[1]]$edge_lengths
       if (inherits(first_el, "units")) {
         units(edge_lengths_combined) <- units(first_el)
       }

       # `self$vertices` is a list of `metric_graph_vertex` objects
       # (each is a numeric length-2 coordinate vector with degree
       # attributes). Concatenate the per-component lists and update
       # the per-vertex `id` attribute to use combined indices.
       vertices_combined <- NULL
       if (!is.null(graphs[[1]]$vertices)) {
         pieces <- lapply(seq_len(n_comp), function(k) {
           v <- graphs[[k]]$vertices
           if (is.null(v)) return(NULL)
           off <- vert_off[k]
           lapply(seq_along(v), function(i) {
             vi <- v[[i]]
             id <- attr(vi, "id")
             attr(vi, "id") <- (if (is.null(id)) i else id) + off
             vi
           })
         })
         vertices_combined <- do.call(c, pieces)
         class(vertices_combined) <- "metric_graph_vertices"
       }

       # Combined mesh
       mesh_combined <- NULL
       all_have_mesh <- all(vapply(graphs,
                                   function(g) !is.null(g$mesh),
                                   logical(1)))
       if (all_have_mesh) {
         # Per-component mesh-V offsets (so mesh-edge indices line up)
         nMeshV_per <- vapply(graphs,
                              function(g) nrow(g$mesh$V), integer(1))
         mv_off     <- c(0L, cumsum(nMeshV_per)[-n_comp])

         mesh_combined <- list()
         mesh_combined$V <- do.call(rbind,
                                    lapply(graphs, function(g) g$mesh$V))
         mesh_combined$E <- do.call(rbind, lapply(seq_len(n_comp),
                                                  function(k) {
                                                    graphs[[k]]$mesh$E + mv_off[k]
                                                  }))
         mesh_combined$h_e <- unlist(
           lapply(graphs, function(g) g$mesh$h_e), use.names = FALSE
         )
         # n_e: per-edge mesh-node counts, just concatenated (lengths
         # match the combined edge ordering).
         mesh_combined$n_e <- unlist(
           lapply(graphs, function(g) g$mesh$n_e), use.names = FALSE
         )

         # PtE: (edge, position) for mesh-only nodes; shift edge index
         ptE_pieces <- lapply(seq_len(n_comp), function(k) {
           p <- graphs[[k]]$mesh$PtE
           if (is.null(p) || nrow(p) == 0L) return(NULL)
           cbind(p[, 1] + edge_off[k], p[, 2])
         })
         mesh_combined$PtE <- do.call(rbind, ptE_pieces)
         vte_pieces <- lapply(seq_len(n_comp), function(k) {
           v <- graphs[[k]]$mesh$VtE
           if (is.null(v) || nrow(v) == 0L) return(NULL)
           cbind(v[, 1] + edge_off[k], v[, 2])
         })
         mesh_combined$VtE <- do.call(rbind, vte_pieces)
         # `ind`: in continuous mode, mesh$ind <- seq_len(self$nV).
         # When stacked, graph vertices are at rows
         # mv_off[k] + 1:nV_per[k] for component k, so we accumulate.
         ind_pieces <- lapply(seq_len(n_comp), function(k) {
           indk <- graphs[[k]]$mesh$ind
           if (is.null(indk)) return(NULL)
           if (length(indk) == 1L && indk == 0L) return(0L)
           indk + mv_off[k]
         })
         mesh_combined$ind <- if (any(vapply(ind_pieces,
                                             function(x) is.null(x) || identical(x, 0L),
                                             logical(1)))) {
           0L
         } else {
           unlist(ind_pieces, use.names = FALSE)
         }
         attr(mesh_combined, "continuous") <-
           attr(graphs[[1]]$mesh, "continuous")

         # FEM matrices: block-diagonal of the per-component matrices
         fem_keys <- c("C", "G", "B", "Cpet", "Gpet", "weights")
         for (key in fem_keys) {
           per <- lapply(graphs, function(g) g$mesh[[key]])
           if (any(vapply(per, is.null, logical(1)))) next
           if (key == "weights") {
             mesh_combined[[key]] <- unlist(per, use.names = FALSE)
           } else {
             mesh_combined[[key]] <- Matrix::bdiag(per)
           }
         }
       }

       # Combined data (with edge-number offsets)
       data_combined <- NULL
       has_data <- any(vapply(graphs,
                              function(g) !is.null(g$.__enclos_env__$private$data),
                              logical(1)))
       if (has_data) {
         pieces <- list()
         # Maintain a global offset for `.loc_idx` so that observations
         # from different components don't collide on the same unique
         # location index — graph_lme uses `.loc_idx` to deduplicate
         # observations and would otherwise treat the i-th observation
         # in two different components as the same location.
         loc_idx_off <- 0L
         for (k in seq_len(n_comp)) {
           g_data <- graphs[[k]]$.__enclos_env__$private$data
           if (is.null(g_data)) next
           d_copy <- g_data
           d_copy[[".edge_number"]] <- d_copy[[".edge_number"]] +
             edge_off[k]
           if (!is.null(d_copy[[".loc_idx"]])) {
             d_copy[[".loc_idx"]] <- d_copy[[".loc_idx"]] + loc_idx_off
             loc_idx_off <- max(d_copy[[".loc_idx"]])
           }
           pieces[[length(pieces) + 1L]] <- d_copy
         }
         all_names <- unique(unlist(lapply(pieces, names)))
         data_combined <- vector("list", length(all_names))
         names(data_combined) <- all_names
         for (nm in all_names) {
           data_combined[[nm]] <- unlist(
             lapply(pieces, function(d) {
               if (nm %in% names(d)) d[[nm]]
               else rep(NA, length(d[[".edge_number"]]))
             }),
             use.names = FALSE
           )
         }
         class(data_combined) <- c("metric_graph_data", "list")
         # Preserve group_variables attribute (default ".none" if any
         # component lacks it)
         gv_first <- attr(graphs[[1]]$.__enclos_env__$private$data,
                          "group_variables")
         attr(data_combined, "group_variables") <- if (is.null(gv_first)) ".none" else gv_first
       }

       # Combined edge weights
       # Stack per-component edge_weights, preserving the per-component
       # format (vector or data.frame). Fall back to ones for
       # components without explicit weights.
       ew_per <- lapply(graphs, function(g) {
         g$.__enclos_env__$private$edge_weights
       })
       any_df <- any(vapply(ew_per,
                            function(x) is.data.frame(x),
                            logical(1)))
       if (any_df) {
         ew_pieces <- lapply(seq_len(n_comp), function(k) {
           ew <- ew_per[[k]]
           if (is.null(ew) || !is.data.frame(ew)) {
             return(data.frame(.weights = rep(1, nE_per[k])))
           }
           ew
         })
         edge_weights_combined <- do.call(rbind, ew_pieces)
       } else {
         # Vector form (the default for metric_graph$new() output)
         ew_pieces <- lapply(seq_len(n_comp), function(k) {
           ew <- ew_per[[k]]
           if (is.null(ew)) rep(1, nE_per[k]) else as.numeric(ew)
         })
         edge_weights_combined <- unlist(ew_pieces, use.names = FALSE)
       }

       # Compute light-weight fields from the assembled V/E
       # ref_edges: vertex -> (edge, end-position). Just two match()
       # calls — cheap on the combined V/E.
       ref_edges_combined <- {
         nV_total <- sum(nV_per)
         idx_pos_0 <- match(seq_len(nV_total), E_combined[, 1], nomatch = 0)
         idx_pos_1 <- match(seq_len(nV_total), E_combined[, 2], nomatch = 0)
         cbind(ifelse(idx_pos_0 > 0, idx_pos_0, idx_pos_1),
               ifelse(idx_pos_0 > 0, 0, 1))
       }
       # bounding_box from combined V (NA-safe)
       bbox_combined <- list(
         min_x = min(V_combined[, 1], na.rm = TRUE),
         max_x = max(V_combined[, 1], na.rm = TRUE),
         min_y = min(V_combined[, 2], na.rm = TRUE),
         max_y = max(V_combined[, 2], na.rm = TRUE)
       )

       # Stack precomputed distance / Laplacian matrices
       # Distance matrices (geo_dist, res_dist) are between vertices in
       # different components are infinite, so the combined matrix
       # has Inf on the cross-component blocks. Laplacian is
       # block-diagonal (zeros off-diagonal).
       stack_distance_blocks <- function(matrices) {
         sizes <- vapply(matrices, function(m) nrow(as.matrix(m)),
                         integer(1))
         n <- sum(sizes)
         out <- matrix(Inf, n, n)
         off <- 0L
         for (m in matrices) {
           ni <- nrow(as.matrix(m))
           if (ni == 0L) next
           out[off + seq_len(ni), off + seq_len(ni)] <- as.matrix(m)
           off <- off + ni
         }
         diag(out) <- 0
         out
       }
       stack_block_diag_keys <- function(per_comp_lists, kind) {
         if (all(vapply(per_comp_lists, is.null, logical(1)))) {
           return(NULL)
         }
         # Use only components that have a list set; absent ones
         # contribute nothing — we still require *all* components to
         # have a given key to combine it (otherwise we can't form a
         # well-shaped matrix).
         all_keys <- unique(unlist(lapply(per_comp_lists, names),
                                   use.names = FALSE))
         if (length(all_keys) == 0L) return(NULL)
         out <- list()
         for (key in all_keys) {
           per_key <- lapply(per_comp_lists,
                             function(L) if (is.null(L)) NULL else L[[key]])
           if (any(vapply(per_key, is.null, logical(1)))) next
           if (kind == "distance") {
             out[[key]] <- stack_distance_blocks(per_key)
           } else { # laplacian
             out[[key]] <- Matrix::bdiag(per_key)
           }
         }
         if (length(out) == 0L) NULL else out
       }
       geo_dist_combined <- stack_block_diag_keys(
         lapply(graphs, function(g) g$geo_dist), "distance"
       )
       res_dist_combined <- stack_block_diag_keys(
         lapply(graphs, function(g) g$res_dist), "distance"
       )
       laplacian_combined <- stack_block_diag_keys(
         lapply(graphs, function(g) g$Laplacian), "laplacian"
       )

       # Build assembly list and instantiate fast
       public_fields <- list(
         V             = V_combined,
         nV            = sum(nV_per),
         E             = E_combined,
         nE            = sum(nE_per),
         edges         = edges_combined,
         edge_lengths  = edge_lengths_combined,
         vertices      = vertices_combined,
         mesh          = mesh_combined,
         geo_dist      = geo_dist_combined,
         res_dist      = res_dist_combined,
         Laplacian     = laplacian_combined
       )
       private_fields <- list(
         data          = data_combined,
         edge_weights  = edge_weights_combined,
         longlat       = if (!is.null(priv_ref$longlat)) priv_ref$longlat else FALSE,
         crs           = priv_ref$crs,
         proj4string   = priv_ref$proj4string,
         vertex_unit   = priv_ref$vertex_unit,
         length_unit   = priv_ref$length_unit,
         which_longlat = priv_ref$which_longlat,
         transform     = if (!is.null(priv_ref$transform)) priv_ref$transform else FALSE,
         connected     = FALSE,
         perform_merges = FALSE,
         ref_edges     = ref_edges_combined,
         bounding_box  = bbox_combined,
         kirchhoff_weights   = priv_ref$kirchhoff_weights,
         directional_weights = priv_ref$directional_weights,
         tolerance     = priv_ref$tolerance
       )

       metric_graph$new(.assemble = list(public  = public_fields,
                                         private = private_fields))
     },

     #' @description For each spatial point, determine which component it
     #' belongs to. The component is the one whose nearest network location
     #' is closest in Euclidean distance to the point.
     #' @param XY An `n x 2` matrix of spatial coordinates.
     #' @return An integer vector of length `n` with the component index
     #' for each point.
     which_component = function(XY) {
       if (is.vector(XY)) {
         if (length(XY) != 2) {
           stop("XY is a vector but does not have length 2")
         }
         XY <- matrix(XY, 1, 2)
       }
       if (ncol(XY) != 2) {
         stop("XY must have two columns!")
       }
       if (self$n == 1L) {
         return(rep(1L, nrow(XY)))
       }
       dists <- matrix(NA_real_, nrow = nrow(XY), ncol = self$n)
       for (k in seq_len(self$n)) {
         res <- snapPointsToLines(
           points = XY,
           lines = self$graphs[[k]]$edges,
           longlat = FALSE,
           crs = NULL
         )
         dists[, k] <- as.numeric(res$df$snap_dist)
       }
       apply(dists, 1, which.min)
     },

     #' @description Add observations to the components. Mirrors
     #' `metric_graph$add_observations`. For `data_coords = "spatial"`,
     #' each observation is routed to the component whose nearest network
     #' location is closest to it. For `data_coords = "PtE"`, the data
     #' specifies a global `edge_number` (in `1:self$nE_total`) and the
     #' component is inferred automatically from it.
     #' @param data A `data.frame`, list, `sf` object, or
     #' `SpatialPointsDataFrame` with the observations.
     #' @param edge_number Name of the (global) edge-number column.
     #' Default is `"edge_number"`.
     #' @param distance_on_edge Name of the distance-on-edge column.
     #' Default is `"distance_on_edge"`.
     #' @param coord_x Name of the x-coordinate column. Default is
     #' `"coord_x"`.
     #' @param coord_y Name of the y-coordinate column. Default is
     #' `"coord_y"`.
     #' @param data_coords Either `"PtE"` (the convention is then
     #' `(edge_number, distance_on_edge)` with global edge numbering) or
     #' `"spatial"`.
     #' @param group Optional grouping variable, see
     #' `metric_graph$add_observations`.
     #' @param normalized If TRUE, distances are assumed normalized to
     #' (0,1).
     #' @param clear_obs If TRUE, all existing observations are cleared
     #' first.
     #' @param verbose Verbosity level.
     #' @param suppress_warnings If TRUE, warnings from the per-component
     #' add_observations are suppressed.
     #' @param ... Additional arguments forwarded to each component's
     #' `add_observations` method.
     #' @return No return value. Called for its side effects.
     add_observations = function(data = NULL,
                                 edge_number = "edge_number",
                                 distance_on_edge = "distance_on_edge",
                                 coord_x = "coord_x",
                                 coord_y = "coord_y",
                                 data_coords = c("PtE", "spatial"),
                                 group = NULL,
                                 normalized = FALSE,
                                 clear_obs = FALSE,
                                 verbose = 1,
                                 suppress_warnings = FALSE,
                                 ...) {
       data_coords <- match.arg(data_coords, c("PtE", "spatial"))

       if (is.null(data)) stop("No data provided!")

       if (clear_obs) {
         self$clear_observations()
       }

       ref_graph <- self$graphs[[1]]

       if (inherits(data, "sf")) {
         data_coords <- "spatial"
         crs_ref <- ref_graph$.__enclos_env__$private$crs
         if (!is.null(crs_ref)) {
           if (!is.na(sf::st_crs(data))) {
             data <- sf::st_transform(data, crs = crs_ref)
           } else {
             data <- sf::st_set_crs(data, crs_ref)
           }
         }
         coord_tmp <- sf::st_coordinates(data)
         data <- sf::st_drop_geometry(data)
         data[[".coord_x"]] <- coord_tmp[, 1]
         data[[".coord_y"]] <- coord_tmp[, 2]
         coord_x <- ".coord_x"
         coord_y <- ".coord_y"
       }

       if ("SpatialPointsDataFrame" %in% is(data)) {
         data_coords <- "spatial"
         proj4ref <- ref_graph$.__enclos_env__$private$proj4string
         if (!is.null(proj4ref)) {
           if (!is.na(sp::proj4string(data))) {
             data <- sp::spTransform(data, sp::CRS(proj4ref))
           } else {
             sp::proj4string(data) <- proj4ref
           }
         }
         coord_tmp <- data@coords
         data <- data@data
         data[[".coord_x"]] <- coord_tmp[, 1]
         data[[".coord_y"]] <- coord_tmp[, 2]
         coord_x <- ".coord_x"
         coord_y <- ".coord_y"
       }

       if (!is.list(data) && !is.data.frame(data)) {
         stop("'data' must be either a list or a data.frame!")
       }

       data_as_list <- as.list(data)

       if (data_coords == "spatial") {
         if (is.null(data_as_list[[coord_x]]) ||
             is.null(data_as_list[[coord_y]])) {
           stop(sprintf("Data does not contain the columns '%s' and '%s'.",
                        coord_x, coord_y))
         }
         XY <- cbind(data_as_list[[coord_x]], data_as_list[[coord_y]])
         component_idx <- self$which_component(XY)
         local_edges <- NULL
       } else {
         if (is.null(data_as_list[[edge_number]]) ||
             is.null(data_as_list[[distance_on_edge]])) {
           stop(sprintf(
             "Data does not contain the columns '%s' and/or '%s'.",
             edge_number, distance_on_edge))
         }
         global_edges <- as.integer(data_as_list[[edge_number]])
         component_idx <- self$edge_to_component(global_edges)
         local_edges <- global_edges - self$edge_offsets[component_idx]
       }

       if (verbose > 0) {
         message(sprintf("Routing %d observations to %d components...",
                         length(component_idx), self$n))
       }

       for (k in seq_len(self$n)) {
         sel <- which(component_idx == k)
         if (length(sel) == 0L) next

         sub_data <- lapply(data_as_list, function(col) col[sel])

         if (data_coords == "PtE") {
           sub_data[[edge_number]] <- local_edges[sel]
           self$graphs[[k]]$add_observations(
             data = sub_data,
             edge_number = edge_number,
             distance_on_edge = distance_on_edge,
             data_coords = "PtE",
             group = group,
             normalized = normalized,
             clear_obs = FALSE,
             verbose = verbose,
             suppress_warnings = suppress_warnings,
             ...
           )
         } else {
           self$graphs[[k]]$add_observations(
             data = sub_data,
             coord_x = coord_x,
             coord_y = coord_y,
             data_coords = "spatial",
             group = group,
             clear_obs = FALSE,
             verbose = verbose,
             suppress_warnings = suppress_warnings,
             ...
           )
         }
       }

       invisible(NULL)
     },

     #' @description Clear observations from all components.
     #' @return No return value. Called for its side effects.
     clear_observations = function() {
       for (k in seq_len(self$n)) {
         self$graphs[[k]]$clear_observations()
       }
       invisible(NULL)
     },

     #' @description Combine observations from all components.
     #' Returns a single data structure with an additional `.component`
     #' column indicating which component each row belongs to.
     #' @param group Optional group filter passed to each component's
     #' `get_data`.
     #' @param format One of `"tibble"`, `"sf"`, `"sp"`, `"list"`.
     #' @param drop_na See `metric_graph$get_data`.
     #' @param drop_all_na See `metric_graph$get_data`.
     #' @return A combined data object with a `.component` column.
     get_data = function(group = NULL,
                         format = c("tibble", "sf", "sp", "list"),
                         drop_na = FALSE, drop_all_na = TRUE) {
       format <- match.arg(format, c("tibble", "sf", "sp", "list"))
       out <- vector("list", self$n)
       any_data <- FALSE
       for (k in seq_len(self$n)) {
         g <- self$graphs[[k]]
         if (is.null(g$.__enclos_env__$private$data)) next
         d <- g$get_data(group = group, format = "list",
                         drop_na = drop_na, drop_all_na = drop_all_na)
         d[[".edge_number"]] <- d[[".edge_number"]] + self$edge_offsets[k]
         d[[".component"]] <- rep(k, length(d[[".edge_number"]]))
         out[[k]] <- d
         any_data <- TRUE
       }
       if (!any_data) {
         stop("No data found in any component.")
       }
       out <- out[!vapply(out, is.null, logical(1))]

       all_names <- unique(unlist(lapply(out, names)))
       combined <- vector("list", length(all_names))
       names(combined) <- all_names
       for (nm in all_names) {
         combined[[nm]] <- unlist(lapply(out, function(d) {
           if (nm %in% names(d)) d[[nm]] else rep(NA, length(d[[".edge_number"]]))
         }), use.names = FALSE)
       }

       if (format == "list") {
         class(combined) <- c("metric_graph_data", "list")
         return(combined)
       }
       combined_df <- as.data.frame(combined, stringsAsFactors = FALSE)
       if (format == "tibble") {
         combined_df <- tidyr::as_tibble(combined_df)
       }
       if (format == "sf") {
         ref_crs <- self$graphs[[1]]$.__enclos_env__$private$crs
         geoms <- lapply(seq_len(nrow(combined_df)), function(i) {
           sf::st_point(as.numeric(combined_df[i,
                                               c(".coord_x", ".coord_y")]))
         })
         combined_df <- sf::st_sf(combined_df,
                                  geometry = sf::st_sfc(geoms),
                                  crs = if (!is.null(ref_crs)) ref_crs else NULL)
       }
       if (format == "sp") {
         combined_df <- as.data.frame(combined_df)
         sp::coordinates(combined_df) <- ~ .coord_x + .coord_y
         return(combined_df)
       }
       if (!inherits(combined_df, "metric_graph_data")) {
         class(combined_df) <- c("metric_graph_data", class(combined_df))
       }
       combined_df
     },

     #' @description Get the unique groups across all components.
     #' @return Character vector of unique group identifiers.
     get_groups = function() {
       groups <- character(0)
       any_data <- FALSE
       for (k in seq_len(self$n)) {
         g_data <- self$graphs[[k]]$.__enclos_env__$private$data
         if (!is.null(g_data)) {
           groups <- union(groups, unique(g_data[[".group"]]))
           any_data <- TRUE
         }
       }
       if (!any_data) {
         warning("There is no data!")
         return(invisible(NULL))
       }
       groups
     },

     #' @description Get the (edge_number, distance_on_edge) pairs for
     #' the observations across all components, using global edge
     #' numbering.
     #' @return A matrix with two columns: `edge_number`,
     #' `distance_on_edge`.
     get_PtE = function() {
       out_list <- list()
       for (k in seq_len(self$n)) {
         g_data <- self$graphs[[k]]$.__enclos_env__$private$data
         if (!is.null(g_data)) {
           pte_k <- self$graphs[[k]]$get_PtE()
           pte_k[, 1] <- pte_k[, 1] + self$edge_offsets[k]
           out_list[[length(out_list) + 1L]] <- pte_k
         }
       }
       if (length(out_list) == 0L) {
         warning("There is no data!")
         return(invisible(NULL))
       }
       do.call(rbind, out_list)
     },

     #' @description Convert between graph coordinates (global
     #' `edge_number`, `distance_on_edge`) and spatial coordinates.
     #' @param PtE A matrix or vector with two columns/entries:
     #' `(edge_number, distance_on_edge)`. Edge numbers are global
     #' across components.
     #' @param XY An `n x 2` matrix of spatial coordinates.
     #' @param normalized If TRUE, the distances are normalized to (0,1).
     #' @return If `PtE` is supplied, an `n x 2` matrix of spatial
     #' coordinates. If `XY` is supplied, an `n x 2` matrix with
     #' `(edge_number, distance_on_edge)` using global edge numbering.
     coordinates = function(PtE = NULL, XY = NULL, normalized = TRUE) {
       if (is.null(PtE) && is.null(XY)) {
         stop("PtE or XY must be provided")
       } else if (!is.null(PtE) && !is.null(XY)) {
         stop("Either PtE or XY must be provided, not both")
       }

       if (!is.null(PtE)) {
         if (is.vector(PtE)) {
           if (length(PtE) != 2) {
             stop("PtE is a vector but does not have length 2 (edge, distance)")
           }
           PtE <- matrix(PtE, 1, 2)
         }
         if (ncol(PtE) != 2) {
           stop("PtE must have two columns: edge_number, distance_on_edge.")
         }
         comp <- self$edge_to_component(PtE[, 1])
         local_edges <- as.integer(PtE[, 1]) - self$edge_offsets[comp]
         Points <- matrix(NA_real_, nrow = nrow(PtE), ncol = 2L)
         for (k in unique(comp)) {
           idx <- which(comp == k)
           Points[idx, ] <- self$graphs[[k]]$coordinates(
             PtE = cbind(local_edges[idx], PtE[idx, 2]),
             normalized = normalized
           )
         }
         return(Points)
       } else {
         if (is.vector(XY)) {
           if (length(XY) != 2) {
             stop("XY is a vector but does not have length 2")
           }
           XY <- matrix(XY, 1, 2)
         }
         if (ncol(XY) != 2) {
           stop("XY must have two columns!")
         }
         component_idx <- self$which_component(XY)
         out <- matrix(NA_real_, nrow = nrow(XY), ncol = 2L)
         for (k in unique(component_idx)) {
           idx <- which(component_idx == k)
           pte_k <- self$graphs[[k]]$coordinates(
             XY = XY[idx, , drop = FALSE], normalized = normalized
           )
           out[idx, 1] <- pte_k[, 1] + self$edge_offsets[k]
           out[idx, 2] <- pte_k[, 2]
         }
         return(out)
       }
     },

     #' @description Build a mesh on each component.
     #' @param ... Arguments forwarded to `metric_graph$build_mesh`.
     #' @return No return value. Called for its side effects.
     build_mesh = function(...) {
       for (k in seq_len(self$n)) {
         self$graphs[[k]]$build_mesh(...)
       }
       invisible(NULL)
     },

     #' @description Compute finite-element matrices on each component
     #' that has a mesh.
     #' @param ... Arguments forwarded to `metric_graph$compute_fem`.
     #' @return No return value. Called for its side effects.
     compute_fem = function(...) {
       for (k in seq_len(self$n)) {
         if (!is.null(self$graphs[[k]]$mesh)) {
           self$graphs[[k]]$compute_fem(...)
         }
       }
       invisible(NULL)
     },

     #' @description Compute geodesic distances per component. The
     #' geodesic distance between vertices in different components is
     #' infinite, so the per-component computation gives well-defined
     #' results. The distances are stored on each component's
     #' `metric_graph` object (in `geo_dist`); when `as_metric_graph()`
     #' is subsequently called, they are stacked into a combined
     #' distance matrix with `Inf` cross-component entries.
     #' @param ... Arguments forwarded to `metric_graph$compute_geodist`.
     #' @return No return value. Called for its side effects.
     compute_geodist = function(...) {
       for (k in seq_len(self$n)) {
         self$graphs[[k]]$compute_geodist(...)
       }
       invisible(NULL)
     },

     #' @description Compute resistance distances per component.
     #' Resistance distance is undefined (infinite) between vertices in
     #' different components, so per-component computation is the only
     #' way to get well-defined values on a disconnected graph.
     #' Stacking via `as_metric_graph()` produces a combined matrix
     #' with `Inf` cross-component entries, making
     #' `graph_lme(model = "isoExp")` and similar isotropic models
     #' fit per-component automatically.
     #' @param ... Arguments forwarded to `metric_graph$compute_resdist`.
     #' @return No return value. Called for its side effects.
     compute_resdist = function(...) {
       for (k in seq_len(self$n)) {
         self$graphs[[k]]$compute_resdist(...)
       }
       invisible(NULL)
     },

     #' @description Compute the (weighted) graph Laplacian per
     #' component. The combined Laplacian on the disjoint union is
     #' block-diagonal in the per-component Laplacians, so
     #' `as_metric_graph()` simply stacks them via `Matrix::bdiag()`.
     #' @param ... Arguments forwarded to `metric_graph$compute_laplacian`.
     #' @return No return value. Called for its side effects.
     compute_laplacian = function(...) {
       for (k in seq_len(self$n)) {
         self$graphs[[k]]$compute_laplacian(...)
       }
       invisible(NULL)
     },

     #' @description Plot a function on the components. Mirrors
     #' `metric_graph$plot_function`. When supplying `newdata`, it must
     #' include `.edge_number` (global edge numbering) and
     #' `.distance_on_edge`; the component for each row is inferred from
     #' the global edge number.
     #' @param data Column name of the stored observations to plot.
     #' @param newdata Optional `metric_graph_data` with `.edge_number`
     #' (global), `.distance_on_edge`, and the value column.
     #' @param group Group identifier passed to each component.
     #' @param type Plot type: `"ggplot"` or `"plotly"`.
     #' @param continuous See `metric_graph$plot_function`.
     #' @param p Optional existing plot to add to.
     #' @param ... Additional arguments forwarded to each component's
     #' `plot_function`.
     #' @return A plot object.
     plot_function = function(data = NULL, newdata = NULL, group = 1,
                              type = c("ggplot", "plotly"),
                              continuous = TRUE, p = NULL, ...) {
       type <- match.arg(type, c("ggplot", "plotly"))

       if (is.null(data) && is.null(newdata)) {
         stop("You must provide either 'data' or 'newdata'.")
       }

       newdata_comp <- NULL
       if (!is.null(newdata)) {
         if (!".edge_number" %in% names(newdata)) {
           stop("'newdata' must include a '.edge_number' column.")
         }
         newdata_comp <- self$edge_to_component(newdata[[".edge_number"]])
       }

       for (k in seq_len(self$n)) {
         g <- self$graphs[[k]]
         sub_newdata <- NULL
         if (!is.null(newdata)) {
           sel <- (newdata_comp == k)
           if (!any(sel)) next
           sub_newdata <- as.data.frame(newdata)[sel, , drop = FALSE]
           sub_newdata[[".edge_number"]] <- sub_newdata[[".edge_number"]] -
             self$edge_offsets[k]
           if (!inherits(sub_newdata, "metric_graph_data")) {
             class(sub_newdata) <- c("metric_graph_data", class(sub_newdata))
           }
         } else if (is.character(data)) {
           g_data <- g$.__enclos_env__$private$data
           if (is.null(g_data) || !(data %in% names(g_data))) next
         }
         p <- suppressMessages(
           g$plot_function(data = data, newdata = sub_newdata,
                           group = group, type = type,
                           continuous = continuous, p = p, ...)
         )
       }
       p
     },

     #' @description Plots all components.
     #' @param edge_colors A 3 x nc matrix with RGB values for the edge colors to
     #' be used when plotting each graph.
     #' @param vertex_colors A 3 x nc matrix with RGB values for the edge colors to
     #' be used when plotting each graph.
     #' @param data Optional column name of the stored observations to
     #' plot. Components without that column are skipped.
     #' @param ... Additional arguments for plotting the individual graphs.
     #' @return A `ggplot` object.
     plot = function(edge_colors = NULL, vertex_colors = NULL,
                     data = NULL, ...) {

       if (is.null(edge_colors)) {
         edge_colors <- matrix(0, nrow = self$n, ncol = 3)
         if(self$n > 1) {
           for(i in 2:self$n) {
             edge_colors[i, ] = runif(3)
           }
         }
       } else {
         if (ncol(edge_colors)!= 3) {
           stop("edge_colors must have three columns!")
         }
         if (nrow(edge_colors)!= self$n) {
           stop("edge_colors must have the same number of rows as there are components!")
         }
       }
       if (is.null(vertex_colors)) {
         vertex_colors <- edge_colors
       } else {
         if (ncol(vertex_colors)!= 3) {
           stop("vertex_colors must have three columns!")
         }
         if (nrow(vertex_colors)!= self$n) {
           stop("vertex_colors must have the same number of rows as there are components!")
         }
       }

       .has_data_for <- function(g) {
         if (is.null(data) || !is.character(data)) return(TRUE)
         g_data <- g$.__enclos_env__$private$data
         !is.null(g_data) && (data %in% names(g_data))
       }

       p <- NULL
       for (i in seq_len(self$n)) {
         g <- self$graphs[[i]]
         if (.has_data_for(g)) {
           args <- list(
             edge_color = rgb(edge_colors[i, 1], edge_colors[i, 2],
                              edge_colors[i, 3]),
             vertex_color = rgb(vertex_colors[i, 1], vertex_colors[i, 2],
                                vertex_colors[i, 3]),
             data = data,
             p = p,
             ...
           )
         } else {
           args <- list(
             edge_color = rgb(edge_colors[i, 1], edge_colors[i, 2],
                              edge_colors[i, 3]),
             vertex_color = rgb(vertex_colors[i, 1], vertex_colors[i, 2],
                                vertex_colors[i, 3]),
             p = p,
             ...
           )
         }
         p <- suppressMessages(do.call(g$plot, args))
       }
       return(p)
     }))

