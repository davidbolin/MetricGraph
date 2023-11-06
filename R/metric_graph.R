#' @title Metric graph
#' @description Class representing a general metric graph.
#' @details A graph object created from vertex and edge matrices, or from an
#' `sp::SpatialLines` object where each line is representing and edge. For more details,
#'  see the vignette:
#' \code{vignette("metric_graph", package = "MetricGraph")}
#' @return Object of \code{\link[R6]{R6Class}} for creating metric graphs.
#' @examples
#' edge1 <- rbind(c(0, 0), c(2, 0))
#' edge2 <- rbind(c(2, 0), c(1, 1))
#' edge3 <- rbind(c(1, 1), c(0, 0))
#' edges <- list(edge1, edge2, edge3)
#' graph <- metric_graph$new(edges)
#' graph$plot()
#'
#' @export
metric_graph <-  R6Class("metric_graph",
  public = list(
  #' @field V Matrix with positions in Euclidean space of the vertices of the
  #' graph.
  V = NULL,

  #' @field nV The number of vertices.
  nV = 0,

  #' @field E Matrix with the edges of the graph, where each row represents an
  #' edge, `E[i,1]` is the vertex at the start of the ith edge and `E[i,2]` is
  #' the vertex at the end of the edge.
  E = NULL,

  #' @field nE The number of edges.
  nE= 0,

  #' @field edge_lengths Vector with the lengths of the edges in the graph.
  edge_lengths = NULL,

  #' @field C Constraint matrix used to set Kirchhoff constraints.
  C = NULL,

  #' @field CoB Change-of-basis object used for Kirchhoff constraints.
  CoB = NULL,

  #' @field PtV Vector with the indices of the vertices which are observation
  #' locations.
  PtV  = NULL,

  #' @field mesh Mesh object used for plotting.
  mesh = NULL,

  #' @field edges The coordinates of the edges in the graph.
  edges = NULL,

  #' @field vertices The coordinates of the vertices in the graph, along with several attributes.
  vertices = NULL,

  #' @field geo_dist Geodesic distances between the vertices in the graph.
  geo_dist = NULL,

  #' @field res_dist Resistance distances between the observation locations.
  res_dist = NULL,

  #' @field Laplacian The weighted graph Laplacian of the vertices in the
  #' graph. The weights are given by the edge lengths.
  Laplacian = NULL,

  #' @field characteristics List with various characteristics of the graph.
  characteristics = NULL,

  #' @description Create a new `metric_graph` object.
  #' @param edges A list containing coordinates as `m x 2` matrices (that is, of `matrix` type) or m x 2 data frames (`data.frame` type) of sequence of points connected by straightlines. Alternatively, you can also prove an object of type `SpatialLinesDataFrame` or `SpatialLines` (from `sp` package) or `MULTILINESTRING` (from `sf` package).
  #' @param V n x 2 matrix with Euclidean coordinates of the n vertices.
  #' @param E m x 2 matrix where each row represents one of the m edges.
  #' @param vertex_unit The unit in which the vertices are specified. The options are 'degrees' (the great circle distance in km), 'km', 'm' and 'miles'. The default is `NULL`, which means no unit. However, if you set `length_unit`, you need to set `vertex_unit`.
  #' @param length_unit The unit in which the lengths will be computed. The options are 'km', 'm' and 'miles'. The default is `vertex_unit`. Observe that if `vertex_unit` is `NULL`, `length_unit` can only be `NULL`.
  #' If `vertex_unit` is 'degrees', then the default value for `length_unit` is 'km'.
  #' @param edge_weights Either a number, a numerical vector with length given by the number of edges, providing the edge weights, or a `data.frame` with the number of rows being equal to the number of edges, where
  #' each row gives a vector of weights to its corresponding edge. Can be changed by using the `set_edge_weights()` method.
  #' @param longlat If `TRUE`, then it is assumed that the coordinates are given.
  #' in Longitude/Latitude and that distances should be computed in meters. If `TRUE` it takes precedence over
  #' `vertex_unit` and `length_unit`, and is equivalent to `vertex_unit = 'degrees'` and `length_unit = 'm'`.
  #' @param crs Coordinate reference system to be used in case `longlat` is set to `TRUE` and `which_longlat` is `sf`. Object of class crs. The default is `sf::st_crs(4326)`.
  #' @param proj4string Projection string of class CRS-class to be used in case `longlat` is set to `TRUE` and `which_longlat` is `sp`. The default is `sp::CRS("+proj=longlat +datum=WGS84")`.
  #' @param which_longlat Compute the distance using which package? The options are `sp` and `sf`. The default is `sp`.
  #' @param project If `longlat` is `TRUE` should a projection be used to compute the distances to be used for the tolerances (see `tolerance` below)? The default is `TRUE`. When `TRUE`, the construction of the graph is faster.
  #' @param project_data If `longlat` is `TRUE` should the vertices be project to planar coordinates? The default is `FALSE`. When `TRUE`, the construction of the graph is faster.
  #' @param which_projection Which projection should be used in case `project` is `TRUE`? The options are `Robinson`, `Winkel tripel` or a proj4string. The default is `Winkel tripel`.
  #' @param tolerance List that provides tolerances during the construction of the graph:
  #' - `vertex_vertex` Vertices that are closer than this number are merged (default = 1e-7).
  #' - `vertex_edge` If a vertex at the end of one edge is closer than this
  #' number to another edge, this vertex is connected to that edge
  #' (default = 1e-7). Previously `vertex_line`, which is now deprecated.
  #' - `edge_edge` If two edges at some point are closer than this number, a new
  #' vertex is added at that point and the two edges are connected (default = 0).
  #' - `vertex_line`, Deprecated. Use `vertex_edge` instead.
  #' - `line_line`, Deprecated. Use `edge_edge` instead.
  #'
  #' In case `longlat = TRUE`, the tolerances are given in `length_unit`.
  #' @param check_connected If `TRUE`, it is checked whether the graph is
  #' connected and a warning is given if this is not the case.
  #' @param remove_deg2 Set to `TRUE` to remove all vertices of degree 2 in the
  #' initialization. Default is `FALSE`.
  #' @param merge_close_vertices should an additional step to merge close vertices be done?
  #' @param factor_merge_close_vertices Which factor to be multiplied by tolerance `vertex_vertex` when merging close vertices at the additional step?
  #' @param remove_circles All circlular edges with a length smaller than this number
  #' are removed. If `TRUE`, the `vertex_vertex` tolerance will be used. If `FALSE`, no circles will be removed.
  #' @param verbose Print progress of graph creation.
  #' @param lines `r lifecycle::badge("deprecated")` Use `edges` instead.
  #' @details A graph object can be initialized in two ways. The first method
  #' is to specify V and E. In this case, all edges are assumed to be straight
  #' lines. The second option is to specify the graph via the `lines` input.
  #' In this case, the vertices are set by the end points of the lines.
  #' Thus, if two lines are intersecting somewhere else, this will not be
  #' viewed as a vertex.
  #' @return A `metric_graph` object.
  initialize = function(edges = NULL,
                        V = NULL,
                        E = NULL,
                        vertex_unit = NULL,
                        length_unit = vertex_unit,
                        edge_weights = 1,
                        longlat = FALSE,
                        crs = NULL,
                        proj4string = NULL,
                        which_longlat = "sp",
                        project = TRUE,
                        project_data = FALSE,
                        which_projection = "Winkel tripel",
                        tolerance = list(vertex_vertex = 1e-7,
                                         vertex_edge = 1e-7,
                                         edge_edge = 0),
                        check_connected = TRUE,
                        remove_deg2 = FALSE,
                        merge_close_vertices = TRUE,
                        factor_merge_close_vertices = 1,
                        remove_circles = TRUE,
                        verbose = FALSE,
                        lines = deprecated()) {

      start_construction_time <- Sys.time()
      if (lifecycle::is_present(lines)) {
         if (is.null(edges)) {
           lifecycle::deprecate_warn("1.2.0", "metric_graph$new(lines)", "metric_graph$new(edges)",
             details = c("`lines` was provided but not `edges`. Setting `edges <- lines`.")
           )
           edges <- lines
         } else {
           lifecycle::deprecate_warn("1.2.0", "metric_graph$new(lines)", "metric_graph$new(edges)",
             details = c("Both `edges` and `lines` were provided. Only `edges` will be considered.")
           )
         }
         lines <- NULL
       }

      if (is.null(tolerance$vertex_edge) && !is.null(tolerance$vertex_line)) {
           lifecycle::deprecate_warn("1.2.0", "metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
             details = c("`tolerance$vertex_line` was provided but not `tolerance$vertex_edge`. Setting `tolerance$vertex_edge <- tolerance$vertex_line`.")
           )
           tolerance$vertex_edge <- tolerance$vertex_line
           tolerance$vertex_line <- NULL
         } else if(!is.null(tolerance$vertex_edge) && !is.null(tolerance$vertex_line)) {
           lifecycle::deprecate_warn("1.2.0","metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
             details = c("Both `tolerance$vertex_edge` and `tolerance$vertex_line` were provided. Only `tolerance$vertex_edge` will be considered.")
           )
            tolerance$vertex_line <- NULL
         }

        if (is.null(tolerance$edge_edge) && !is.null(tolerance$line_line)) {
           lifecycle::deprecate_warn("1.2.0", "metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
             details = c("`tolerance$line_line` was provided but not `tolerance$edge_edge`. Setting `tolerance$edge_edge <- tolerance$line_line`.")
           )
           tolerance$edge_edge <- tolerance$line_line
           tolerance$line_line <- NULL
         } else if(!is.null(tolerance$edge_edge) && !is.null(tolerance$line_line)) {
           lifecycle::deprecate_warn("1.2.0","metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
             details = c("Both `tolerance$edge_edge` and `tolerance$line_line` were provided. Only `tolerance$edge_edge` will be considered.")
           )
          tolerance$line_line <- NULL
         }

      valid_units_vertex <- c("m", "km", "miles", "degrees")
      valid_units_length <- c("m", "km", "miles")

      if(!(which_longlat %in% c("sp", "sf"))){
        stop("The options for 'which_longlat' are 'sp' and 'sf'!")
      }

      if(longlat){
        private$longlat <- TRUE
        private$which_longlat <- which_longlat
      }

      if(longlat && (which_longlat == "sp") && is.null(proj4string)){
        proj4string <- sp::CRS("+proj=longlat +datum=WGS84")
        private$crs <- sf::st_crs(proj4string)
        private$proj4string <- proj4string
      }

      if(longlat && (which_longlat == "sf") && is.null(crs)){
        crs <- sf::st_crs(4326)
        private$crs <- crs
      }

    # private$longlat <- longlat

    if((is.null(vertex_unit) && !is.null(length_unit)) || (is.null(length_unit) && !is.null(vertex_unit))){
      stop("If one of 'vertex_unit' or 'length_unit' is NULL, then the other must also be NULL.")
    }

    if(!is.null(vertex_unit)){
      vertex_unit <- vertex_unit[[1]]
      if(!is.character(vertex_unit)){
        stop("'vertex_unit' must be a string!")
      }
      if(!(vertex_unit %in% valid_units_vertex)){
        stop(paste("The possible options for 'vertex_unit' are ", toString(valid_units_vertex)))
      }
      private$vertex_unit <- vertex_unit
    }

    if(!is.null(length_unit)){
      length_unit <- length_unit[[1]]
      if(!is.character(length_unit)){
        stop("'length_unit' must be a string!")
      }
      if(length_unit == "degrees"){
        length_unit <- "km"
      }
      if(!(length_unit %in% valid_units_length)){
        stop(paste("The possible options for 'length_unit' are ", toString(valid_units_length)))
      }
      private$length_unit <- length_unit
    }

    if(longlat){
      private$vertex_unit <- "degrees"
      private$length_unit <- "km"
    } else if(!is.null(vertex_unit)){
        if(private$vertex_unit == "degrees"){
          longlat <- TRUE
        }
    }

    factor_unit <- process_factor_unit(private$vertex_unit, private$length_unit)



    tolerance_default = list(vertex_vertex = 1e-7,
                             vertex_edge = 1e-7,
                              edge_edge = 0)



    for(i in 1:length(tolerance_default)){
      if(!(names(tolerance_default)[i] %in% names(tolerance))){
        tolerance[names(tolerance_default)[i]] <- tolerance_default[i]
      }
    }

    if(is.null(tolerance$buffer_edge_edge)){
      tolerance$buffer_edge_edge <- max(tolerance$edge_edge/2 - 1e-10,0)
    }
    max_tol <- max(c(tolerance$vertex_vertex,
                     tolerance$vertex_edge,
                     tolerance$edge_edge))

    private$tolerance <- tolerance

    PtE_tmp_edge_edge <- NULL
    PtE_tmp_edge_vertex <- NULL

    if(is.null(edges) && is.null(V) && is.null(E)) {
      edges <- logo_lines()
    }
    if(!is.null(edges)){
      if(!is.null(V) || !is.null(E)){
        warning("object initialized from edges, then E and V are ignored")
      }
      if (inherits(edges,"SpatialLinesDataFrame")) {
        tmp_lines = SpatialLines(edges@lines)
        self$edges <- lapply(1:length(tmp_lines), function(i){tmp_lines@lines[[i]]@Lines[[1]]@coords})
      } else if (inherits(edges,"SpatialLines")) {
        self$edges = lapply(1:length(edges), function(i){edges@lines[[i]]@Lines[[1]]@coords})
      } else if(inherits(edges, "MULTILINESTRING")) {
        coords_multilinestring <- sf::st_coordinates(edges)
        lines_ids <- unique(coords_multilinestring[,"L1"])
        self$edges <- lapply(1:length(lines_ids), function(i){coords_multilinestring[coords_multilinestring[,"L1"]==i ,1:2]})
      } else if(is.list(edges)){
        self$edges <- check_lines_input(edges)
      } else {
        stop("edges should either be a list, or of class MULTILINESTRING, SpatialLines or SpatialLinesDataFrame")
      }

    } else {
      if(is.null(V) || is.null(E)){
        stop("You must supply edges or V and E")
      }
      if(ncol(V)!=2 || ncol(E)!=2){
        stop("V and E must have two columns!")
      }
      edges <- list()
      for(i in 1:dim(E)[1]) {
        edges[[i]] <- rbind(V[E[i,1], ], V[E[i,2], ])
      }
      self$edges <- edges
    }


    if(verbose){
      message("Setup edges and merge close vertices")
    }

    t <- system.time(
      private$line_to_vertex(tolerance = tolerance$vertex_vertex,
                           longlat = private$longlat, factor_unit, verbose=verbose,
                           private$crs, private$proj4string, which_longlat, private$length_unit, private$vertex_unit,
                           project, which_projection, project_data)
      )


    if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
    }


    if(length(self$edges) > 1){

    if (tolerance$edge_edge > 0) {
    private$addinfo <- TRUE

    if(verbose){
      message("Find edge-edge intersections")
    }

    t <- system.time(
      points_add <- private$find_edge_edge_points(tol = tolerance$edge_edge, verbose=verbose,
      crs=private$crs, proj4string = private$proj4string, longlat=private$longlat, fact = factor_unit, which_longlat = which_longlat)
      )

    if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
    }

    PtE <- points_add$PtE

    PtE[,2] <- PtE[,2]/self$edge_lengths[PtE[,1]]

    filter_tol <- ((PtE[,2] > max_tol/self$edge_lengths[PtE[,1]]) &
                     (PtE[,2] < 1- max_tol/self$edge_lengths[PtE[,1]]))

    PtE <- PtE[filter_tol,]

    if(!is.null(PtE)){
      if(nrow(PtE) == 0){
        PtE <- NULL
      }
    }

    if(!is.null(PtE)){
      if(verbose){
        message(sprintf("Add %d new vertices", nrow(PtE)))
      }

      PtE <- na.omit(PtE)

      t <- system.time(
      private$add_vertices(PtE, tolerance = tolerance$edge_edge, verbose = verbose)
      )

      if(verbose){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    }

    private$clear_initial_info()
    }

    if(tolerance$vertex_edge > 0){
      private$addinfo <- TRUE
      if(verbose){
        message("Snap vertices to close edges")
      }

      t <- system.time(
        PtE_tmp <- private$coordinates_multiple_snaps(XY = self$V,
                                              tolerance = tolerance$vertex_edge, verbose = verbose,
      crs=private$crs, proj4string = private$proj4string, longlat=private$longlat, fact = factor_unit, which_longlat = which_longlat)
        )

      if(verbose){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
      edge_length_filter <- self$edge_lengths[PtE_tmp[,1]]

      filter_tol <- ((PtE_tmp[,2] > max_tol/edge_length_filter) &
                       (PtE_tmp[,2] < 1- max_tol/edge_length_filter))

      PtE_tmp <- PtE_tmp[filter_tol,,drop = FALSE]
      PtE_tmp <- unique(PtE_tmp)
      PtE_tmp <- PtE_tmp[order(PtE_tmp[,1], PtE_tmp[,2]),,drop = FALSE]

      if(!is.null(PtE_tmp)){
        if(nrow(PtE_tmp) == 0){
          PtE_tmp <- NULL
        }
      }

      if(!is.null(PtE_tmp)){
        if(verbose){
          message(sprintf("Add %d new vertices", nrow(PtE_tmp)))
        }

        PtE_tmp <- na.omit(PtE_tmp)

        t <- system.time(
          private$add_vertices(PtE_tmp, tolerance = tolerance$vertex_edge, verbose=verbose)
          )

        if(verbose){
          message(sprintf("time: %.3f s", t[["elapsed"]]))
        }
      }
      private$clear_initial_info()
    }


    if(merge_close_vertices){
      private$merge_close_vertices(factor_merge_close_vertices * tolerance$vertex_vertex, factor_unit)
    }

    if(is.logical(remove_circles)){
      if(remove_circles){
        private$remove_circles(tolerance$vertex_vertex, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
      }
    } else {
        private$remove_circles(remove_circles, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
        remove_circles <- TRUE
    }

    if(merge_close_vertices || remove_circles){
      if(verbose){
        message("Recomputing edge lengths")
      }
      t <- system.time({
        self$edge_lengths <- private$compute_lengths(private$longlat, private$length_unit, private$crs, private$proj4string, private$which_longlat, private$vertex_unit, project_data)
      })
       if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
       }
    }
    # End of cond of having more than 1 edge
    }

    if (remove_deg2) {
      if (verbose) {
        message("Remove degree 2 vertices")
      }
      t <- system.time(
        self$prune_vertices(verbose = verbose)
      )
      if(verbose){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    }

    # Cleaning the edges

    if(verbose){
      message("Post-processing the edges")
    }

    t <- system.time(
          self$edges <- lapply(self$edges, function(edge){
            tmp_edge <- edge[1:(nrow(edge)-1),]
            tmp_edge <- unique(tmp_edge)
            tmp_edge <- rbind(tmp_edge, edge[nrow(edge),,drop=FALSE])
            if(nrow(tmp_edge)>2){
              tmp_edge <- tmp_edge[2:nrow(tmp_edge),]
              tmp_edge <- unique(tmp_edge)
              tmp_edge <- rbind(edge[1,,drop=FALSE], tmp_edge)
            }
            rownames(tmp_edge) <- NULL
            return(tmp_edge)
          }
            )
    )

    if(verbose){
          message(sprintf("time: %.3f s", t[["elapsed"]]))
    }

    # Checking if there is some edge with infinite length
    if(any(!is.finite(self$edge_lengths))){
      warning("There is at least one edge of infinite length. Please, consider redefining the graph.")
    }

    # Checking if there is some edge with zero length
    if(any(self$edge_lengths == 0)){
      warning("There is at least one edge of length zero. Please, consider redefining the graph.")
    }

    end_construction_time <- Sys.time()
    construction_time <- end_construction_time - start_construction_time

    if(verbose){
      message(sprintf('Total construction time: %.2f %s', construction_time, units(construction_time)))

    }

    # Checking if graph is connected
    if (check_connected) {
      g <- graph(edges = c(t(self$E)), directed = FALSE)
      components <- igraph::clusters(g, mode="weak")
      nc <- components$no
      if(nc>1){
        message("The graph is disconnected. You can use the function 'graph_components' to obtain the different connected components.")
        private$connected = FALSE
      }
    }
    self$set_edge_weights(weights = edge_weights)
    private$create_update_vertices()

    # Adding IDs to edges and setting up their class

    # for(i in 1:length(self$edges)){
    #   attr(self$edges[[i]], "id") <- i
    #   class(self$edges[[i]]) <- "metric_graph_edge"
    # }

    # Cloning the initial graph

    private$initial_graph <- self$clone()

    # Cloning again to add the initial graph to the initial graph
    private$initial_graph <- self$clone()

  },

  #' @description Sets the edge weights
  #' @param weights Either a number, a numerical vector with length given by the number of edges, providing the edge weights, or a `data.frame` with the number of rows being equal to the number of edges, where
  #' each row gives a vector of weights to its corresponding edge.
  #' @return No return value. Called for its side effects.

  set_edge_weights = function(weights = rep(1, self$nE)){
    if(!is.vector(weights) && !is.data.frame(weights)){
      stop("'weights' must be either a vector or a data.frame!")
    }

    edge_lengths_ <- self$get_edge_lengths()

    if(is.vector(weights)){
      if ( (length(weights) != 1) && (length(weights) != self$nE)){
        stop(paste0("The length of 'weights' must be either 1 or ", self$nE))
      }
      if(length(weights)==1){
        private$edge_weights <- rep(weights, self$nE)
      } else{
        private$edge_weights <- weights
      }
    } else{
      if(nrow(weights) != self$nE){
        stop("The number of rows of weights must be equal to the number of edges!")
      }
      private$edge_weights <- weights
    }
    self$edges <- lapply(1:self$nE, function(i){
      edge <- self$edges[[i]]
      if(is.vector(private$edge_weights)){
        attr(edge,"weight") <- private$edge_weights[i]
      } else{
        attr(edge,"weight") <- private$edge_weights[i,]
      }
      attr(edge, "longlat") <- private$longlat
      attr(edge, "crs") <- private$crs$input
      attr(edge, "length") <- edge_lengths_[i]
      attr(edge, "id") <- i
      class(edge) <- "metric_graph_edge"
      return(edge)
    })
    class(self$edges) <- "metric_graph_edges"

  },


  #' @description Gets the edge weights
  #' @return A vector containing the edge weights.

  get_edge_weights = function(){
    return(private$edge_weights)
  },



  #' @description Gets vertices with incompatible directions
  #' @return A vector containing the vertices with incompatible directions.

  get_vertices_incomp_dir = function(){
    start.deg <- end.deg <- rep(0,self$nV)
    for(i in 1:self$nV) {
      start.deg[i] <- sum(self$E[,1]==i)
      end.deg[i] <- sum(self$E[,2]==i)
    }

    degrees <- self$get_degrees()

    # Finding problematic vertices, that is, vertices with incompatible directions
    # They will not be pruned.
    problematic <- (degrees > 1) & (start.deg == 0 | end.deg == 0)
    return(which(problematic))
  },

  #' @description Prints a summary of various informations of the graph
  #' @param messages Should message explaining how to build the results be given for missing quantities?
  #' @param compute_characteristics Should the characteristics of the graph be computed?
  #' @param check_euclidean Check if the graph has Euclidean edges?
  #' @param check_distance_consistency Check the distance consistency assumption?
  #' @return No return value. Called for its side effects.

  summary = function(messages = FALSE, compute_characteristics = TRUE, check_euclidean = TRUE, check_distance_consistency = TRUE){
    if(compute_characteristics){
      self$compute_characteristics()
    }
    if(check_distance_consistency){
      self$check_distance_consistency()
    }
    if(check_euclidean){
      self$check_euclidean()
    }
    cat("A metric graph object with:\n\n")
    cat("Vertices:\n")
    cat("\t Total:", self$nV,"\n")
    degrees <- self$get_degrees()
    cat("\t")
    degrees_u <- sort(unique(degrees))
    for(i in 1:length(degrees_u)){
      if((i>1) && (i%%5 == 1)){
        cat("\n\t")
      }
      cat(paste0(" Degree ", degrees_u[i],": ",sum(degrees == degrees_u[i]), "; "))
    }
    cat("\n")

    cat("\t With incompatible directions: ", length(self$get_vertices_incomp_dir()), "\n\n")
    cat("Edges: \n")
    cat("\t Total:", self$nE,"\n")
    cat("\t Lengths: \n")
    cat("\t\t Min:", min(self$get_edge_lengths()), " ; Max:", max(self$get_edge_lengths()), " ; Total:", sum(self$get_edge_lengths()), "\n")
    cat("\t Weights: \n")
    if(is.vector(private$edge_weights)){
      cat("\t\t Min:", min(private$edge_weights), " ; Max:", max(private$edge_weights), "\n")
    } else{
      if(!is.null(colnames(private$edge_weights))){
        cat("\t\t Columns:", colnames(private$edge_weights),"\n")
      } else{
        cat("\t\t Number of columns:", ncol(private$edge_weights), "\n")
      }
    }
    cat("\t That are circles: ", sum(self$E[,1] == self$E[,2]), "\n\n")
    cat("Graph units: \n")
    cat("\t Vertices unit: ", ifelse(is.null(private$vertex_unit), "None", private$vertex_unit), " ; Lengths unit: ", ifelse(is.null(private$length_unit), "None", private$length_unit), "\n\n")
    cat("Longitude and Latitude coordinates: ", private$longlat)
    if(private$longlat){
      cat("\n\t Which spatial package: ", private$which_longlat, "\n")
      cat("\t CRS: ", private$crs$input)
    }
    cat("\n\n")
    if(!is.null(self$characteristics)) {
      cat("Some characteristics of the graph:\n")
      if(self$characteristics$connected){
        cat("\t Connected: TRUE\n")
      } else {
        cat("\t Connected: FALSE\n")
      }
      if(self$characteristics$has_loops){
        cat("\t Has loops: TRUE\n")
      } else {
        cat("\t Has loops: FALSE\n")
      }
      if(self$characteristics$has_multiple_edges){
        cat("\t Has multiple edges: TRUE\n")
      } else {
        cat("\t Has multiple edges: FALSE\n")
      }
      if(self$characteristics$is_tree){
        cat("\t Is a tree: TRUE\n")
      } else {
        cat("\t Is a tree: FALSE\n")
      }
      if(!is.null(self$characteristics$distance_consistency)){
        if(self$characteristics$distance_consistency){
          cat("\t Distance consistent: TRUE\n")
        } else {
          cat("\t Distance consistent: FALSE\n")
        }
      } else{
        cat("\t Distance consistent: unknown\n")
        if(messages){
          message("To check if the graph satisfies the distance consistency, run the `check_distance_consistency()` method.")
        }
      }
      if(!is.null(self$characteristics$euclidean)){
        if(self$characteristics$euclidean){
          cat("\t Has Euclidean edges: TRUE\n")
        } else {
          cat("\t Has Euclidean edges: FALSE\n")
        }
      } else{
        cat("\t Has Euclidean edges: unknown\n")
        if(messages){
          message("To check if the graph has Euclidean edges, run the `check_euclidean()` method.")
        }
      }
    } else{
      cat("Some characteristics of the graph: Not computed.\n")
      if(messages){
        message("To compute the characteristics, run the `compute_characteristics()` method.")
      }
    }
    cat("\n")
    cat("Computed quantities inside the graph: \n")
    cat("\t Laplacian: ", !is.null(self$Laplacian), " ; Geodesic distances: ", !is.null(self$geo_dist), "\n")
    if(is.null(self$Laplacian)){
      if(messages){
        message("To compute the Laplacian, run the 'compute_laplacian()' method.")
      }
    }
    if(is.null(self$geo_dist)){
      if(messages){
        message("To compute the geodesic distances, run the 'compute_geodist()' method.")
      }
    }
    cat("\t Resistance distances: ", !is.null(self$res_dist), " ; Finite element matrices: ", !is.null(self$mesh$C), "\n")
    if(is.null(self$res_dist)){
      if(messages){
        message("To compute the resistance distances, run the 'compute_resdist()' method.")
      }
    }
    if(is.null(self$mesh$C)){
      if(messages){
        message("To compute the finite element matrices, run the 'compute_fem()' method.")
      }
    }
    cat("\n")
    if(is.null(self$mesh)){
      cat("Mesh: The graph has no mesh! \n")
      if(messages){
        message("To build the mesh, run the 'build_mesh()' method.")
      }
    } else{
      cat("Mesh: \n")
      cat("\t Max h_e: ", max(self$mesh$h_e), " ; Min n_e: ", min(self$mesh$n_e), "\n")
    }
    cat("\n")
    if(is.null(private$data)){
      cat("Data: The graph has no data!\n")
      if(messages){
        message("To add observations, use the 'add_observations()' method.")
      }
    } else{
      cat("Data: \n")
      col_names_valid <- grepl("^[^.]+$", names(private$data))
      cat("\t Columns: ", names(private$data)[col_names_valid], "\n")
      cat("\t Groups: ", ifelse(is.null(private$group_col), "None", private$group_col), "\n")
    }
    cat("\n")
    cat("Tolerances: \n")
    cat("\t vertex-vertex: ", private$tolerance$vertex_vertex, "\n")
    cat("\t vertex-edge: ", private$tolerance$vertex_edge, "\n")
    cat("\t edge-edge: ", private$tolerance$edge_edge, "\n")
  },


  #' @description Prints various characteristics of the graph
  #' @return No return value. Called for its side effects.

  print = function() {
    cat("A metric graph with ", self$nV, " vertices and ", self$nE, " edges.\n\n")
    cat("Vertices:\n")
    degrees <- self$get_degrees()
    cat("\t")
    degrees_u <- sort(unique(degrees))
    for(i in 1:length(degrees_u)){
      if((i>1) && (i%%5 == 1)){
        cat("\n\t")
      }
      cat(paste0(" Degree ", degrees_u[i],": ",sum(degrees == degrees_u[i]), "; "))
    }
    cat("\n")

    cat("\t With incompatible directions: ", length(self$get_vertices_incomp_dir()), "\n\n")
    cat("Edges: \n")
    cat("\t Lengths: \n")
    cat("\t\t Min:", min(self$get_edge_lengths()), " ; Max:", max(self$get_edge_lengths()), " ; Total:", sum(self$get_edge_lengths()), "\n")
    cat("\t Weights: \n")
    if(is.vector(private$edge_weights)){
      cat("\t\t Min:", min(private$edge_weights), " ; Max:", max(private$edge_weights), "\n")
    } else{
      if(!is.null(colnames(private$edge_weights))){
        cat("\t\t Columns:", colnames(private$edge_weights),"\n")
      } else{
        cat("\t\t Number of columns:", ncol(private$edge_weights), "\n")
      }
    }
    cat("\t That are circles: ", sum(self$E[,1] == self$E[,2]), "\n\n")
    cat("Graph units: \n")
    cat("\t Vertices unit: ", ifelse(is.null(private$vertex_unit), "None", private$vertex_unit), " ; Lengths unit: ", ifelse(is.null(private$length_unit), "None", private$length_unit), "\n\n")
    cat("Longitude and Latitude coordinates: ", private$longlat)
    if(private$longlat){
      cat("\n\t Which spatial package: ", private$which_longlat, "\n")
      cat("\t CRS: ", private$crs$input)
    }
    cat("\n\n")

    if(!is.null(self$characteristics)) {
      cat("Some characteristics of the graph:\n")
      if(self$characteristics$connected){
        cat("  Connected: TRUE\n")
      } else {
        cat("  Connected: FALSE\n")
      }
      if(self$characteristics$has_loops){
        cat("  Has loops: TRUE\n")
      } else {
        cat("  Has loops: FALSE\n")
      }
      if(self$characteristics$has_multiple_edges){
        cat("  Has multiple edges: TRUE\n")
      } else {
        cat("  Has multiple edges: FALSE\n")
      }
      if(self$characteristics$is_tree){
        cat("  Is a tree: TRUE\n")
      } else {
        cat("  Is a tree: FALSE\n")
      }
      if(!is.null(self$characteristics$distance_consistency)){
        if(self$characteristics$distance_consistency){
          cat("  Distance consistent: TRUE\n")
        } else {
          cat("  Distance consistent: FALSE\n")
        }
      } else{
        cat("  Distance consistent: unknown\n")
        message("To check if the graph satisfies the distance consistency, run the `check_distance_consistency()` method.")
      }
      if(!is.null(self$characteristics$euclidean)){
        if(self$characteristics$euclidean){
          cat("  Has Euclidean edges: TRUE\n")
        } else {
          cat("  Has Euclidean edges: FALSE\n")
        }
      } else{
        cat("  Has Euclidean edges: unknown\n")
        message("To check if the graph has Euclidean edges, run the `check_euclidean()` method.")
      }
    }
    invisible(self)
  },
  #' @description Computes various characteristics of the graph
  #' @param check_euclidean Also check if the graph has Euclidean edges? This essentially means that the distance consistency check will also be perfomed. If the graph does not have Euclidean edges due to another reason rather than the distance consistency, then it will already be indicated that the graph does not have Euclidean edges.
  #' @return No return value. Called for its side effects. The computed characteristics
  #' are stored in the `characteristics` element of the `metric_graph` object.
  compute_characteristics = function(check_euclidean = FALSE) {
    if(is.null(self$characteristics)){
      self$characteristics <- list()
    }

    #check for loops
    if(is.null(self$characteristics$has_loops)){
      if(sum(self$E[,1]==self$E[,2])>0) {
        self$characteristics$has_loops <- TRUE
      } else {
        self$characteristics$has_loops <- FALSE
      }
    }

    self$characteristics$connected <- private$connected

    #check for multiple edges
    if(is.null(self$characteristics$has_multiple_edges)){
      self$characteristics$has_multiple_edges <- FALSE
      k <- 1
      while(k < self$nV && self$characteristics$has_multiple_edges == FALSE) {
        ind <- which(self$E[,1]==k | self$E[,2]==k) #edges starting or ending in k
        if(length(ind) > length(unique(rowSums(self$E[ind,,drop=FALSE])))) {
          self$characteristics$has_multiple_edges <- TRUE
        } else {
          k <- k + 1
        }
      }
    }


    #check for tree structure
    if(!self$characteristics$has_loops && !self$characteristics$has_multiple_edges){
      self$characteristics$is_tree <- self$is_tree()
    } else {
      self$characteristics$is_tree <- FALSE
    }

    if(!self$characteristics$connected || self$characteristics$has_loops || self$characteristics$has_multiple_edges){
      self$characteristics$euclidean <- FALSE
    } else if(self$characteristics$is_tree){
      self$characteristics$euclidean <- TRUE
    }

  },

  #' @description Check if the graph has Euclidean edges.
  #' @return Returns `TRUE` if the graph has Euclidean edges, or `FALSE` otherwise.
  #' The result is stored in the `characteristics` element of the `metric_graph` object.
  #' The result is displayed when the graph is printed.

  check_euclidean = function(){
    self$compute_characteristics()
    if(!is.null(self$characteristics$euclidean)){
      return(invisible(NULL))
    }

    if(is.null(self$characteristics$distance_consistency)){
      self$check_distance_consistency()
      if(self$characteristics$distance_consistency){
        self$characteristics$euclidean <- TRUE
      } else{
        self$characteristics$euclidean <- FALSE
      }
    } else{
      if(self$characteristics$distance_consistency){
        self$characteristics$euclidean <- TRUE
      } else{
        self$characteristics$euclidean <- FALSE
      }
    }
  },

  #' @description Checks distance consistency of the graph.
  #' @return No return value.
  #' The result is stored in the `characteristics` element of the `metric_graph` object.
  #' The result is displayed when the graph is printed.

  check_distance_consistency = function(){
    self$compute_characteristics()
    if(is.null(self$geo_dist)){
      self$geo_dist <- list()
    }

    if(is.null(self$geo_dist[[".vertices"]])){
      g <- graph(edges = c(t(self$E)), directed = FALSE)
      E(g)$weight <- self$edge_lengths
      self$geo_dist[[".vertices"]] <- distances(g)
    }

    geo_dist_edges <- self$geo_dist[[".vertices"]][self$E]
    if(any(abs(geo_dist_edges - self$edge_lengths) > 1e-8)){
      self$characteristics$distance_consistency <- FALSE
    } else{
      self$characteristics$distance_consistency <- TRUE
    }
  },

  #' @description Computes shortest path distances between the vertices in the
  #' graph
  #' @param full Should the geodesic distances be computed for all
  #' the available locations? If `FALSE`, it will be computed
  #' separately for the locations of each group.
  #' @param obs Should the geodesic distances be computed at the observation
  #' locations?
  #' @param group Vector or list containing which groups to compute the distance
  #' for. If `NULL`, it will be computed for all groups.
  #' @return No return value. Called for its side effects. The computed geodesic
  #' distances are stored in the `geo_dist` element of the `metric_graph` object.
  compute_geodist = function(full = FALSE, obs = TRUE, group = NULL) {
    if(is.null(self$geo_dist)){
      self$geo_dist <- list()
    }

    if(is.null(private$data)){
      obs <- FALSE
    }
    if(!obs){
      g <- graph(edges = c(t(self$E)), directed = FALSE)
      E(g)$weight <- self$edge_lengths
      self$geo_dist[[".vertices"]] <- distances(g)
    } else if(full){
      PtE_full <- self$get_PtE()
      self$geo_dist[[".complete"]] <- self$compute_geodist_PtE(PtE = PtE_full,
                                                              normalized = TRUE)
    } else{
      if(is.null(group)){
          group <- unique(private$data[[".group"]])
      }
      for(grp in group){
          data_grp <- select_group(private$data, grp)
          idx_notna <- idx_not_all_NA(data_grp)
          PtE_group <- cbind(data_grp[[".edge_number"]][idx_notna],
                     data_grp[[".distance_on_edge"]][idx_notna])
          self$geo_dist[[grp]] <- self$compute_geodist_PtE(PtE = PtE_group,
                                                              normalized = TRUE)
      }
    }
  },
  #' @description Computes shortest path distances between the vertices in the
  #' graph.
  #' @param PtE Points to compute the metric for.
  #' @param normalized are the locations in PtE in normalized distance?
  #' @param include_vertices Should the original vertices be included in the
  #' distance matrix?
  #' @return A matrix containing the geodesic distances.
  compute_geodist_PtE = function(PtE,
                                 normalized = TRUE,
                                 include_vertices = TRUE){
      graph.temp <- self$clone()
      graph.temp$clear_observations()
      df_temp <- data.frame(y = rep(0, dim(PtE)[1]),
                            edge_number = PtE[,1],
                            distance_on_edge = PtE[,2])
      if(sum(duplicated(df_temp))>0){
        warning("Duplicated locations were found when computing geodist. The returned values are given for unique locations.")
        df_temp <- unique(df_temp)
      }

      graph.temp$build_mesh(h = 1000)

      df_temp2 <- data.frame(y = 0, edge_number = graph.temp$mesh$VtE[1:nrow(self$V),1],
                                  distance_on_edge = graph.temp$mesh$VtE[1:nrow(self$V),2])

      df_temp$included <- TRUE
      temp_merge <- merge(df_temp, df_temp2, all = TRUE)

      df_temp$included <- NULL

      df_temp2 <- temp_merge[is.na(temp_merge["included"]),]

      nV_new <- sum(is.na(temp_merge["included"]))

      df_temp2$included <- NULL

      df_temp <- rbind(df_temp2, df_temp)

      df_temp[["__dummy"]] <- 1:nrow(df_temp)

      graph.temp$add_observations(data = df_temp,
                                     normalized = normalized)
      graph.temp$observation_to_vertex(mesh_warning = FALSE)
      g <- graph(edges = c(t(graph.temp$E)), directed = FALSE)
      E(g)$weight <- graph.temp$edge_lengths
      geodist_temp <- distances(g)
      geodist_temp <- geodist_temp[graph.temp$PtV, graph.temp$PtV]
      #Ordering back in the input order
      geodist_temp[graph.temp$.__enclos_env__$private$data[["__dummy"]],graph.temp$.__enclos_env__$private$data[["__dummy"]]] <- geodist_temp
      if(!include_vertices){
        geodist_temp <- geodist_temp[(nV_new+1):nrow(geodist_temp), (nV_new+1):nrow(geodist_temp)]
      }
      return(geodist_temp)
  },

  #' @description Computes shortest path distances between the vertices in the
  #' mesh.
  #' @return No return value. Called for its side effects. The geodesic distances
  #' on the mesh are stored in `mesh$geo_dist` in the `metric_graph` object.
  compute_geodist_mesh = function() {
    g <- graph(edges = c(t(self$mesh$E)), directed = FALSE)
    E(g)$weight <- self$mesh$h_e
    self$mesh$geo_dist <- distances(g)
  },

  #' @description Computes the resistance distance between the observation
  #' locations.
  #' @param full Should the resistance distances be computed for all
  #' the available locations. If `FALSE`, it will be computed
  #' separately for the locations of each group.
  #' @param obs Should the resistance distances be computed at the observation
  #' locations?
  #' @param group Vector or list containing which groups to compute the distance
  #' for. If `NULL`, it will be computed for all groups.
  #' @return No return value. Called for its side effects. The geodesic distances
  #' are stored in the `res_dist` element of the `metric_graph` object.
  compute_resdist = function(full = FALSE, obs = TRUE, group = NULL) {
    self$res_dist <- list()
    if(is.null(private$data)){
      obs <- FALSE
    }
    if(!obs){
      graph.temp <- self$clone()
      graph.temp$build_mesh(h=1000)

      PtE <- graph.temp$mesh$VtE[1:nrow(self$V),]
      rm(graph.temp)
      self$res_dist[[".vertices"]] <- self$compute_resdist_PtE(PtE,
                                                                normalized=TRUE)
    } else if(full){
      PtE <- self$get_PtE()
      self$res_dist[[".complete"]] <- self$compute_resdist_PtE(PtE,
                                                                normalized=TRUE)
    } else{
      if(is.null(group)){
          group <- unique(private$data[[".group"]])
      }
      for(grp in group){
        data_grp <- select_group(private$data, grp)
        idx_notna <- idx_not_all_NA(data_grp)
        if(sum(idx_notna) == 0){
          stop("There are no non-NA observations.")
        }
        PtE <- cbind(data_grp[[".edge_number"]][idx_notna],
                     data_grp[[".distance_on_edge"]][idx_notna])
        self$res_dist[[as.character(grp)]] <- self$compute_resdist_PtE(PtE,
                                                                       normalized=TRUE)
      }
    }
  },

  #' @description Computes the resistance distance between the observation
  #' locations.
  #' @param PtE Points to compute the metric for.
  #' @param normalized Are the locations in PtE in normalized distance?
  #' @param include_vertices Should the original vertices be included in the
  #' Laplacian matrix?
  #' @return A matrix containing the resistance distances.
  compute_resdist_PtE = function(PtE,
                                 normalized = TRUE,
                                 include_vertices = FALSE) {
      graph.temp <- self$clone()
      graph.temp$clear_observations()
      df_temp <- data.frame(y = rep(0, dim(PtE)[1]),
                            edge_number = PtE[,1],
                            distance_on_edge = PtE[,2])
      if(sum(duplicated(df_temp))>0){
        warning("Duplicated locations were found when computing geodist. The returned values are given for unique locations.")
        df_temp <- unique(df_temp)
      }

      graph.temp$build_mesh(h = 1000)

      df_temp2 <- data.frame(y = 0,
                             edge_number = graph.temp$mesh$VtE[1:nrow(self$V), 1],
                             distance_on_edge = graph.temp$mesh$VtE[1:nrow(self$V), 2])

      df_temp$included <- TRUE
      temp_merge <- merge(df_temp, df_temp2, all = TRUE)

      df_temp$included <- NULL

      df_temp2 <- temp_merge[is.na(temp_merge["included"]),]

      nV_new <- sum(is.na(temp_merge["included"]))

      df_temp2$included <- NULL

      df_temp <- rbind(df_temp2, df_temp)

      df_temp[["__dummy"]] <- 1:nrow(df_temp)

      graph.temp$add_observations(data = df_temp,
                                     normalized = normalized)

        graph.temp$observation_to_vertex(mesh_warning=FALSE)
        graph.temp$compute_geodist(full=TRUE)
        geodist_temp <- graph.temp$geo_dist[[".complete"]]
        geodist_temp[graph.temp$PtV, graph.temp$PtV] <- geodist_temp

      L <- Matrix(0, graph.temp$nV, graph.temp$nV)

      for (i in 1:graph.temp$nE) {
        tmp <- -1 / geodist_temp[graph.temp$E[i, 1],
                                                         graph.temp$E[i, 2]]
        L[graph.temp$E[i, 2], graph.temp$E[i, 1]] <- tmp
        L[graph.temp$E[i, 1], graph.temp$E[i, 2]] <- tmp
      }
      for(i in 1:graph.temp$nV){
        L[i, i] <- -sum(L[i, -i])
      }
      L[1, 1] <- L[1, 1] + 1

      Li <- solve(L)
      R <- -2*Li + t(diag(Li)) %x% rep(1, graph.temp$nV) +
        t(rep(1, graph.temp$nV)) %x% diag(Li)

      R <- R[graph.temp$PtV, graph.temp$PtV]
      R[graph.temp$.__enclos_env__$private$data[["__dummy"]],graph.temp$.__enclos_env__$private$data[["__dummy"]]] <- R

      if(!include_vertices){
        R <- R[(nV_new+1):nrow(R), (nV_new+1):nrow(R)]
      }

      return(R)
  },

  #' @description Returns the degrees of the vertices in the metric graph.
  #' @param which If "degree", returns the degree of the vertex. If "indegree", returns the indegree,
  #' and if "outdegree", it returns the outdegree.
  #' @return A vector containing the degrees of the vertices.
  get_degrees = function(which = "degree"){
    which <- which[[1]]
    if(!(which %in% c("degree", "indegree", "outdegree"))){
      stop("'which' must be either 'degree', 'indegree' or 'outdegree'!")
    }
    if(which == "degree"){
      degrees <- sapply(self$vertices, function(vert){attr(vert, "degree")})
    } else if(which == "indegree"){
       degrees <- sapply(self$vertices, function(vert){attr(vert, "indegree")})     
    } else{
        degrees <- sapply(self$vertices, function(vert){attr(vert, "outdegree")})      
    }
    return(degrees)
  },

  #' @description Computes the relative positions of the coordinates of the edges and save it as an attribute to each edge. This improves the quality of plots obtained by the `plot_function()` method, however it might be costly to compute.
  #' @return No return value, called for its side effects.
  compute_PtE_edges = function(){
    edges_PtE <- lapply(self$edges, function(edge){self$coordinates(XY = edge)})
     for(j in 1:length(self$edges)){
        attr(self$edges[[j]], "PtE") <- edges_PtE[[j]]
    }
    return(invisible(NULL))
  },

  #' @description Computes the resistance metric between the vertices in the
  #' mesh.
  #' @return No return value. Called for its side effects. The geodesic distances
  #' on the mesh are stored in the `mesh$res_dist` element in the `metric_graph`
  #' object.
  compute_resdist_mesh = function() {
    if (is.null(self$mesh)) {
      stop("no mesh provided")
    }
    if(is.null(self$mesh$geo_dist)){
      self$compute_geodist_mesh()
    }
    L <- Matrix(0, dim(self$mesh$V)[1], dim(self$mesh$V)[1])
    for (i in 1:dim(self$mesh$E)[1]) {
      tmp <- -1 / self$mesh$geo_dist[self$mesh$E[i, 1], self$mesh$E[i, 2]]
      L[self$mesh$E[i, 2], self$mesh$E[i, 1]] <- tmp
      L[self$mesh$E[i, 1], self$mesh$E[i, 2]] <- tmp
    }
    for (i in 1:dim(self$mesh$V)[1]) {
      L[i, i] <- -sum(L[i, -i])
    }
    L[1, 1] <- L[1, 1] + 1

    Li <- solve(L)
    self$mesh$res_dist <- -2*Li + t(diag(Li)) %x% rep(1,dim(self$mesh$V)[1]) +
      t(rep(1,dim(self$mesh$V)[1])) %x% diag(Li)
  },


  #' @description Computes the weigthed graph Laplacian for the graph.
  #' @param full Should the resistance distances be computed for all
  #' the available locations. If `FALSE`, it will be computed
  #' separately for the locations of each group.
  #' @param obs Should the resistance distances be computed at the observation
  #' locations? It will only compute for locations in which there is at least one observations that is not NA.
  #' @param group Vector or list containing which groups to compute the
  #' Laplacian for. If `NULL`, it will be computed for all groups.
  #' @return No reutrn value. Called for its side effects. The Laplacian is stored
  #' in the `Laplacian` element in the `metric_graph` object.
  compute_laplacian = function(full = FALSE, obs = TRUE, group = NULL) {
    self$Laplacian <- list()
    if(is.null(private$data)){
      obs <- FALSE
    }
    if(!obs){
      graph.temp <- self$clone()
      graph.temp$build_mesh(h=1000)

      PtE <- graph.temp$mesh$VtE[1:nrow(self$V),]
      rm(graph.temp)
      self$Laplacian[[".vertices"]] <- private$compute_laplacian_PtE(PtE,
                                                            normalized = TRUE)
    } else if(full){
      PtE <- self$get_PtE()
      self$Laplacian[[".complete"]] <- private$compute_laplacian_PtE(PtE,
                                                            normalized = TRUE)
    } else{
      if(is.null(group)){
          group <- unique(private$data[[".group"]])
      }
      for(grp in group){
          data_grp <- select_group(private$data, grp)
          idx_notna <- idx_not_all_NA(data_grp)
          PtE <- cbind(data_grp[[".edge_number"]][idx_notna],
                       data_grp[[".distance_on_edge"]][idx_notna])
          if(nrow(PtE) == 0){
            stop("All the observations are NA.")
          }
          self$Laplacian[[grp]] <- private$compute_laplacian_PtE(PtE,
                                                              normalized = TRUE)
      }
    }
  },

  #' @description Removes vertices of degree 2 from the metric graph.
  #' @return No return value. Called for its side effects.
  #' @param verbose Show progress? Default is `FALSE`.
  #' @details
    #' Vertices of degree 2 are removed as long as the corresponding edges that
    #' would be merged are compatible in terms of direction.
    #'
  prune_vertices = function(verbose = FALSE){
    t <- system.time({
    degrees <- private$compute_degrees()$degrees
    
    # Finding problematic vertices, that is, vertices with incompatible directions
    # They will not be pruned.

    if(is.null(self$vertices)){
      start.deg <- end.deg <- rep(0,self$nV)
      for(i in 1:self$nV) {
        start.deg[i] <- sum(self$E[,1]==i)
        end.deg[i] <- sum(self$E[,2]==i)
      }
      problematic <- (degrees > 1) & (start.deg == 0 | end.deg == 0)
    } else{
      problematic <- sapply(self$vertices, function(vert){attr(vert,"problematic")})
    }

    res <- list(degrees = degrees, problematic = problematic)
    if(verbose){
      to.prune <- sum(res$degrees==2 & !res$problematic)
      k <- 1
      message(sprintf("removing %d vertices", to.prune))
      if(to.prune > 0) {
        # pb = txtProgressBar(min = 1, max = to.prune, initial = 1, style = 3)
        bar_prune <- msg_progress_bar(to.prune)
      }

    }

   while(sum(res$degrees==2 & !res$problematic)>0) {
     if(verbose && to.prune > 0){
      #  setTxtProgressBar(pb,k)
      bar_prune$increment()
       #message(sprintf("removing vertex %d of %d.", k, to.prune))
       k <- k + 1
     }
     res <- private$remove.first.deg2(res)
   }
   })
    if(verbose){
          message(sprintf("time: %.3f s", t[["elapsed"]]))
    }
    # if(verbose && to.prune > 0){
    #   close(pb)
    # }

   if(verbose){
    message("Updating attributes of the edges and vertices")
   }

   t <- system.time({
      private$create_update_vertices()
      for(i in 1:length(self$edges)){
         attr(self$edges[[i]], "id") <- i
         attr(self$edges[[i]], "longlat") <- private$longlat
         attr(self$edges[[i]], "crs") <- private$crs$input
         attr(self$edges[[i]], "length") <- self$edge_lengths[i]
         class(self$edges[[i]]) <- "metric_graph_edge"
        if(!is.null(private$length_unit)){
          units(attr(self$edges[[i]], "length")) <- private$length_unit
        }
        if(is.vector(private$edge_weights)){
          attr(self$edges[[i]], "weight") <- private$edge_weights[i]
        } else{
          attr(self$edges[[i]], "weight") <- private$edge_weights[i,]
        }
      }
   })
   if(verbose){
            message(sprintf("time: %.3f s", t[["elapsed"]]))
   }


   if(!is.null(private$data)){
    if(verbose){
      message("Updating data locations.")
    }
      t <- system.time({
      x_coord <- private$data[[".coord_x"]]
      y_coord <- private$data[[".coord_y"]]
      new_PtE <- self$coordinates(XY = cbind(x_coord, y_coord))
      group_vec <- private$data[[".group"]]
      private$data[[".edge_number"]] <- new_PtE[,1]
      private$data[[".distance_on_edge"]] <- new_PtE[,2]
      order_idx <- order(group_vec, new_PtE[,1], new_PtE[,2])
      private$data <- lapply(private$data, function(dat){dat[order_idx]})
      })
      if(verbose){
            message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
   }

   for(j in 1:length(self$edges)){
        attr(self$edges[[j]], "PtE") <- NULL
    }


   if(!is.null(self$mesh)){
    max_h <- max(self$mesh$h_e)
    self$mesh <- NULL
    self$build_mesh(h = max_h)
   }
   private$pruned <- TRUE
  },

  #' @description Gets the groups from the data.
  #' @param get_cols Should the names of the columns that created the group variable be returned?
  #' @return A vector containing the available groups in the internal data.

  get_groups = function(get_cols = FALSE){
    if(is.null(private$data)){
      warning("There is no data!")
      return(invisible(NULL))
    }
    if(get_cols){
      return(private$group_col)
    }
    return(unique(private$data[[".group"]]))
  },

  #' @description Gets PtE from the data.
  #' @param group For which group, should the PtE be returned? `NULL` means that all PtEs available will be returned.
  #' @param include_group Should the group be included as a column? If `TRUE`, the PtEs for each group will be concatenated, otherwise a single matrix containing the unique PtEs will be returned.
  #' @return A matrix with two columns, where the first column contains the edge
  #' number and the second column contains the distance on edge of the
  #' observation locations.
  get_PtE = function() {
    if(is.null(private$data)){
      warning("There is no data!")
      return(invisible(NULL))
    }
    group <- private$data[[".group"]]
    group <- which(group == group[1])
    PtE <- cbind(private$data[[".edge_number"]][group],
                 private$data[[".distance_on_edge"]][group])
    return(PtE)
  },

  #' @description Gets the edge lengths with the corresponding unit.
  #' @param unit If non-NULL, changes from `length_unit` from the graph construction to `unit`.
  #' @return a vector with the length unit (if the graph was constructed with a length unit).

  get_edge_lengths = function(unit = NULL){
    el <- self$edge_lengths
    units(el) <- private$length_unit
    if(!is.null(unit)){
      units(el) <- unit
    }
    return(el)
  },

  #' @description Gets the spatial locations from the data.
  #' @return A `data.frame` object with observation locations. If `longlat = TRUE`, the column names are lon and lat, otherwise the column names are x and y.
  get_locations = function(){
     if(is.null(private$data)){
      warning("There is no data!")
      return(invisible(NULL))
    }
    group <- private$data[[".group"]]
    group <- which(group == group[1])
    Spoints <- data.frame(x = private$data[[".coord_x"]][group], y = private$data[[".coord_y"]][group])
    if(private$longlat){
      colnames(Spoints) <- c("lon", "lat")
    }
    return(Spoints)
  },

  #' @description Adds observation locations as vertices in the graph.
  #' @param tolerance Observations locations are merged to a single vertex if
  #' they are closer than this number (given in relative edge distance between
  #' 0 and 1). The default is `1e-15`.
  #' @param share_weights Should the same weight be shared among the split edges? If `FALSE`, the weights will be removed, and a common weight given by 1 will be given.
  #' @param mesh_warning Display a warning if the graph structure change and the metric graph has a mesh object.
  #' @return No return value. Called for its side effects.
  observation_to_vertex = function(tolerance = 1e-15, mesh_warning = TRUE) {
    if(tolerance <= 0 || tolerance >=1){
      stop("tolerance should be between 0 and 1.")
    }
    private$temp_PtE <- self$get_PtE()
    n_group <- length(unique(private$data[[".group"]]))
    l <- length(private$temp_PtE[, 1])
    self$PtV <- rep(NA, l)
    self$nE <- nrow(self$E)
    for (i in 1:l) {
      e <- as.vector(private$temp_PtE[i, 1])
      t <- as.vector(private$temp_PtE[i, 2])
      l_e <- self$edge_lengths[e]
      if (abs(t) < tolerance) {
        private$temp_PtE[i, 2] <- 0
        self$PtV[i] <- self$E[e, 1]
      } else if (t > 1 - tolerance) {
        private$temp_PtE[i, 2] <- 1
        self$PtV[i] <- self$E[e, 2]
      } else {
        PtV_tmp <- private$split_edge(e, t, tolerance)
        if(!is.null(PtV_tmp)){
          self$PtV[i] <- PtV_tmp
        }
      }
    }
    self$PtV <- self$PtV[!is.na(self$PtV)]

    private$data[[".edge_number"]] <- rep(private$temp_PtE[,1],
                                        times = n_group)
    private$data[[".distance_on_edge"]] <- rep(private$temp_PtE[,2],
                                             times = n_group)

    tmp_df <- data.frame(PtE1 = private$data[[".edge_number"]],
              PtE2 = private$data[[".distance_on_edge"]],
              group = private$data[[".group"]])
    index_order <- order(tmp_df$group, tmp_df$PtE1, tmp_df$PtE2)
    private$data <- lapply(private$data, function(dat){ dat[index_order]})

    self$PtV <- self$PtV[index_order[1:length(self$PtV)]]

    private$temp_PtE <- NULL

    if (!is.null(self$geo_dist)) {
      self$geo_dist <- NULL
    }
    if (!is.null(self$res_dist)) {
      self$res_dist <- NULL
    }
    if (!is.null(self$CoB)) {
      self$buildC(2)
    }

    if (!is.null(self$mesh)) {
      self$mesh <- NULL
      if(mesh_warning){
        warning("Removing the existing mesh due to the change in the graph structure, please create a new mesh if needed.")
      }
    }

    private$create_update_vertices()
      
    # Updating the edge attributes
    self$set_edge_weights(weights = private$edge_weights)

    for(j in 1:length(self$edges)){
        attr(self$edges[[j]], "PtE") <- NULL
    }

  },

  #' @description Returns a list or a matrix with the mesh locations.
  #' @param bru Should an 'inlabru'-friendly list be returned?
  #' @param loc If `bru` is set to `TRUE`, the name of the location variable.
  #' The default name is 'loc'.
  #' @param normalized If TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default TRUE.
  #'
  #' @return A list or a matrix containing the mesh locations.
  get_mesh_locations = function(bru = FALSE, loc = NULL, normalized = TRUE) {
    if(is.null(self$mesh)){
      warning("There is no mesh!")
      return(invisible(NULL))
    }

    if(!bru){
      return(self$mesh$VtE)
    } else{
      if(is.null(loc)){
        stop("If bru is TRUE, then the loc argument must be provided!")
      }
      data_list <- list()
      tmp_VtE <- self$mesh$VtE
      if(!normalized){
        tmp_VtE[,2] <- tmp_VtE[,2] * self$edge_lengths[tmp_VtE[, 1]]
      }
      data_list[[loc]] <- tmp_VtE
      return(data_list)
    }
  },

  #' @description Clear all observations from the `metric_graph` object.
  #' @return No return value. Called for its side effects.
  clear_observations = function() {
    private$data <- NULL
    self$geo_dist <- NULL
    self$res_dist <- NULL
    self$PtV <- NULL
  },


  #' @description Process data to the metric graph data format.
  #' @param Spoints `SpatialPoints` or `SpatialPointsDataFrame` containing the
  #' observations. It may include the coordinates of the observations only, or
  #' the coordinates as well as the observations.
  #' @param data A `data.frame` or named list containing the observations. In
  #' case of groups, the data.frames for the groups should be stacked vertically,
  #' with a column indicating the index of the group. If `data` is not `NULL`,
  #' it takes priority over any eventual data in `Spoints`.
  #' @param edge_number Column (or entry on the list) of the `data` that
  #' contains the edge numbers. If not supplied, the column with name
  #' "edge_number" will be chosen. Will not be used if `Spoints` is not `NULL`.
  #' @param distance_on_edge Column (or entry on the list) of the `data` that
  #' contains the edge numbers. If not supplied, the column with name
  #' "distance_on_edge" will be chosen.  Will not be used if `Spoints` is not
  #' `NULL`.
  #' @param coord_x Column (or entry on the list) of the `data` that contains
  #' the x coordinate. If not supplied, the column with name "coord_x" will be
  #' chosen. Will not be used if `Spoints` is not `NULL` or if `data_coords` is
  #' `PtE`.
  #' @param coord_y Column (or entry on the list) of the `data` that contains
  #' the y coordinate. If not supplied, the column with name "coord_x" will be
  #' chosen. Will not be used if `Spoints` is not `NULL` or if `data_coords` is
  #' `PtE`.
  #' @param data_coords To be used only if `Spoints` is `NULL`. It decides which
  #' coordinate system to use. If `PtE`, the user must provide `edge_number` and
  #' `distance_on_edge`, otherwise if `spatial`, the user must provide
  #' `coord_x` and `coord_y`. The option `euclidean` is `r lifecycle::badge("deprecated")`. Use `spatial` instead.
  #' @param group Vector. If the data is grouped (for example measured at different time
  #' points), this argument specifies the columns (or entries on the list) in
  #' which the group variables are stored. It will be stored as a single column `.group` with the combined entries.
  #' @param group_sep separator character for creating the new group variable when grouping two or more variables.
  #' @param normalized if TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default FALSE. Will not be used if
  #' `Spoints` is not `NULL`.
  #' @param tibble Should the data be returned as a `tidyr::tibble`?
  #' @param tolerance Parameter to control a warning when adding observations.
  #' If the distance of some location and the closest point on the graph is
  #' greater than the tolerance, the function will display a warning.
  #' This helps detecting mistakes on the input locations when adding new data.
  #' @param verbose If `TRUE`, report steps and times.
  #' @return No return value. Called for its side effects. The observations are
  #' stored in the `data` element of the `metric_graph` object.

  process_data = function(Spoints = NULL,
                              data = NULL,
                              edge_number = "edge_number",
                              distance_on_edge = "distance_on_edge",
                              coord_x = "coord_x",
                              coord_y = "coord_y",
                              data_coords = c("PtE", "spatial"),
                              group = NULL,
                              group_sep = ".",
                              normalized = FALSE,
                              tibble = TRUE,
                              tolerance = max(self$edge_lengths)/2,
                              verbose = FALSE) {


        data_coords <- data_coords[[1]]
        if(!is.null(group)){
          group <- unique(group)
        }

    if(length(tolerance)>1){
      tolerance <- tolerance[[1]]
      warning("'tolerance' had more than one element, only the first one will be used.")
    }

        if(is.null(data)){
          if(is.null(Spoints)){
            stop("No data provided!")
          }
          if("SpatialPointsDataFrame"%in%is(Spoints)){
            data <- Spoints@data
          } else{
            stop("No data provided!")
          }
        }

        if(!is.null(data)){
          if(!is.list(data) && !is.data.frame(data)){
            stop("'data' must be either a list or a data.frame!")
          }
        }

        data <- as.list(data)


        if(is.null(Spoints)){
        if(data_coords == "PtE"){
          if(any( !(c(edge_number, distance_on_edge) %in% names(data)))){
            stop(paste("The data does not contain either the column", edge_number,"or the column",distance_on_edge))
          }
        } else{
          if(any( !(c(coord_x, coord_y) %in% names(data)))){
            stop(paste("The data does not contain either the column", coord_x,"or the column",coord_y))
          }
        }
        }

        if(!is.null(group)){
          if(!all(group%in%names(data))){
            stop("There were group variables that are not columns of 'data'!")
          }
          data_group_tmp <- lapply(group, function(lab){data[[lab]]})
          ord_tmp <- do.call(order, data_group_tmp)
          rm(data_group_tmp)
          data <- lapply(data, function(dat){dat[ord_tmp]})
          rm(ord_tmp)
          data[[".dummy_var"]] <- as.character(data[[group[1]]])
          if(length(group)>1){
            for(j in 2:length(group)){
              data[[".dummy_var"]] <- sapply(1:length(data[[".dummy_var"]]), function(i){paste0(data[[".dummy_var"]][i], group_sep,data[[group[j]]][i])})
            }
          }
          data[[".group"]] <- data[[".dummy_var"]]
          data[[".dummy_var"]] <- NULL
        }

        ## convert everything to PtE
        if(verbose){
          if(data_coords == "spatial" || !is.null(Spoints)){
          message("Converting data to PtE")
          if(private$longlat){
            message("This step may take long. If this step is taking too long consider pruning the vertices to possibly obtain some speed up.")
          }
          }
        }

        ## Check data for repeated observations
        if (!is.null(Spoints)){
            if(is.null(group)){
            data_tmp <- Spoints@coords
          } else{
            data_tmp <- cbind(Spoints@coords, data[[".group"]])
          }
        } else if(data_coords == "spatial"){
          if(is.null(group)){
            data_tmp <- cbind(data[[coord_x]], data[[coord_y]])
          } else{
            data_tmp <- cbind(data[[coord_x]], data[[coord_y]], data[[".group"]])
          }
        } else{
          if(is.null(group)){
            data_tmp <- cbind(data[[edge_number]], data[[distance_on_edge]])
          } else{
            data_tmp <- cbind(data[[edge_number]], data[[distance_on_edge]], data[[".group"]])
          }
        }

        if(nrow(unique(data_tmp)) != nrow(data_tmp)){
          warning("There is at least one 'column' of the data with repeated (possibly different) values at the same location for the same group variable. Only one of these values will be used. Consider using the group variable to differentiate between these values or provide different names for such variables.")
          if(data_coords == "spatial" || !is.null(Spoints)){
            warning("It is also possible that two different points were projected to the same location on the metric graph.")
          }
        }

        t <- system.time({
          if(!is.null(Spoints)){
            PtE <- self$coordinates(XY = Spoints@coords)
            XY_new <- self$coordinates(PtE = PtE, normalized = TRUE)
            # norm_XY <- max(sqrt(rowSums( (Spoints@coords-XY_new)^2 )))
            fact <- process_factor_unit(private$vertex_unit, private$length_unit)
            norm_XY <- compute_aux_distances(lines = Spoints@coords, points = XY_new, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit)
            # norm_XY <- max(norm_XY)
            # if(norm_XY > tolerance){
            #   warning("There was at least one point whose location is far from the graph,
            #   please consider checking the input.")
            # }
            rm(Spoints)
            far_points <- (norm_XY > tolerance)
            rm(norm_XY)
            data <- lapply(data, function(dat){dat[!far_points]})
            if(any(far_points)){
              warning("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance.")
            }
            PtE <- PtE[!far_points,,drop=FALSE]
            rm(far_points)
          } else{
            if(data_coords == "PtE"){
                PtE <- cbind(data[[edge_number]], data[[distance_on_edge]])
                if(!normalized){
                  PtE[, 2] <- PtE[,2] / self$edge_lengths[PtE[, 1]]
                }
              } else if(data_coords == "spatial"){
                point_coords <- cbind(data[[coord_x]], data[[coord_y]])
                PtE <- self$coordinates(XY = point_coords)
                XY_new <- self$coordinates(PtE = PtE, normalized = TRUE)
                # norm_XY <- max(sqrt(rowSums( (point_coords-XY_new)^2 )))
                fact <- process_factor_unit(private$vertex_unit, private$length_unit)
                norm_XY <- compute_aux_distances(lines = point_coords, points = XY_new, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit)
                # norm_XY <- max(norm_XY)
                far_points <- (norm_XY > tolerance)
                rm(norm_XY)
                data <- lapply(data, function(dat){dat[!far_points]})
                PtE <- PtE[!far_points,,drop=FALSE]
                # if(norm_XY > tolerance){
                #   warning("There was at least one point whose location is far from the graph,
                #     please consider checking the input.")
                #   }
                if(any(far_points)){
                  warning("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance.")
                }
                rm(far_points)
            } else{
                stop("The options for 'data_coords' are 'PtE' and 'spatial'.")
            }
          }
        })

      if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))

      message("Processing data")
    }

    t <- system.time({
     if(!is.null(group)){
       group_vector <- data[[".group"]]
     } else{
       group <- ".group"
       group_vector <- NULL
     }

    lapply(data, function(dat){if(nrow(matrix(PtE, ncol=2)) != length(dat)){
        stop(paste(dat,"has a different number of elements than the number of
                   coordinates!"))
       }})


    group_vals <- unique(group_vector)


    # n_group <- length(unique(group_vector))
    n_group <- length(group_vals)
    n_group <- ifelse(n_group == 0, 1, n_group)

    data[[edge_number]] <- NULL
    data[[distance_on_edge]] <- NULL
    data[[coord_x]] <- NULL
    data[[coord_y]] <- NULL
    data[[".group"]] <- NULL
    data[[".coord_x"]] <- NULL
    data[[".coord_y"]] <- NULL

    # Process the data (find all the different coordinates
    # across the different replicates, and also merge the new data to the old data)
    data <- process_data_add_obs(PtE, new_data = data, old_data = NULL,
                                        group_vector)
    ## convert to Spoints and add
    group_1 <- data[[".group"]]
    group_1 <- which(group_1 == group_1[1])
    PtE <- cbind(data[[".edge_number"]][group_1],
                 data[[".distance_on_edge"]][group_1])
    spatial_points <- self$coordinates(PtE = PtE, normalized = TRUE)
    data[[".coord_x"]] <- rep(spatial_points[,1], times = n_group)
    data[[".coord_y"]] <- rep(spatial_points[,2], times = n_group)
    if(tibble){
      data <- tidyr::as_tibble(data)
    }
    class(data) <- c("metric_graph_data", class(data))
    })
          if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
          }
          return(data)

                              },

  #' @description Add observations to the metric graph.
  #' @param Spoints `SpatialPoints` or `SpatialPointsDataFrame` containing the
  #' observations. It may include the coordinates of the observations only, or
  #' the coordinates as well as the observations.
  #' @param data A `data.frame` or named list containing the observations. In
  #' case of groups, the data.frames for the groups should be stacked vertically,
  #' with a column indicating the index of the group. If `data` is not `NULL`,
  #' it takes priority over any eventual data in `Spoints`.
  #' @param edge_number Column (or entry on the list) of the `data` that
  #' contains the edge numbers. If not supplied, the column with name
  #' "edge_number" will be chosen. Will not be used if `Spoints` is not `NULL`.
  #' @param distance_on_edge Column (or entry on the list) of the `data` that
  #' contains the edge numbers. If not supplied, the column with name
  #' "distance_on_edge" will be chosen.  Will not be used if `Spoints` is not
  #' `NULL`.
  #' @param coord_x Column (or entry on the list) of the `data` that contains
  #' the x coordinate. If not supplied, the column with name "coord_x" will be
  #' chosen. Will not be used if `Spoints` is not `NULL` or if `data_coords` is
  #' `PtE`.
  #' @param coord_y Column (or entry on the list) of the `data` that contains
  #' the y coordinate. If not supplied, the column with name "coord_x" will be
  #' chosen. Will not be used if `Spoints` is not `NULL` or if `data_coords` is
  #' `PtE`.
  #' @param data_coords To be used only if `Spoints` is `NULL`. It decides which
  #' coordinate system to use. If `PtE`, the user must provide `edge_number` and
  #' `distance_on_edge`, otherwise if `spatial`, the user must provide
  #' `coord_x` and `coord_y`. The option `euclidean` is `r lifecycle::badge("deprecated")`. Use `spatial` instead.
  #' @param group Vector. If the data is grouped (for example measured at different time
  #' points), this argument specifies the columns (or entries on the list) in
  #' which the group variables are stored. It will be stored as a single column `.group` with the combined entries.
  #' @param group_sep separator character for creating the new group variable when grouping two or more variables.
  #' @param normalized if TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default FALSE. Will not be used if
  #' `Spoints` is not `NULL`.
  #' @param clear_obs Should the existing observations be removed before adding the data?
  #' @param tibble Should the data be returned as a `tidyr::tibble`?
  #' @param tolerance Parameter to control a warning when adding observations.
  #' If the distance of some location and the closest point on the graph is
  #' greater than the tolerance, the function will display a warning.
  #' This helps detecting mistakes on the input locations when adding new data.
  #' @param verbose If `TRUE`, report steps and times.
  #' @return No return value. Called for its side effects. The observations are
  #' stored in the `data` element of the `metric_graph` object.
  add_observations = function(Spoints = NULL,
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
                              verbose = FALSE) {

    if(clear_obs){
      df_temp <- data
      self$clear_observations()
      data <- df_temp
    }

    if(length(tolerance)>1){
      tolerance <- tolerance[[1]]
      warning("'tolerance' had more than one element, only the first one will be used.")
    }

    if(!is.null(group)){
      group <- unique(group)
    }

    if(inherits(data, "metric_graph_data")){
      if(!any(c(".edge_number", ".distance_on_edge", ".group", ".coord_x", ".coord_y") %in% names(data))){
        warning("The data is of class 'metric_graph_data', but it is not a proper 'metric_graph_data' object. The data will be added as a regular data.")
        class(data) <- setdiff(class(data), "metric_graph_data")
      } else{
        data_coords <- "PtE"
        edge_number <- ".edge_number"
        distance_on_edge <- ".distance_on_edge"
        group <- ".group"
        normalized <- TRUE
      }
    }

        data_coords <- data_coords[[1]]
        if(data_coords == "euclidean"){
          lifecycle::deprecate_warn("1.2.0", "add_observations(data_coords = 'must be either PtE or spatial')")
          data_coords <- "spatial"
        }
        if(is.null(data)){
          if(is.null(Spoints)){
            stop("No data provided!")
          }
          if("SpatialPointsDataFrame"%in%is(Spoints)){
            data <- Spoints@data
          } else{
            stop("No data provided!")
          }
        }

        if(!is.null(data)){
          if(!is.list(data) && !is.data.frame(data)){
            stop("'data' must be either a list or a data.frame!")
          }
        }


        data <- as.list(data)

        if(is.null(Spoints)){
        if(data_coords == "PtE"){
          if(any( !(c(edge_number, distance_on_edge) %in% names(data)))){
            stop(paste("The data does not contain either the colum", edge_number,"or the column",distance_on_edge))
          }
        } else{
          if(any( !(c(coord_x, coord_y) %in% names(data)))){
            stop(paste("The data does not contain either the colum", coord_x,"or the column",coord_y))
          }
        }
        }

        if(!is.null(group)){
          if(!all(group%in%names(data))){
            stop("There were group variables that are not columns of 'data'!")
          }
          data_group_tmp <- lapply(group, function(lab){data[[lab]]})
          ord_tmp <- do.call(order, data_group_tmp)
          rm(data_group_tmp)
          data <- lapply(data, function(dat){dat[ord_tmp]})
          if(!is.null(Spoints)){
            Spoints@coords <- Spoints@coords[ord_tmp,]
          }
          rm(ord_tmp)
          data[[".dummy_var"]] <- as.character(data[[group[1]]])
          if(length(group)>1){
            for(j in 2:length(group)){
              data[[".dummy_var"]] <- sapply(1:length(data[[".dummy_var"]]), function(i){paste0(data[[".dummy_var"]][i], group_sep,data[[group[j]]][i])})
            }
          }
          data[[".group"]] <- data[[".dummy_var"]]
          data[[".dummy_var"]] <- NULL
        }


        ## convert everything to PtE
        if(verbose){
          if(data_coords == "spatial" || !is.null(Spoints)){
          message("Converting data to PtE")
          if(private$longlat){
            message("This step may take long. If this step is taking too long consider pruning the vertices to possibly obtain some speed up.")
          }
          }
        }



        ## Check data for repeated observations
        if (!is.null(Spoints)){
            if(is.null(group)){
            data_tmp <- Spoints@coords
          } else{
            data_tmp <- Spoints@coords
            data_tmp <- cbind(Spoints@coords, data[[".group"]])
          }
        } else if(data_coords == "spatial"){
          if(is.null(group)){
            data_tmp <- cbind(data[[coord_x]], data[[coord_y]])
          } else{
            data_tmp <- cbind(data[[coord_x]], data[[coord_y]], data[[".group"]])
          }
        } else{
          if(is.null(group)){
            data_tmp <- cbind(data[[edge_number]], data[[distance_on_edge]])
          } else{
            data_tmp <- cbind(data[[edge_number]], data[[distance_on_edge]], data[[".group"]])
          }
        }

        if(nrow(unique(data_tmp)) != nrow(data_tmp)){
          warning("There is at least one 'column' of the data with repeated (possibly different) values at the same location for the same group variable. Only one of these values will be used. Consider using the group variable to differentiate between these values or provide different names for such variables.")
        }

        t <- system.time({
          if(!is.null(Spoints)){
            PtE <- self$coordinates(XY = Spoints@coords)
            XY_new <- self$coordinates(PtE = PtE, normalized = TRUE)
            # norm_XY <- max(sqrt(rowSums( (Spoints@coords-XY_new)^2 )))
            fact <- process_factor_unit(private$vertex_unit, private$length_unit)
            norm_XY <- compute_aux_distances(lines = Spoints@coords, points = XY_new, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit)
            rm(Spoints)
            # norm_XY <- max(norm_XY)
            # if(norm_XY > tolerance){
            #   warning("There was at least one point whose location is far from the graph,
            #   please consider checking the input.")
            # }
            far_points <- (norm_XY > tolerance)
            rm(norm_XY)
            data <- lapply(data, function(dat){dat[!far_points]})
            if(any(far_points)){
              warning("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance.")
            }
            PtE <- PtE[!far_points, ,drop=FALSE]
            rm(far_points)
          } else{
            if(data_coords == "PtE"){
                PtE <- cbind(data[[edge_number]], data[[distance_on_edge]])
                if(!normalized){
                  PtE[, 2] <- PtE[,2] / self$edge_lengths[PtE[, 1]]
                }
              } else if(data_coords == "spatial"){
                point_coords <- cbind(data[[coord_x]], data[[coord_y]])
                PtE <- self$coordinates(XY = point_coords)
                XY_new <- self$coordinates(PtE = PtE, normalized = TRUE)
                # norm_XY <- max(sqrt(rowSums( (point_coords-XY_new)^2 )))
                fact <- process_factor_unit(private$vertex_unit, private$length_unit)
                norm_XY <- compute_aux_distances(lines = point_coords, points = XY_new, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit)
                # norm_XY <- max(norm_XY)
                # if(norm_XY > tolerance){
                #   warning("There was at least one point whose location is far from the graph,
                #     please consider checking the input.")
                #   }
                far_points <- (norm_XY > tolerance)
                rm(norm_XY)
                data <- lapply(data, function(dat){dat[!far_points]})
                if(any(far_points)){
                  warning("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance.")
                }
                PtE <- PtE[!far_points,,drop=FALSE]
                rm(far_points)
            } else{
                stop("The options for 'data_coords' are 'PtE' and 'spatial'.")
            }
          }
        })

      if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))

      message("Processing data")
    }

    t <- system.time({
     if(!is.null(group)){
       group_vector <- data[[".group"]]
     } else{
       group_vector <- NULL
     }

    lapply(data, function(dat){if(nrow(matrix(PtE, ncol=2)) != length(dat)){
        stop(paste(dat,"has a different number of elements than the number of
                   coordinates!"))
       }})

    if(!is.null(private$data[[".group"]])){
      group_vals <- unique(private$data[[".group"]])
      group_vals <- unique(union(group_vals, group_vector))
    } else{
      group_vals <- unique(group_vector)
    }

    # n_group <- length(unique(group_vector))
    n_group <- length(group_vals)
    n_group <- ifelse(n_group == 0, 1, n_group)

    data[[edge_number]] <- NULL
    data[[distance_on_edge]] <- NULL
    data[[coord_x]] <- NULL
    data[[coord_y]] <- NULL
    data[[".group"]] <- NULL
    private$data[[".coord_x"]] <- NULL
    private$data[[".coord_y"]] <- NULL

    # Process the data (find all the different coordinates
    # across the different replicates, and also merge the new data to the old data)
    private$data <- process_data_add_obs(PtE, new_data = data, private$data,
                                        group_vector)

    ## convert to Spoints and add
    PtE <- self$get_PtE()
    spatial_points <- self$coordinates(PtE = PtE, normalized = TRUE)
    private$data[[".coord_x"]] <- rep(spatial_points[,1], times = n_group)
    private$data[[".coord_y"]] <- rep(spatial_points[,2], times = n_group)
    if(tibble){
      private$data <- tidyr::as_tibble(private$data)
    }
    private$group_col <- group
    class(private$data) <- c("metric_graph_data", class(private$data))
    })
          if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
          }
  },


  #' @description Use `dplyr::mutate` function on the internal metric graph data object.
  #' @param ... Arguments to be passed to `dplyr::mutate()`.
   #' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
  #' @details A wrapper to use `dplyr::mutate()` within the internal metric graph data object.
  #' @return A `tidyr::tibble` object containing the resulting data list after the mutate.
  mutate = function(..., .drop_na = FALSE, .drop_all_na = TRUE) {
    if(!inherits(private$data, "tbl_df")){
      data_res <- tidyr::as_tibble(private$data)
    } else{
      data_res <- private$data
    }

    if(.drop_all_na){
      is_tbl <- inherits(data_res, "tbl_df")
        idx_temp <- idx_not_all_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
        if(is_tbl){
          data_res <- tidyr::as_tibble(data_res)
        }
    }

    if(.drop_na){
      if(!inherits(data_res, "tbl_df")){
        idx_temp <- idx_not_any_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
      } else{
        data_res <- tidyr::drop_na(data_res)
      }
    }


    data_res <- dplyr::mutate(.data = data_res, ...)

    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }

    return(data_res)
  },


  #' @description Use `tidyr::drop_na()` function on the internal metric graph data object.
  #' @param ... Arguments to be passed to `tidyr::drop_na()`.
  #' @details A wrapper to use `dplyr::drop_na()` within the internal metric graph data object.
  #' @return A `tidyr::tibble` object containing the resulting data list after the drop_na.
  drop_na = function(...) {
    if(!inherits(private$data, "tbl_df")){
      data_res <- tidyr::as_tibble(private$data)
    } else{
      data_res <- private$data
    }

    data_res <- tidyr::drop_na(data = data_res, ...)

    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }

    return(data_res)
  },



  #' @description Use `dplyr::select` function on the internal metric graph data object.
  #' @param ... Arguments to be passed to `dplyr::select()`.
 #' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
  #' @details A wrapper to use `dplyr::select()` within the internal metric graph data object. Observe that it is a bit different from directly using `dplyr::select()` since it does not allow to remove the internal positions that are needed for the metric_graph methods to work.
  #' @return A `tidyr::tibble` object containing the resulting data list after the selection.
  select = function(..., .drop_na = FALSE, .drop_all_na = TRUE) {
    if(!inherits(private$data, "tbl_df")){
      data_res <- tidyr::as_tibble(private$data)
    } else{
      data_res <- private$data
    }
    data_res <- dplyr::select(.data = data_res, ...)
    data_res[[".group"]] <- private$data[[".group"]]
    data_res[[".edge_number"]] <- private$data[[".edge_number"]]
    data_res[[".distance_on_edge"]] <- private$data[[".distance_on_edge"]]
    data_res[[".coord_x"]] <- private$data[[".coord_x"]]
    data_res[[".coord_y"]] <- private$data[[".coord_y"]]


    if(.drop_all_na){
      is_tbl <- inherits(data_res, "tbl_df")
        idx_temp <- idx_not_all_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
        if(is_tbl){
          data_res <- tidyr::as_tibble(data_res)
        }
    }

    if(.drop_na){
      if(!inherits(data_res, "tbl_df")){
        idx_temp <- idx_not_any_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
      } else{
        data_res <- tidyr::drop_na(data_res)
      }
    }

    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }

    return(data_res)
  },

    #' @description Use `dplyr::filter` function on the internal metric graph data object.
  #' @param ... Arguments to be passed to `dplyr::filter()`.
   #' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
  #' @details A wrapper to use `dplyr::filter()` within the internal metric graph data object.
  #' @return A `tidyr::tibble` object containing the resulting data list after the filter.
  filter = function(..., .drop_na = FALSE, .drop_all_na = TRUE) {
    if(!inherits(private$data, "tbl_df")){
      data_res <- tidyr::as_tibble(private$data)
    } else{
      data_res <- private$data
    }


    if(.drop_all_na){
      is_tbl <- inherits(data_res, "tbl_df")
        idx_temp <- idx_not_all_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
        if(is_tbl){
          data_res <- tidyr::as_tibble(data_res)
        }
    }

    if(.drop_na){
      if(!inherits(data_res, "tbl_df")){
        idx_temp <- idx_not_any_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
      } else{
        data_res <- tidyr::drop_na(data_res)
      }
    }

    data_res <- dplyr::filter(.data = data_res, ...)
    data_res <- dplyr::arrange(.data = data_res, `.group`, `.edge_number`, `.distance_on_edge`)

    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }

    return(data_res)
  },


  #' @description Use `dplyr::summarise` function on the internal metric graph data object grouped by the spatial locations and the internal group variable.
  #' @param ... Arguments to be passed to `dplyr::summarise()`.
  #' @param .include_graph_groups Should the internal graph groups be included in the grouping variables? The default is `FALSE`. This means that, when summarising, the data will be grouped by the internal group variable together with the spatial locations.
  #' @param .groups A vector of strings containing the names of the columns to be additionally grouped, when computing the summaries. The default is `NULL`.
   #' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
  #' @details A wrapper to use `dplyr::summarise()` within the internal metric graph data object grouped by manually inserted groups (optional), the internal group variable (optional) and the spatial locations. Observe that if the integral group variable was not used as a grouping variable for the summarise, a new column, called `.group`, will be added, with the same value 1 for all rows.
  #' @return A `tidyr::tibble` object containing the resulting data list after the summarise.
  summarise = function(..., .include_graph_groups = FALSE, .groups = NULL, .drop_na = FALSE, .drop_all_na = TRUE) {
    if(!inherits(private$data, "tbl_df")){
      data_res <- tidyr::as_tibble(private$data)
    } else{
      data_res <- private$data
    }


    if(.drop_all_na){
      is_tbl <- inherits(data_res, "tbl_df")
        idx_temp <- idx_not_all_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
        if(is_tbl){
          data_res <- tidyr::as_tibble(data_res)
        }
    }

    if(.drop_na){
      if(!inherits(data_res, "tbl_df")){
        idx_temp <- idx_not_any_NA(data_res)
        data_res <- lapply(data_res, function(dat){dat[idx_temp]})
      } else{
        data_res <- tidyr::drop_na(data_res)
      }
    }


    group_vars <- c(".edge_number", ".distance_on_edge", ".coord_x", ".coord_y")
    if(.include_graph_groups){
      group_vars <- c(".group", group_vars)
    }
    group_vars <- c(.groups, group_vars)
    data_res <- dplyr::group_by_at(.tbl = data_res, .vars = group_vars)
    data_res <- dplyr::summarise(.data = data_res, ...)
    data_res <- dplyr::ungroup(data_res)
    if(is.null(data_res[[".group"]])){
      data_res[[".group"]] <- 1
    }

    data_res <- dplyr::arrange(.data = data_res, `.group`, `.edge_number`, `.distance_on_edge`)

    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }

    return(data_res)

  },

 #' @description Return the internal data with the option to filter by groups.
 #' @param group A vector contaning which groups should be returned? The default is `NULL`, which gives the result for the all groups.
 #' @param tibble Should the data be returned as a `tidyr::tibble`?
 #' @param drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.

  get_data = function(group = NULL, tibble = TRUE, drop_na = FALSE, drop_all_na = TRUE){
    if(is.null(private$data)){
      stop("The graph does not contain data.")
    }
    if(!is.null(group)){
      total_groups <- self$get_groups()
      if(!(all(group %in% total_groups))){
        if(all(group%%1 == 0)){
          group <- total_groups[group]
        } else{
          stop("At least one entry of 'group' is a not an existing group or an existing index.")
        }
      }

      data_temp <- select_group(private$data, group)
    } else{
      data_temp <- private$data
    }
    if(tibble){
      data_temp <- tidyr::as_tibble(data_temp)
    }

    if(drop_all_na){
      is_tbl <- inherits(data_temp, "tbl_df")
        idx_temp <- idx_not_all_NA(data_temp)
        data_temp <- lapply(data_temp, function(dat){dat[idx_temp]})
        if(is_tbl){
          data_temp <- tidyr::as_tibble(data_temp)
        }
    }

    if(drop_na){
      if(!inherits(data_temp, "tbl_df")){
        idx_temp <- idx_not_any_NA(data_temp)
        data_temp <- lapply(data_temp, function(dat){dat[idx_temp]})
      } else{
        data_temp <- tidyr::drop_na(data_temp)
      }
    }

    if(!inherits(data_temp, "metric_graph_data")){
      class(data_temp) <- c("metric_graph_data", class(data_temp))
    }
    return(data_temp)
  },

  #' @description Build Kirchoff constraint matrix from edges.
  #' @param alpha the type of constraint (currently only supports 2)
  #' @param edge_constraint if TRUE, add constraints on vertices of degree 1
  #' @details Currently not implemented for circles (edges that start and end
  #' in the same vertex)
  #' @return No return value. Called for its side effects.
  buildC = function(alpha = 2, edge_constraint = FALSE) {

    if(alpha==2){
      i_  =  rep(0, 2 * self$nE)
      j_  =  rep(0, 2 * self$nE)
      x_  =  rep(0, 2 * self$nE)

      count_constraint <- 0
      count <- 0
      for (v in 1:self$nV) {
        lower_edges  <- which(self$E[, 1] %in% v)
        upper_edges  <- which(self$E[, 2] %in% v)
        n_e <- length(lower_edges) + length(upper_edges)

        #derivative constraint
        if ((edge_constraint & n_e ==1) | n_e > 1) {
          i_[count + 1:n_e] <- count_constraint + 1
          j_[count + 1:n_e] <- c(4 * (lower_edges-1) + 2, 4 * (upper_edges-1) + 4)
          x_[count + 1:n_e] <- c(rep(1,length(lower_edges)),
                                 rep(-1,length(upper_edges)))
          count <- count + n_e
          count_constraint <- count_constraint + 1
        }
        if (n_e > 1) {
          if (length(upper_edges) == 0) {
            edges <- cbind(lower_edges, 1)
          } else if(length(lower_edges) == 0){
            edges <- cbind(upper_edges, 3)
          }else{
            edges <- rbind(cbind(lower_edges, 1),
                           cbind(upper_edges, 3))
          }
          for (i in 2:n_e) {
            i_[count + 1:2] <- count_constraint + 1
            j_[count + 1:2] <- c(4 * (edges[i-1,1] - 1) + edges[i-1, 2],
                                 4 * (edges[i,1]   - 1) + edges[i,   2])
            x_[count + 1:2] <- c(1,-1)
            count <- count + 2
            count_constraint <- count_constraint + 1
          }
        }
      }
      C <- Matrix::sparseMatrix(i = i_[1:count],
                                j = j_[1:count],
                                x = x_[1:count],
                                dims = c(count_constraint, 4*self$nE))
      self$C = C

      self$CoB <- c_basis2(self$C)
      self$CoB$T <- t(self$CoB$T)
    }else{
      error("only alpha=2 implemented")
    }
  },

  #' @description Builds mesh object for graph.
  #' @param h Maximum distance between mesh nodes (should be provided if n is
  #' not provided).
  #' @param n Maximum number of nodes per edge (should be provided if h is not
  #' provided).
  #' @param  continuous If `TRUE` (default), the mesh contains only one node per vertex.
  #' If `FALSE`, each vertex v is split into deg(v) disconnected nodes to allow
  #' for the creation of discontinuities at the vertices.
  #' @param continuous.outs If `continuous = FALSE` and `continuous.outs = TRUE`, continuity is
  #' assumed for the outgoing edges from each vertex.
  #' @param continuous.deg2 If `TRUE`, continuity is assumed at degree 2 vertices.
  #' @details The mesh is a list with the objects:
  #' - `PtE` The mesh locations excluding the original vertices;
  #' - `V` The verties of the mesh;
  #' - `E` The edges of the mesh;
  #' - `n_e` The number of vertices in the mesh per original edge in the graph;
  #' - `h_e` The mesh width per edge in the graph;
  #' - `ind` The indices of the vertices in the mesh;
  #' - `VtE` All mesh locations including the original vertices.
  #' @return No return value. Called for its side effects. The mesh is stored in
  #' the `mesh` element of the `metric_graph` object.
 build_mesh = function(h=NULL, n=NULL, continuous = TRUE,
                       continuous.outs = FALSE, continuous.deg2 = FALSE) {

   if(is.null(h) && is.null(n)){
     stop("You should specify either h or n!")
   }

   if(!is.null(h)){
     if(length(h)>1 || (!is.numeric(h))){
       stop("h should be a single number")
     }

     if(h<=0){
       stop("h must be positive!")
     }
   }

   if(!is.null(n)){
     if(length(n)>1 || (!is.numeric(n))){
       stop("n should be a single number")
     }

     if(n<=0){
       stop("n must be positive!")
     }

     if(n%%1!=0){
       warning("A noninteger n was given, we are rounding it to an integer.")
       n <- round(n)
     }
   }

   if(continuous) {
     self$mesh <- list(PtE = NULL,
                       V = NULL,
                       E = NULL,
                       n_e = rep(0, self$nV),
                       h_e = NULL,
                       ind = 1:self$nV,
                       VtE = NULL)
     attr(self$mesh, 'continuous') <- TRUE
     self$mesh$V <- self$V

     for (i in 1:length(self$edges)) {
       if (is.null(n)) {
         #remove boundary points
         self$mesh$n_e[i] <- ceiling(self$edge_lengths[i] / h) + 1 - 2
       } else {
         self$mesh$n_e[i] <- n
       }
       if (self$mesh$n_e[i] > 0) {
         d.e <- seq(from = 0, to = 1, length.out = self$mesh$n_e[i] + 2)
         d.e <- d.e[2:(1+self$mesh$n_e[i])]

         self$mesh$PtE <- rbind(self$mesh$PtE, cbind(rep(i, self$mesh$n_e[i]),
                                                     d.e))

         self$mesh$h_e <- c(self$mesh$h_e,
                            rep(self$edge_lengths[i] * d.e[1],
                                self$mesh$n_e[i] + 1))

         V.int <- (max(self$mesh$ind) + 1):(max(self$mesh$ind) + self$mesh$n_e[i])
         self$mesh$ind <- c(self$mesh$ind, V.int)
         self$mesh$E <- rbind(self$mesh$E, cbind(c(self$E[i, 1], V.int),
                                                 c(V.int, self$E[i, 2])))
       } else {
         self$mesh$E <- rbind(self$mesh$E, self$E[i, ])
         self$mesh$h_e <- c(self$mesh$h_e,self$edge_lengths[i])
       }
     }

     self$mesh$VtE <- rbind(self$VtEfirst(), self$mesh$PtE)
     if(!is.null(self$mesh$PtE)) {
       self$mesh$V <- rbind(self$mesh$V, self$coordinates(PtE = self$mesh$PtE))
     } else {
       self$mesh$V <- rbind(self$mesh$V)
     }
   } else {
     self$mesh <- list(PtE = NULL,
                       V = NULL,
                       E = NULL,
                       n_e = NULL,
                       h_e = NULL,
                       ind = 0,
                       VtE = NULL)
     attr(self$mesh, 'continuous') <- FALSE
     for (i in 1:length(self$edges)) {
       if (is.null(n)) {
         #remove boundary points
         self$mesh$n_e[i] <- ceiling(self$edge_lengths[i] / h) + 1
       } else {
         self$mesh$n_e[i] <- n + 2
       }
       if (self$mesh$n_e[i] > 0) {
         d.e <- seq(from = 0, to = 1, length.out = self$mesh$n_e[i])

         self$mesh$PtE <- rbind(self$mesh$PtE, cbind(rep(i, self$mesh$n_e[i]),
                                                     d.e))

         self$mesh$h_e <- c(self$mesh$h_e,
                            rep(self$edge_lengths[i] * d.e[2],
                                self$mesh$n_e[i] - 1))

         V.int <- (max(self$mesh$ind) + 1):(max(self$mesh$ind) + self$mesh$n_e[i])
         self$mesh$ind <- c(self$mesh$ind, V.int)

         self$mesh$E <- rbind(self$mesh$E, cbind(V.int[-length(V.int)], V.int[-1]))

       } else {
         self$mesh$E <- rbind(self$mesh$E, self$E[i, ])
         self$mesh$h_e <- c(self$mesh$h_e,self$edge_lengths[i])
       }
     }
     self$mesh$VtE <- self$mesh$PtE
     self$mesh$V <- self$coordinates(PtE = self$mesh$PtE)
     if(continuous.outs) {
       private$mesh_merge_outs()
     }
     private$move_V_first()
     if(continuous.deg2) {
       private$mesh_merge_deg2()
     }

   }
 },

  #' @description Build mass and stiffness matrices for given mesh object.
  #' @details The function builds: The matrix `C` which is the mass matrix with
  #' elements \eqn{C_{ij} = <\phi_i, \phi_j>}, the matrix `G` which is the stiffness
  #' matrix with elements \eqn{G_{ij} = <d\phi_i, d\phi_j>}, the matrix `B` with
  #' elements \eqn{B_{ij} = <d\phi_i, \phi_j>}, the matrix `D` with elements
  #' \eqn{D_{ij} = \sum_{v\in V}\phi_i(v)\phi_j(v)}, and the vector with weights
  #' \eqn{<\phi_i, 1>}.
  #' @param petrov Compute Petrov-Galerkin matrices? (default `FALSE`). These
  #' are defined as \eqn{Cpet_{ij} = <\phi_i, \psi_j>} and \eqn{Gpet_{ij} = <d\phi_i, \psi_j>},
  #' where \eqn{\psi_{i}} are piecewise constant basis functions on the edges of
  #' the mesh.
  #' @return No return value. Called for its side effects. The finite element
  #' matrices `C`, `G` and `B` are stored in the `mesh` element in the
  #' `metric_graph` object. If `petrov=TRUE`, the corresponding Petrov-Galerkin
  #' matrices are stored in `Cpet` and `Gpet`.
  compute_fem = function(petrov = FALSE) {
    if (is.null(self$mesh)) {
      stop("no mesh provided")
    }
    nV <- dim(self$mesh$V)[1]
    fem_temp <- assemble_fem(E = self$mesh$E, h_e = self$mesh$h_e, nV = nV, petrov = petrov)
    self$mesh$C <- fem_temp$C
    self$mesh$G <- fem_temp$G
    self$mesh$B <- fem_temp$B
    self$mesh$D <- Diagonal(dim(self$mesh$C)[1],
                            c(rep(1, self$nV), rep(0, dim(self$mesh$C)[1] - self$nV)))
    #set weighted Krichhoff matrix
    self$mesh$K <- Diagonal(dim(self$mesh$C)[1],
                            c(rep(0, self$nV), rep(0, dim(self$mesh$C)[1] - self$nV)))

    if(!all(self$get_edge_weights()==1)){
      for(i in 1:self$nV) {
        if(attr(self$vertices[[i]],"degree") > 1) {
          edges.i <- which(rowSums(self$E==i)>0)
          edges.mesh <- which(rowSums(self$mesh$E==i)>0)
          w <- rep(0,length(edges.mesh))
          h <- rep(0,length(edges.mesh))
          for(j in 1:length(edges.mesh)) {
            V.e <- self$mesh$E[edges.mesh[j],which(self$mesh$E[edges.mesh[j],]!=i)]
            E.e <- self$mesh$VtE[V.e,1] #the edge the mesh node is on
            w[j] <- attr(self$edges[[E.e]],"weight")
            h[j] <- self$mesh$h_e[edges.mesh[j]]
          }
          for(j in 2:attr(self$vertices[[i]],"degree")){
            self$mesh$K[i,i] <- self$mesh$K[i,i] +  (w[j]/w[1] - 1)/h[j]
          }
        }
      }
    }

    if(petrov) {
      self$mesh$Cpet <- fem_temp$Cpet
      self$mesh$Gpet <- fem_temp$Gpet
      private$set_petrov_matrices()
    }

    self$mesh$weights <- rowSums(self$mesh$C)
  },

  #' @description Deprecated - Computes observation matrix for mesh.
  #'
  #'  `r lifecycle::badge("deprecated")` in favour of `metric_graph$fem_basis()`.
  #' @param PtE Locations given as (edge number in graph, normalized location on
  #' edge)
  #' @details For n locations and a mesh with m nodes, `A` is an n x m matrix with
  #' elements \eqn{A_{ij} = \phi_j(s_i)}{A_{ij} = \phi_j(s_i)}.
  #' @return The observation matrix.
  #'

  mesh_A = function(PtE){
    lifecycle::deprecate_warn("1.2.0", "metric_graph$mesh_A()", "metric_graph$fem_basis()")
    self$fem_basis(PtE)
  },

  #' @description Computes observation matrix for mesh.
  #' @param PtE Locations given as (edge number in graph, normalized location on
  #' edge)
  #' @details For n locations and a mesh with m nodes, `A` is an n x m matrix with
  #' elements \eqn{A_{ij} = \phi_j(s_i)}{A_{ij} = \phi_j(s_i)}.
  #' @return The observation matrix.
  fem_basis = function(PtE) {
    if(ncol(PtE)!= 2){
      stop("PtE must have two columns!")
    }

    if (min(PtE[,2]) < 0) {
      stop("PtE[, 2] has negative values")
    }
    if ((max(PtE[,2]) > 1)) {
      stop("For normalized distances, the values in PtE[, 2] should not be
             larger than 1")
    }
    if (is.null(self$mesh)) {
      stop("no mesh given")
    }

    x <- private$PtE_to_mesh(PtE)
    n <- dim(x)[1]

    A <- sparseMatrix(i = c(1:n, 1:n),
                      j = c(self$mesh$E[x[, 1], 1], self$mesh$E[x[, 1], 2]),
                      x = c(1 - x[, 2], x[, 2]),
                      dims = c(n, dim(self$mesh$V)[1]))
    return(A)
  },

  #' @description Find one edge corresponding to each vertex.
  #' @return A nV x 2 matrix the first element of the `i`th row is the edge
  #' number corresponding to the `i`th vertex and the second value is 0
  #' if the vertex is at the start of the edge and 1 if the vertex
  #' is at the end of the edge.
  VtEfirst = function() {
    n.V <- dim(self$V)[1]
    VtE <- matrix(0, n.V, 2)

    for (i in 1:n.V) {
      Ei <- which(self$E[, 1] == i)[1]
      pos <- 0
      if (is.na(Ei) == 1) {
        pos <- 1
        Ei <- which(self$E[, 2] == i)[1]
      }
      VtE[i,] <- c(Ei, pos)
    }
    return(VtE)
  },

  #' @description Plots the metric graph.
  #' @param data Which column of the data to plot? If `NULL`, no data will be
  #' plotted.
  #' @param newdata A dataset of class `metric_graph_data`, obtained by any `get_data()`, `mutate()`, `filter()`, `summarise()`, `drop_na()` methods of metric graphs, see the vignette on data manipulation for more details.
  #' @param group If there are groups, which group to plot? If `group` is a
  #' number, it will be the index of the group as stored internally. If `group`
  #' is a character, then the group will be chosen by its name.
  #' @param plotly Use plot_ly for 3D plot (default `FALSE`). This option
  #' requires the 'plotly' package.
  #' @param vertex_size Size of the vertices.
  #' @param vertex_color Color of vertices.
  #' @param edge_width Line width for edges.
  #' @param edge_color Color of edges.
  #' @param data_size Size of markers for data.
  #' @param support_width For 3D plot, width of support lines.
  #' @param support_color For 3D plot, color of support lines.
  #' @param mesh Plot the mesh locations?
  #' @param X Additional values to plot.
  #' @param X_loc Locations of the additional values in the format
  #' (edge, normalized distance on edge).
  #' @param p Existing objects obtained from 'ggplot2' or 'plotly' to add the graph to
  #' @param degree Show the degrees of the vertices?
  #' @param direction Show the direction of the edges?
##  # ' @param mutate A string containing the commands to be passed to `dplyr::mutate` function in order to obtain new variables as functions of the existing variables.
##  # ' @param filter A string containing the commands to be passed to `dplyr::filter` function in order to obtain new filtered data frame.
##  # ' @param summarise A string containing the commands to be passed to `dplyr::summarise` function in order to obtain new  data frame containing the summarised variable.
##  # ' @param summarise_group_by A vector of strings containing the names of the columns to be additionally grouped, when computing the summaries. The default is `NULL`.
##  # ' @param summarise_by_graph_group Should the internal graph groups be included in the grouping variables? The default is `FALSE`. This means that, when summarising, the data will be grouped by the internal group variable together with the spatial locations.
  #' @param ... Additional arguments to pass to `ggplot()` or `plot_ly()`
  #' @return A `plot_ly` (if `plotly = TRUE`) or `ggplot` object.
  plot = function(data = NULL,
                  newdata = NULL,
                  group = 1,
                  plotly = FALSE,
                  vertex_size = 3,
                  vertex_color = 'black',
                  edge_width = 0.3,
                  edge_color = 'black',
                  data_size = 1,
                  support_width = 0.5,
                  support_color = "gray",
                  mesh = FALSE,
                  X = NULL,
                  X_loc = NULL,
                  p = NULL,
                  degree = FALSE,
                  direction = FALSE,
                  # mutate = NULL,
                  # filter = NULL,
                  # summarise = NULL,
                  # summarise_group_by = NULL,
                  # summarise_by_graph_group = FALSE,
                  ...) {
    if(!is.null(data) && is.null(private$data)) {
      stop("The graph does not contain data.")
    }
    if(is.numeric(group) && !is.null(data)) {
      unique_group <- unique(private$data[[".group"]])
      group <- unique_group[group]
    }
    if(!is.null(newdata)){
      if(!inherits(newdata, "metric_graph_data")){
        stop("'newdata' must be of class 'metric_graph_data'!")
      }
    }

    if(!plotly) {
      p <- private$plot_2d(line_width = edge_width,
                           marker_size = vertex_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           newdata = newdata,
                           data_size = data_size,
                           group = group,
                           mesh = mesh,
                           X = X,
                           X_loc = X_loc,
                           p = p,
                           degree = degree,
                           direction = direction,
                           ...)
      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degrees"){
          p <- p + labs(x = "Longitude",  y = "Latitude")
        } else{
          p <- p + labs(x = paste0("x (in ",private$vertex_unit, ")"),  y = paste0("y (in ",private$vertex_unit, ")"))
        }
      }
    } else {
      requireNamespace("plotly")
      p <- private$plot_3d(line_width = edge_width,
                           marker_size = vertex_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           newdata = newdata,
                           data_size = data_size,
                           group = group,
                           mesh = mesh,
                           X = X,
                           X_loc = X_loc,
                           support_color = support_color,
                           support_width = support_width,
                           p = p,
                           ...)
      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degrees"){
          p <- plotly::layout(p, scene = list(xaxis = list(title = "Longitude"), yaxis = list(title = "Latitude")))
        } else{
          p <- plotly::layout(p, scene = list(xaxis = list(title = paste0("x (in ",private$vertex_unit, ")")), yaxis = list(title = paste0("y (in ",private$vertex_unit, ")"))))
        }
      }
    }
    return(p)
  },

  #' @description Plots the connections in the graph
  #' @return No return value. Called for its side effects.
  plot_connections = function(){
        g <- graph(edges = c(t(self$E)), directed = FALSE)
        plot(g)
  },

  #' @description Checks if the graph is a tree (without considering directions)
  #' @return TRUE if the graph is a tree and FALSE otherwise.
  is_tree = function(){
        g <- graph(edges = c(t(self$E)), directed = FALSE)
        return(igraph::is_tree(g, mode = "all"))
  },

  #' @description Plots continuous function on the graph.
  #' @param data Which column of the data to plot? If `NULL`, no data will be
  #' plotted.
  #' @param newdata A dataset of class `metric_graph_data`, obtained by any `get_data()`, `mutate()`, `filter()`, `summarise()`, `drop_na()` methods of metric graphs, see the vignette on data manipulation for more details.
  #' @param group If there are groups, which group to plot? If `group` is a
  #' number, it will be the index of the group as stored internally. If `group`
  #' is a character, then the group will be chosen by its name.
  #' @param X A vector with values for the function
  #' evaluated at the mesh in the graph
  #' @param plotly If `TRUE`, then the plot is shown in 3D. This option requires
  #' the package 'plotly'.
  #' @param improve_plot Should the original edge coordinates be added to the data with linearly interpolated values to improve the plot?
  #' @param continuous Should continuity be assumed when the plot uses `newdata`?
  #' @param vertex_size Size of the vertices.
  #' @param vertex_color Color of vertices.
  #' @param edge_width Width for edges.
  #' @param edge_color For 3D plot, color of edges.
  #' @param line_width For 3D plot, line width of the function curve.
  #' @param line_color Color of the function curve.
  #' @param support_width For 3D plot, width of support lines.
  #' @param support_color For 3D plot, color of support lines.
  #' @param p Previous plot to which the new plot should be added.
  #' @param ... Additional arguments for `ggplot()` or `plot_ly()`
  #' @return Either a `ggplot` (if `plotly = FALSE`) or a `plot_ly` object.
  plot_function = function(data = NULL,
                           newdata = NULL,
                           group = 1,
                           X = NULL,
                           plotly = FALSE,
                           improve_plot = FALSE,
                           continuous = TRUE,
                           vertex_size = 5,
                           vertex_color = "black",
                           edge_width = 1,
                           edge_color = 'black',
                           line_width = NULL,
                           line_color = 'rgb(0,0,200)',
                           support_width = 0.5,
                           support_color = "gray",
                           p = NULL,
                           ...){
    if (is.null(line_width)) {
      line_width = edge_width
    }

    mesh <- FALSE

    if(is.null(data) && is.null(X)){
      stop("You should provide either 'data' or 'X'.")
    }

    if(!is.null(data) && !is.null(X)){
      warning("Both 'data' and 'X' were provided. Only 'data' will be considered.")
      X <- NULL
    }

    if(!is.character(data) && (is.vector(data) || !is.null(dim(data)))){
      X <- data
    }

    if(!is.null(X)){
      mesh <- TRUE
      if(!is.vector(X) && is.null(dim(X))){
        stop("'X' should be a vector, or a row-matrix or a column-matrix!")
      }
      if(!is.null(dim(X)) && min(dim(X)) > 1){
        stop("If 'X' is a matrix, it needs to be either a row matrix or a column matrix!")
      }
      X <- as.vector(X)
    }

    if (mesh) {
      if (is.null(self$mesh)) {
        stop("X is a vector but no mesh provided")
      }

      if (is.null(self$mesh$PtE)) {
        PtE_dim <- 0
      } else {
        PtE_dim <- dim(self$mesh$PtE)[1]
      }

      if (length(X) == PtE_dim && attr(self$mesh, "continuous")) {
        X <- c(rep(NA, dim(self$V)[1]), X)
      }


      if (length(X) != dim(unique(self$mesh$V))[1] && length(X) != dim(self$mesh$V)[1]) {
        stop(paste0("X does not have the correct size (the possible sizes are ",
                    PtE_dim,", ", dim(self$mesh$V)[1], " and ",
                    dim(unique(self$mesh$V))[1],")"))
      }

      if(dim(unique(self$mesh$V))[1] != dim(self$mesh$V)[1]){
        if(length(X) == dim(unique(self$mesh$V))[1]){
          X_temp <-  rep(NA, dim(self$mesh$V)[1])
          X_temp[which(!duplicated(self$mesh$V))] <- X
          X <- X_temp
        }
      }

      n.v <- dim(self$V)[1]
      XV <- X[1:n.v]
    } else{
      if(is.null(newdata)){
        X <- self$get_data(group = group)
        X <- X[,c(".edge_number", ".distance_on_edge", data)]
      } else{
        if(!inherits(newdata, "metric_graph_data")){
          stop("'newdata' must be of class 'metric_graph_data'!")
        }
        X <- newdata[,c(".edge_number", ".distance_on_edge", data)]
      }
    }


      if(improve_plot){
        if(is.null(attr(self$edges[[1]], "PtE"))){
          self$compute_PtE_edges()
        }
        PtE_edges <- lapply(1:length(self$edges), function(i){attr(self$edges[[i]], "PtE")})
      }


    x.loc <- y.loc <- z.loc <- i.loc <- NULL
    kk = 1
    for (i in 1:self$nE) {
      Vs <- self$E[i, 1]
      Ve <- self$E[i, 2]
      if (mesh) {
        ind <- self$mesh$PtE[, 1] == i

        if (sum(ind)==0) {
          vals <- rbind(c(0, XV[Vs]),
                        c(1, XV[Ve]))

        } else {
          if(attr(self$mesh,"continuous")) {
            vals <- rbind(c(0, XV[Vs]),
                          cbind(self$mesh$PtE[ind, 2], X[n.v + which(ind)]),
                          c(1, XV[Ve]))
          } else {
            if(min(self$mesh$PtE[ind,2])==0 && max(self$mesh$PtE[ind,2])==1) {
              vals <- cbind(self$mesh$PtE[ind, 2], X[which(ind)])
            } else if (min(self$mesh$PtE[ind,2])==0) {
                vals <- rbind(cbind(self$mesh$PtE[ind, 2], X[which(ind)]),
                              c(1, XV[Ve]))
            } else if (max(self$mesh$PtE[ind,2])==1){
              vals <- rbind(c(0, XV[Vs]),
                            cbind(self$mesh$PtE[ind, 2], X[which(ind)]))
            } else {
              vals <- rbind(c(0, XV[Vs]),
                            cbind(self$mesh$PtE[ind, 2], X[which(ind)]),
                            c(1, XV[Ve]))
            }
          }


        }

        if(improve_plot){
          PtE_tmp <- PtE_edges[[i]]
          PtE_tmp <- PtE_tmp[PtE_tmp[,1] == i,, drop=FALSE]
          PtE_tmp <- PtE_tmp[,2, drop=TRUE]
          PtE_tmp <- setdiff(PtE_tmp, vals[,1])
          if(length(PtE_tmp)>0){
                PtE_tmp <- cbind(PtE_tmp, NA)
                vals <- rbind(vals,PtE_tmp)
          }
            ord_idx <- order(vals[,1])
            vals <- vals[ord_idx,]
            if(vals[1,1] > 0){
              vals <- rbind(c(0,NA), vals)
            }
            if(vals[nrow(vals),1] < 1){
              vals <- rbind(vals, c(1,NA))
            }
            max_val <- max(vals[,2], na.rm=TRUE)
            min_val <- min(vals[,2], na.rm=TRUE)
            vals[,2] <- na.const(pmax(pmin(object = zoo::na.approx(object = vals[,2],
                                                  x = vals[,1],
                                                      na.rm=FALSE, ties = "mean"),
                                               max_val), min_val))
            vals <- vals[(vals[,1] >= 0) & (vals[,1]<=1),]
        }


      } else {
        X <- as.data.frame(X)
        vals <- X[X[, 1]==i, 2:3, drop = FALSE]

    if(continuous){
      if(nrow(vals)>0){

        if(!improve_plot){
            if (max(vals[, 1]) < 1) {
              # #check if we can add end value from other edge
              Ei <- self$E[, 1] == Ve #edges that start in Ve
              Ei <- which(Ei)
              if (sum(Ei) > 0) {
                ind <- which(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE] == 0)
                if(sum(ind)>0){
                  ind <- which.min(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE])
                  min.val <- X[X[,1,drop=TRUE] %in% Ei, 3,drop=TRUE][ind]
                } else {
                ind <- NULL
                ind.val <- which.min(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE])
                min.val <- X[X[,1,drop=TRUE] %in% Ei, 3,drop=TRUE][ind.val]
              }} else{
                ind <- NULL
                ind.val <- integer(0)
              }
              if (length(ind) > 0) {
                # vals <- rbind(vals, c(1, X[ind, 3,drop=TRUE]))
                vals <- rbind(vals, c(1, min.val))
              }
              else {
                Ei <- self$E[, 2] == Ve #edges that end in Ve
                Ei <- which(Ei)
                if (sum(Ei)  > 0) {
                  ind <- which(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE] == 1)
                  if(sum(ind)>0){
                    ind <- which.max(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE])
                    max.val <- X[X[,1,drop=TRUE] %in% Ei, 3,drop=TRUE][ind]
                  } else {
                  ind.val.max <- which.max(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE])
                  max.val <- X[X[,1,drop=TRUE] %in% Ei, 3,drop=TRUE][ind.val.max]
                  if(length(ind.val) == 0){
                    ind <- ind.val.max
                  } else if (length(ind.val.max) == 0){
                    ind <- ind.val
                  } else{
                    ind <- ifelse(1-max.val < min.val, ind.val.max, ind.val)
                  }
                } } else{
                  if(length(ind.val)>0){
                    ind <- ind.val
                  } else{
                    ind <- NULL
                  }
                }
                if (length(ind) > 0){
                  # vals <- rbind(vals, c(1, X[ind, 3, drop=TRUE]))
                  vals <- rbind(vals, c(1, max.val))
                }
              }
            }

            if (min(vals[, 1] > 0)) {
              #check if we can add start value from other edge
              Ei <- self$E[, 1] == Vs #edges that start in Vs
              Ei <- which(Ei)
              if (sum(Ei) > 0) {
                  ind <- which(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE] == 0)
                if(sum(ind)>0){
                    ind <- ind[1]
                  } else {
                ind <- NULL
                ind.val <- which.min(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE])
                min.val <- X[ind.val, 2,drop=TRUE]
              }} else{
                ind <- NULL
                ind.val <- integer(0)
              }
              if (length(ind) > 0) {
                vals <- rbind(c(0, X[ind, 3, drop=TRUE]), vals)
              } else {
                Ei <- self$E[, 2] == Vs #edges that end in Vs
                Ei <- which(Ei)
                if (sum(Ei) > 0) {
                  ind <- which(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE] == 1)
                  if(sum(ind)>0){
                    ind <- ind[1]
                  } else {
                  ind.val.max <- which.max(X[X[,1,drop=TRUE] %in% Ei, 2,drop=TRUE])
                  max.val <- X[ind.val.max, 2,drop=TRUE]
                  if(length(ind.val) == 0){
                    ind <- ind.val.max
                  } else if (length(ind.val.max) == 0){
                    ind <- ind.val
                  } else{
                    ind <- ifelse(1-max.val < min.val, ind.val.max, ind.val)
                  }
                } } else{
                  if(length(ind.val)>0){
                    ind <- ind.val
                  } else{
                    ind <- NULL
                  }
                }
                if (length(ind) > 0) {
                  vals <- rbind(c(0, X[ind, 3, drop=TRUE]), vals)
                }
              }
            }
        } else {
            PtE_tmp <- PtE_edges[[i]]
            if(PtE_tmp[1,1] != i){
              edge_new <- PtE_tmp[1,1]
              idx_new <- which(X[,1] == edge_new)
              new_val <- X[idx_new, 2:3, drop=FALSE]
              if(nrow(new_val)>0){
                sum_fact <- ifelse(max(vals[,1]==1),  -1e-6,  min(new_val[,1]))
                sub_fact <- ifelse(min(vals[,1]==0), 1+1e-6,  max(new_val[,1]))
                pos_edge <- which(self$E[edge_new,] == self$E[i,1])
                if(pos_edge == 2){
                  new_val[,1] <- new_val[,1] - sub_fact
                } else {
                  new_val[,1] <- -new_val[,1] + sum_fact
                }
                vals <- rbind(vals, new_val)
              }
            } else{
              if(any(self$E[,2] == self$E[i,1])){
                edge_new <- which(self$E[,2] == self$E[i,1])[1]
                idx_new <- which(X[,1] == edge_new)
                new_val <- X[idx_new, 2:3, drop=FALSE]
                if(nrow(new_val)>0){
                  sub_fact <- ifelse(min(vals[,1]==0), 1+1e-6,  max(new_val[,1]))
                  new_val[,1] <- new_val[,1] - sub_fact
                  vals <- rbind(vals, new_val)
                }
              } else if(any(self$E[-i,1] == self$E[i,1])){
                edge_new <- which(self$E[-i,1] == self$E[i,1])[1]
                idx_new <- which(X[,1] == edge_new)
                new_val <- X[idx_new, 2:3, drop=FALSE]
                if(nrow(new_val)>0){
                  sum_fact <- ifelse(max(vals[,1]==1),  -1e-6,  min(new_val[,1]))
                  new_val[,1] <- -new_val[,1] + sum_fact
                  vals <- rbind(vals, new_val)
                }
              }
            }

            if (PtE_tmp[nrow(PtE_tmp), 1]!=i){
              edge_new <- PtE_tmp[nrow(PtE_tmp), 1]
              idx_new <- which(X[,1] == edge_new)
              new_val <- X[idx_new, 2:3, drop=FALSE]
              if(nrow(new_val)>0){
                sub_fact <- ifelse(min(vals[,1]==0), 1+1e-6,  -max(new_val[,1]) - 1)
                sum_fact <- ifelse(max(vals[,1]==1),  1+1e-6,  1-min(new_val[,1]))
                pos_edge <- which(self$E[edge_new,] == self$E[i,2])
                if(pos_edge == 2){
                  new_val[,1] <- -new_val[,1] - sub_fact
                } else {
                  new_val[,1] <- new_val[,1] + sum_fact
                }
                vals <- rbind(vals, new_val)
              }
            } else {
                if (any(self$E[,1] == self$E[i,2])){
                  edge_new <- which(self$E[,1] == self$E[i,2])[1]
                  idx_new <- which(X[,1] == edge_new)
                  new_val <- X[idx_new, 2:3, drop=FALSE]
                  if(nrow(new_val)>0){
                    sum_fact <- ifelse(max(vals[,1]==1),  1+1e-6,  1-min(new_val[,1]))
                    new_val[,1] <- new_val[,1] + sum_fact
                    vals <- rbind(vals, new_val)
                  }
              } else if (any(self$E[,2] == self$E[i,2])){
                  edge_new <- which(self$E[,2] == self$E[i,2])[1]
                  idx_new <- which(X[,1] == edge_new)
                  new_val <- X[idx_new, 2:3, drop=FALSE]
                  if(nrow(new_val)>0){
                    sub_fact <- ifelse(min(vals[,1]==0), 1+1e-6,  -max(new_val[,1]) - 1)
                    new_val[,1] <- -new_val[,1] - sum_fact
                    vals <- rbind(vals, new_val)
                  }
              }
            }

            if(length(PtE_tmp[,1]==i) > 0){
              PtE_tmp <- PtE_tmp[PtE_tmp[,1]==i, ,drop=FALSE]
              PtE_tmp <- PtE_tmp[,2,drop = TRUE]
              PtE_tmp <- setdiff(PtE_tmp, vals[,1])
              if(length(PtE_tmp)>0){
                PtE_tmp <- cbind(PtE_tmp, NA)
                PtE_tmp <- as.data.frame(PtE_tmp)
                colnames(PtE_tmp) <- c(".distance_on_edge", data)
                vals <- rbind(vals,PtE_tmp)
              }
            }

            ord_idx <- order(vals[,1])
            vals <- vals[ord_idx,]
            max_val <- max(vals[,2], na.rm=TRUE)
            min_val <- min(vals[,2], na.rm=TRUE)
            vals[,2] <- na.const(pmax(pmin(object = zoo::na.approx(object = vals[,2],
                                                  x = vals[,1],
                                                      na.rm=FALSE, ties = "mean"),
                                               max_val), min_val))
            vals <- vals[(vals[,1] >= 0) & (vals[,1]<=1),]
        }
      }
      } else if(improve_plot){
       
          PtE_tmp <- PtE_edges[[i]]
          PtE_tmp <- PtE_tmp[PtE_tmp[,1] == i,, drop=FALSE]
          PtE_tmp <- PtE_tmp[,2, drop=TRUE]
          PtE_tmp <- setdiff(PtE_tmp, vals[,1])
          if(length(PtE_tmp)>0){
                PtE_tmp <- cbind(PtE_tmp, NA)
                PtE_tmp <- as.data.frame(PtE_tmp)
                colnames(PtE_tmp) <- c(".distance_on_edge", data)
                vals <- rbind(vals,PtE_tmp)
          }
            ord_idx <- order(vals[,1])
            vals <- vals[ord_idx,]
            if(vals[1,1] > 0){
              vals <- rbind(c(0,NA), vals)
            }
            if(vals[nrow(vals),1] < 1){
              vals <- rbind(vals, c(1,NA))
            }
            max_val <- max(vals[,2], na.rm=TRUE)
            min_val <- min(vals[,2], na.rm=TRUE)
            vals[,2] <- na.const(pmax(pmin(object = zoo::na.approx(object = vals[,2],
                                                  x = vals[,1],
                                                      na.rm=FALSE, ties = "mean"),
                                               max_val), min_val))
            vals <- vals[(vals[,1] >= 0) & (vals[,1]<=1),]

      }
      }

        data.to.plot <- vals

        data.to.plot.order <- data.to.plot[order(vals[, 1,drop=TRUE]), ,
                                           drop = FALSE]

        coords <- interpolate2(self$edges[[i]],
                               pos = data.to.plot.order[, 1, drop = TRUE],
                               normalized = TRUE)

        x.loc <- c(x.loc, coords[, 1])
        y.loc <- c(y.loc, coords[, 2])
        z.loc <- c(z.loc, data.to.plot.order[, 2, drop=TRUE])
        i.loc <- c(i.loc, rep(kk, length(coords[, 1])))
        kk = kk+1
    }

    data <- data.frame(x = x.loc, y = y.loc, z = z.loc, i = i.loc)

    if(plotly){
      requireNamespace("plotly")
      if(is.null(p)){
        p <- self$plot(plotly = TRUE,
                       vertex_color = vertex_color,
                       vertex_size = vertex_size,
                       edge_width = edge_width,
                       edge_color = edge_color)
      }
      p <- plotly::add_trace(p, data = data, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = line_color),
                           split = ~i, showlegend = FALSE, ...)
      if(support_width > 0) {
        data.mesh <- data.frame(x = c(x.loc, x.loc), y = c(y.loc, y.loc),
                                z = c(rep(0, length(z.loc)), z.loc),
                                i = rep(1:length(z.loc), 2))
        p <- plotly::add_trace(p, data = data.mesh, x = ~y, y = ~x, z = ~z,
                             mode = "lines", type = "scatter3d",
                             line = list(width = support_width,
                                         color = support_color),
                             split = ~i, showlegend = FALSE)
      }

      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degrees"){
          p <- plotly::layout(p, scene = list(xaxis = list(title = "Longitude"), yaxis = list(title = "Latitude")))
        } else{
          p <- plotly::layout(p, scene = list(xaxis = list(title = paste0("x (in ",private$vertex_unit, ")")), yaxis = list(title = paste0("y (in ",private$vertex_unit, ")"))))
        }
      }

    } else {
      if(is.null(p)) {
          p <- ggplot(data = data, aes(x = x, y = y,
                                     group = i,
                                     colour = z)) +
          geom_path(linewidth = line_width) + scale_color_viridis() +
          labs(colour = "")
      } else {
        p <- p + geom_path(data = data,
                           aes(x = x, y = y,
                               group = i, colour = z),
                           linewidth = line_width) +
          scale_color_viridis() + labs(colour = "")
      }
          p <- self$plot(edge_width = 0, vertex_size = vertex_size,
                     vertex_color = vertex_color, p = p)
      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degrees"){
          p <- p + labs(x = "Longitude",  y = "Latitude")
        } else{
          p <- p + labs(x = paste0("x (in ",private$vertex_unit, ")"),  y = paste0("y (in ",private$vertex_unit, ")"))
        }
      }
    }
    return(p)
  },

  #' @description Plots a movie of a continuous function evolving on the graph.
  #' @param X A m x T matrix where the ith column represents the function at the
  #' ith time, evaluated at the mesh locations.
  #' @param plotly If `TRUE`, then plot is shown in 3D. This option requires the
  #' package 'plotly'.
  #' @param vertex_size Size of the vertices.
  #' @param vertex_color Color of vertices.
  #' @param edge_width Width for edges.
  #' @param edge_color For 3D plot, color of edges.
  #' @param line_width For 3D plot, line width of the function curve.
  #' @param line_color Color of the function curve.
  #' @param ... Additional arguments for ggplot or plot_ly.
  #' @return Either a `ggplot` (if `plotly=FALSE`) or a `plot_ly` object.
  plot_movie = function(X,
                        plotly = TRUE,
                        vertex_size = 5,
                        vertex_color = "black",
                        edge_width = 1,
                        edge_color = 'black',
                        line_width = NULL,
                        line_color = 'rgb(0,0,200)',
                           ...){
    if (is.null(line_width)) {
      line_width = edge_width
    }


    if (is.null(self$mesh)) {
      stop("X is a vector but no mesh provided")
    }

    if (is.null(self$mesh$PtE)) {
      PtE_dim <- 0
    } else {
      PtE_dim <- dim(self$mesh$PtE)[1]
    }

    if (length(X) == PtE_dim) {
      X <- c(rep(NA, dim(self$V)[1]), X)
    }

    if (dim(X)[1] != dim(unique(self$mesh$V))[1]) {
      stop(paste0("X does not have the correct size (the possible sizes are ", PtE_dim, " and ", dim(unique(self$mesh$V))[1],")"))
    }

    n.v <- dim(self$V)[1]
    XV <- X[1:n.v,]

    x.loc <- y.loc <- z.loc <- i.loc <- f.loc <-  NULL
    kk = 1
    frames <- dim(X)[2]
    for (i in 1:self$nE) {
      Vs <- self$E[i, 1]
      Ve <- self$E[i, 2]

      ind <- self$mesh$PtE[, 1] == i

      if (sum(ind)==0) {
        vals <- rbind(c(0, XV[Vs,]),
                      c(1, XV[Ve,]))

      } else {
        vals <- rbind(c(0, XV[Vs,]),
                      cbind(self$mesh$PtE[ind, 2], X[n.v + which(ind),]),
                      c(1, XV[Ve,]))

      }

        data.to.plot <- vals
        data.to.plot.order <- data.to.plot[order(data.to.plot[, 1]), ,
                                           drop = FALSE]

        coords <- interpolate2(self$edges[[i]],
                               pos = data.to.plot.order[, 1, drop = TRUE],
                               normalized = TRUE)
        x.loc <- c(x.loc, rep(coords[, 1], frames))
        y.loc <- c(y.loc, rep(coords[, 2], frames))
        z.loc <- c(z.loc, c(data.to.plot.order[, 2:(frames+1)]))
        i.loc <- c(i.loc, rep(rep(kk, length(coords[, 1])), frames))
        f.loc <- c(f.loc, rep(1:frames, each = length(coords[, 1])))
        kk = kk+1
    }
    data <- data.frame(x = x.loc, y = y.loc, z = z.loc, i = i.loc, f = f.loc)

    if(plotly){
      requireNamespace("plotly")
      x <- y <- ei <- NULL
      for (i in 1:self$nE) {
        xi <- self$edges[[i]][, 1]
        yi <- self$edges[[i]][, 2]
        ii <- rep(i,length(xi))
        x <- c(x, xi)
        y <- c(y, yi)
        ei <- c(ei, ii)
      }
      frames <- dim(X)[2]
      data.graph <- data.frame(x = rep(x, frames),
                               y = rep(y, frames),
                               z = rep(rep(0,length(x)), frames),
                               i = rep(ei, frames),
                               f = rep(1:frames, each = length(x)))

      p <- plotly::plot_ly(data=data.graph, x = ~y, y = ~x, z = ~z,frame = ~f)
      p <- plotly::add_trace(p, data = data.graph, x = ~y, y = ~x, z = ~z,
                             frame = ~f,
                             mode = "lines", type = "scatter3d",
                             line = list(width = line_width,
                                         color = edge_color),
                             split = ~i, showlegend = FALSE)

      p <- plotly::add_trace(p, data = data, x = ~y, y = ~x, z = ~z,
                             frame = ~f,
                             mode = "lines", type = "scatter3d",
                             line = list(width = line_width,
                                         color = line_color),
                             split = ~i, showlegend = FALSE, ...)

      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degrees"){
          p <- plotly::layout(p, scene = list(xaxis = list(title = "Longitude"), yaxis = list(title = "Latitude")))
        } else{
          p <- plotly::layout(p, scene = list(xaxis = list(title = paste0("x (in ",private$vertex_unit, ")")), yaxis = list(title = paste0("y (in ",private$vertex_unit, ")"))))
        }
      }

    } else {
      stop("not implemented")
    }
    return(p)
  },

  #' @description Add observations on mesh to the object.
  #' @param data A `data.frame` or named list containing the observations.
  #' In case of groups, the data.frames for the groups should be stacked vertically,
  #' with a column indicating the index of the group. If `data_frame` is not
  #' `NULL`, it takes priority over any eventual data in `Spoints`.
  #' @param group If the data_frame contains groups, one must provide the column
  #' in which the group indices are stored.
  #' @return No return value. Called for its side effects. The observations are
  #' stored in the `data` element in the `metric_graph` object.
  add_mesh_observations = function(data = NULL, group = NULL) {

    if(is.null(self$mesh)){
      stop("You should have a mesh!")
    }
    PtE_mesh <- self$mesh$PtE
    data[[".edge_number"]] <- PtE_mesh[,1]
    data[[".distance_on_edge"]] <- PtE_mesh[,2]
    self$add_observations(data = data, group = group,
                          edge_number = ".edge_number",
                          distance_on_edge = ".distance_on_edge",
                          normalized = TRUE)
  },

  #' @description Returns a copy of the initial metric graph.
  #' @return A `metric_graph` object.
  get_initial_graph = function() {
    tmp_graph <- private$initial_graph$clone()
    if(private$pruned){
      tmp_graph$prune_vertices()
   }
    return(tmp_graph)
  },

  #' @description Convert between locations on the graph and Euclidean
  #' coordinates.
  #' @param PtE Matrix with locations on the graph (edge number and normalized
  #' position on the edge).
  #' @param XY Matrix with locations in Euclidean space
  #' @param normalized If `TRUE`, it is assumed that the positions in `PtE` are
  #' normalized to (0,1), and the object returned if `XY` is specified contains
  #' normalized locations.
  #' @return If `PtE` is specified, then a matrix with Euclidean coordinates of
  #' the locations is returned. If `XY` is provided, then a matrix with the
  #' closest locations on the graph is returned.
  coordinates = function(PtE = NULL, XY = NULL, normalized = TRUE) {
    if(is.null(PtE) && is.null(XY)) {
      stop("PtE or XY must be provided")
    } else if(!is.null(PtE) && !is.null(XY)) {
      stop("Either PtE or XY must be provided, not both")
    }
    x <- y <- NULL
    if (!is.null(PtE)) {
      if(is.vector(PtE)) {
        if(length(PtE) != 2)
          stop("PtE is a vector but does not have length 2")

        PtE <- matrix(PtE,1,2)
      }
      if(ncol(PtE)!= 2){
        stop("PtE must have two columns!")
      }

      if (min(PtE[,2]) < 0) {
        stop("PtE[, 2] has negative values")
      }
      if ((max(PtE[,2]) > 1) && normalized) {
        stop("For normalized distances, the values in PtE[, 2] should not be
             larger than 1")
      }
      if(max(PtE[,2] - self$edge_lengths[PtE[, 1]]) > 0 && !normalized) {
        stop("PtE[, 2] contains values which are larger than the edge lengths")
      }

      if(!normalized) {
        PtE <- cbind(PtE[, 1], PtE[, 2] / self$edge_lengths[PtE[, 1]])
      }

      Points <- matrix(NA, nrow=nrow(PtE), ncol=ncol(PtE))

      for (i in 1:dim(PtE)[1]) {
        Points[i,] <- interpolate2(self$edges[[PtE[i, 1]]] ,
                                   pos = PtE[i, 2], normalized = TRUE)
      }
      return(Points)
    } else {
      SP <- snapPointsToLines(XY, self$edges, longlat = private$longlat, crs = private$crs)
      # coords.old <- XY
      # colnames(coords.old) <- paste(colnames(coords.old), '_old', sep="")
      XY = t(SP[["coords"]])
      PtE <- cbind(match(SP[["df"]][["nearest_line_index"]], 1:length(self$edges)), 0)

      for (ind in unique(PtE[, 1])) {
        index.p <- PtE[, 1] == ind
        PtE[index.p,2]=projectVecLine2(self$edges[[ind]], XY[index.p, , drop=FALSE],
                                       normalized=TRUE)

      }

      if (!normalized) {
        PtE <- cbind(PtE[, 1], PtE[, 2] * self$edge_lengths[PtE[, 1]])
      }
      return(PtE)
    }
  }
  ),

  private = list(
  #function for creating Vertex and Edges from self$edges
  line_to_vertex = function(tolerance = 0, longlat = FALSE, fact, verbose, crs, proj4string, which_longlat, length_unit, vertex_unit, project, which_projection, project_data) {



    if(project_data && longlat){
     if(verbose){
      message("Projecting edges")
      bar_edges_proj <- msg_progress_bar(length(self$edges))
    }
      private$vertex_unit <- length_unit
      if(which_longlat == "sf"){
        if(which_projection == "Robinson"){
          str_proj <- "+proj=robin +datum=WGS84 +no_defs +over"
        } else if (which_projection == "Winkel tripel") {
          str_proj <- "+proj=wintri +datum=WGS84 +no_defs +over"
        } else{
          str_proj <- which_projection
        }
        fact <- process_factor_unit("m", length_unit)
        for(i in 1:length(self$edges)){
          if(verbose){
            bar_edges_proj$increment()
          }
          sf_points <- sf::st_as_sf(as.data.frame(self$edges[[i]]), coords = 1:2, crs = crs)
          sf_points_eucl <- sf::st_transform(sf_points, crs=sf::st_crs(str_proj))
          self$edges[[i]] <- sf::st_coordinates(sf_points_eucl) * fact
        }
      } else{
        if(which_projection == "Robinson"){
          str_proj <- "+proj=robin +datum=WGS84 +no_defs +over"
        } else if (which_projection == "Winkel tripel") {
          str_proj <- "+proj=wintri +datum=WGS84 +no_defs +over"
        } else{
          str_proj <- which_projection
        }
        fact <- process_factor_unit("km", length_unit)
        for(i in 1:length(self$edges)){
          if(verbose){
            bar_edges_proj$increment()
          }
          sp_points <- sp::SpatialPoints(coords = self$edges[[i]], proj4string = proj4string)
          sp_points_eucl <- sp::spTransform(sp_points,CRSobj=sp::CRS(str_proj))
          self$edges[[i]] <- sp::coordinates(sp_points_eucl) * fact
        }
      }
      private$longlat <- FALSE
      private$crs <- NULL
      private$proj4string <- NULL
    }

    if(verbose){
      message("Part 1/2")
      bar_line_vertex <- msg_progress_bar(length(self$edges))
    }


    lines <- matrix(nrow = 2*length(self$edges), ncol = 3)
    for(i in 1:length(self$edges)){
      if(verbose){
        bar_line_vertex$increment()
      }
      points <- self$edges[[i]]
      n <- dim(points)[1]
      lines[2*i-1,] <- c(i, points[1,])
      lines[2*i, ] <- c(i, points[n,])
    }

    #save all vertices that are more than tolerance distance apart
    if(verbose){
      message("Computing auxiliary distances")
    }

    if(!project_data || !longlat){
    if(!project || !longlat){
        fact <- process_factor_unit(vertex_unit, length_unit)
          dists <- compute_aux_distances(lines = lines[,2:3,drop=FALSE], crs = crs, longlat = longlat, proj4string = proj4string, fact = fact, which_longlat = which_longlat, length_unit = private$length_unit)
    } else if (which_longlat == "sf"){
        if(which_projection == "Robinson"){
          str_proj <- "+proj=robin +datum=WGS84 +no_defs +over"
        } else if (which_projection == "Winkel tripel") {
          str_proj <- "+proj=wintri +datum=WGS84 +no_defs +over"
        } else{
          str_proj <- which_projection
        }
        sf_points <- sf::st_as_sf(as.data.frame(lines), coords = 2:3, crs = crs)
        sf_points_eucl <- sf::st_transform(sf_points, crs=sf::st_crs(str_proj))
        fact <- process_factor_unit("m", length_unit)
        dists <- dist(sf::st_coordinates(sf_points_eucl)) * fact
    } else{
        if(which_projection == "Robinson"){
          str_proj <- "+proj=robin +datum=WGS84 +no_defs +over"
        } else if (which_projection == "Winkel tripel") {
          str_proj <- "+proj=wintri +datum=WGS84 +no_defs +over"
        } else{
          str_proj <- which_projection
        }
        sp_points <- sp::SpatialPoints(coords = lines[,2:3], proj4string = proj4string)
        sp_points_eucl <- sp::spTransform(sp_points,CRSobj=sp::CRS(str_proj))
        fact <- process_factor_unit("m", length_unit)
        dists <- dist(sp_points_eucl@coords) * fact
    }
    } else{
      fact <- process_factor_unit(vertex_unit, length_unit)
      dists <- dist(lines[,2:3,drop=FALSE]) * fact
    }


    if(verbose){
      message("Done!")
    }


    if(!inherits(dists,"dist")){
        idx_keep <- sapply(1:nrow(lines), function(i){ifelse(i==1,TRUE,all(dists[i, 1:(i-1)] > tolerance))})
      vertex <- lines[idx_keep,, drop=FALSE]
    } else{
        idx_keep <- sapply(1:nrow(lines), function(i){ifelse(i==1,TRUE,all(dists[ nrow(lines)*(1:(i-1)-1) - (1:(i-1))*(1:(i-1) -1)/2 + i -1:(i-1)] > tolerance))})
      vertex <- lines[idx_keep,, drop=FALSE]
    }
      # if(inherits(dists,"dist")) dists <- dist2mat(dists,256)
      # idx_keep <- sapply(1:nrow(lines), function(i){ifelse(i==1,TRUE,all(dists[i, 1:(i-1)] > tolerance))})
      # vertex <- lines[idx_keep,]

    if(verbose){
      message("Part 2/2")
      bar_line_vertex <- msg_progress_bar(max(lines[, 1]))
    }

    lvl <- matrix(0, nrow = max(lines[,1]), 4)
    k=1
    lines_keep_id <- NULL

    for (i in 1:max(lines[, 1])) {

      if(verbose){
        bar_line_vertex$increment()
      }

      which.line <- sort(which(lines[, 1] == i))
      line <- lines[which.line, , drop=FALSE]

      #index of vertex corresponding to the start of the line
      ind1 <- which.min((vertex[, 2] - line[1, 2])^2 +
                          (vertex[, 3] - line[1, 3])^2)
      #index of vertex corresponding to the end of the line
      ind2 <- which.min((vertex[, 2] - line[2, 2])^2 +
                          (vertex[, 3] - line[2, 3])^2)
      if(length(ind1)>0){
        self$edges[[i]][1,] <- vertex[ind1, 2:3, drop=FALSE]
        i.e <- dim(self$edges[[i]])[1]
      } else{
        i.e <- 1
      }
      if(length(ind2)>0){
        self$edges[[i]][i.e,] <- vertex[ind2, 2:3, drop=FALSE]
      }
      ll <- compute_line_lengths(self$edges[[i]], longlat = longlat, unit = length_unit, crs = crs, proj4string, which_longlat, vertex_unit, project_data)
      if(ll > tolerance) {
        lvl[k,] <- c(i, ind1, ind2, ll)
        k=k+1
        lines_keep_id <- c(lines_keep_id, i)
      }
    }

    lvl <- lvl[1:(k-1),,drop = FALSE]
    self$edges <- self$edges[lines_keep_id]
    self$V <- vertex[, 2:3, drop = FALSE]
    self$E <- lvl[, 2:3, drop = FALSE]
    self$edge_lengths <- lvl[,4]
    # units(self$edge_lengths) <- length_unit
    self$nV <- dim(self$V)[1]
    self$nE <- dim(self$E)[1]
  },

  #Compute PtE for mesh given PtE for graph
  PtE_to_mesh = function(PtE){
    Vertexes <- self$VtEfirst()
    VtE <- rbind(Vertexes, self$mesh$PtE)
    PtE_update <- matrix(0, dim(PtE)[1], 2)
    for (i in 1:dim(PtE)[1]) {
      ei <- PtE[i, 1] #extract edge

      #extract distances of all mesh nodes on edge including end points
      ind.nodes <- which(self$mesh$PtE[,1] == ei)
      dist.nodes <- self$mesh$PtE[ind.nodes,2]
      if(length(ind.nodes)>0) {
        ind <- c(self$E[ei,1], ind.nodes + self$nV, self$E[ei,2])
        dists <- c(0,dist.nodes,1)
      } else {
        ind <- c(self$E[ei,1], self$E[ei,2])
        dists <- c(0,1)
      }

      ind2 <- sort(sort(abs(dists - PtE[i, 2]), index.return = TRUE)$ix[1:2])
      v1 <- ind[ind2[1]] #vertex "before" the point
      v2 <- ind[ind2[2]] #vertex "after" the point
      d1 <- dists[ind2[1]] #distance on the edge of the point before
      d2 <- dists[ind2[2]] #distance on the edge of the point after

      #find the "edge" in the mesh on which the point is
      e <- which(rowSums((self$mesh$E == v1) + (self$mesh$E == v2)) == 2)

      # Handle the case of multiple edges
      # In the case the edge lengths are different, we can identify
      # In the case they are equal, we cannot identify, but it does not matter, as there is no difference then.
      if(length(e)>1){
       ind <- which(self$edge_lengths[ei] == self$mesh$h_e[e])
       e <- e[ind]
       e <- e[1]
      }

    for(ei in e){
      if (self$mesh$E[ei, 1] == v1) { #edge starts in the vertex before
        d <- (PtE[i, 2] - d1)/(d2 - d1)
      } else {
        d <- 1- (PtE[i, 2] - d1)/(d2 - d1)
      }
        PtE_update[i, ] <- c(e, d)
      }
    }
    return(PtE_update)
  },

  ## Adding column_y argument which tells which column of y to get

  plot_2d = function(line_width = 0.1,
                     marker_size = 1,
                     vertex_color = 'black',
                     edge_color = 'black',
                     data,
                     newdata,
                     data_size = 1,
                     group = 1,
                     mesh = FALSE,
                     X = NULL,
                     X_loc = NULL,
                     p = NULL,
                     degree = FALSE,
                     direction = FALSE,
                     ...){
    xyl <- c()

    nc <- do.call(rbind,lapply(self$edges, function(x) dim(x)[1]))
    xyl <- cbind(do.call(rbind,self$edges), rep(1:length(nc), times = nc))

    if(is.null(p)){
        p <- ggplot() + geom_path(data = data.frame(x = xyl[, 1],
                                                    y = xyl[, 2],
                                                    grp = xyl[, 3]),
                                  mapping = aes(x = x, y = y, group = grp),
                                  linewidth = line_width,
                                  colour = edge_color,
                                  ...)
    } else {
        p <- p + geom_path(data = data.frame(x = xyl[, 1],
                                             y = xyl[,2],
                                             grp = xyl[,3]),
                           mapping = aes(x = x, y = y, group = grp),
                           linewidth = line_width,
                           colour = edge_color, ...)
    }
    if(direction) {
      mid.l <- self$coordinates(PtE = cbind(1:self$nE, rep(0.49,self$nE)))
      mid.u <- self$coordinates(PtE = cbind(1:self$nE, rep(0.5,self$nE)))
      p <- p + geom_path(data = data.frame(x = c(mid.l[,1], mid.u[,1]),
                                           y = c(mid.l[, 2], mid.u[,2]),
                                           edge = c(1:self$nE,1:self$nE)),
                          mapping = aes(x = x, y = y, group = edge),
                         arrow = ggplot2::arrow(),
                          size= marker_size/2, ...)
    }
    if (marker_size > 0) {
      if(degree) {
        x <- self$V[,1]
        y <- self$V[,2]
        degrees <- self$get_degrees()
        p <- p + geom_point(data = data.frame(x = self$V[, 1],
                                              y = self$V[, 2],
                                              degree = degrees),
                            mapping = aes(x, y, colour = factor(degree)),
                            size= marker_size, ...) +
    scale_color_viridis(discrete = TRUE, guide_legend(title = ""))
      } else if (direction) {
        degrees <- self$get_degrees()
        start.deg <- end.deg <- rep(0,self$nV)
        for(i in 1:self$nV) {
          start.deg[i] <- sum(self$E[,1]==i)
          end.deg[i] <- sum(self$E[,2]==i)
        }
        problematic <- (degrees > 1) & (start.deg == 0 | end.deg == 0)
        p <- p + geom_point(data = data.frame(x = self$V[problematic, 1],
                                              y = self$V[problematic, 2]),
                            mapping = aes(x, y),
                            colour = "red",
                            size= marker_size, ...)
        p <- p + geom_point(data = data.frame(x = self$V[!problematic, 1],
                                              y = self$V[!problematic, 2]),
                            mapping = aes(x, y),
                            colour = "green",
                            size= marker_size, ...)

      } else {
        p <- p + geom_point(data = data.frame(x = self$V[, 1],
                                              y = self$V[, 2]),
                            mapping = aes(x, y),
                            colour = vertex_color,
                            size= marker_size, ...)
      }
    }
    if (!is.null(data)) {
      x <- y <- NULL

      if(!is.null(newdata)){
        data_group <- select_group(newdata, group)
      } else{
        data_group <- select_group(private$data, group)
      }

      if(!(data%in%names(data_group))){
        stop(paste(data,"is not an existing column name of the dataset."))
      }

      y_plot <-data_group[[data]]

      x <- data_group[[".coord_x"]]
      y <- data_group[[".coord_y"]]

      p <- p + geom_point(data = data.frame(x = x[!is.na(as.vector(y_plot))],
                                            y = y[!is.na(as.vector(y_plot))],
                                            val = as.vector(y_plot[!is.na(as.vector(y_plot))])),
                          mapping = aes(x, y, color = val),
                          size = data_size, ...) +
        scale_colour_gradientn(colours = viridis(100), guide_legend(title = ""))

    }
    if (mesh) {
      p <- p + geom_point(data=data.frame(x = self$mesh$V[, 1],
                                          y = self$mesh$V[, 2]),
                          mapping = aes(x, y), size = marker_size * 0.5,
                          pch = 21,
                          colour = "black",
                          fill = "gray")
    }
    if(!is.null(X)){
      if(is.null(X_loc)){
        stop("X supplied but not X_loc")
      }
      x <- y <- NULL
      if(length(X) != nrow(X_loc)){
        stop("The number of observations does not match the number of observations!")
      }
      points_xy <- self$coordinates(PtE = X_loc)
      x <- points_xy[,1]
      y <- points_xy[,2]
      p <- p + geom_point(data = data.frame(x = x, y = y,
                                            val = as.vector(X)),
                          mapping = aes(x, y, color = val),
                          size = data_size) +
        scale_color_viridis() + labs(colour = "")
    }
    p <- p + coord_fixed()
    return(p)
  },

  ## Adding column_y argument which tells which column of y to get

  plot_3d = function(line_width = 1,
                     marker_size = 1,
                     vertex_color = 'rgb(0,0,0)',
                     edge_color = 'rgb(0,0,0)',
                     data,
                     newdata,
                     data_size = 1,
                     group = 1,
                     mesh = FALSE,
                     p = NULL,
                     support_width = 0.5,
                    support_color = "gray",
                     ...){
      x <- y <- ei <- NULL
      for (i in 1:self$nE) {
        xi <- self$edges[[i]][, 1]
        yi <- self$edges[[i]][, 2]
        ii <- rep(i,length(xi))
        x <- c(x, xi)
        y <- c(y, yi)
        ei <- c(ei, ii)
      }
      data.plot <- data.frame(x = x, y = y, z = rep(0,length(x)), i = ei)

    if(is.null(p)) {
      p <- plotly::plot_ly(data=data.plot, x = ~y, y = ~x, z = ~z,...)
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = edge_color),
                           split = ~i, showlegend = FALSE)
    } else {
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = edge_color),
                           split = ~i, showlegend = FALSE)
    }


    if(marker_size > 0) {
      data.plot2 <- data.frame(x = self$V[, 1], y = self$V[, 2],
                               z = rep(0, self$nV))
      p <- plotly::add_trace(p, data = data.plot2, x = ~y, y = ~x, z = ~z,
                           type = "scatter3d", mode = "markers",
                           marker = list(size = marker_size,
                                         color = vertex_color))
    }

    if (!is.null(data)) {
      x <- y <- NULL
      if(!is.null(newdata)){
        data_group <- select_group(newdata, group)
      } else{
        data_group <- select_group(private$data, group)
      }
      y_plot <- data_group[[data]]

      x <- data_group[[".coord_x"]]
      y <- data_group[[".coord_y"]]

      data.plot <- data.frame(x = x[!is.na(as.vector(y_plot))],
                              y = y[!is.na(as.vector(y_plot))],
                              z = rep(0,length(x[!is.na(as.vector(y_plot))])),
                              val = as.vector(y_plot[!is.na(as.vector(y_plot))]),
                              i =rep(1:length(y_plot), 2))

      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~val,
                           type = "scatter3d", mode = "markers",
                           marker = list(size = marker_size,
                                         color = ~val,
                                         colorbar=list(title='', len = 0.5),
                                         colorscale='Viridis'),
                           showlegend=FALSE)
      if(support_width > 0) {
          data.support <- data.frame(x = c(data.plot$x, data.plot$x), y = c(data.plot$y, data.plot$y),
                                z = c(rep(0, length(data.plot$z)), data.plot$val),
                                i = rep(1:length(data.plot$val), 2))
        p <- plotly::add_trace(p, data = data.support, x = ~y, y = ~x, z = ~z,
                             mode = "lines", type = "scatter3d",
                             line = list(width = support_width,
                                         color = support_color),
                             split = ~i, showlegend = FALSE)
      }
    }
    if (mesh) {
      data.plot <- data.frame(x = self$mesh$V[, 1],
                              y = self$mesh$V[, 2],
                              z = rep(0, dim(self$mesh$V)[1]))
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           type = "scatter3d", mode = "markers",
                           marker = list(size = marker_size/2,
                                         color = 'rgb(100,100,100)'),
                           showlegend = FALSE)
    }
    return(p)
  },

  ## Coordinates function to return all the lines intersecting within a tolerance

    coordinates_multiple_snaps = function(XY, tolerance, verbose = verbose, crs, proj4string, longlat, fact, which_longlat) {

      coords_line <- c()
      coords_tmp <- c()

      if(!is.null(XY)){
        class(XY) <- setdiff(class(XY), "metric_graph_edge")
      }


      if(!private$longlat){
      lines_sf <- sf::st_sfc(lapply(self$edges, function(i){sf::st_linestring(i)}))
          points_sf <- sf::st_as_sf(as.data.frame(XY), coords = 1:2)
      crs <- NULL
    } else if (private$which_longlat == "sf"){
      points_sf <- sf::st_as_sf(as.data.frame(XY), coords = 1:2, crs = private$crs)
      lines_sf <- sf::st_sfc(lapply(self$edges, function(i){sf::st_linestring(i)}), crs = private$crs)
    } else{
      lines_sf <- sf::st_sfc(lapply(self$edges, function(i){sf::st_linestring(i)}), crs = sf::st_crs(private$proj4string))
      crs <- sf::st_crs(private$proj4string)
      points_sf <- sf::st_as_sf(as.data.frame(XY), coords = 1:2, crs = private$crs)
    }

      if(verbose){
        message("Computing auxiliary distances")
      }

      within_dist <- t(as.matrix(sf::st_is_within_distance(points_sf, lines_sf, dist = tolerance)))

      if(verbose){
        message("Done!")
      }

      if(verbose){
        message("Snapping vertices")
        bar_multiple_snaps <- msg_progress_bar(length(self$edges))
      }

      for(i in 1:length(self$edges)){
        if(verbose){
          bar_multiple_snaps$increment()
        }
        select_points <- matrix(XY[within_dist[i,],], ncol=2)
        if(nrow(select_points) > 0){
          SP <- snapPointsToLines(select_points, self$edges[i], longlat, crs, i)
          idx_tol <- (SP[["df"]][["snap_dist"]] <= tolerance)
          coords_line <- c(coords_line, SP[["df"]][["nearest_line_index"]][idx_tol])
          coords_tmp <- rbind(coords_tmp, (t(SP[["coords"]]))[idx_tol,])
        }
      }

      XY <- coords_tmp

      PtE = cbind(match(coords_line, 1:length(self$edges)), 0)

      for (ind in unique(PtE[, 1])) {
        if(verbose){
          bar_multiple_snaps$increment()
        }
        index.p <- PtE[, 1] == ind

        PtE[index.p,2] <- projectVecLine2(self$edges[[ind]], XY[index.p,,drop=FALSE],
                                       normalized=TRUE)
      }
      return(PtE)
      },

  #utility function to merge close vertices
  merge_close_vertices = function(tolerance, fact) {
    if(tolerance > 0) {
      dists <- compute_aux_distances(lines = self$V, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit)
      v.merge <- NULL
      k <- 0
      for (i in 2:self$nV) {
            if(!inherits(dists,"dist")){
                i.min <- which.min(dists[i, 1:(i-1)])
                cond_tmp <- (dists[i, i.min] < tolerance)
            } else{
                i.min <- which.min(dists[ nrow(self$V)*(1:(i-1)-1) - (1:(i-1))*(1:(i-1) -1)/2 + i -1:(i-1)])
                cond_tmp <- (dists[ nrow(self$V)*(i.min-1) - (i.min)*(i.min -1)/2 + i -i.min] < tolerance)
            }

        # i.min <- which.min(dists[i, 1:(i-1)])
        if (cond_tmp) {
          k <- k + 1
          v.merge <- rbind(v.merge, sort(c(i, i.min)))
        }
      }
      if(k>0){
        for( j in 1:k) {
          v.keep <- v.merge[1,1]
          v.rem <- v.merge[1,2]
          if(j < k) {
            v.merge <- v.merge[-1,,drop = FALSE]
            v.merge[v.merge == v.rem] <- v.keep
            v.merge[v.merge > v.rem] <- v.merge[v.merge > v.rem] - 1
          }

          idx_E1 <- which(self$E[,1] == v.rem)
          idx_E2 <- which(self$E[,2] == v.rem)

          for(k1 in idx_E1){
            self$edges[[k1]][1,] <- self$V[v.keep,]
          }
          for(k2 in idx_E2){
            self$edges[[k2]][nrow(self$edges[[k2]]),] <- self$V[v.keep,]
          }

          self$V <- self$V[-v.rem,]
          self$nV <- self$nV - 1
          self$E[self$E == v.rem] <- v.keep
          self$E[self$E > v.rem] <- self$E[self$E > v.rem] - 1
        }
      }
    }
  },

  # utility function to remove small circles
  remove_circles = function(threshold, verbose,longlat, unit, crs, proj4string, which_longlat, vertex_unit, project_data) {
    if(verbose){
      message("Small circles found!")
      message("Removing small circles")
    }

    if(threshold > 0) {
      loop.ind <- which(self$E[,1] == self$E[,2])
      if(length(loop.ind)>0) {
        loop.size <- self$edge_lengths[loop.ind]
        ind <- loop.ind[loop.size < threshold]
        if(length(ind)>0) {
          v.loop <- self$E[ind,1]
          v.degrees <- private$compute_degrees()$degrees[v.loop]

          ind.rem <- v.loop[v.degrees == 2]
          ind.keep <- v.loop[v.degrees > 2]

          if(length(ind.rem)>0) {
            self$V <- self$V[-ind.rem,]
            self$nV <- self$nV - length(ind.rem)
            i.sort <- sort(ind.rem, decreasing = TRUE)
            for(i in i.sort) {
              self$E[self$E >= i] <- self$E[self$E >= i] - 1
            }
          }

          self$edges <- self$edges[-ind]
          self$E <- self$E[-ind,]
          self$nE <- self$nE - length(ind)
        }
      }
    }
  },

  # utility function to merge lines connected by degree 2 vertices
  merge.all.deg2 = function() {
   while(sum(private$compute_degrees()$degrees==2)>0) {
     private$remove.first.deg2()
   }
  },


  remove.first.deg2 = function(res) {
    ind <- which(res$degrees==2 & !res$problematic)
    res.out <- res
    if(length(ind)>0) {

      ind <- ind[1]
      e1 <- which(self$E[,2]==ind)
      e2 <- which(self$E[,1]==ind)
      order_edges <- order(c(e1,e2))
      e_rem <- c(e1,e2)[order_edges]

      # Finding the right order, so the edges can be merged.

      which_line_starts <- which(order_edges == 1)

      v1 <- setdiff(self$E[e_rem[1],],ind)
      v2 <- setdiff(self$E[e_rem[2],],ind)

      if(v1 > ind) {
        v1 <- v1-1
      }
      if(v2 > ind) {
        v2 <- v2 - 1
      }

      # Making sure it is not a single circle

      if( e_rem[1] != e_rem[2]){

         coords <- self$edges[[e_rem[1]]] #line from v1 to v.rem
         tmp <- self$edges[[e_rem[2]]] #line from v.rem to v2

        if(which_line_starts == 1){
           coords <- rbind(coords, tmp)
           E_new <- matrix(c(v1,v2),1,2)
        } else{
           coords <- rbind(tmp,coords)
           E_new <- matrix(c(v2,v1),1,2)
        }


      # } else{
      #   E_new1 <- self$E[e1,1]
      #   E_new2 <- self$E[e2,2]

      #   if(E_new1 > ind){
      #     E_new1 <- E_new1 - 1
      #   }
      #   if(E_new2 > ind){
      #     E_new2 <- E_new2 -1
      #   }
      #   E_new <- matrix(c(E_new1, E_new2),1,2)
      # }

      # Updating the merged graph

      res.out$degrees <- res$degrees[-ind]
      res.out$problematic <- res$problematic[-ind]

      #update vertices
      self$V <- self$V[-ind,]
      self$vertices[[ind]] <- NULL
      self$nV <- self$nV - 1
      # for(i in ind:length(self$vertices)){
      #   attr(self$vertices[[i]], "id") <- attr(self$vertices[[i]], "id") - 1
      # }

      #update edges
      self$E[self$E >= ind] <- self$E[self$E >= ind] - 1
      self$E <- self$E[-e_rem[2],,drop=FALSE]
      self$E[e_rem[1],] <- E_new
      self$edges[[e_rem[1]]] <- coords
      self$edges <- self$edges[-e_rem[2]]
      # for(i in e_rem[2]:length(self$edges)){
      #     attr(self$edges[[i]], "id") <- attr(self$edges[[i]], "id") - 1
      # }
      # attr(self$edges[[e_rem[1]]], "id") <- e_rem[1]

      self$edge_lengths[e_rem[1]] <- self$edge_lengths[e_rem[1]] + self$edge_lengths[e_rem[2]]
      self$edge_lengths <- self$edge_lengths[-e_rem[2]]

      self$nE <- self$nE - 1


      if(is.vector(private$edge_weights)){
        private$edge_weights <- private$edge_weights[-e_rem[2]]
      } else{
        private$edge_weights <- private$edge_weights[-e_rem[2],]
      }

      # attr(self$edges[[e_rem[1]]], "longlat") <- private$longlat
      # class(self$edges[[e_rem[1]]]) <- "metric_graph_edge"
      # attr(self$edges[[e_rem[1]]], "crs") <- private$crs$input
      # attr(self$edges[[e_rem[1]]], "length") <- self$edge_lengths[e_rem[1]]
      # if(!is.null(private$length_unit)){
      #   units(attr(self$edges[[e_rem[1]]], "length")) <- private$length_unit
      # }
      # if(is.vector(private$edge_weights)){
      #   attr(self$edges[[e_rem[1]]], "weight") <- private$edge_weights[e_rem[1]]
      # } else{
      #   attr(self$edges[[e_rem[1]]], "weight") <- private$edge_weights[e_rem[1],]
      # }
      }
    }
    class(self$edges) <- "metric_graph_edges"
    return(res.out)
  },

  # Compute lengths

  compute_lengths = function(longlat, unit, crs, proj4string, which_longlat, vertex_unit, project_data){
          ll <- sapply(self$edges,
          function(edge){compute_line_lengths(edge, longlat = longlat, unit = unit, crs = crs, proj4string, which_longlat, vertex_unit, project_data)})
          return(ll)
  },

 ## @description Get the observation/prediction matrix A
 ## @param group A vector. If `NULL`, the A matrix for the first group will be
 ## returned. One can use all groups by simply setting the `group` variable
 ## to `.all`. Otherwise, the A matrix for the groups in the vector will be
 ## returned.
 ## @param obs_to_vert Should the observations be turned into vertices?
 ## @param include_NA Should the locations for which all observations are NA be
 ## included?
 ## @return The observation or prediction matrix.
  A = function(group = NULL,
               obs_to_vert = FALSE,
               drop_na = FALSE,
               drop_all_na = FALSE){

    if(is.null(self$PtV) && !obs_to_vert){
        stop("The A matrix was not computed. If you want to compute rerun this
             method with 'obs_to_vertex=TRUE', in which the observations will be
             turned to vertices and the A matrix will then be computed")
    } else if(is.null(self$PtV)){
      self$observation_to_vertex(mesh_warning=FALSE)
    }

    if(is.null(group)){
      group <- unique(private$data[[".group"]])
      group <- group[1]
    } else if (group[1] == ".all"){
      group <- unique(private$data[[".group"]])
    }
    n_group <- length(unique(group))

    if(!drop_na && !drop_all_na){
      A <- Matrix::Diagonal(self$nV)[self$PtV, ]
      return(Matrix::kronecker(Diagonal(n_group),A))
    } else {
        data_group <- select_group(private$data, group[1])
        if(drop_na){
          idx_notna <- idx_not_any_NA(data_group)
        } else if(drop_all_na){
          idx_notna <- idx_not_all_NA(data_group)
        }
        # nV_tmp <- sum(idx_notna)
        # A <- Matrix::Diagonal(nV_tmp)[self$PtV[idx_notna], ]
        A <- Matrix::Diagonal(self$nV)[self$PtV[idx_notna], ]
        if(n_group > 1){
          for (i in 2:length(group)) {
            data_group <- select_group(private$data, group[i])
            if(drop_na){
              idx_notna <- idx_not_any_NA(data_group)
            } else if(drop_all_na){
              idx_notna <- idx_not_all_NA(data_group)
            }
            # nV_tmp <- sum(idx_notna)
            A <- Matrix::bdiag(A, Matrix::Diagonal(self$nV)[self$PtV[idx_notna], ])
          }
        }
        return(A)
      }
  },


  #' data List containing data on the metric graph.

  data = NULL,

  # Initial edges added

  initial_edges_added = NULL,

  # Initial graph

  initial_graph = NULL,

  # pruned

  pruned = FALSE,

  # should the information be saved when splitting edges?

  addinfo = FALSE,

  # vertex unit

  vertex_unit = NULL,

  # Length unit

  length_unit = NULL,

  # longlat

  longlat = FALSE,

  # crs

  crs = NULL,

  # proj4string

  proj4string = NULL,

  # which_longlat

  which_longlat = NULL,

  # connected

  connected = TRUE,

  # group columns

  group_col = NULL,

  clear_initial_info = function(){
    private$addinfo = FALSE
    private$initial_edges_added = NULL
  },

  # @description function for splitting lines in the graph
  # @param Ei index of line to split
  # @param t  position on line to split (normalized)
  # @param tolerance tolerance for merging overlapping vertices
  split_edge = function(Ei, t, tolerance = 0) {
    edge <- self$edges[[Ei]]

     val_line <- interpolate2(edge, pos = t, normalized = TRUE, get_idx = TRUE)

     idx_pos <- val_line[["idx"]]
     val_line <- val_line[["coords"]]

    closest_vertex <- which.min(sapply(1:nrow(self$V), function(i){
      (self$V[i,1]-val_line[1])^2 + (self$V[i,2] - val_line[2])^2
    }))

    min_dist <- sqrt(sum((val_line - self$V[closest_vertex,])^2))
    add_V <- FALSE
    if(min_dist <= tolerance){
      newV <- closest_vertex
    } else{
      newV <- self$nV + 1
      add_V <- TRUE
    }

    if((newV != self$E[Ei, 1]) && newV != self$E[Ei,2]){

        # Val_line is in the current edge and closest_vertex is elsewhere
        # So, if we need to add V, we will consider val_line in the edge split,
        # but if we have already added V, then we consider the closest_vertex in the
        # edge split.
        if(add_V){
          self$V <- rbind(self$V, c(val_line))
          coords1 <- rbind(matrix(edge[1:idx_pos,],ncol=2),
                         matrix(val_line,ncol=2))

          coords2 <- rbind(matrix(val_line, ncol=2),
                         matrix(edge[(idx_pos+1):nrow(edge),],
                                ncol=2))
        } else{
          coords1 <- rbind(matrix(edge[1:idx_pos,],ncol=2),
                         matrix(self$V[closest_vertex,],ncol=2))

          coords2 <- rbind(matrix(self$V[closest_vertex,], ncol=2),
                         matrix(edge[(idx_pos+1):nrow(edge),],
                                ncol=2))
        }
        self$nE <- self$nE + 1
        self$E <- rbind(self$E, c(newV, self$E[Ei, 2]))
        self$E[Ei, 2] <- newV
        self$nV <- dim(self$V)[1]

        self$edges[[Ei]] <- coords1
        length(self$edges) <- length(self$edges)+1
        self$edges[[length(self$edges)]] <- coords2
        l_e <- self$edge_lengths[Ei]
        self$edge_lengths[Ei] <- t * l_e
        self$edge_lengths <- c(self$edge_lengths, (1 - t) * l_e)
        if(private$addinfo){
            private$initial_edges_added <- rbind(private$initial_edges_added,cbind(Ei, nrow(self$E)))
        }

        if(is.vector(private$edge_weights)){
          private$edge_weights <- c(private$edge_weights, private$edge_weights[Ei])
        } else{
          private$edge_weights <- rbind(private$edge_weights, private$edge_weights[Ei,])
        }          

        # attr(self$edges[[Ei]], "length") <- t * l_e
        # attr(self$edges[[length(self$edges)]], "length") <- (1 - t) * l_e
        # if(!is.null(private$length_unit)){
        #   units(attr(self$edges[[Ei]], "length")) <- private$length_unit
        #   units(attr(self$edges[[length(self$edges)]], "length")) <- private$length_unit
        # }
        # attr(self$edges[[Ei]], "longlat") <- private$longlat
        # attr(self$edges[[length(self$edges)]], "longlat") <- private$longlat
        # attr(self$edges[[Ei]], "crs") <- private$crs$input
        # attr(self$edges[[length(self$edges)]], "crs") <- private$crs$input
        # if(is.vector(private$edge_weights)){
        #   attr(self$edges[[Ei]],"weight") <- private$edge_weights[Ei]
        #   attr(self$edges[[length(self$edges)]],"weight") <- private$edge_weights[Ei]
        # } else{
        #   attr(self$edges[[Ei]],"weight") <- private$edge_weights[Ei,]
        #   attr(self$edges[[length(self$edges)]],"weight") <- private$edge_weights[Ei,]
        # }
        # attr(self$edges[[Ei]], "id") <- Ei
        # attr(self$edges[[length(self$edges)]], "id") <- length(self$edges)
        # class(self$edges[[Ei]]) <- "metric_graph_edge"
        # class(self$edges[[length(self$edges)]]) <- "metric_graph_edge"

        # class(self$edges) <- "metric_graph_edges"

        if(!is.null(private$data)){
          ind <- which(private$temp_PtE[, 1] %in% Ei)
          for (i in ind) {
            if (private$temp_PtE[i, 2] >= t - tolerance) {
              private$temp_PtE[i, 1] <- self$nE
              private$temp_PtE[i, 2] <- abs(private$temp_PtE[i, 2] - t) / (1 - t)
            }
          }
        }
        return(newV)
    } else{
      return(NULL)
    }
  },

  compute_laplacian_PtE = function(PtE, normalized = TRUE) {

    graph.temp <- self$clone()
    graph.temp$clear_observations()
    df_temp <- data.frame(y = rep(0, dim(PtE)[1]),
                         edge_number = PtE[,1],
                         distance_on_edge = PtE[,2])
    if(sum(duplicated(df_temp))>0){
      warning("Duplicated locations were found when computing the laplacian. The returned values are given for unique locations.")
      df_temp <- unique(df_temp)
    }

    graph.temp$build_mesh(h = 1000)

    df_temp2 <- data.frame(y = 0, edge_number = graph.temp$mesh$VtE[1:nrow(self$V),1],
                                  distance_on_edge = graph.temp$mesh$VtE[1:nrow(self$V),2])
    df_temp$included <- TRUE
    temp_merge <- merge(df_temp, df_temp2, all = TRUE)

    df_temp$included <- NULL
    df_temp2 <- temp_merge[is.na(temp_merge["included"]),]
    df_temp2$included <- NULL
    df_temp <- rbind(df_temp2, df_temp)

    nV_new <- sum(is.na(temp_merge["included"]))

    df_temp[["__dummy"]] <- 1:nrow(df_temp)

    graph.temp$add_observations(data = df_temp, normalized = normalized)
    graph.temp$observation_to_vertex(mesh_warning = FALSE)
    Wmat <- Matrix(0,graph.temp$nV,graph.temp$nV)

    for (i in 1:graph.temp$nE) {
      Wmat[graph.temp$E[i, 1], graph.temp$E[i, 2]] <- 1 / graph.temp$edge_lengths[i]
      Wmat[graph.temp$E[i, 2], graph.temp$E[i, 1]] <- 1 / graph.temp$edge_lengths[i]
    }
    Laplacian <- Matrix::Diagonal(graph.temp$nV,
                                  as.vector(Matrix::rowSums(Wmat))) - Wmat

    # Reordering from vertices to points
    Laplacian <- Laplacian[graph.temp$PtV, graph.temp$PtV]
    # Order back to the input order
    Laplacian[graph.temp$.__enclos_env__$private$data[["__dummy"]], graph.temp$.__enclos_env__$private$data[["__dummy"]]] <- Laplacian

    attr(Laplacian, "nV_idx") <- nV_new

    return(Laplacian)
  },

  find_edge_edge_points = function(tol,verbose, crs, proj4string, longlat, fact, which_longlat) {

    if(!private$longlat){
      lines_sf <- sf::st_sfc(lapply(self$edges, function(i){sf::st_linestring(i)}))
      crs <- NULL
    } else if (private$which_longlat == "sf"){
      lines_sf <- sf::st_sfc(lapply(self$edges, function(i){sf::st_linestring(i)}), crs = private$crs)
    } else{
      lines_sf <- sf::st_sfc(lapply(self$edges, function(i){sf::st_linestring(i)}), crs = sf::st_crs(private$proj4string))
      crs <- sf::st_crs(proj4string)
    }

  dists <- t(as.matrix(sf::st_is_within_distance(lines_sf, dist = tol)))
  points_add <- NULL
  points_add_PtE <- NULL

  if(verbose){
    bar_line_line <- msg_progress_bar(length(self$edges)-1)
  }
  for(i in 1:(length(self$edges)-1)) {
    if(verbose){
      bar_line_line$increment()
    }
    #lines within tol of line i
    inds <- i+which(as.vector(dists[i, (i+1):length(self$edges)]))
    if(length(inds)>0) {
      for(j in inds) {
        #first check if there are intersections

        intersect_tmp <- intersection3(lines_sf[i], lines_sf[j])

        if( any(c("GEOMETRYCOLLECTION","GEOMETRY") %in% sf::st_geometry_type(intersect_tmp)) || (length(sf::st_geometry_type(intersect_tmp)) > 1)){
          intersect_tmp1 <- sf::st_collection_extract(intersect_tmp, type = "LINESTRING")
          intersect_tmp2 <- sf::st_collection_extract(intersect_tmp, type = "POINT")
          if(nrow(sf::st_coordinates(intersect_tmp1))>0){
            intersect_tmp <- intersect_tmp1
          } else{
            intersect_tmp <- intersect_tmp2
          }
        }

        p_cur <- NULL
        # if(!is.null(intersect_tmp)) {
        if(nrow(sf::st_coordinates(intersect_tmp))>0){
          # if("SpatialPoints"%in%is(intersect_tmp)){
          if( ("POINT"%in%sf::st_geometry_type(intersect_tmp))||("MULTIPOINT"%in%sf::st_geometry_type(intersect_tmp))){
            # coord_tmp <- coordinates(intersect_tmp)
            coord_tmp <- sf::st_coordinates(intersect_tmp)
            coord_tmp <- coord_tmp[,1:2]
          # } else if ("SpatialLines"%in%is(intersect_tmp)){
          } else if ( ("LINESTRING"%in%sf::st_geometry_type(intersect_tmp)) || ("MULTILINESTRING"%in%sf::st_geometry_type(intersect_tmp))){
            intersect_tmp <- sf::st_coordinates(intersect_tmp)
            intersect_tmp <- intersect_tmp[,1:2]
            # coord_tmp <-gInterpolate(intersect_tmp, d=0.5, normalized = TRUE)
            coord_tmp <- interpolate2(intersect_tmp, pos=0.5, normalized = TRUE)
            coord_tmp <- matrix(coordinates(coord_tmp),1,2)
          }

            # for(k in 1:length(intersect_tmp)) {
            coord_tmp <- matrix(coord_tmp, ncol=2)

            for(k in 1:nrow(coord_tmp)){
              p <- matrix(coord_tmp[k,],1,2)
              #add points if they are not close to V or previous points
              if(!is.matrix(self$V)){
                self$V <- matrix(self$V,ncol=2)
              }
              if(min(compute_aux_distances(lines = self$V, crs=private$crs, longlat=private$longlat, proj4string = private$proj4string, points = p, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit))>tol) {
                p_cur <- rbind(p_cur,p)
                p2 <- snapPointsToLines(p,self$edges[i], longlat, crs)
                p2 <- t(p2[["coords"]])
                # p2 <- cbind(p2["X",], p2["Y",])

                points_add <- rbind(points_add, p, p2)
                points_add_PtE <- rbind(points_add_PtE,
                                        c(i,projectVecLine2(self$edges[[i]],
                                                     p)),
                                        c(j,projectVecLine2(self$edges[[j]],p)))
              }
            }
        }
        #now check if there are intersections with buffer

        lines_tmp_sf <- lines_sf[i]
        lines2_tmp_sf <- lines_sf[j]
        tmp_line <- sf::st_buffer(lines_tmp_sf, dist = tol)

        intersect_tmp <- intersection3(tmp_line, lines2_tmp_sf)

        if( "GEOMETRYCOLLECTION" %in% sf::st_geometry_type(intersect_tmp)){
          intersect_tmp <- sf::st_collection_extract(intersect_tmp, type = "LINESTRING")
        }
        # if(!is.null(intersect_tmp)) {
        if(nrow(sf::st_coordinates(intersect_tmp))>0){
            for(k in 1:length(intersect_tmp)) {
              if ( ("LINESTRING"%in%sf::st_geometry_type(intersect_tmp)) || ("MULTILINESTRING"%in%sf::st_geometry_type(intersect_tmp))){
              # coord_tmp <-gInterpolate(intersect_tmp[k], d=0.5, normalized = TRUE)
              coord_int <- sf::st_coordinates(intersect_tmp[k])
              coord_int <- coord_int[,c("X","Y")]
              coord_tmp <- interpolate2(coord_int, pos=0.5, normalized = TRUE)
              p <- matrix(coord_tmp,1,2)
            } else {
              coord_int <- sf::st_coordinates(intersect_tmp[k])
              p <- matrix(coord_int[,c("X","Y")],1,2)
              # p <- matrix(sf::st_coordinates(intersect_tmp[k]),1,2)
            }

            #add points if they are not close to V or previous points
            if(min(compute_aux_distances(lines = self$V, crs=private$crs, longlat=private$longlat, proj4string = private$proj4string, points = p, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit))>tol) {
              # if(is.null(p_cur) || gDistance(SpatialPoints(p_cur), intersect_tmp[k])>tol) {
                if(!private$longlat && !is.null(p_cur)){
                  dist_tmp <- sf::st_distance(sf::st_as_sf(as.data.frame(p_cur), coords = 1:2), intersect_tmp[k])
                } else if (!is.null(p_cur)) {
                  intersect_tmp_sfc <- sf::st_sfc(intersect_tmp[k], crs = private$crs)
                  dist_tmp <- sf::st_distance(sf::st_as_sf(as.data.frame(p_cur), coords = 1:2, crs = private$crs), intersect_tmp_sfc)
                  units(dist_tmp) <- private$length_unit
                  units(dist_tmp) <- NULL
                }
              if(is.null(p_cur) || min(dist_tmp) >tol) {
                p2 <- snapPointsToLines(p,self$edges[i], private$longlat, private$crs)
                p2 <- t(p2[["coords"]])
                points_add <- rbind(points_add, p, p2)
                points_add_PtE <- rbind(points_add_PtE,
                                        c(i,projectVecLine2(self$edges[[i]],p)),
                                        c(j,projectVecLine2(self$edges[[j]], p)))

              }
            }
          }
        }
      }
    }
  }
  # ret_tmp <- list(points = points_add, PtE = points_add_PtE)
  return(list(points = points_add, PtE = points_add_PtE))
},

add_vertices = function(PtE, tolerance = 1e-10, verbose) {
  e.u <- unique(PtE[,1])
  if(verbose){
    bar_eu <- msg_progress_bar(length(e.u))
  }
  for (i in 1:length(e.u)) {
    if(verbose){
      bar_eu$increment()
    }
    dists <- sort(PtE[which(PtE[,1]==e.u[i]),2])
    if(length(dists) > 0){
    private$split_edge(e.u[i], dists[1], tolerance)
    }
    if(length(dists)>1) {
      dists_up <- dists
      for(j in 2:length(dists)){
        dists_up[j] <- (dists[j] - dists[j-1])/(1 - dists[j-1])
        private$split_edge(self$nE, dists_up[j], tolerance)
      }
    }
  }
  return(self)
},

 # Function to compute the degrees of the vertices,

  compute_degrees = function(){
    degrees_in <- rep(0,self$nV)
    degrees_out <- rep(0,self$nV)
    for(i in 1:self$nV) {
          degrees_out[i] <- sum(self$E[,1]==i)
          degrees_in[i] <- sum(self$E[,2]==i)
    }
    degrees <- degrees_in + degrees_out
    return(list(degrees = degrees, indegrees = degrees_in,
              outdegrees = degrees_out))
  },

  #  Creates/updates the vertices element of the metric graph list

  create_update_vertices = function(){
    degrees <- private$compute_degrees()
    self$vertices <- lapply(1:nrow(self$V),
        function(i){
          vert <- self$V[i,]
          attr(vert, "degree") <- degrees$degrees[i]
          attr(vert, "indegree") <- degrees$indegrees[i]
          attr(vert, "outdegree") <- degrees$outdegrees[i]
          attr(vert, "problematic") <- ifelse((degrees$degrees[i]>1) && ((degrees$indegrees[i] == 0) || (degrees$outdegrees[i] == 0)), TRUE, FALSE)
          attr(vert, "longlat") <- private$longlat
          attr(vert, "crs") <- private$crs$input
          attr(vert, "id") <- i
          class(vert) <- "metric_graph_vertex"
          return(vert)
        })
    class(self$vertices) <- "metric_graph_vertices"
  },

  # Merge degree 2 vertices in mesh
  mesh_merge_deg2 = function() {
    outs <- rep(0,self$nV)
    ins <- rep(0, self$nV)
    for(i in 1:self$nV) {
      outs[i] <- sum(self$E[,1] == i)
      ins[i] <- sum(self$E[,2] == i)
    }
    ind <- which(outs == 1 & ins == 1)
    if(length(ind)>0) {
      for(i in 1:length(ind)) {
        V.keep <- self$mesh$V[ind[i],]
        V.rem <- which(self$mesh$V[,1] == V.keep[1] & self$mesh$V[,2] == V.keep[2])[2]
        self$mesh$PtE <- self$mesh$PtE[-V.rem,]
        self$mesh$V <- self$mesh$V[-V.rem,]
        self$mesh$E[self$mesh$E == V.rem] <- ind[i]
        self$mesh$E[self$mesh$E>V.rem] <- self$mesh$E[self$mesh$E>V.rem] - 1
      }
    }
  },

  mesh_merge_outs = function() {
    outs <- rep(0,self$nV)
    for(i in 1:self$nV) {
      outs[i] <- sum(self$E[,1] == i)
    }
    ind <- which(outs > 1)
    while(length(ind)>0) {
      #find edges going out
      e.ind <- which(self$E[,1]==ind[1])
      V.keep <- which(self$mesh$PtE[,1]==e.ind[1])[1]
      V.rem <- which(self$mesh$PtE[,1]==e.ind[2])[1]
      self$mesh$PtE <- self$mesh$PtE[-V.rem,]
      self$mesh$V <- self$mesh$V[-V.rem,]
      self$mesh$E[self$mesh$E == V.rem] <- V.keep
      self$mesh$E[self$mesh$E>V.rem] <- self$mesh$E[self$mesh$E>V.rem] - 1
      outs[ind[1]] <- outs[ind[1]] - 1
      ind <- which(outs > 1)
    }
    #now merge ins with outs if there is only one in
    ins <- rep(0,self$nV)
    outs <- rep(0,self$nV)
    for(i in 1:self$nV) {
      ins[i] <- sum(self$E[,2] == i)
      outs[i] <- sum(self$E[,1] == i)
    }
    ind <- which(ins == 1 & outs > 1)
    while(length(ind)>0) {
      #merge in vertex with out
      e.ind.in <- which(self$E[,2]==ind[1])
      e.ind <- which(self$E[,1]==ind[1])
      V.keep <- which(self$mesh$PtE[,1]==e.ind[1])[1]
      V.rem <- which(self$mesh$PtE[,1]==e.ind.in)
      V.rem <- V.rem[length(V.rem)]
      self$mesh$PtE <- self$mesh$PtE[-V.rem,]
      self$mesh$V <- self$mesh$V[-V.rem,]
      self$mesh$E[self$mesh$E == V.rem] <- V.keep
      self$mesh$E[self$mesh$E>V.rem] <- self$mesh$E[self$mesh$E>V.rem] - 1
      ins[ind[1]] <- ins[ind[1]] - 1
      ind <- which(ins == 1 & outs > 1)
    }
  },

  #find one mesh node corresponding to each vertex and move it first
  move_V_first = function() {
    nv <- dim(self$mesh$V)[1]
    for(i in 1:self$nV) {
      ind <- which(self$mesh$V[,1] == self$V[i,1] & self$mesh$V[,2] == self$V[i,2])[1]

      if(ind > i && i < nv) {
        if (i == 1) {
          reo <- c(ind, setdiff(i:nv,ind))
        } else {
          reo <- c(1:(i-1), ind, setdiff(i:nv,ind))
        }
        self$mesh$V <- self$mesh$V[reo,]
        self$mesh$PtE <- self$mesh$PtE[reo,]
        self$mesh$VtE <- self$mesh$VtE[reo,]
        Etmp <- self$mesh$E
        ind1 <- Etmp == ind
        ind2 <- Etmp >= i & Etmp < ind
        self$mesh$E[ind1] = i
        self$mesh$E[ind2] = self$mesh$E[ind2] + 1
      }
    }
  },

  # find vertices and connections for Petrov-Galerkin matrices
  find_mesh_bc = function() {
    if(attr(self$mesh,"continuous")) {
      stop("mesh discontinuous")
    }
    starts <- which(self$mesh$PtE[,2]==0)
    edges <- self$mesh$PtE[starts,1]

    ind.keep <- rep(TRUE,length(starts))
    for(i in 1:length(starts)) {
      vert <- self$E[edges[i],1]
      if(length(which(self$E[,2] == vert)) == 1) {
        ind.keep[i] <- FALSE
      }
    }
    starts <- starts[ind.keep]
    edges <- edges[ind.keep]
    connections <-list()
    for(i in 1:length(starts)) {
      vert <- self$E[edges[i],1] #the vertex of the start point
      n.in <- sum(self$E[,2] == vert)
      if(n.in == 0) {
        connections[[i]] <- 0 #true starting value
      } else {
        #find edges ending in the vertex
        edge.in <- which(self$E[,2] == vert)
        #find edges starting in the vertex
        edge.out <- which(self$E[,1] == vert)

        mesh.in <- NULL
        for(j in 1:length(edge.in)){
          tmp <- which(self$mesh$PtE[,1] == edge.in[j]) #mesh nodes on the edge
          mesh.in <- c(mesh.in, tmp[length(tmp)])
        }
        connections[[i]] <- mesh.in
      }
    }
    return(connections)
  },

  # Function to add vertex conditions to Petrov-Galerkin matrices
  set_petrov_matrices = function() {
    if(attr(self$mesh,"continuous")) {
      V <- self$mesh$V[1:self$nV,]
      deg1 <- which(self$get_degrees() == 1)
      starts <- deg1[deg1 %in% self$E[,1]]

      if(length(starts)>1){

        G <- rbind(sparseMatrix(i = 1:length(starts), j = starts,
                                x = rep(0, length(starts)),
                                dims = c(length(starts), dim(self$mesh$Cpet)[1])),
                   t(self$mesh$Gpet[,-starts[-1]]))
      } else {
        C <- rbind(sparseMatrix(i = 1, j = starts, x = 1,
                                dims = c(1,dim(self$mesh$Cpet)[1])),
                   t(self$mesh$Cpet))
        G <- rbind(sparseMatrix(i = 1, j = starts, x= 0,
                                dims = c(1, dim(self$mesh$Cpet)[1])),
                   t(self$mesh$Gpet))
      }
    } else {
      bc <- private$find_mesh_bc()
      starts <- which(self$mesh$PtE[,2]==0)
      edges <- self$mesh$PtE[starts,1]

      ind.keep <- rep(TRUE,length(starts))
      for(i in 1:length(starts)) {
        vert <- self$E[edges[i],1]
        if(length(which(self$E[,2] == vert)) == 1) {
          ind.keep[i] <- FALSE
        }
      }
      starts <- starts[ind.keep]
      edges <- edges[ind.keep]
      C <- t(self$mesh$Cpet)
      G <- t(self$mesh$Gpet)

      for(i in length(bc):1) {
        G <- rbind(sparseMatrix(i = 1, j = 1, x = 0,
                                dims = c(1, dim(self$mesh$Cpet)[1])), G)
        if(bc[[i]][1] == 0) {
          C <- rbind(sparseMatrix(i = 1,j = starts[i], x = 1,
                                  dims = c(1, dim(self$mesh$Cpet)[1])), C)
        } else {
          C <- rbind(sparseMatrix(i = rep(1,length(bc[[i]])+1),
                                  j = c(starts[i], bc[[i]]),
                                  x = c(1, rep(-1/length(bc[[i]]), length(bc[[i]]))),
                                  dims = c(1, dim(self$mesh$Cpet)[1])), C)
        }

      }
    }
    self$mesh$Cpet = C
    self$mesh$Gpet = G
    self$mesh$n.bc = length(starts)
    self$mesh$h0 = which(unlist(lapply(bc,length))>1)
  },

  # Temp PtE

  temp_PtE = NULL,

  # # longlat
  # longlat = NULL,

  # tolerance

  tolerance = NULL,

  # edge_weights

  edge_weights = NULL

))

#' @title Connected components of metric graph
#' @description Class representing connected components of a metric graph.
#' @details A list of `metric_graph` objects (representing the different
#' connected components in the full graph) created from vertex and edge matrices,
#' or from an sp::SpatialLines object where each line is representing and edge.
#' For more details, see the vignette:
#' \code{vignette("metric_graph", package = "MetricGraph")}
#' @return Object of \code{\link[R6]{R6Class}} for creating metric graph components.
#' @examples
#' library(sp)
#' edge1 <- rbind(c(0, 0), c(1, 0))
#' edge2 <- rbind(c(1, 0), c(2, 0))
#' edge3 <- rbind(c(1, 1), c(2, 1))
#' edges <- list(edge1, edge2, edge3)
#'
#' graphs <- graph_components$new(edges)
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

   #' Create metric graphs for connected components
   #'
  #' @param edges A list containing coordinates as `m x 2` matrices (that is, of `matrix` type) or m x 2 data frames (`data.frame` type) of sequence of points connected by straightlines. Alternatively, you can also prove an object of type `SpatialLinesDataFrame` or `SpatialLines` (from `sp` package) or `MULTILINESTRING` (from `sf` package).
   #' @param V n x 2 matrix with Euclidean coordinates of the n vertices.
   #' @param E m x 2 matrix where each row represents an edge.
  #' @param vertex_unit The unit in which the vertices are specified. The options are 'degrees' (the great circle distance in km), 'km', 'm' and 'miles'. The default is `NULL`, which means no unit. However, if you set `length_unit`, you need to set `vertex_unit`.
  #' @param length_unit The unit in which the lengths will be computed. The options are 'km', 'm' and 'miles'. The default is `vertex_unit`. Observe that if `vertex_unit` is `NULL`, `length_unit` can only be `NULL`.
  #' If `vertex_unit` is 'degrees', then the default value for `length_unit` is 'km'.
  #' @param longlat If TRUE, then it is assumed that the coordinates are given.
  #' in Longitude/Latitude and that distances should be computed in meters. It takes precedence over
  #' `vertex_unit` and `length_unit`, and is equivalent to `vertex_unit = 'degrees'` and `length_unit = 'm'`.
   #' @param tolerance Vertices that are closer than this number are merged when
   #' constructing the graph (default = 1e-10). If `longlat = TRUE`, the
   #' tolerance is given in km.
   #' @param by_length Sort the components by total edge length? If `FALSE`,
   #' the components are sorted by the number of vertices.
   #' @param ... Additional arguments used when specifying the graphs
   #' @param lines `r lifecycle::badge("deprecated")` Use `edges` instead.
   #' @return A `graph_components` object.
   initialize = function(edges = NULL,
                         V = NULL,
                         E = NULL,
                         by_length = TRUE,
                         ...,
                         lines = deprecated()) {


      if (lifecycle::is_present(lines)) {
         if (is.null(edges)) {
           lifecycle::deprecate_warn("1.2.0", "graph_components$new(lines)", "graph_components$new(edges)",
             details = c("`lines` was provided but not `edges`. Setting `edges <- lines`.")
           )
           edges <- lines
         } else {
           lifecycle::deprecate_warn("1.2.0", "graph_components$new(lines)", "graph_components$new(edges)",
             details = c("Both `edges` and `lines` were provided. Only `edges` will be considered.")
           )
         }
         lines <- NULL
       }


      dots_args <- list(...)
      dots_list <- as.list(dots_args)
      if(!is.null(dots_list[["project_data"]])){
        warning("The argument project_data is not compatible with graph_components. Setting project_data to FALSE.")
        dots_list[["project_data"]] <- FALSE
        dots_list[["edges"]] <- edges
        dots_list[["V"]] <- V
        dots_list[["E"]] <- E
        dots_list[["check_connected"]] <- FALSE
        graph <- do.call(metric_graph$new, dots_list)
      } else{
            graph <- metric_graph$new(edges = edges, V = V, E = E,
                               check_connected = FALSE, ...)
      }


     g <- graph(edges = c(t(graph$E)), directed = FALSE)
     igraph::E(g)$weight <- graph$edge_lengths
     components <- igraph::clusters(g, mode="weak")

     self$n <- components$no
     if(self$n > 1) {
       self$graphs <- vector(mode = "list", length = self$n)
       for(k in 1:self$n) {
         vert_ids <- igraph::V(g)[components$membership == k]
         edge_rem <- NULL
         for (i in 1:graph$nE) {
           if(!(graph$E[i, 1] %in% vert_ids) && !(graph$E[i, 2] %in% vert_ids))
             edge_rem <- c(edge_rem, i)
         }
         edge_keep <- setdiff(1:graph$nE, edge_rem)
         ind_keep <- rep(0,graph$nE)
         ind_keep[edge_keep] <- 1
         if(length(graph$edges[which(ind_keep!=0)]) > 0){
          self$graphs[[k]] = metric_graph$new(edges = graph$edges[which(ind_keep!=0)],
                                             check_connected = FALSE, ...)
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
   },

   #' @description Returns the largest component in the graph.
   #' @return A `metric_graph` object.
   get_largest = function() {
     return(self$graphs[[1]])
   },

   #' @description Plots all components.
   #' @param edge_colors A 3 x nc matrix with RGB values for the edge colors to
   #' be used when plotting each graph.
   #' @param vertex_colors A 3 x nc matrix with RGB values for the edge colors to
   #' be used when plotting each graph.
   #' @param ... Additional arguments for plotting the individual graphs.
   #' @return A `ggplot` object.
   plot = function(edge_colors = NULL, vertex_colors = NULL, ...) {

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
     p <- self$graphs[[1]]$plot(edge_color = rgb(edge_colors[1, 1],
                                                 edge_colors[1, 2],
                                                 edge_colors[1, 3]),
                                vertex_color = rgb(vertex_colors[1, 1],
                                                   vertex_colors[1, 2],
                                                   vertex_colors[1, 3]), ...)
     if (self$n > 1) {
       for(i in 2:self$n){
         suppressMessages(p <- self$graphs[[i]]$plot(edge_color = rgb(edge_colors[i, 1],
                                                     edge_colors[i, 2],
                                                     edge_colors[i, 3]),
                                    vertex_color = rgb(vertex_colors[i, 1],
                                                       vertex_colors[i, 2],
                                                       vertex_colors[i, 3]),
                                    p = p, ...))
       }
     }
     return(p)
   }))
