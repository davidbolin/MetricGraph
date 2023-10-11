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

  #' @field data List containing data on the metric graph.
  data = NULL,

  #' @field PtV Vector with the indices of the vertices which are observation
  #' locations.
  PtV  = NULL,

  #' @field mesh Mesh object used for plotting.
  mesh = NULL,

  #' @field edges The coordinates of the edges in the graph.
  edges = NULL,

  #' @field geo_dist Geodesic distances between the vertices in the graph.
  geo_dist = NULL,

  #' @field res_dist Resistance distances between the observation locations.
  res_dist = NULL,

  #' @field Laplacian The weighted graph Laplacian of the vertices in the
  #' graph. The weights are given by the edge lengths.
  Laplacian = NULL,

  #' @description Create a new `metric_graph` object.
  #' @param edges A list containing coordinates as `m x 2` matrices (that is, of `matrix` type) or m x 2 data frames (`data.frame` type) of sequence of points connected by straightlines. Alternatively, you can also prove an object of type `SpatialLinesDataFrame` or `SpatialLines` (from `sp` package) or `MULTILINESTRING` (from `sf` package).
  #' @param V n x 2 matrix with Euclidean coordinates of the n vertices.
  #' @param E m x 2 matrix where each row represents one of the m edges.
  #' @param vertex_unit The unit in which the vertices are specified. The options are 'degrees' (the great circle distance in km), 'km', 'm' and 'miles'. The default is `NULL`, which means no unit. However, if you set `length_unit`, you need to set `vertex_unit`.
  #' @param length_unit The unit in which the lengths will be computed. The options are 'km', 'm' and 'miles'. The default is `vertex_unit`. Observe that if `vertex_unit` is `NULL`, `length_unit` can only be `NULL`.
  #' If `vertex_unit` is 'degrees', then the default value for `length_unit` is 'km'.
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
                        remove_circles = TRUE,
                        verbose = FALSE,
                        lines = deprecated()) {

      start_construction_time <- Sys.time()
      if (lifecycle::is_present(lines)) {
         if (is.null(edges)) {
           lifecycle::deprecate_warn("1.1.2.9000", "metric_graph$new(lines)", "metric_graph$new(edges)",
             details = c("`lines` was provided but not `edges`. Setting `edges <- lines`.")
           )
           edges <- lines
         } else {
           lifecycle::deprecate_warn("1.1.2.9000", "metric_graph$new(lines)", "metric_graph$new(edges)",
             details = c("Both `edges` and `lines` were provided. Only `edges` will be considered.")
           )
         }
         lines <- NULL
       }            

      if (is.null(tolerance$vertex_edge) && !is.null(tolerance$vertex_line)) {
           lifecycle::deprecate_warn("1.1.2.9000", "metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
             details = c("`tolerance$vertex_line` was provided but not `tolerance$vertex_edge`. Setting `tolerance$vertex_edge <- tolerance$vertex_line`.")
           )
           tolerance$vertex_edge <- tolerance$vertex_line
           tolerance$vertex_line <- NULL
         } else if(!is.null(tolerance$vertex_edge) && !is.null(tolerance$vertex_line)) {
           lifecycle::deprecate_warn("1.1.2.9000","metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
             details = c("Both `tolerance$vertex_edge` and `tolerance$vertex_line` were provided. Only `tolerance$vertex_edge` will be considered.")
           )
            tolerance$vertex_line <- NULL
         }           

        if (is.null(tolerance$edge_edge) && !is.null(tolerance$line_line)) {
           lifecycle::deprecate_warn("1.1.2.9000", "metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
             details = c("`tolerance$line_line` was provided but not `tolerance$edge_edge`. Setting `tolerance$edge_edge <- tolerance$line_line`.")
           )
           tolerance$edge_edge <- tolerance$line_line
           tolerance$line_line <- NULL
         } else if(!is.null(tolerance$edge_edge) && !is.null(tolerance$line_line)) {
           lifecycle::deprecate_warn("1.1.2.9000","metric_graph$new(tolerance = 'must contain either vertex_vertex, vertex_edge or edge_edge')",
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
        warning("object initialized from lines, then E and V are ignored")
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
        message("Snap vertices to close lines")
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
      private$merge_close_vertices(tolerance$vertex_vertex, factor_unit)
    }
    
    if(is.logical(remove_circles)){
      if(remove_circles){
        private$remove_circles(tolerance$vertex_vertex, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
      }
    } else {
        private$remove_circles(remove_circles, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
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

    private$initial_graph <- self$clone()

    # Cloning again to add the initial graph to the initial graph
    private$initial_graph <- self$clone()



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
      }
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
    self$geo_dist <- list()

    if(is.null(self$data)){
      obs <- FALSE
    }
    if(!obs){
      g <- graph(edges = c(t(self$E)), directed = FALSE)
      E(g)$weight <- self$edge_lengths
      self$geo_dist[["__vertices"]] <- distances(g)
    } else if(full){
      PtE_full <- self$get_PtE()
      self$geo_dist[["__complete"]] <- self$compute_geodist_PtE(PtE = PtE_full,
                                                              normalized = TRUE)
    } else{
      if(is.null(group)){
          group <- unique(self$data[["__group"]])
      }
      for(grp in group){
          data_grp <- select_group(self$data, grp)
          idx_notna <- idx_not_all_NA(data_grp)
          PtE_group <- cbind(data_grp[["__edge_number"]][idx_notna],
                     data_grp[["__distance_on_edge"]][idx_notna])
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
      geodist_temp[graph.temp$data[["__dummy"]],graph.temp$data[["__dummy"]]] <- geodist_temp
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
    if(is.null(self$data)){
      obs <- FALSE
    }
    if(!obs){
      graph.temp <- self$clone()
      graph.temp$build_mesh(h=1000)

      PtE <- graph.temp$mesh$VtE[1:nrow(self$V),]
      rm(graph.temp)
      self$res_dist[["__vertices"]] <- self$compute_resdist_PtE(PtE,
                                                                normalized=TRUE)
    } else if(full){
      PtE <- self$get_PtE()
      self$res_dist[["__complete"]] <- self$compute_resdist_PtE(PtE,
                                                                normalized=TRUE)
    } else{
      if(is.null(group)){
          group <- unique(self$data[["__group"]])
      }
      for(grp in group){
        data_grp <- select_group(self$data, grp)
        idx_notna <- idx_not_all_NA(data_grp)
        if(sum(idx_notna) == 0){
          stop("There are no non-NA observations.")
        }
        PtE <- cbind(data_grp[["__edge_number"]][idx_notna],
                     data_grp[["__distance_on_edge"]][idx_notna])
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
        geodist_temp <- graph.temp$geo_dist[["__complete"]]
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
      R[graph.temp$data[["__dummy"]],graph.temp$data[["__dummy"]]] <- R

      if(!include_vertices){
        R <- R[(nV_new+1):nrow(R), (nV_new+1):nrow(R)]
      }

      return(R)
  },

  #' @description Returns the degrees of the vertices in the metric graph.
  #' @return A vector containing the degrees of the vertices.
  get_degrees = function(){
    degrees <- rep(0,self$nV)
    for(i in 1:self$nV) {
          degrees[i] <- sum(self$E[,1]==i) + sum(self$E[,2]==i)
    }
    return(degrees)
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
  #' locations?
  #' @param group Vector or list containing which groups to compute the
  #' Laplacian for. If `NULL`, it will be computed for all groups.
  #' @return No reutrn value. Called for its side effects. The Laplacian is stored
  #' in the `Laplacian` element in the `metric_graph` object.
  compute_laplacian = function(full = FALSE, obs = TRUE, group = NULL) {
    self$Laplacian <- list()
    if(is.null(self$data)){
      obs <- FALSE
    }
    if(!obs){
      graph.temp <- self$clone()
      graph.temp$build_mesh(h=1000)

      PtE <- graph.temp$mesh$VtE[1:nrow(self$V),]
      rm(graph.temp)
      self$Laplacian[["__vertices"]] <- private$compute_laplacian_PtE(PtE,
                                                            normalized = TRUE)
    } else if(full){
      PtE <- self$get_PtE()
      self$Laplacian[["__complete"]] <- private$compute_laplacian_PtE(PtE,
                                                            normalized = TRUE)
    } else{
      if(is.null(group)){
          group <- unique(self$data[["__group"]])
      }
      for(grp in group){
          data_grp <- select_group(self$data, grp)
          idx_notna <- idx_not_all_NA(data_grp)
          PtE <- cbind(data_grp[["__edge_number"]][idx_notna],
                       data_grp[["__distance_on_edge"]][idx_notna])
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
    degrees <- self$get_degrees()
    start.deg <- end.deg <- rep(0,self$nV)
    for(i in 1:self$nV) {
      start.deg[i] <- sum(self$E[,1]==i)
      end.deg[i] <- sum(self$E[,2]==i)
    }

    # Finding problematic vertices, that is, vertices with incompatible directions
    # They will not be pruned.
    problematic <- (degrees > 1) & (start.deg == 0 | end.deg == 0)
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

   if(!is.null(self$data)){
    if(verbose){
      message("Updating data locations.")
    }    
      t <- system.time({
      x_coord <- self$data[["__coord_x"]]
      y_coord <- self$data[["__coord_y"]]
      new_PtE <- self$coordinates(XY = cbind(x_coord, y_coord))
      group_vec <- self$data[["__group"]]
      self$data[["__edge_number"]] <- new_PtE[,1]
      self$data[["__distance_on_edge"]] <- new_PtE[,2]
      order_idx <- order(group_vec, new_PtE[,1], new_PtE[,2])
      self$data <- lapply(self$data, function(dat){dat[order_idx]})
      })
      if(verbose){
            message(sprintf("time: %.3f s", t[["elapsed"]]))
      }      
   }

   if(!is.null(self$mesh)){
    max_h <- max(self$mesh$h_e)
    self$mesh <- NULL
    self$build_mesh(h = max_h)
   }
   private$pruned <- TRUE
  },

  #' @description Gets PtE from the data.
  #' @return A matrix with two columns, where the first column contains the edge
  #' number and the second column contains the distance on edge of the
  #' observation locations.
  get_PtE = function() {
    if(is.null(self$data)){
      stop("There is no data!")
    }
    group <- self$data[["__group"]]
    group <- which(group == group[1])
    PtE <- cbind(self$data[["__edge_number"]][group],
                 self$data[["__distance_on_edge"]][group])

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
     if(is.null(self$data)){
      stop("There is no data!")
    }
    group <- self$data[["__group"]]
    group <- which(group == group[1])
    Spoints <- data.frame(x = self$data[["__coord_x"]][group], y = self$data[["__coord_y"]][group])
    if(private$longlat){
      colnames(Spoints) <- c("lon", "lat")
    }
    return(Spoints)
  },

  #' @description Adds observation locations as vertices in the graph.
  #' @param tolerance Observations locations are merged to a single vertex if
  #' they are closer than this number (given in relative edge distance between
  #' 0 and 1). The default is `1e-15`.
  #' @param mesh_warning Display a warning if the graph structure change and the metric graph has a mesh object.
  #' @return No return value. Called for its side effects.
  observation_to_vertex = function(tolerance = 1e-15, mesh_warning = TRUE) {
    if(tolerance <= 0 || tolerance >=1){
      stop("tolerance should be between 0 and 1.")
    }
    private$temp_PtE <- self$get_PtE()
    n_group <- length(unique(self$data[["__group"]]))
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

    self$data[["__edge_number"]] <- rep(private$temp_PtE[,1],
                                        times = n_group)
    self$data[["__distance_on_edge"]] <- rep(private$temp_PtE[,2],
                                             times = n_group)

    tmp_df <- data.frame(PtE1 = self$data[["__edge_number"]],
              PtE2 = self$data[["__distance_on_edge"]],
              group = self$data[["__group"]])
    index_order <- order(tmp_df$group, tmp_df$PtE1, tmp_df$PtE2)
    self$data <- lapply(self$data, function(dat){ dat[index_order]})

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
      stop("There is no mesh!")
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
    self$data <- NULL
    self$geo_dist <- NULL
    self$res_dist <- NULL
    self$PtV <- NULL
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
  #' @param group If the data is grouped (for example measured at different time
  #' points), this argument specifies the the column (or entry on the list) in
  #' which the group varialbe is stored.
  #' @param normalized if TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default FALSE. Will not be used if
  #' `Spoints` is not `NULL`.
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
                              normalized = FALSE,
                              tolerance = max(self$edge_lengths)/2,
                              verbose = FALSE) {
    data_coords <- data_coords[[1]]
    if(data_coords == "euclidean"){
      lifecycle::deprecate_warn("1.1.2.9000", "add_observations(data_coords = 'must be either PtE or spatial')")
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

    data <- as.list(data)

    ## convert everything to PtE
    if(verbose){
      if(data_coords == "spatial" || !is.null(Spoints)){
      message("Converting data to PtE")
      if(private$longlat){
        message("This step may take long. If this step is taking too long consider pruning the vertices, and if it still takes long, consider setting 'project_data' to 'TRUE' to project the coordinates in the plane and obtain a significant speed up.")
      } 
      }
    }

    t <- system.time({
      if(!is.null(Spoints)){
        PtE <- self$coordinates(XY = Spoints@coords)
        XY_new <- self$coordinates(PtE = PtE, normalized = TRUE)
        norm_XY <- max(sqrt(rowSums( (Spoints@coords-XY_new)^2 )))
        if(norm_XY > tolerance){
          warning("There was at least one point whose location is far from the graph,
          please consider checking the input.")
        }
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
            norm_XY <- max(sqrt(rowSums( (point_coords-XY_new)^2 )))
            if(norm_XY > tolerance){
              warning("There was at least one point whose location is far from the graph,
                please consider checking the input.")
              }
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
       group_vector <- data[[group]]
     } else{
       group <- "__group"
       group_vector <- NULL
     }

    lapply(data, function(dat){if(nrow(matrix(PtE, ncol=2)) != length(dat)){
        stop(paste(dat,"has a different number of elements than the number of
                   coordinates!"))
       }})

    if(!is.null(self$data[["__group"]])){
      group_vals <- unique(self$data[["__group"]])
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
    data[[group]] <- NULL
    self$data[["__coord_x"]] <- NULL
    self$data[["__coord_y"]] <- NULL

    # Process the data (find all the different coordinates
    # across the different replicates, and also merge the new data to the old data)
    self$data <- process_data_add_obs(PtE, new_data = data, self$data,
                                        group_vector)

    ## convert to Spoints and add
    PtE <- self$get_PtE()
    spatial_points <- self$coordinates(PtE = PtE, normalized = TRUE)
    self$data[["__coord_x"]] <- rep(spatial_points[,1], times = n_group)
    self$data[["__coord_y"]] <- rep(spatial_points[,2], times = n_group)
    })
          if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
          }
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
  build_mesh = function(h=NULL,n=NULL) {

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

    self$mesh <- list(PtE = NULL,
                      V = NULL,
                      E = NULL,
                      n_e = rep(0, self$nV),
                      h_e = NULL,
                      ind = 1:self$nV,
                      VtE = NULL)

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
  },

  #' @description Build mass and stiffness matrices for given mesh object.
  #' @details The function builds: The matrix `C` which is the mass matrix with
  #' elements \eqn{C_{ij} = <\phi_i, \phi_j>}, the matrix `G` which is the stiffness
  #' matrix with elements \eqn{G_{ij} = <d\phi_i, d\phi_j>}, the matrix `B` with
  #' elements \eqn{B_{ij} = <d\phi_i, \phi_j>}, the matrix `D` with elements
  #' \eqn{D_{ij} = \sum_{v\in V}\phi_i(v)\phi_j(v)}, and the vector with weights
  #' \eqn{<\phi_i, 1>}.
  #' @return No return value. Called for its side effects. The finite element
  #' matrices `C`, `G` and `B` are stored in the `mesh` element in the
  #' `metric_graph` object.
  compute_fem = function() {
    if (is.null(self$mesh)) {
      stop("no mesh provided")
    }
    nV <- dim(self$mesh$V)[1]
    fem_temp <- assemble_fem(E = self$mesh$E, h_e = self$mesh$h_e, nV = nV)
    self$mesh$C <- fem_temp$C
    self$mesh$G <- fem_temp$G
    self$mesh$B <- fem_temp$B
    self$mesh$D <- Diagonal(dim(self$mesh$C)[1],
                            c(rep(1, self$nV), rep(0, dim(self$mesh$C)[1] - self$nV)))

    self$mesh$weights <- rowSums(self$mesh$C)
  },

  #' @description Computes observation matrix for mesh.
  #' @param PtE Locations given as (edge number in graph, normalized location on
  #' edge)
  #' @details For n locations and a mesh with m nodes, `A` is an n x m matrix with
  #' elements \eqn{A_{ij} = \phi_j(s_i)}{A_{ij} = \phi_j(s_i)}.
  #' @return The observation matrix.
  mesh_A = function(PtE) {
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
  #' @param mesh Plot the mesh locations?
  #' @param X Additional values to plot.
  #' @param X_loc Locations of the additional values in the format
  #' (edge, normalized distance on edge).
  #' @param p Existing objects obtained from 'ggplot2' or 'plotly' to add the graph to
  #' @param degree Show the degrees of the vertices?
  #' @param direction Show the direction of the edges?
  #' @param ... Additional arguments to pass to `ggplot()` or `plot_ly()`
  #' @return A `plot_ly` (if `plotly = TRUE`) or `ggplot` object.
  plot = function(data = NULL,
                  group = 1,
                  plotly = FALSE,
                  vertex_size = 3,
                  vertex_color = 'black',
                  edge_width = 0.3,
                  edge_color = 'black',
                  data_size = 1,
                  mesh = FALSE,
                  X = NULL,
                  X_loc = NULL,
                  p = NULL,
                  degree = FALSE,
                  direction = FALSE,
                  ...) {
    if(!is.null(data) && is.null(self$data)) {
      stop("The graph does not contain data.")
    }
    if(is.numeric(group) && !is.null(data)) {
      unique_group <- unique(self$data[["__group"]])
      group <- unique_group[group]
    }
    if(!plotly) {
      p <- private$plot_2d(line_width = edge_width,
                           marker_size = vertex_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
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
                           data_size = data_size,
                           group = group,
                           mesh = mesh,
                           X = X,
                           X_loc = X_loc,
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
        return(is_tree(g, mode = "all"))
  },

  #' @description Plots continuous function on the graph.
  #' @param X Either an m x 3 matrix with (edge number, position on
  #' curve (in length), value) or a vector with values for the function
  #' evaluated at the mesh in the graph
  #' @param plotly If `TRUE`, then the plot is shown in 3D. This option requires
  #' the package 'plotly'.
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
  plot_function = function(X,
                           plotly = FALSE,
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

    if(!is.matrix(X) || is.matrix(X) && dim(X)[1] == 1) {
      mesh <- TRUE
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

      if (length(X) == PtE_dim) {
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
    }

    if (mesh) {
      n.v <- dim(self$V)[1]
      XV <- X[1:n.v]
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
          vals <- rbind(c(0, XV[Vs]),
                        cbind(self$mesh$PtE[ind, 2], X[n.v + which(ind)]),
                        c(1, XV[Ve]))

        }
      } else {
        vals <- X[X[, 1]==i, 2:3, drop = FALSE]
        if (max(vals[, 1]) < 1) {
          #check if we can add end value from other edge
          Ei <- self$E[, 1] == Ve #edges that start in Ve
          if (sum(Ei) > 0) {
            ind <- which(X[Ei, 2] == 0)[1]
          } else {
            ind <- NULL
          }
          if (length(ind) > 0) {
            vals <- rbind(vals, c(1, X[ind, 3]))
          } else {
            Ei <- self$E[, 2] == Ve #edges that end in Ve
            if (sum(Ei)  > 0) {
              ind <- which(X[Ei, 2] == 1)[1]
            } else {
              ind <- NULL
            }
            if (length(ind) > 0){
              vals <- rbind(vals, c(1, X[ind, 3]))
            }
          }
        }
        if (min(vals[, 1] > 0)) {
          #check if we can add start value from other edge
          Ei <- self$E[, 1] == Vs #edges that start in Vs
          if (sum(Ei) > 0) {
            ind <- which(X[Ei, 2] == 0)[1]
          } else {
            ind <- NULL
          }
          if (length(ind) > 0) {
            vals <- rbind(c(0, X[ind, 3]), vals)
          } else {
            Ei <- self$E[, 2] == Vs #edges that end in Vs
            if (sum(Ei) > 0) {
              ind <- which(X[Ei, 2] == 1)[1]
            } else {
              ind <- NULL
            }
            if (length(ind) > 0) {
              vals <- rbind(c(0, X[ind, 3]), vals)
            }
          }
        }
      }

        data.to.plot <- vals

        data.to.plot.order <- data.to.plot[order(vals[, 1]), ,
                                           drop = FALSE]

        coords <- interpolate2(self$edges[[i]],
                               pos = data.to.plot.order[, 1, drop = TRUE],
                               normalized = TRUE)

        x.loc <- c(x.loc, coords[, 1])
        y.loc <- c(y.loc, coords[, 2])
        z.loc <- c(z.loc, data.to.plot.order[, 2])
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
    data[["__edge_number"]] <- PtE_mesh[,1]
    data[["__distance_on_edge"]] <- PtE_mesh[,2]
    self$add_observations(data = data, group = group,
                          edge_number = "__edge_number",
                          distance_on_edge = "__distance_on_edge",
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

  #' @description Get the observation/prediction matrix A
  #' @param group A vector. If `NULL`, the A matrix for the first group will be
  #' returned. One can use all groups by simply setting the `group` variable
  #' to `__all`. Otherwise, the A matrix for the groups in the vector will be
  #' returned.
  #' @param obs_to_vert Should the observations be turned into vertices?
  #' @param include_NA Should the locations for which all observations are NA be
  #' included?
  #' @return The observation or prediction matrix.
  A = function(group = NULL,
               obs_to_vert = FALSE,
               include_NA = TRUE){

    if(is.null(self$PtV) && !obs_to_vert){
        stop("The A matrix was not computed. If you want to compute rerun this
             method with 'obs_to_vertex=TRUE', in which the observations will be
             turned to vertices and the A matrix will then be computed")
    } else if(is.null(self$PtV)){
      self$observation_to_vertex(mesh_warning=FALSE)
    }

    if(is.null(group)){
      group <- unique(self$data[["__group"]])
      group <- group[1]
    } else if (group[1] == "__all"){
      group <- unique(self$data[["__group"]])
    }
    n_group <- length(unique(group))

    if(include_NA){
      A <- Matrix::Diagonal(self$nV)[self$PtV, ]
      return(Matrix::kronecker(Diagonal(n_group),A))
    } else{
      if(length(group) == 1){
        A <- Matrix::Diagonal(self$nV)[self$PtV, ]
        return(A)
      } else {
        data_group <- select_group(self$data, group[1])
        idx_notna <- idx_not_all_NA(data_group)
        nV_tmp <- sum(idx_notna)
        A <- Matrix::Diagonal(nV_tmp)[self$PtV[idx_notna], ]
        for (i in 2:length(group)) {
          data_group <- select_group(self$data, group[i])
          idx_notna <- idx_not_all_NA(data_group)
          nV_tmp <- sum(idx_notna)
          A <- bdiag(A, Matrix::Diagonal(nV_tmp)[self$PtV[idx_notna], ])
        }
        return(A)
      }
    }
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

      if (self$mesh$E[e, 1] == v1) { #edge starts in the vertex before
        d <- (PtE[i, 2] - d1)/(d2 - d1)
      } else {
        d <- 1- (PtE[i, 2] - d1)/(d2 - d1)
      }
      PtE_update[i, ] <- c(e, d)
    }
    return(PtE_update)
  },

  ## Adding column_y argument which tells which column of y to get

  plot_2d = function(line_width = 0.1,
                     marker_size = 1,
                     vertex_color = 'black',
                     edge_color = 'black',
                     data,
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
      data_group <- select_group(self$data, group)
      y_plot <-data_group[[data]]
      points_xy <- self$coordinates(PtE = self$get_PtE())
      x <- points_xy[,1]
      y <- points_xy[,2]

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
                     data_size = 1,
                     group = 1,
                     mesh = FALSE,
                     p = NULL,
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
      data_group <- select_group(self$data, group)
      y_plot <- self$data_group[[data]]
      PtE <- self$get_PtE()
      points_xy <- self$coordinates(PtE = PtE)

      x <- points_xy[,1]
      y <- points_xy[,2]

      data.plot <- data.frame(x = x[!is.na(as.vector(y_plot))],
                              y = y[!is.na(as.vector(y_plot))],
                              z = rep(0,length(x[!is.na(as.vector(y_plot))])),
                              val = as.vector(y_plot[!is.na(as.vector(y_plot))]))
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           type = "scatter3d", mode = "markers",
                           marker = list(size = marker_size,
                                         color = ~val,
                                         colorbar=list(title='', len = 0.5),
                                         colorscale='Viridis'),
                           showlegend=FALSE)
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
            self$edges[[k1]] <- rbind(self$V[v.rem,], self$edges[[k1]])
          }
          for(k2 in idx_E2){
            self$edges[[k2]] <- rbind(self$edges[[k2]],self$V[v.rem,])
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
          v.degrees <- self$get_degrees()[v.loop]

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
      if(verbose){
        message("Recomputing edge lengths")
      }
      t <- system.time({
        self$edge_lengths <- private$compute_lengths(longlat, private$length_unit, crs, proj4string, which_longlat, private$vertex_unit, project_data)
      })
       if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
       }     
    }
  },

  # utility function to merge lines connected by degree 2 vertices
  merge.all.deg2 = function() {
   while(sum(self$get_degrees()==2)>0) {
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

        self$edges[[e_rem[1]]] <- coords

        self$edges <- self$edges[-e_rem[2]]
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
      self$nV <- self$nV - 1

      #update edges
      self$E[self$E >= ind] <- self$E[self$E >= ind] - 1
      self$E <- self$E[-e_rem[2],,drop=FALSE]
      self$E[e_rem[1],] <- E_new

      self$edge_lengths[e_rem[1]] <- self$edge_lengths[e_rem[1]] + self$edge_lengths[e_rem[2]]
      self$edge_lengths <- self$edge_lengths[-e_rem[2]]

      self$nE <- self$nE - 1
      } 
    }
    return(res.out)
  },

  # Compute lengths

  compute_lengths = function(longlat, unit, crs, proj4string, which_longlat, vertex_unit, project_data){
          ll <- sapply(self$edges, 
          function(edge){compute_line_lengths(edge, longlat = longlat, unit = unit, crs = crs, proj4string, which_longlat, vertex_unit, project_data)})
          return(ll)
  },


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

        if(!is.null(self$data)){
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
    Laplacian[graph.temp$data[["__dummy"]], graph.temp$data[["__dummy"]]] <- Laplacian

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

  # Temp PtE

  temp_PtE = NULL,

  # # longlat
  # longlat = NULL,

  # tolerance

  tolerance = NULL

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
           lifecycle::deprecate_warn("1.1.2.9000", "mgraph_components$new(lines)", "graph_components$new(edges)",
             details = c("`lines` was provided but not `edges`. Setting `edges <- lines`.")
           )
           edges <- lines
         } else {
           lifecycle::deprecate_warn("1.1.2.9000", "graph_components$new(lines)", "graph_components$new(edges)",
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

       units(self$lengths) <- units(self$graphs[[1]]$get_edge_lengths())

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
