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

  #' @field DirectionalWeightFunction_in Function for inwards weights in directional models
 DirectionalWeightFunction_in  =NULL,
  #' @field DirectionalWeightFunction_out Function for outwards weights in directional models
 DirectionalWeightFunction_out  = NULL,

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
  #' @param edges A list containing coordinates as `m x 2` matrices (that is, of `matrix` type) or m x 2 data frames (`data.frame` type) of sequence of points connected by straightlines. Alternatively, you can also prove an object of type `SSN`, `osmdata_sp`, `osmdata_sf`, `SpatialLinesDataFrame` or `SpatialLines` (from `sp` package) or `MULTILINESTRING` (from `sf` package).
  #' @param V n x 2 matrix with Euclidean coordinates of the n vertices. If non-NULL, no merges will be performed.
  #' @param E m x 2 matrix where each row represents one of the m edges. If non-NULL, no merges will be performed.
  #' @param vertex_unit The unit in which the vertices are specified. The options are 'degree' (the great circle distance in km), 'km', 'm' and 'miles'. The default is `NULL`, which means no unit. However, if you set `length_unit`, you need to set `vertex_unit`.
  #' @param length_unit The unit in which the lengths will be computed. The options are 'km', 'm' and 'miles'. The default, when longlat is `TRUE`, or an `sf` or `sp` objects are provided, is 'km'.
  #' @param edge_weights Either a number, a numerical vector with length given by the number of edges, providing the edge weights, or a `data.frame` with the number of rows being equal to the number of edges, where
  #' each row gives a vector of weights to its corresponding edge. Can be changed by using the `set_edge_weights()` method.
  #' @param kirchhoff_weights If non-null, the name (or number) of the column of `edge_weights` that contain the Kirchhoff weights. Must be equal to 1 (or `TRUE`) in case `edge_weights` is a single number and those are the Kirchhoff weights.
  #' @param directional_weights If non-null, the name (or number) of the column of `edge_weights` that contain the directional weights. The default is the first column of the edge weights.
  #' @param longlat There are three options: `NULL`, `TRUE` or `FALSE`. If `NULL` (the default option), the `edges` argument will be checked to see if there is a CRS or proj4string available, if so, `longlat` will be set to `TRUE`, otherwise, it will be set to `FALSE`. If `TRUE`, then it is assumed that the coordinates are given.
  #' in Longitude/Latitude and that distances should be computed in meters. If `TRUE` it takes precedence over
  #' `vertex_unit` and `length_unit`, and is equivalent to `vertex_unit = 'degree'` and `length_unit = 'm'`.
  #' @param include_obs If the object is of class `SSN`, should the observations be added? If `NULL` and the edges are of class `SSN`, the data will be automatically added. If `FALSE`, the data will not be added. Alternatively, one can set this argument to the numbers or names of the columns of the observations to be added as observations.
  #' @param add_obs_options List containing additional options to be passed to the `add_observations()` method when adding observations from `SSN` data?
  #' @param include_edge_weights If the object is of class `SSN`, `osmdata_sp`, `osmdata_sf`, `SpatialLinesDataFrame`, `MULTILINESTRING`, `LINESTRING`, `sfc_LINESTRING`, `sfc_MULTILINESTRING`, should the edge data (if any) be added as edge weights? If `NULL`, the edge data will be added as edge weights, if `FALSE` they will not be added. Alternatively, one can set this argument to the numbers or names of the columns of the edge data to be added as edge weights.
  #' @param crs Coordinate reference system to be used in case `longlat` is set to `TRUE` and `which_longlat` is `sf`. Object of class crs. The default choice, if the `edges` object does not have CRS nor proj4string, is `sf::st_crs(4326)`.
  #' @param proj4string Projection string of class CRS-class to be used in case `longlat` is set to `TRUE` and `which_longlat` is `sp`. The default choice, if the `edges` object does not have CRS nor proj4string, is `sp::CRS("+proj=longlat +datum=WGS84")`.
  #' @param which_longlat Compute the distance using which package? The options are `sp` and `sf`. The default is `sp`.
  #' @param project If `longlat` is `TRUE` should a projection be used to compute the distances to be used for the tolerances (see `tolerance` below)? The default is `FALSE`. When `TRUE`, the construction of the graph is faster.
  #' @param project_data If `longlat` is `TRUE` should the vertices be project to planar coordinates? The default is `FALSE`. When `TRUE`, the construction of the graph is faster.
  #' @param which_projection Which projection should be used in case `project` is `TRUE`? The options are `Robinson`, `Winkel tripel` or a proj4string. The default is `Winkel tripel`.
  #' @param manual_edge_lengths If non-NULL, a vector containing the edges lengths, and all the quantities related to edge lengths will be computed in terms of these. If merges are performed, it is likely that the merges will override the manual edge lengths. In such a case, to provide manual edge lengths, one should either set the `perform_merges` argument to `FALSE` or use the `set_manual_edge_lengths()` method.
  #' @param perform_merges There are three options, `NULL`, `TRUE` or `FALSE`. The default option is `NULL`. If `NULL`, it will be set to `FALSE` unless 'edges', 'V' and 'E' are `NULL`, in which case it will be set to `TRUE`. If FALSE, this will take priority over the other arguments, and no merges (except the optional `merge_close_vertices` below) will be performed. Note that the merge on the additional `merge_close_vertices` might still be performed, if it is set to `TRUE`.
  #' @param approx_edge_PtE Should the relative positions on the edges be approximated? The default is `TRUE`. If `FALSE`, the speed can be considerably slower, especially for large metric graphs.
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
  #' @param merge_close_vertices Should an additional step to merge close vertices be done? The options are `NULL` (the default), `TRUE` or `FALSE`. If `NULL`, it will be determined automatically. If `TRUE` this step will be performed even if `perfom_merges` is set to `FALSE`.
  #' @param factor_merge_close_vertices Which factor to be multiplied by tolerance `vertex_vertex` when merging close vertices at the additional step?
  #' @param remove_circles All circlular edges with a length smaller than this number
  #' are removed. If `TRUE`, the `vertex_vertex` tolerance will be used. If `FALSE`, no circles will be removed.
  #' @param auto_remove_point_edges Should edges of length zero, that is, edges that are actually points, be automatically removed? 
  #' @param verbose Print progress of graph creation. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
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
                        tolerance = list(vertex_vertex = 1e-3,
                                         vertex_edge = 1e-3,
                                         edge_edge = 0),
                        check_connected = TRUE,
                        remove_deg2 = FALSE,
                        merge_close_vertices = NULL,
                        factor_merge_close_vertices = 1,
                        remove_circles = FALSE,
                        auto_remove_point_edges = TRUE,                        
                        verbose = 1,
                        add_obs_options = list(return_removed = FALSE,
                                                verbose = verbose),
                        lines = deprecated()) {

      start_construction_time <- Sys.time()

      add_data_tmp <- FALSE

      private$project_data <- project_data

      if(is.null(merge_close_vertices)){
          merge_close_vertices <- TRUE
      }

      if(inherits(edges, "SSN")){
        if(is.null(include_obs) || include_obs[[1]] == TRUE){
          dataset_tmp <- edges$obs
          add_data_tmp <- TRUE
        } else if(is.numeric(include_obs) || is.character(include_obs)){
          dataset_tmp <- edges$obs[,include_obs]
          add_data_tmp <- TRUE
        } else if(!is.logical(include_obs[[1]])){
          stop("invalid option passed to include_obs.")
        }
        edges <- edges$edges        
        if(is.null(include_edge_weights) || include_edge_weights[[1]] == TRUE){
          edge_weights <- sf::st_drop_geometry(edges)
        } else if(is.numeric(include_edge_weights) || is.character(include_edge_weights)){
          edge_weights <- sf::st_drop_geometry(edges)          
          edge_weights <- edge_weights[,include_edge_weights]
        } else if(!is.logical(include_edge_weights[[1]])){
          stop("invalid option passed to include_edge_weights.")
        }
      }

      if(inherits(edges, c("osmdata_sp", "osmdata_sf"))){
        edges <- edges$osm_lines
        if(is.null(include_edge_weights) || include_edge_weights[[1]] == TRUE){
          if(inherits(edges, "osmdata_sp")){
            edge_weights <- edges@data
          } else{
            edge_weights <- sf::st_drop_geometry(edges)
          }
        } else if(is.numeric(include_edge_weights) || is.character(include_edge_weights)){
          if(inherits(edges, "osmdata_sp")){
            edge_weights <- edges@data
          } else{
            edge_weights <- sf::st_drop_geometry(edges)
          }          
          edge_weights <- edge_weights[,include_edge_weights]
        } else if(!is.logical(include_edge_weights[[1]])){
          stop("invalid option passed to include_edge_weights.")
        }        
      }

      if (inherits(edges,"SpatialLines") || inherits(edges,"SpatialLinesDataFrame")) {
        if(is.null(longlat) || longlat){
          if(!is.na(sp::proj4string(edges))){
            longlat <- TRUE
            if(is.null(proj4string) && is.null(crs)){
              proj4string <- sp::proj4string(edges)
              crs_tmp <- sf::st_crs(proj4string, parameters = TRUE)
              if(is.null(vertex_unit)){
                vertex_unit <- crs_tmp$units_gdal
                if(vertex_unit == "metre"){
                  vertex_unit <- "m"
                } else if(vertex_unit == "kilometre"){
                  vertex_unit <- "km"
                }
              }
            }
            if(is.null(length_unit)){
              length_unit <- "km"          
            }
          } else{
            longlat <- FALSE
          }
        }
        if(inherits(edges,"SpatialLinesDataFrame")){
          if(is.null(include_edge_weights) || include_edge_weights[[1]] == TRUE){
            edge_weights <- edges@data
        } else if(is.numeric(include_edge_weights) || is.character(include_edge_weights)){
            edge_weights <- edges@data      
            edge_weights <- edge_weights[,include_edge_weights]
        } else if(!is.logical(include_edge_weights[[1]])){
          stop("invalid option passed to include_edge_weights.")
        }        
        }
      } else if(inherits(edges, c("MULTILINESTRING", "LINESTRING", "sfc_LINESTRING", "sfc_MULTILINESTRING", "sf"))){
        if(is.null(longlat) || longlat){
          if(!is.na(sf::st_crs(edges))){
            longlat <- TRUE
            if(is.null(proj4string) && is.null(crs)){
              crs <- sf::st_crs(edges)
              crs_tmp <- sf::st_crs(edges, parameters = TRUE)
              if(is.null(vertex_unit)){
                vertex_unit <- crs_tmp$units_gdal
                if(vertex_unit == "metre"){
                  vertex_unit <- "m"
                } else if(vertex_unit == "kilometre"){
                  vertex_unit <- "km"
                }
              }
            }
            if(is.null(length_unit)){
              length_unit <- "km"              
            }
          } else{
            longlat <- FALSE
          }
        }

        if(inherits(edges, "data.frame")){
        if(is.null(include_edge_weights) || include_edge_weights[[1]] == TRUE){
            edge_weights <- sf::st_drop_geometry(edges)
        } else if(is.numeric(include_edge_weights) || is.character(include_edge_weights)){
            edge_weights <- sf::st_drop_geometry(edges)    
          edge_weights <- edge_weights[,include_edge_weights]
        } else if(!is.logical(include_edge_weights[[1]])){
          stop("invalid option passed to include_edge_weights.")
        }        

        }


      } else{
        if(is.null(longlat)){
          longlat <- FALSE
        }
      }

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

       if(is.null(edge_weights)){
          edge_weights <- 1
       }

       if(!is.null(kirchhoff_weights)){
        if(length(kirchhoff_weights)>1){
          warning("Only the first entry of 'kirchhoff_weights' was used.")
          kirchhoff_weights <- kirchhoff_weights[[1]]
        }
        if(!is.numeric(kirchhoff_weights)){
          if(!is.character(kirchhoff_weights)){
            stop("'kirchhoff_weights' must be either a number of a string.")
          }
          if(!(kirchhoff_weights%in%colnames(edge_weights))){
            stop(paste(kirchhoff_weights, "is not a column of 'edge_weights'!"))
          }
        } else{
          if(!is.data.frame(edge_weights)){
            if(kirchhoff_weights != 1){
              stop("Since 'edge_weights' is not a data.frame, 'kirchhoff_weights' must be either NULL or 1.")
            }
          } else{
            if(kirchhoff_weights %%1 != 0){
              stop("'kirchhoff_weights' must be an integer.")
            }
            if((kirchhoff_weights < 1) || (kirchhoff_weights > ncol(edge_weights))){
              stop("'kirchhoff_weights' must be a positive integer number smaller or equal to the number of columns of 'edge_weights'.")
            }
          }
        }
        private$kirchhoff_weights <- kirchhoff_weights
       } else{
        if(is.vector(edge_weights)){
          private$kirchhoff_weights <- 1
        } else{
          if(".weights" %in% colnames(edge_weights)){
            private$kirchhoff_weights <- ".weights"
          } else{
            edge_weights[, ".weights"] <- 1
            private$kirchhoff_weights <- ".weights"
          }
        }        
       }



      if(!is.null(directional_weights)){
        if(!is.numeric(directional_weights)){
          if(!is.character(directional_weights)){
            stop("'directional_weights' must be either a numerical vector of a string vector.")
          }
          if(!(directional_weights%in%colnames(edge_weights))){
            stop(paste(directional_weights, "is not a column of 'weights'!"))
          }
        } else{
          if(!is.data.frame(edge_weights)){
            if(directional_weights != 1){
              stop("Since 'weights' is not a data.frame, 'directional_weights' must be either NULL or 1.")
            }
          } else{
            if(directional_weights %%1 != 0){
              stop("'directional_weights' must be an integer.")
            }
            if((directional_weights < 1) || (directional_weights > ncol(edge_weights))){
              stop("'directional_weights' must be a positive integer number smaller or equal to the number of columns of 'weights'.")
            }
          }
        }
        private$directional_weights <- directional_weights
       } else{
        if(is.vector(edge_weights)){
          private$directional_weights <- 1
        } else{
          if(".weights" %in% colnames(edge_weights)){
            private$directional_weights <- ".weights"
          } else{
            edge_weights[, ".weights"] <- 1
            private$directional_weights <- ".weights"
          }
        }
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

      valid_units_vertex <- c("m", "km", "miles", "degree")
      valid_units_length <- c("m", "km", "miles")

      if(!(which_longlat %in% c("sp", "sf"))){
        stop("The options for 'which_longlat' are 'sp' and 'sf'!")
      }

      if(longlat){
        private$longlat <- TRUE
        private$which_longlat <- which_longlat
      }

      if(!is.null(proj4string)){
        if(!longlat){
          warning("proj4string was passed, so setting longlat to TRUE")
          longlat <- TRUE
          private$longlat <- TRUE
          private$which_longlat <- which_longlat
        }
        private$crs <- sf::st_crs(proj4string)
        private$proj4string <- proj4string
        crs <- private$crs
        private$transform <- !(sf::st_is_longlat(private$crs))
      }

      if(!is.null(crs)){
        if(!longlat){
          warning("crs was passed, so setting longlat to TRUE")
          longlat <- TRUE
          private$longlat <- TRUE
          private$which_longlat <- which_longlat
        }
        private$crs <- sf::st_crs(crs)
        private$proj4string <- sp::CRS(private$crs$proj4string)
        proj4string <- private$proj4string
        private$transform <- !(sf::st_is_longlat(private$crs))
      }

      if(longlat && (which_longlat == "sp") && is.null(proj4string)){
        proj4string <- sp::CRS("+proj=longlat +datum=WGS84")
        private$crs <- sf::st_crs(proj4string)
        private$proj4string <- proj4string
        private$transform <- !(sf::st_is_longlat(private$crs))
      }

      if(longlat && (which_longlat == "sf") && is.null(crs)){
        crs <- sf::st_crs(4326)
        private$crs <- crs
        private$transform <- !(sf::st_is_longlat(private$crs))
      }

    # private$longlat <- longlat

    if((is.null(vertex_unit) && !is.null(length_unit)) || (is.null(length_unit) && !is.null(vertex_unit))){
      stop("If one of 'vertex_unit' or 'length_unit' is NULL, and the edges are not sf nor sp objects, then the other must also be NULL.")
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
      if(length_unit == "degree"){
        length_unit <- "km"
      }
      if(!(length_unit %in% valid_units_length)){
        stop(paste("The possible options for 'length_unit' are ", toString(valid_units_length)))
      }
      private$length_unit <- length_unit
    }

    if(longlat){
      private$vertex_unit <- vertex_unit
      if(is.null(vertex_unit)){
        private$vertex_unit <- "degree"
      }
      if(!is.null(length_unit)){
        private$length_unit <- length_unit
      } else{
        private$length_unit <- "km"
      }
    } else if(!is.null(vertex_unit)){
        if(private$vertex_unit == "degree"){
          longlat <- TRUE
          private$longlat <- TRUE
        }
    }

    if(is.null(edges) && is.null(V) && is.null(E)) {
      edges <- logo_lines()
        if(is.null(perform_merges)){
          perform_merges <- TRUE
          remove_circles <- TRUE
        }      
    }

    if(!is.null(manual_edge_lengths)){
        if(is.null(perform_merges)){
          perform_merges <- FALSE
        }
      } else{
        if(is.null(perform_merges)){
          perform_merges <- FALSE
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

    if(verbose > 0){
      message("Starting graph creation...")
      message(paste("LongLat is set to",longlat))
      if(longlat){
        message(paste("The unit for edge lengths is", private$length_unit))
        if(perform_merges){
          message(paste0("The current tolerances (in ",private$length_unit,") are:"))
          message(paste("\t Vertex-Vertex", tolerance$vertex_vertex))
          message(paste("\t Vertex-Edge", tolerance$vertex_edge))
          message(paste("\t Edge-Edge", tolerance$edge_edge))
        }
      } else if (perform_merges){
        message("The current tolerances are:")
        message(paste("\t Vertex-Vertex", tolerance$vertex_vertex))
        message(paste("\t Vertex-Edge", tolerance$vertex_edge))
        message(paste("\t Edge-Edge", tolerance$edge_edge))
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


    if(!(perform_merges %in% c(TRUE,FALSE))){
      stop("perform_merges should be either TRUE or FALSE.")
    }

    private$perform_merges <- perform_merges

    if(verbose > 0){
      message("Creating edges...")
    }

    if(!is.null(edges)){
      if(!is.null(V) || !is.null(E)){
        warning("object initialized from edges, then E and V are ignored")
      }
      if (inherits(edges,"SpatialLinesDataFrame")) {
        tmp_lines = sp::SpatialLines(edges@lines)
        self$edges <- lapply(1:length(tmp_lines), function(i){tmp_lines@lines[[i]]@Lines[[1]]@coords})
      } else if (inherits(edges,"SpatialLines")) {
        self$edges = lapply(1:length(edges), function(i){edges@lines[[i]]@Lines[[1]]@coords})
      } else if (inherits(edges, c("MULTILINESTRING", "LINESTRING", "sfc_LINESTRING", "sfc_MULTILINESTRING", "sf"))) {
        # Ensure 'edges' is an 'sf' object if it is not already
        if (!inherits(edges, "sf")) {
          # Convert to 'sfc' (simple feature geometry list column)
          edges <- sf::st_sfc(edges)
          edges <- sf::st_sf(geometry = edges)  # Wrap in an 'sf' data frame
        }

        # Define valid geometry types for filtering
        valid_types <- c("LINESTRING", "MULTILINESTRING")

        # Filter for valid geometry types in 'edges'
        valid_indices <- sf::st_geometry_type(edges) %in% valid_types
        valid_edges <- edges[valid_indices, , drop = FALSE]  # Filter only valid geometries

        # Extract coordinates for the valid edges
        coords_multilinestring <- sf::st_coordinates(sf::st_geometry(valid_edges))

        split_coords <- split(coords_multilinestring[, 1:2, drop = FALSE], coords_multilinestring[, "L1"])
        self$edges <- lapply(split_coords, function(coords) matrix(coords, ncol=2, byrow=FALSE))        
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
      self$E <- E
      self$V <- V
      private$perform_merges <- FALSE
    }

    self$edges <- lapply(self$edges, function(x) if(nrow(x) < 2) NULL else x)

    self$nE <- length(self$edges)

    if(verbose > 0){
      message("Setting edge weights...")
    }

    private$set_first_weights(weights = edge_weights)


    if(verbose > 0){
      message("Computing bounding box...")
    }
    private$compute_bounding_box()

    if(!is.null(manual_edge_lengths)){
      self$set_manual_edge_lengths(edge_lengths = manual_edge_lengths, unit = length_unit)
    }

    if(private$perform_merges){

    if(verbose > 0){
      message("Setup edges and merge close vertices")
    }      

      if(!is.null(manual_edge_lengths)){
        warning("Since 'perform_merges' is TRUE, the manual edge lengths will not be used. Either set 'perform_merges' to FALSE or use the 'set_manual_edge_lengths()' method on the graph after the graph construction.")
      }

      t <- system.time(
        private$line_to_vertex(tolerance = tolerance$vertex_vertex,
                           longlat = private$longlat, factor_unit, verbose=verbose,
                           private$crs, private$proj4string, which_longlat, private$length_unit, private$vertex_unit,
                           project, which_projection, project_data)
        )

      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    
      self$compute_PtE_edges(approx = approx_edge_PtE, verbose=verbose)

      if(length(self$edges) > 1){

      if (tolerance$edge_edge > 0) {
      private$addinfo <- TRUE

      if(verbose > 0){
        message("Find edge-edge intersections")
      }

      t <- system.time(
        points_add <- private$find_edge_edge_points(tol = tolerance$edge_edge, verbose=verbose,
        crs=private$crs, proj4string = private$proj4string, longlat=private$longlat, fact = factor_unit, which_longlat = which_longlat)
        )

      if(verbose == 2){
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
        if(verbose == 2){
          message(sprintf("Add %d new vertices", nrow(PtE)))
        }

        PtE <- na.omit(PtE)

        t <- system.time(
        private$add_vertices(PtE, tolerance = tolerance$edge_edge, verbose = verbose)
        )

        if(verbose == 2){
          message(sprintf("time: %.3f s", t[["elapsed"]]))
        }
      }

      private$clear_initial_info()
      }

      if(tolerance$vertex_edge > 0){
        private$addinfo <- TRUE
        if(verbose > 0){
          message("Snap vertices to close edges")
        }

        t <- system.time(
          PtE_tmp <- private$coordinates_multiple_snaps(XY = self$V,
                                              tolerance = tolerance$vertex_edge, verbose = verbose,
          crs=private$crs, proj4string = private$proj4string, longlat=private$longlat, fact = factor_unit, which_longlat = which_longlat)
          )

        if(verbose == 2){
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
          if(verbose == 2){
            message(sprintf("Add %d new vertices", nrow(PtE_tmp)))
          }

          PtE_tmp <- na.omit(PtE_tmp)

          t <- system.time(
            private$add_vertices(PtE_tmp, tolerance = tolerance$vertex_edge, verbose=verbose)
            )

          if(verbose == 2){
            message(sprintf("time: %.3f s", t[["elapsed"]]))
          }
        }
        private$clear_initial_info()
      }

      if(merge_close_vertices){
        private$merge_close_vertices(factor_merge_close_vertices * tolerance$vertex_vertex, factor_unit)
      }

      if(auto_remove_point_edges){
          private$remove_circles(1e-15, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
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
        if(verbose == 2){
          message("Recomputing edge lengths")
        }
        t <- system.time({
          self$edge_lengths <- private$compute_lengths(private$longlat, private$length_unit, private$crs, private$proj4string, private$which_longlat, private$vertex_unit, project_data,private$transform)
        })
       if(verbose == 2){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
         }
      }
      # End of cond of having more than 1 edge
      }

      # Cleaning the edges

      if(verbose == 2){
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

      if(verbose == 2){
            message(sprintf("time: %.3f s", t[["elapsed"]]))
      }

      # Checking if there is some edge with infinite length
      if(any(!is.finite(self$edge_lengths))){
        warning("There is at least one edge of infinite length. Please, consider redefining the graph.")
      }

        # end of if do merges
    } else{

      if (verbose > 0) {
        message("Setting up edges")
      }

      # Extract the start and end vertices of all edges at once
      edges_vertices <- do.call(rbind, lapply(self$edges, function(edge) {
        edge[c(1, nrow(edge)), ]
      }))

      # Round the vertex coordinates and ensure uniqueness
      self$V <- unique(round(edges_vertices * 10^15) / 10^15)
      self$nV <- nrow(self$V)

      lvl <- matrix(0, nrow = length(self$edges), ncol = 2)

      # Collect all points and split them into start and end points
      all_points <- do.call(rbind, self$edges)
      n_points <- sapply(self$edges, nrow)

      # Split all points into start and end points
      start_indices <- cumsum(c(1, n_points[-length(n_points)]))
      end_indices <- cumsum(n_points)
      start_points <- all_points[start_indices, , drop = FALSE]
      end_points <- all_points[end_indices, , drop = FALSE]

      # Use RANN for efficient nearest-neighbor search
      nn_start <- nn2(self$V, start_points, k = 1)
      nn_end <- nn2(self$V, end_points, k = 1)

      # Extract the indices of the closest vertices
      lvl[, 1] <- nn_start$nn.idx
      lvl[, 2] <- nn_end$nn.idx

      self$E <- lvl
  
      if (merge_close_vertices) {
        if (verbose > 0) {
          message("Merging close vertices")
        }
        private$merge_close_vertices(factor_merge_close_vertices * tolerance$vertex_vertex, factor_unit)
      }
    }

    if(!private$perform_merges && is.null(manual_edge_lengths)){
        if(verbose == 2){
          message("Computing edge lengths")
        }
        t <- system.time({
          self$edge_lengths <- private$compute_lengths(private$longlat, private$length_unit, private$crs, private$proj4string, private$which_longlat, private$vertex_unit, project_data,private$transform)
        })
       if(verbose == 2){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
         }
    }


      if(auto_remove_point_edges){
          private$remove_circles(1e-15, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
      }

    if(!private$perform_merges){
      if(is.logical(remove_circles)){
        if(remove_circles){
          private$remove_circles(tolerance$vertex_vertex, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
        }
      } else {
        private$remove_circles(remove_circles, verbose=verbose,longlat = private$longlat, unit=length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=which_longlat, vertex_unit=vertex_unit, project_data)
        remove_circles <- TRUE
      }
    }

    end_construction_time <- Sys.time()
    construction_time <- end_construction_time - start_construction_time

    if(verbose > 0){
      message(sprintf('Total construction time: %.2f %s', construction_time, units(construction_time)))
    }

    # Checking if graph is connected
    if (check_connected) {
      g <- make_graph(edges = c(t(self$E)), directed = FALSE)
      # components <- igraph::clusters(g, mode="weak")
      components <- igraph::components(g, mode="weak")
      nc <- components$no
      if(nc>1){
        message("The graph is disconnected. You can either use the function 'graph_components' to obtain the different connected components or set 'perform_merges' to 'TRUE' and adjust the 'tolerances' to create a single connected graph.")
        private$connected = FALSE
      }
    }
    # creating/updating reference edges
    private$ref_edges <- map_into_reference_edge(self, verbose=verbose)

    # Cloning the initial graph

    if(verbose > 0){
      message("Creating and updating vertices...")
    }

   private$create_update_vertices(verbose=verbose)

    if(verbose > 0){
      message("Storing the initial graph...")
    }
    
    self$compute_PtE_edges(approx = approx_edge_PtE, verbose=verbose)

    private$initial_graph <- self$clone()

    # Cloning again to add the initial graph to the initial graph
    private$initial_graph <- self$clone()
    
    self$set_edge_weights(weights = private$edge_weights, kirchhoff_weights = private$kirchhoff_weights, directional_weights = private$directional_weights, verbose=verbose)
   
    self$setDirectionalWeightFunction()

    if (remove_deg2) {
      if (verbose > 0) {
        message("Remove degree 2 vertices")
      }
      t <- system.time(
        self$prune_vertices(verbose = verbose)
      )
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    }
    # Adding IDs to edges and setting up their class

    # for(i in 1:length(self$edges)){
    #   attr(self$edges[[i]], "id") <- i
    #   class(self$edges[[i]]) <- "metric_graph_edge"
    # }

    if(add_data_tmp){
      if(verbose > 0){
        message("Adding observations...")
      }
      add_obs_options[["data"]] <- dataset_tmp
      do.call(self$add_observations, add_obs_options)
    }

    # Checking if there is some edge with zero length
    if(any(self$edge_lengths == 0)){
        warning("There is at least one edge of length zero. Please, consider redefining the graph.")
    }

  },

  #' @description Sets the edge weights
  #' @param tolerance Tolerance at which circles with length less than this will be removed.
#' @param verbose Print progress of graph creation. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @return No return value. Called for its side effects.

  remove_small_circles = function(tolerance, verbose = 1){
    private$remove_circles(tolerance, verbose=verbose,longlat = private$longlat, unit=private$length_unit, crs=private$crs, proj4string=private$proj4string, which_longlat=private$which_longlat, vertex_unit=private$vertex_unit, project_data = private$project_data)
    private$create_update_vertices(verbose=verbose)
  },

  #' @description Exports the edges of the MetricGraph object as an `sf` or `sp`.
  #' @param format The format for the exported object. The options are `sf` (default), `sp` and `list`.
  #' @return 
#' For `format == "sf"`, the function returns an `sf` object of `LINESTRING` geometries, where the associated data frame includes edge weights. 
#' 
#' For `format == "sp"`, the function returns a `SpatialLinesDataFrame` where the data frame includes edge weights. 

  get_edges = function(format = c("sf", "sp", "list")){
    format <- format[[1]]
    if(!(format%in%c("sf", "sp", "list"))){
      stop("'format' must be one of the following: 'sf', 'sp' or 'list'.")
    }
      if(format == "sf"){
        edges_geometries <- lapply(self$edges, sf::st_linestring) 
      
        if(is.vector(private$edge_weights)){
          ew_tmp <- data.frame(.weights = private$edge_weights)
        } else{
          ew_tmp <- as.data.frame(private$edge_weights)
        }

        ew_tmp[[".edge_lengths"]] <- self$edge_lengths

        edges_sf <- sf::st_sf(ew_tmp, geometry = sf::st_sfc(edges_geometries), crs = if(!is.null(private$crs)) private$crs else sf::NA_crs_)
        return(edges_sf)
      } else if(format == "sp"){
             edges_list <- lapply(1:length(self$edges), function(i) {
        sp::Line(coords = matrix(self$edges[[i]], nrow = dim(self$edges[[i]])[1], ncol = dim(self$edges[[i]])[2]))
      })
      sp_edges <- sp::SpatialLines(
          lapply(1:length(edges_list), function(i) sp::Lines(list(edges_list[[i]]), ID = as.character(i))),
          proj4string = if (!is.null(private$crs)) private$proj4string else sp::CRS(NA_character_)
      )
      if(is.vector(private$edge_weights)){
        ew_tmp <- data.frame(.weights = private$edge_weights)
      } else{
        ew_tmp <- as.data.frame(private$edge_weights)
      }
      ew_tmp[[".ID"]] <- as.character(1:length(sp_edges))
      ew_tmp[[".edge_lengths"]] <- self$edge_lengths
      edges_sldf <- sp::SpatialLinesDataFrame(sp_edges, data = ew_tmp, match.ID = ".ID")
      return(edges_sldf)
      } else{
        return(self$edges)
      }
  },

  #' @description Bounding box of the metric graph
  #' @param format If the metric graph has a coordinate reference system, the format for the exported object. The options are `sf` (default), `sp` and `matrix`.
  #' @return A bounding box of the metric graph

  get_bounding_box = function(format = "sf"){
    bounding_box <- private$bounding_box
      if (private$longlat) {
      if (format == "sf") {
        bounding_box <- sf::st_bbox(
          c(xmin = bounding_box$min_x, 
            ymin = bounding_box$min_y,
            xmax = bounding_box$max_x, 
            ymax = bounding_box$max_y),
          crs = private$crs 
        )
      } else {
        sp_bbox_matrix <- matrix(
          c(bounding_box$min_x, bounding_box$max_x, 
            bounding_box$min_y, bounding_box$max_y), 
          ncol = 2, 
          dimnames = list(c("x", "y"), c("min", "max"))
        )

        bounding_box <- sp::SpatialPolygons(
          list(sp::Polygons(list(sp::Polygon(cbind(
            c(sp_bbox_matrix[1, "min"], sp_bbox_matrix[1, "max"], sp_bbox_matrix[1, "max"], sp_bbox_matrix[1, "min"], sp_bbox_matrix[1, "min"]),
            c(sp_bbox_matrix[2, "min"], sp_bbox_matrix[2, "min"], sp_bbox_matrix[2, "max"], sp_bbox_matrix[2, "max"], sp_bbox_matrix[2, "min"])
          ))), ID = "1")), 
          proj4string = private$proj4string  
        )
      }
    }
    return(bounding_box)
  },

  #' @description Exports the vertices of the MetricGraph object as an `sf`, `sp` or as a matrix.
  #' @param format The format for the exported object. The options are `sf` (default), `sp` and `matrix`.
  #' @return 
#' For `which_format == "sf"`, the function returns an `sf` object of `POINT` geometries. 
#' 
#' For `which_format == "sp"`, the function returns a `SpatialPointsDataFrame` object.

  get_vertices = function(format = c("sf", "sp", "list")){
    format <- format[[1]]
    if(!(format%in%c("sf", "sp", "list"))){
      stop("'format' must be one of the following: 'sf', 'sp' or 'list'.")
    }
      vertices_df <- do.call(rbind, lapply(self$vertices, function(v) {
        data.frame(
          X = v[1],
          Y = v[2],
          degree = attr(v, "degree"),
          indegree = attr(v, "indegree"),
          outdegree = attr(v, "outdegree"),
          problematic = attr(v, "problematic"),
          id = attr(v, "id"),
          longlat = attr(v, "longlat"),
          stringsAsFactors = FALSE
          )}))    
    if(format == "sf"){
      vertices_df <- sf::st_as_sf(vertices_df, coords = c("X", "Y"))
      if(!is.null(private$crs)){
        sf::st_crs(vertices_df) <- private$crs
      } else{
        sf::st_crs(vertices_df) <- sf::NA_crs_
      }
    } else if(format == "sp"){
      sp::coordinates(vertices_df) <- ~ X + Y
      if(!is.null(private$proj4string)){
        sp::proj4string(vertices_df) <- private$proj4string
      }
    } 
    return(vertices_df)
  }, 

  #' @description Exports the MetricGraph object as an `sf` or `sp` object.
  #' @param format The format for the exported object. The options are `sf` (default) and `sp`.
  #' @return Returns a list with three elements: `edges`, `vertices`, and `data`.
  #' 
  #' For `format == "sf"`, `edges` is an `sf` object of `LINESTRING` geometries with edge weights, and `vertices` and `data` are `sf` objects with `POINT` geometries.
  #' 
  #' For `format == "sp"`, `edges` is a `SpatialLinesDataFrame` with edge weights, and `vertices` and `data` are `SpatialPointsDataFrame`.

  export = function(format = "sf"){
    edges <- self$get_edges(format = format)
    vertices <- self$get_vertices(format = format)
    if(!is.null(private$data)){
      data = self$get_data(format = format, drop_all_na=FALSE, drop_na = FALSE, group = NULL)
    } else{
      data = NULL
    }
    exported_metric_graph = list(edges = edges, vertices = vertices, data = data)
    # if(format == "sf"){
    #   edges_geometries <- lapply(self$edges, sf::st_linestring) 
      
    #   if(is.vector(private$edge_weights)){
    #     ew_tmp <- data.frame(.weights = private$edge_weights)
    #   } else{
    #     ew_tmp <- as.data.frame(private$edge_weights)
    #   }

    #   edges_sf <- sf::st_sf(ew_tmp, geometry = sf::st_sfc(edges_geometries), crs = if(!is.null(private$crs)) private$crs else NULL)

    #   if(!is.null(private$data)){
    #     data_tmp <- as.data.frame(private$data)
    #     data_geometries <- lapply(1:nrow(data_tmp), function(i) sf::st_point(as.numeric(data_tmp[i, c('.coord_x', '.coord_y')])))
    #     data_sf <- sf::st_sf(data_tmp, geometry = sf::st_sfc(data_geometries), crs = if(!is.null(private$crs)) private$crs else NULL)
    #     class(data_sf) <- c("metric_graph_data", class(data_sf))
    #   } else{
    #     data_sf <- NULL
    #   }
    #   exported_metric_graph <- list(edges = edges_sf, data = data_sf)
    # } else if(format == "sp"){
    #   edges_list <- lapply(1:length(self$edges), function(i) {
    #     sp::Line(coords = matrix(self$edges[[i]], nrow = dim(self$edges[[i]])[1], ncol = dim(self$edges[[i]])[2]))
    #   })
    #   sp_edges <- sp::SpatialLines(
    #       lapply(1:length(edges_list), function(i) sp::Lines(list(edges_list[[i]]), ID = as.character(i))),
    #       proj4string = if (!is.null(private$crs)) private$proj4string else NULL
    #   )
    #   if(is.vector(private$edge_weights)){
    #     ew_tmp <- data.frame(.weights = private$edge_weights)
    #   } else{
    #     ew_tmp <- as.data.frame(private$edge_weights)
    #   }
    #   ew_tmp[[".ID"]] <- as.character(1:length(sp_edges))
    #   edges_sldf <- sp::SpatialLinesDataFrame(sp_edges, data = ew_tmp, match.ID = ".ID")
    #   data_sp <- as.data.frame(private$data)
    #   sp::coordinates(data_sp) <- ~ .coord_x + .coord_y
    #   exported_metric_graph <- list(edges = edges_sldf, data = data_sp)
    # } else{
    #   stop(paste(format,"is not currently not a valid format to be exported."))
    # }
    return(exported_metric_graph)
  },

  #' @description Return the metric graph as a `leaflet::leaflet()` object to be built upon.
  #' @param width	the width of the map
  #' @param height the height of the map
  #' @param padding	the padding of the map
  #' @param options	the map options
  #' @param elementId	Use an explicit element ID for the widget (rather than an automatically generated one).
  #' @param sizingPolicy	htmlwidgets sizing policy object. Defaults to `leafletSizingPolicy()`.
  
  leaflet = function(width = NULL, height = NULL, padding = 0, options = leafletOptions(), elementId = NULL,
  sizingPolicy = leafletSizingPolicy(padding = padding)){
    edges_sf <- self$get_edges(format = "sf")
    return(leaflet::leaflet(data = edges_sf, width = width, height = height, padding = padding, options = options, elementId = elementId, sizingPolicy = sizingPolicy))
  },

  #' @description Returns a `mapview::mapview()` object of the metric graph
  #' @param ... Additional arguments to be passed to `mapview::mapview()`. The `x` argument of mapview, containing the metric graph is already passed internally.
  
  mapview = function(...){
    edges_sf <- self$get_edges(format = "sf")
    return(mapview::mapview(x = edges_sf, ...))
  },

  #' @description Sets the edge weights
  #' @param weights Either a number, a numerical vector with length given by the number of edges, providing the edge weights, or a `data.frame` with the number of rows being equal to the number of edges, where
  #' each row gives a vector of weights to its corresponding edge.
  #' @param kirchhoff_weights If non-null, the name (or number) of the column of `weights` that contain the Kirchhoff weights. Must be equal to 1 (or `TRUE`) in case `weights` is a single number and those are the Kirchhoff weights.
  #' @param directional_weights If non-null, the name (or number) of the column of `weights` that contain the directional weights.
  #' @param verbose There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @return No return value. Called for its side effects.

  set_edge_weights = function(weights = NULL, kirchhoff_weights = NULL,
                            directional_weights = NULL, verbose = 0) {
  # Validate input types early
    if(!is.vector(weights) && !is.data.frame(weights) && !is.null(weights)){
      stop("'weights' must be either a vector or a data.frame!")
    }

    if(verbose == 2){
      message("Setting edge weights...")
    }

    if(!is.null(weights)){
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
          if(!(".weights" %in% colnames(private$edge_weights))){
            private$edge_weights[, ".weights"] <- 1
          }
        }
    } else{
      weights <- private$edge_weights
    }

    if(!is.null(kirchhoff_weights)){
        if(length(kirchhoff_weights)>1){
          warning("Only the first entry of 'kirchhoff_weights' was used.")
          kirchhoff_weights <- kirchhoff_weights[[1]]
        }
        if(!is.numeric(kirchhoff_weights)){
          if(!is.character(kirchhoff_weights)){
            stop("'kirchhoff_weights' must be either a number of a string.")
          }
          if(!(kirchhoff_weights%in%colnames(weights))){
            stop(paste(kirchhoff_weights, "is not a column of 'weights'!"))
          }
        } else{
          if(!is.data.frame(weights)){
            if(kirchhoff_weights != 1){
              stop("Since 'weights' is not a data.frame, 'kirchhoff_weights' must be either NULL or 1.")
            }
          } else{
            if(kirchhoff_weights %%1 != 0){
              stop("'kirchhoff_weights' must be an integer.")
            }
            if((kirchhoff_weights < 1) || (kirchhoff_weights > ncol(weights))){
              stop("'kirchhoff_weights' must be a positive integer number smaller or equal to the number of columns of 'weights'.")
            }
          }
        }
        private$kirchhoff_weights <- kirchhoff_weights
       } else{
        if(is.vector(weights)){
          private$kirchhoff_weights <- 1
        } else{
          private$kirchhoff_weights <- ".weights"
        }        
       }

      if(!is.null(directional_weights)){
        if(!is.numeric(directional_weights)){
          if(!is.character(directional_weights)){
            stop("'directional_weights' must be either a numerical vector of a string vector.")
          }
          if(!(directional_weights%in%colnames(weights))){
            stop(paste(directional_weights, "is not a column of 'weights'!"))
          }
        } else{
          if(!is.data.frame(weights)){
            if(directional_weights != 1){
              stop("Since 'weights' is not a data.frame, 'directional_weights' must be either NULL or 1.")
            }
          } else{
            if(directional_weights %%1 != 0){
              stop("'directional_weights' must be an integer.")
            }
            if((directional_weights < 1) || (directional_weights > ncol(weights))){
              stop("'directional_weights' must be a positive integer number smaller or equal to the number of columns of 'weights'.")
            }
          }
        }
        private$directional_weights <- directional_weights
       } else{
        if(is.vector(weights)){
          private$directional_weights <- 1
        } else{
          private$directional_weights <- ".weights"
        }
       }

  # Get edge lengths and preallocate attributes
  edge_lengths <- self$get_edge_lengths()
  edges_PtE <- lapply(self$edges, function(edge) attr(edge, "PtE"))
  nE <- self$nE

  # Use `lapply` with indexing for varying attributes
  self$edges <- lapply(seq_along(self$edges), function(i) {
    edge <- self$edges[[i]]
    if(is.data.frame(private$edge_weights)){
      attr(edge, "weight") <- private$edge_weights[i, , drop=FALSE]
    } else{
      attr(edge, "weight") <- private$edge_weights[i]
    }
    attr(edge, "length") <- edge_lengths[i]
    attr(edge, "id") <- i
    attr(edge, "PtE") <- edges_PtE[[i]]
    attr(edge, "longlat") <- private$longlat
    attr(edge, "crs") <- private$crs$input
    attr(edge, "kirchhoff_weight") <- private$kirchhoff_weights
    attr(edge, "directional_weights") <- private$directional_weights

    # Set the class for each edge
    class(edge) <- "metric_graph_edge"
    edge
  })

  # Set the class of `self$edges` to `metric_graph_edges`
  class(self$edges) <- "metric_graph_edges"
  },


  #' @description Gets the edge weights
  #' @param data.frame If the edge weights are given as vectors, should the result be returned as a data.frame?
 #' @param format Which format should the data be returned? The options are `tibble` for `tidyr::tibble`, `sf` for `POINT`, `sp` for `SpatialPointsDataFrame` and `list` for the internal list format.
 #' @param tibble `r lifecycle::badge("deprecated")` Use `format` instead.
  #' @return A vector or `data.frame` containing the edge weights.

  get_edge_weights = function(data.frame = FALSE, format = c("tibble", "sf", "sp", "list"), tibble = deprecated()){

  if (lifecycle::is_present(tibble)) {
      lifecycle::deprecate_warn("1.3.0.9000", "get_edge_weights(tibble)", "get_edge_weights(format)",
        details = c("The argument `tibble` was deprecated in favor of the argument `format`.")
      )
    if(tibble){
      format <- "tibble"
    }
  }   

    format <- format[[1]]

    format <- tolower(format)

    if(!(format %in% c("tibble", "sf", "sp", "list"))){
      stop("The possible formats are 'tibble', 'sf', 'sp' and 'list'.")
    }  
    
    tmp <- private$edge_weights
    row.names(tmp) <- NULL
    if(!is.data.frame(tmp) && data.frame){
      tmp <- data.frame(.weights = tmp)
    }
    if(format == "tibble"){
      if(!is.data.frame(tmp)){
        tmp <- data.frame(.weights = tmp)
      }
      tmp <- dplyr::as_tibble(tmp)
    } else{
      return(self$get_edges(format = format))
    }
    class(tmp) <- c("metric_graph_weights", class(tmp))
    return(tmp)
  },

  #' @description Gets vertices with incompatible directions
  #' @return A vector containing the vertices with incompatible directions.

  get_vertices_incomp_dir = function(){
    if(is.null(self$vertices)){
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
    } else{
      problematic <- sapply(self$vertices, function(vert){attr(vert, "problematic")})
      return(which(problematic))
    }

  },

  #' @description Prints a summary of various informations of the graph
  #' @param messages Should message explaining how to build the results be given for missing quantities?
  #' @param compute_characteristics Should the characteristics of the graph be computed? If `NULL` it will be determined based on the size of the graph.
  #' @param check_euclidean Check if the graph has Euclidean edges? If `NULL` it will be determined based on the size of the graph.
  #' @param check_distance_consistency Check the distance consistency assumption? If `NULL` it will be determined based on the size of the graph.
  #' @return No return value. Called for its side effects.

  summary = function(messages = FALSE, compute_characteristics = NULL, check_euclidean = NULL, check_distance_consistency = NULL){

    if(self$nV > 10000){
      if(is.null(compute_characteristics)){
        compute_characteristics <- FALSE
      }
      if(is.null(check_euclidean)){
        check_euclidean <- FALSE
      }
      if(is.null(check_distance_consistency)){
        check_distance_consistency <- FALSE
      }
    } else{
      if(is.null(compute_characteristics)){
        compute_characteristics <- TRUE
      }
      if(is.null(check_euclidean)){
        check_euclidean <- TRUE
      }
      if(is.null(check_distance_consistency)){
        check_distance_consistency <- TRUE
      }
    }

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
      g <- make_graph(edges = c(t(self$E)), directed = FALSE)
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
  #' @param verbose Print progress of the computation of the geodesic distances. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @return No return value. Called for its side effects. The computed geodesic
  #' distances are stored in the `geo_dist` element of the `metric_graph` object.
  compute_geodist = function(full = FALSE, obs = TRUE, group = NULL, verbose = 0) {
    if(is.null(self$geo_dist)){
      self$geo_dist <- list()
    }

    if(is.null(private$data)){
      obs <- FALSE
    }
    if(!obs){
      if(verbose == 2){
        message("Creating auxiliary graph...")
      }
      t <- system.time({
        g <- make_graph(edges = c(t(self$E)), directed = FALSE)
        E(g)$weight <- self$edge_lengths
      })
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
      if(verbose > 0){
        message("Computing geodesic distances...")
      }
      t <- system.time(
        self$geo_dist[[".vertices"]] <- distances(g)
      )
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    } else if(full){
      PtE_full <- self$get_PtE()
      self$geo_dist[[".complete"]] <- self$compute_geodist_PtE(PtE = PtE_full,
                                                              normalized = TRUE, verbose = verbose)
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
                                                              normalized = TRUE, verbose = verbose)
      }
    }
  },
  #' @description Computes shortest path distances between the vertices in the
  #' graph.
  #' @param PtE Points to compute the metric for.
  #' @param normalized are the locations in PtE in normalized distance?
  #' @param include_vertices Should the original vertices be included in the
  #' distance matrix?
  #' @param verbose Print progress of the computation of the geodesic distances. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @return A matrix containing the geodesic distances.
  compute_geodist_PtE = function(PtE,
                                 normalized = TRUE,
                                 include_vertices = TRUE, verbose = 0){
      if(verbose == 2){
        message("Processing the graph locations...")
      }
      t <- system.time({
        graph.temp <- self$clone()
        graph.temp$clear_observations()
        df_temp <- data.frame(y = rep(0, dim(PtE)[1]),
                            edge_number = PtE[,1],
                            distance_on_edge = PtE[,2])
        if(sum(duplicated(df_temp))>0){
          warning("Duplicated locations were found when computing geodist. The returned values are given for unique locations.")
          df_temp <- unique(df_temp)
        }

        df_temp <- standardize_df_positions(df_temp, self)

        graph.temp$build_mesh(h = 10000)

        df_temp2 <- data.frame(y = 0, edge_number = graph.temp$mesh$VtE[1:nrow(self$V),1],
                                  distance_on_edge = graph.temp$mesh$VtE[1:nrow(self$V),2])

        df_temp2 <- standardize_df_positions(df_temp2, self)

        df_temp$included <- TRUE
        temp_merge <- merge(df_temp, df_temp2, all = TRUE)

        df_temp$included <- NULL

        df_temp2 <- temp_merge[is.na(temp_merge["included"]),]

        nV_new <- sum(is.na(temp_merge["included"]))

        df_temp2$included <- NULL

        df_temp <- rbind(df_temp2, df_temp)

        df_temp[["__dummy"]] <- 1:nrow(df_temp)

        graph.temp$add_observations(data = df_temp,
                                     normalized = normalized,
                                     verbose=0,
                  suppress_warnings = TRUE)

      })

      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }

      if(verbose == 2){
        message("Turning observations of the auxiliary graph to vertices...")
      }
      t <- system.time(
        graph.temp$observation_to_vertex(mesh_warning = FALSE)
      )
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
       }
      if(verbose == 2){
        message("Creating auxiliary graph...")
      }
      t <- system.time({
        g <- make_graph(edges = c(t(graph.temp$E)), directed = FALSE)
        E(g)$weight <- graph.temp$edge_lengths
      })
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
       }
      if(verbose>0){
        message("Computing geodesic distances...")
      }
      t <- system.time(
        geodist_temp <- distances(g)
      )
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
       }

      if(length(graph.temp$PtV)[1]!=nrow(geodist_temp)){
        un_PtV <- unique(graph.temp$PtV)
        un_coords <- !duplicated(graph.temp$PtV)

        geodist_temp <- geodist_temp[un_PtV, un_PtV]

        tmp_vec <- graph.temp$.__enclos_env__$private$data[["__dummy"]][un_coords]

        un_ord <- order(tmp_vec)

        tmp_vec[un_ord] <- 1:length(tmp_vec)
        #Ordering back in the input order
        geodist_temp[tmp_vec,tmp_vec] <- geodist_temp
      } else{
          geodist_temp <- geodist_temp[graph.temp$PtV, graph.temp$PtV]
          #Ordering back in the input order
          geodist_temp[graph.temp$.__enclos_env__$private$data[["__dummy"]],graph.temp$.__enclos_env__$private$data[["__dummy"]]] <- geodist_temp
      }

      if(!include_vertices){
        geodist_temp <- geodist_temp[(nV_new+1):nrow(geodist_temp), (nV_new+1):nrow(geodist_temp)]
      }

      attr(geodist_temp, "unit") <- private$length_unit

      return(geodist_temp)
  },

  #' @description Computes shortest path distances between the vertices in the
  #' mesh.
  #' @return No return value. Called for its side effects. The geodesic distances
  #' on the mesh are stored in `mesh$geo_dist` in the `metric_graph` object.
  compute_geodist_mesh = function() {
    g <- make_graph(edges = c(t(self$mesh$E)), directed = FALSE)
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
  #' @param check_euclidean Check if the graph used to compute the resistance distance has Euclidean edges? The graph used to compute the resistance distance has the observation locations as vertices.
  #' @param include_vertices Should the vertices of the graph be also included in the resulting matrix when using `FULL=TRUE`?
  #' @param verbose Print progress of the computation of the resistance distances. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @return No return value. Called for its side effects. The geodesic distances
  #' are stored in the `res_dist` element of the `metric_graph` object.
  compute_resdist = function(full = FALSE, obs = TRUE, group = NULL,
                                 check_euclidean = FALSE, include_vertices = FALSE, verbose = 0) {
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
                                                                normalized=TRUE,
                                                                       check_euclidean = check_euclidean, verbose = verbose)
    } else if(full){
      PtE <- self$get_PtE()
      if(!include_vertices){
        self$res_dist[[".complete"]] <- self$compute_resdist_PtE(PtE,
                                                                normalized=TRUE, include_vertices = FALSE,
                                                                       check_euclidean = check_euclidean, verbose = verbose)
      } else{
        self$res_dist[[".complete"]] <- self$compute_resdist_PtE(PtE,
                                                                normalized=TRUE, include_vertices = TRUE,
                                                                       check_euclidean = check_euclidean, verbose = verbose)
      }

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
                                                                       normalized=TRUE,
                                                                       check_euclidean = check_euclidean, verbose = verbose)
      }
    }
  },

  #' @description Computes the resistance distance between the observation
  #' locations.
  #' @param PtE Points to compute the metric for.
  #' @param normalized Are the locations in PtE in normalized distance?
  #' @param include_vertices Should the original vertices be included in the
  #' Laplacian matrix?
  #' @param check_euclidean Check if the graph used to compute the resistance distance has Euclidean edges? The graph used to compute the resistance distance has the observation locations as vertices.
  #' @param verbose Print progress of the computation of the resistance distances. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @return A matrix containing the resistance distances.
  compute_resdist_PtE = function(PtE,
                                 normalized = TRUE,
                                 include_vertices = FALSE,
                                 check_euclidean = FALSE, verbose = 0) {
      if(verbose == 2){
        message("Processing the graph locations...")
      }
      t <- system.time({
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
                                     normalized = normalized, verbose = 0,
                  suppress_warnings = TRUE)
      })
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
      if(verbose == 2){
        message("Turning observations of the auxiliary graph to vertices...")
      }
      t <- system.time(
        graph.temp$observation_to_vertex(mesh_warning=FALSE)
      )
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
       }
      if(verbose > 0){
        message("Computing auxiliary geodesic distances...")
      }
      t <- system.time(
        graph.temp$compute_geodist(full=TRUE,verbose = verbose)
      )
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
       }

        if(check_euclidean){
          graph.temp$check_euclidean()
          is_euclidean <- graph.temp$characteristics$euclidean
        }

      if(verbose > 0){
        message("Computing the resistance distances...")
      }
      t <- system.time({
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
      })

      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
       }

      if(check_euclidean){
        attr(R, "euclidean") <- is_euclidean
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
  #' @param approx Should the computation of the relative positions be approximate? Default is `TRUE`. If `FALSE`, the speed can be considerably slower, especially for large metric graphs.
  #' @param verbose Level of verbosity, 0, 1 or 2. The default is 0.
  #' @return No return value, called for its side effects.
  compute_PtE_edges = function(approx = TRUE, verbose = 0){
    
    if(verbose > 0){
      message("Computing the relative positions of the edges...")
      if(verbose == 2){
        bar_edges_pte <- msg_progress_bar(length(self$edges))
      }
    }
    
    if(approx){
      edges_PtE <- lapply(self$edges, function(edge){
      if(verbose == 2){
        bar_edges_pte$increment()
      }
      private$approx_coordinates(edge = edge)
      })
    } else{
      edges_PtE <- lapply(self$edges, function(edge){
      if(verbose == 2){
        bar_edges_pte$increment()
      }
      private$exact_PtE_coordinates(edge = edge)
      })
    }
    self$edges <- lapply(1:length(self$edges), function(j){
      edge <- self$edges[[j]]
      attr(edge, "PtE") <- edges_PtE[[j]]
      return(edge)
    })
    class(self$edges) <- "metric_graph_edges"
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
  #' @param verbose Print progress of the computation of the Laplacian. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @return No reutrn value. Called for its side effects. The Laplacian is stored
  #' in the `Laplacian` element in the `metric_graph` object.
  compute_laplacian = function(full = FALSE, obs = TRUE, group = NULL, verbose = 0) {
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
                                                            normalized = TRUE, verbose = verbose)
    } else if(full){
      PtE <- self$get_PtE()
      self$Laplacian[[".complete"]] <- private$compute_laplacian_PtE(PtE,
                                                            normalized = TRUE, verbose = verbose)
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
                                                              normalized = TRUE, verbose = verbose)
      }
    }
  },

  #' @description Removes vertices of degree 2 from the metric graph.
  #' @return No return value. Called for its side effects.
  #' @param check_weights If `TRUE` will only prune edges with different weights.
  #' @param check_circles If `TRUE` will not prune a vertex such that the resulting edge is a circle.
 #' @param verbose Print progress of pruning. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @details
    #' Vertices of degree 2 are removed as long as the corresponding edges that
    #' would be merged are compatible in terms of direction.
    #'
  prune_vertices = function(check_weights = TRUE, check_circles = TRUE, verbose = FALSE){
    t <- system.time({
    degrees <- private$degrees$degrees

    # Finding problematic vertices, that is, vertices with incompatible directions
    # They will not be pruned.

    problematic <- sapply(self$vertices, function(vert){attr(vert,"problematic")})

    if((verbose > 0) && (sum(problematic) > 0)){
      message(paste(sum(problematic), "vertices were not pruned due to incompatible directions."))
    }

    if (check_weights) {
      # Identify the vertices with degree 2 that are not problematic
      idx_tmp <- which(degrees == 2 & !problematic)
      problematic_weights <- rep(FALSE, self$nV)

      if (verbose == 2) {
        message("Checking weight compatibility")
      }

      # Vectorized operations
      start_deg <- match(idx_tmp, self$E[, 1], nomatch = 0)
      end_deg <- match(idx_tmp, self$E[, 2], nomatch = 0)

      # Combine indices to get the edges related to each vertex
      edges_tmp <- cbind(start_deg, end_deg)
      valid_edges <- rowSums(edges_tmp == 0) == 0  # Exclude invalid matches

      if (any(valid_edges)) {
        # Only check weights for valid edges
        edges_tmp <- edges_tmp[valid_edges, , drop = FALSE]
        idx_tmp <- idx_tmp[valid_edges]

        # Check weight compatibility in a vectorized way
        if (is.vector(private$edge_weights)) {
          # Create temporary copies with NA replaced by a unique placeholder
          # edge_weights_copy <- private$edge_weights
          # edge_weights_copy[is.na(edge_weights_copy)] <- ".dummy_na_val"

          # Compare the edges normally, treating NA == NA as TRUE via the placeholder
          # cnd_tmp <- edge_weights_copy[edges_tmp[, 1]] != edge_weights_copy[edges_tmp[, 2]]          
          # Compare the edges normally, treating NA == NA as TRUE via the placeholder
          cnd_tmp <- compare_with_na(
            private$edge_weights[edges_tmp[, 1]], 
            private$edge_weights[edges_tmp[, 2]]
          )
        } else {
          # For matrices, replace NA with a unique placeholder in a temporary copy
          # edge_weights_copy <- private$edge_weights
          # edge_weights_copy[is.na(edge_weights_copy)] <- ".dummy_na_val"

          # Perform row-wise comparison, each row must have all columns matching
          # cnd_tmp <- rowSums(
          #   edge_weights_copy[edges_tmp[, 1], , drop = FALSE] != edge_weights_copy[edges_tmp[, 2], , drop = FALSE]
          # ) == ncol(edge_weights_copy)
          cnd_tmp <- compare_with_na(
            private$edge_weights[edges_tmp[, 1], , drop = FALSE],
            private$edge_weights[edges_tmp[, 2], , drop = FALSE],
            is_matrix = TRUE
          )
        }

        # Update problematic_weights vector based on the results
        problematic_weights[idx_tmp] <- cnd_tmp
      }

      # Update problematic vertices
      problematic <- (problematic | problematic_weights)

      if ((verbose > 0) && (sum(problematic_weights) > 0)) {
        message(paste(sum(problematic_weights), "vertices were not pruned due to incompatible weights. Turn 'check_weights' to FALSE to prune these vertices."))
      }
    }

    res <- list(degrees = degrees, problematic = problematic)

    res[["problematic_circles"]] <- rep(FALSE, length(degrees))

    if(verbose > 0){
      to.prune <- sum(res$degrees==2 & !res$problematic)
      k <- 1
      message(sprintf("removing %d vertices", to.prune))
      if(to.prune > 0) {
        # pb = txtProgressBar(min = 1, max = to.prune, initial = 1, style = 3)
        bar_prune <- msg_progress_bar(to.prune)
      }

    }

   while(sum(res$degrees==2 & !res$problematic & !res$problematic_circles)>0) {
     if((verbose == 2) && to.prune > 0){
      #  setTxtProgressBar(pb,k)
      bar_prune$increment()
       #message(sprintf("removing vertex %d of %d.", k, to.prune))
       k <- k + 1
     }
     res <- private$remove.first.deg2(res, check_circles = check_circles)
   }
    if(verbose == 2){
      if(!is.null(res$circles_avoided)){
        message(paste(sum(res$problematic_circles),"vertices were not pruned in order to avoid creating circles. Turn 'check_circles' to FALSE to prune these vertices."))
      }
    }
   })
    if(verbose  == 2){
          message(sprintf("time: %.3f s", t[["elapsed"]]))
    }
    # if(verbose && to.prune > 0){
    #   close(pb)
    # }

   t <- system.time({
      private$create_update_vertices(verbose=verbose)
      # creating/updating reference edges
      private$ref_edges <- map_into_reference_edge(self, verbose=verbose)
      if(verbose==2){
                message("Updating attributes of the edges")
                bar_update_attr_edges <- msg_progress_bar(length(self$edges))
      }
      # for(i in 1:length(self$edges)){
      #    attr(self$edges[[i]], "id") <- i
      #    attr(self$edges[[i]], "longlat") <- private$longlat
      #    attr(self$edges[[i]], "crs") <- private$crs$input
      #    attr(self$edges[[i]], "length") <- self$edge_lengths[i]
      #    class(self$edges[[i]]) <- "metric_graph_edge"
      #   if(!is.null(private$length_unit)){
      #     units(attr(self$edges[[i]], "length")) <- private$length_unit
      #   }
      #   if(is.vector(private$edge_weights)){
      #     attr(self$edges[[i]], "weight") <- private$edge_weights[i]
      #   } else{
      #     attr(self$edges[[i]], "weight") <- private$edge_weights[i,]
      #   }
      #   attr(self$edges[[i]], "kirchhoff_weight") <- private$kirchhoff_weights
      #   if(verbose == 2){
      #     bar_update_attr_edges$increment()
      #   }        
      # }

    edge_lengths_ <- self$get_edge_lengths()

    self$edges <- lapply(1:self$nE, function(i){
      edge <- self$edges[[i]]
      if(is.vector(private$edge_weights)){
        attr(edge,"weight") <- private$edge_weights[i]
      } else{
        attr(edge,"weight") <- private$edge_weights[i, ,drop=FALSE]
      }
      attr(edge, "longlat") <- private$longlat
      attr(edge, "crs") <- private$crs$input
      attr(edge, "length") <- edge_lengths_[i]
      attr(edge, "id") <- i
      attr(edge, "kirchhoff_weight") <- private$kirchhoff_weights
      attr(edge, "directional_weights") <- private$directional_weights
      class(edge) <- "metric_graph_edge"
        if(verbose == 2){
          bar_update_attr_edges$increment()
        }      
      return(edge)
    })



   })
   if(verbose == 2){
            message(sprintf("time: %.3f s", t[["elapsed"]]))
   }


   if(!is.null(private$data)){
    if(verbose > 0){
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
      old_group_variable <- attr(private$data, "group_variable")
      private$data <- lapply(private$data, function(dat){dat[order_idx]})
      attr(private$data, "group_variable") <- old_group_variable
      })
      if(verbose == 2){
            message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
   }

   if(!is.null(self$mesh)){
    max_h <- max(self$mesh$h_e)
    self$mesh <- NULL
    self$build_mesh(h = max_h)
   }

   self$C <- NULL
   self$CoB <- NULL

   private$pruned <- TRUE
   private$degrees <- NULL
   if(private$prune_warning){
    warning("At least two edges with different weights were merged due to pruning. Only one of the weights has been assigned to the merged edge. Please, review carefully.")
    private$prune_warning <- FALSE
   }
  },


  #' @description Gets the groups from the data.
  #' @param edge_lengths edge lengths to be set to the metric graph edges.
  #' @param unit set or override the edge lengths unit.
  #' @return does not return anything. Called for its side effects.

    set_manual_edge_lengths = function(edge_lengths, unit = NULL){
      if(is.null(edge_lengths)){
        warning("edge_lengths is NULL, edge lengths were not set.")
        return(invisible(NULL))
      }
      if(length(edge_lengths)!=length(self$edges)){
        stop("edge_lengths must have length equal to the number of edges.")
      }
      private$manual_edge_lengths <- TRUE
      self$edge_lengths <- edge_lengths
      private$length_unit <- unit
      return(invisible(NULL))
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
    names(el) <- NULL
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
  #' @param share_weights Should the same weight be shared among the split edges? If `FALSE`, the weights will be removed, and a common weight given by 1 will be given.
  #' @param mesh_warning Display a warning if the graph structure change and the metric graph has a mesh object.
  #' @param verbose Print progress of the steps when adding observations. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @param tolerance `r lifecycle::badge("deprecated")`. Not used anymore
  #' @return No return value. Called for its side effects.

  observation_to_vertex = function(mesh_warning = TRUE, verbose = 0, tolerance = deprecated()) {
    if (lifecycle::is_present(tolerance)) {
        lifecycle::deprecate_warn("1.3.0.9000", "observation_to_vertex(tolerance)")
    }
    
    # Initialize temp data
    if (is.null(private$data)) {
        stop("There is no data!")
    }

    private$temp_PtE <- self$get_PtE()
    n <- nrow(private$temp_PtE)
    self$PtV <- rep(NA, n)
    private$temp_PtE <- cbind(private$temp_PtE, seq_len(n))

    # Identify start and end vertices based on the tolerance
    is_start_vertex <- abs(private$temp_PtE[, 2]) < 1e-15
    is_end_vertex <- private$temp_PtE[, 2] > 1 - 1e-15

    # Assign known vertices directly
    self$PtV[is_start_vertex] <- self$E[private$temp_PtE[is_start_vertex, 1], 1]
    self$PtV[is_end_vertex] <- self$E[private$temp_PtE[is_end_vertex, 1], 2]

    # Get remaining indices that need to be split
    remaining_indices <- which(!(is_start_vertex | is_end_vertex))

    # Group PtE by edges
    edge_groups <- split(private$temp_PtE[remaining_indices, c(2, 3)], private$temp_PtE[remaining_indices, 1])
    edge_groups <- lapply(edge_groups, function(coords) matrix(coords, ncol = 2, byrow = FALSE))   

    # Progress bar setup
    if (verbose == 2) {
        bar_otv <- msg_progress_bar(length(edge_groups))
    }

    # Loop over each edge in `edge_groups`
    for (Ei in names(edge_groups)) {
        
        # Extract `t_values` and `indices`
        t_values <- edge_groups[[Ei]][, 1]  # positions on edges
        indices <- edge_groups[[Ei]][, 2]   # original row number
        
        # Progress bar increment (if verbose mode is enabled)
        if (verbose == 2) {
            bar_otv$increment()
        }
        
        # Perform `private$split_edge`
        new_vertices <- private$split_edge(as.numeric(Ei), t_values,  indices = indices)
        
        # Assign new vertices
        self$PtV[indices] <- new_vertices
    }

    # Remove NA values from PtV
    self$PtV <- self$PtV[!is.na(self$PtV)]

    # Update temp_PtE for the known vertices

    # Find the positions in `self$E[,1]` where they match `self$PtV[is_start]` and `self$PtV[is_end]`
    start_positions <- match(self$PtV[is_start_vertex], self$E[, 1])
    end_positions <- match(self$PtV[is_end_vertex], self$E[, 2])

    # Use vectorized assignment for `private$temp_PtE` values
    private$temp_PtE[is_start_vertex, 1] <- start_positions
    private$temp_PtE[is_end_vertex, 1] <- end_positions

    # Assign 1 for `is_start_vertex` rows in the second column and 0 for `is_end_vertex` rows
    private$temp_PtE[is_start_vertex, 2] <- 0
    private$temp_PtE[is_end_vertex, 2] <- 1

    # Replicate edge numbers and distances for the number of groups
    n_group <- length(unique(private$data[[".group"]]))
    private$data[[".edge_number"]] <- rep(private$temp_PtE[, 1], times = n_group)
    private$data[[".distance_on_edge"]] <- rep(private$temp_PtE[, 2], times = n_group)

    # Reorder the data based on group and edge information
    tmp_df <- data.frame(
        PtE1 = private$data[[".edge_number"]],
        PtE2 = private$data[[".distance_on_edge"]],
        group = private$data[[".group"]]
    )
    index_order <- order(tmp_df$group, tmp_df$PtE1, tmp_df$PtE2)
    old_group_variable <- attr(private$data, "group_variable")
    
    private$data <- lapply(private$data, function(dat) dat[index_order])
    attr(private$data, "group_variable") <- old_group_variable

    # Reorder PtV according to index_order
    self$PtV <- self$PtV[index_order[1:length(self$PtV)]]

    # Reset temporary data
    private$temp_PtE <- NULL

    # Invalidate cached distances and recompute if necessary
    self$geo_dist <- NULL
    self$res_dist <- NULL
    if (!is.null(self$CoB)) self$buildC(2)

    if (!is.null(self$mesh)) {
        self$mesh <- NULL
        if (mesh_warning) {
            warning("Removing the existing mesh due to the change in the graph structure, please create a new mesh if needed.")
        }
    }

    # Update vertices and reference edges
    private$ref_edges <- map_into_reference_edge(self, verbose = verbose)
    private$data <- standardize_df_positions(private$data, self, edge_number = ".edge_number", distance_on_edge = ".distance_on_edge")
    private$create_update_vertices(verbose = verbose)
    self$set_edge_weights(
        weights = private$edge_weights,
        kirchhoff_weights = private$kirchhoff_weights,
        directional_weights = private$directional_weights,
        verbose = verbose
    )
  },
  
  #' @description Turns edge weights into data on the metric graph
  #' @param loc A `matrix` or `data.frame` with two columns containing the locations to generate the data from the edge weights. If `data_coords` is 'spatial', the first column must be the x-coordinate of the data, and the second column must be the y-coordinate. If `data_coords` is 'PtE', the first column must be the edge number and the second column must be the distance on edge.
  #' @param data_loc Should the data be generated to the data locations? In this case, the `loc` argument will be ignored. Observe that the metric graph must have data for one to use this option. CAUTION: To add edgeweight to data to both the data locations and mesh locations, please, add at the data locations first, then to mesh locations.
  #' @param mesh Should the data be generated to the mesh locations? In this case, the `loc` argument will be ignored. Observe that the metric graph must have a mesh built for one to use this option. CAUTION: To add edgeweight to data to both the data locations and mesh locations, please, add at the data locations first, then to mesh locations.
  #' @param weight_col Which columns of the edge weights should be turned into data? If `NULL`, all columns will be turned into data.
  #' @param add Should the data generated be added to the metric graph internal data?
  #' @param data_coords To be used only if `mesh` is `FALSE`. It decides which
  #' coordinate system to use. If `PtE`, the user must provide `edge_number` and
  #' `distance_on_edge`, otherwise if `spatial`, the user must provide
  #' `coord_x` and `coord_y`.
  #' @param normalized if TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default FALSE.
  #' @param tibble Should the data be returned in a `tibble` format?
  #' @param format If `return` is `TRUE`, the format of the output: "tibble", "sf", or "sp". Default is "tibble". 
  #' @param verbose Print progress of the steps when adding observations. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @param suppress_warnings Suppress warnings related to duplicated observations?
  #' @param return Should the data be returned? If `return_removed` is `TRUE`, only the removed locations will be return (if there is any).
  edgeweight_to_data = function(loc = NULL, mesh = FALSE,
                data_loc = FALSE,
                weight_col = NULL, add = TRUE,
                data_coords = c("PtE", "spatial"),
                normalized = FALSE,
                tibble = FALSE,
                format = c("tibble", "sf", "sp", "list"),
                verbose = 1,
                suppress_warnings = FALSE,
                return = FALSE){

              if(is.null(loc) && !mesh && !data_loc){
                stop("Either 'loc' must be provided, 'mesh' must be TRUE or 'data_loc' must be TRUE.")
              }
              if(mesh){
                if(is.null(self$mesh)){
                  stop("There is no mesh! Build a mesh using the build_mesh() method.")
                }
                loc <- self$mesh$VtE
                normalized <- TRUE
                data_coords <- "PtE"
              }

              if(data_loc){
                if(is.null(private$data)){
                  stop("The graph has no data!")
                }
                loc <- self$get_PtE()
                normalized <- TRUE
                data_coords <- "PtE"
              }

              data_coords <- data_coords[[1]]

              if(tolower(data_coords) == "pte"){
                df_ew <- data.frame(.edge_number = loc[,1],
                .distance_on_edge = loc[,2])
                if(normalized){
                  if(max(loc[,2])>1){
                    stop("distance_on_edge of normalized locations cannot be greater than 1!")
                    }
                  }
              } else if(tolower(data_coords) == "spatial"){
                loc_tmp <- self$coordinates(XY = loc)
                nr_tmp <- nrow(loc_tmp)
                loc_tmp <- unique(loc_tmp)
                if(nr_tmp != nrow(loc_tmp)){
                  if(!suppress_warnings){
                    warning("Some locations were projected to the same point of the metric graph and were not considered.")
                  }
                }
                df_ew <- data.frame(.edge_number = loc_tmp[,1],
                .distance_on_edge = loc_tmp[,2])
                normalized <- TRUE
              } else{
                stop("'data_coords' must be either 'PtE' or 'spatial'!")
              }

              if(!is.null(private$data[[".group"]])){
                df_ew <- cbind(df_ew, .group = rep(unique(private$data[[".group"]]),
                            each = nrow(df_ew)))
              }

              ew <- private$get_edge_weights_internal()

              if(is.vector(ew)){
                ew <- data.frame(.weight = ew)
              }

              if(!is.null(weight_col)){
                ew <- ew[, weight_col, drop=FALSE]
              }

              ew[[".edge_number"]] <- 1:nrow(ew)

              ew[[".edge_lengths"]] <- self$edge_lengths

              df_ew <- merge(df_ew, ew, by = ".edge_number")

              if(!normalized){
                if(max(df_ew[[".distance_on_edge"]]/df_ew[[".edge_lengths"]]) > 1){
                  stop("There is at least one distance on edge that is greater than the corresponding edge length!")
                }
              }

              df_ew[[".edge_lengths"]] <- NULL

              if(add){
                self$add_observations(data = df_ew,
                                      edge_number = ".edge_number",
                                      distance_on_edge = ".distance_on_edge",
                                      data_coords = "PtE",
                                      group = ".group",
                                      tibble = tibble,
                                      normalized = normalized,
                                      verbose = verbose)
              }

              if(return){
                  return(self$process_data(data = df_ew,
                                      edge_number = ".edge_number",
                                      distance_on_edge = ".distance_on_edge",
                                      data_coords = "PtE",
                                      group = ".group",
                                      normalized = normalized,
                                      format = format,
                                      verbose = verbose))
              } else{
                return(invisible(NULL))
              }
  },

  #' @description Returns a list or a matrix with the mesh locations.
  #' @param bru Should an 'inlabru'-friendly list be returned?
  #' @param loc If `bru` is set to `TRUE`, the column names of the location variables.
  #' The default name is `c('.edge_number', '.distance_on_edge')`.
  #' @param normalized If TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default TRUE.
  #' @param loc_name The name of the location variables. Not needed for `rSPDE` models.
  #'
  #' @return A list or a matrix containing the mesh locations.
  get_mesh_locations = function(bru = FALSE, loc = c(".edge_number", ".distance_on_edge"), loc_name = NULL, normalized = TRUE) {
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
      tmp_VtE <- self$mesh$VtE
      if(!normalized){
        tmp_VtE[,2] <- tmp_VtE[,2] * self$edge_lengths[tmp_VtE[, 1]]
      }
      data_list <- list()
      data_list[[loc[1]]] <- tmp_VtE[,1]
      data_list[[loc[2]]] <- tmp_VtE[,2]
      data_list <- as.data.frame(data_list)
      if(!is.null(loc_name)){
        ret_list <- list()
        ret_list[[loc_name]] <- data_list
        return(ret_list)
      }
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
  #' @param data_coords It decides which
  #' coordinate system to use. If `PtE`, the user must provide `edge_number` and
  #' `distance_on_edge`, otherwise if `spatial`, the user must provide
  #' `coord_x` and `coord_y`. The option `euclidean` is `r lifecycle::badge("deprecated")`. Use `spatial` instead.
  #' @param group Vector. If the data is grouped (for example measured at different time
  #' points), this argument specifies the columns (or entries on the list) in
  #' which the group variables are stored. It will be stored as a single column `.group` with the combined entries.
  #' @param group_sep separator character for creating the new group variable when grouping two or more variables.
  #' @param normalized if TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default FALSE.
  #' @param format Which format should the data be returned? The options are `tibble` for `tidyr::tibble`, `sf` for `POINT`, `sp` for `SpatialPointsDataFrame` and `list` for the internal list format.
  #' @param duplicated_strategy Which strategy to handle observations on the same location on the metric graph (that is, if there are two or more observations projected at the same location).
  #' The options are 'closest' and 'jitter'. If 'closest', only the closest observation will be used. If 'jitter', a small perturbation will be performed on the projected observation location. The default is 'closest'.
  #' @param include_distance_to_graph When `data_coord` is 'spatial', should the distance of the observations to the graph be included as a column?
  #' @param only_return_removed Should the removed data (if it exists) when using 'closest' `duplicated_strategy` be returned instead of the processed data?
  #' @param tolerance Parameter to control a warning when adding observations.
  #' If the distance of some location and the closest point on the graph is
  #' greater than the tolerance, the function will display a warning.
  #' This helps detecting mistakes on the input locations when adding new data.
  #' @param verbose If `TRUE`, report steps and times.
  #' @param suppress_warnings Suppress warnings related to duplicated observations?
  #' @param Spoints `r lifecycle::badge("deprecated")` Use `data` instead.
  #' @param tibble  `r lifecycle::badge("deprecated")` Use `format` instead.
  #' @return No return value. Called for its side effects. The observations are
  #' stored in the `data` element of the `metric_graph` object.

  process_data = function(    data = NULL,
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
                              tibble = lifecycle::deprecated()) {

  if (lifecycle::is_present(tibble)) {
      lifecycle::deprecate_warn("1.3.0.9000", "get_edge_weights(tibble)", "get_edge_weights(format)",
        details = c("The argument `tibble` was deprecated in favor of the argument `format`.")
      )
    if(tibble){
      format <- "tibble"
    }
  }       

    format <- format[[1]]

    format <- tolower(format)

    if(!(format %in% c("tibble", "sf", "sp", "list"))){
      stop("The possible formats are 'tibble', 'sf', 'sp' and 'list'.")
    }

        data_coords <- data_coords[[1]]
        if(!is.null(group)){
          group <- unique(group)
        }

    if(length(tolerance)>1){
      tolerance <- tolerance[[1]]
      warning("'tolerance' had more than one element, only the first one will be used.")
    }

    if(lifecycle::is_present(Spoints)){
           lifecycle::deprecate_warn("1.2.9000", "add_observations(Spoints)", "add_observations(data)",
             details = c("`Spoints` is deprecated, use `data` instead.")
           )
           data <- Spoints
    }

    removed_data <- NULL
    strc_data <- FALSE

    if(data_coords == "spatial" && 
       !inherits(data, "sf") && 
       !inherits(data, "metric_graph_data") && 
       !"SpatialPointsDataFrame" %in% is(data) && 
       !is.null(private$crs)) {
      
      if(!is.null(data[[coord_x]]) && !is.null(data[[coord_y]])) {
        warning("The crs of the data is not set, the data will be assumed to be in the same crs as the graph.")
        
        # Create sf object from the data
        data_df <- as.data.frame(data)
        data <- sf::st_as_sf(data_df, 
                               coords = c(coord_x, coord_y), 
                               crs = private$crs)
        
      } else{
        stop("`coord_x` and `coord_y` must be columns of the data when `data_coords` is 'spatial'.")
      }
    }

    if(inherits(data, "sf")){
      if(!inherits(data, "data.frame")){
        stop("No data was found in 'data'")
      }
      data_coords <- "spatial"
      coord_x = ".coord_x"
      coord_y = ".coord_y"

      if(!is.null(private$crs)){
        if(!is.na((sf::st_crs(data)))){
          data <- sf::st_transform(data, crs = private$crs)
        } else{
          warning("The crs of the data is not set, the data will be assumed to be in the same crs as the graph.")
          data <- sf::st_set_crs(data, private$crs)
        }
      }
      coord_tmp <- sf::st_coordinates(data,geometry)
      data <- sf::st_drop_geometry(data)
      data[[".coord_x"]] <- coord_tmp[,1]
      data[[".coord_y"]] <- coord_tmp[,2]
      strc_data <- TRUE
    }

    if("SpatialPointsDataFrame"%in%is(data)){
      data_coords <- "spatial"
      coord_x = ".coord_x"
      coord_y = ".coord_y"
      if(!is.null(private$proj4string)){
        if(!is.na(sp::proj4string(data))){
          data <- sp::spTransform(data,sp::CRS(private$proj4string))
        } else {
          warning("The crs of the data is not set, the data will be assumed to be in the same crs as the graph.")
          sp::proj4string(data) <- private$proj4string
        }
      }
      coord_tmp <- data@coords
      data <- data@data
      data[[".coord_x"]] <- coord_tmp[,1]
      data[[".coord_y"]] <- coord_tmp[,2]
      strc_data <- TRUE
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

        if(is.null(data)){
            stop("No data provided!")
        }

        if(!is.null(data)){
          if(!is.list(data) && !is.data.frame(data)){
            stop("'data' must be either a list or a data.frame!")
          }
        }

        data <- as.list(data)


        if(!strc_data){
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
        if(verbose > 0){
          if(data_coords == "spatial"){
          message("Converting data to PtE")
          if(private$longlat){
            message("This step may take long. If this step is taking too long consider pruning the vertices to possibly obtain some speed up.")
          }
          }
        }

        ## Check data for repeated observations
        if(data_coords == "spatial"){
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
          if(!suppress_warnings){
            warning("There is at least one 'column' of the data with repeated (possibly different) values at the same location for the same group variable. Only one of these values will be used. Consider using the group variable to differentiate between these values or provide different names for such variables.")
          }
        }


        t <- system.time({
            if(data_coords == "PtE"){
                PtE <- cbind(data[[edge_number]], data[[distance_on_edge]])
                if(!normalized){
                  PtE[, 2] <- PtE[,2] / self$edge_lengths[PtE[, 1]]
                }
              } else if(data_coords == "spatial"){
                point_coords <- cbind(data[[coord_x]], data[[coord_y]])
                PtE <- self$coordinates(XY = point_coords)


                if(tolower(duplicated_strategy) == "closest"){
                  norm_XY <- NULL
                  fact <- process_factor_unit(private$vertex_unit, private$length_unit)
                  PtE_new <- NULL
                  far_points <- NULL
                  dup_points <- NULL
                  closest_points <- NULL
                  grp_dat <- data[[".group"]]
                  if(length(grp_dat) == 0){
                    grp_dat <- rep(1, nrow(PtE))
                  }

                  for(grp in unique(grp_dat)){
                      idx_grp <- which(grp_dat == grp)
                      PtE_grp <- PtE[idx_grp,, drop=FALSE]
                      XY_new_grp <- self$coordinates(PtE = PtE_grp, normalized = TRUE)
                      point_coords_grp <- point_coords[idx_grp,, drop=FALSE]
                    # norm_XY <- max(sqrt(rowSums( (point_coords-XY_new)^2 )))
                      norm_XY_grp <- compute_aux_distances(lines = point_coords_grp, points = XY_new_grp, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit, transform = private$transform)
                    # norm_XY <- max(norm_XY)
                    # if(norm_XY > tolerance){
                    #   warning("There was at least one point whose location is far from the graph,
                    #     please consider checking the input.")
                    #   }
                      far_points_grp <- (norm_XY_grp > tolerance)
                      PtE_grp <- PtE_grp[!far_points_grp,,drop=FALSE]
                      XY_new_grp <- XY_new_grp[!far_points_grp,,drop=FALSE]
                      point_coords_grp <- point_coords_grp[!far_points_grp,,drop=FALSE]
                      dup_points_grp <- duplicated(XY_new_grp) | duplicated(XY_new_grp, fromLast=TRUE)
                      norm_XY_grp <- norm_XY_grp[!far_points_grp]
                      if(any(dup_points_grp)){
                        old_new_coords <- cbind(point_coords_grp[dup_points_grp,], XY_new_grp[dup_points_grp,], norm_XY_grp[dup_points_grp], which(dup_points_grp))
                        old_new_coords <- as.data.frame(old_new_coords)
                        colnames(old_new_coords) <- c("coordx", "coordy", "pcoordx", "pcoordy", "dist", "idx")
                        old_new_coords <- dplyr::as_tibble(old_new_coords)
                        old_new_coords <- old_new_coords %>% dplyr::group_by(pcoordx, pcoordy) %>% dplyr::mutate(min_dist = min(dist)) %>% dplyr::mutate(min_idx = dist == min_dist, min_idx = get_only_first(min_idx)) %>% dplyr::ungroup()
                        min_dist_idx <- old_new_coords[["idx"]][!old_new_coords[["min_idx"]]]
                        closest_points_grp <- rep(FALSE, length(dup_points_grp))
                        closest_points_grp[min_dist_idx] <- TRUE
                        norm_XY_grp <- norm_XY_grp[!closest_points_grp]
                        PtE_grp <- PtE_grp[!closest_points_grp,,drop=FALSE]
                        closest_points <- c(closest_points, closest_points_grp)
                      } else{
                        closest_points_grp <- rep(FALSE, length(norm_XY_grp))
                        closest_points <- c(closest_points, closest_points_grp)
                      }

                      dup_points <- c(dup_points, dup_points_grp)
                      far_points <- c(far_points, far_points_grp)
                      PtE_new <- rbind(PtE_new, PtE_grp)
                      norm_XY <- c(norm_XY, norm_XY_grp)
                    }
                    PtE <- PtE_new

                    if(sum(dup_points)>0){
                      if(!suppress_warnings){
                        warning("There were points projected at the same location. Only the closest point was kept. To keep all the observations change 'duplicated_strategy'   to 'jitter'.")
                      }
                    }

                    data <- lapply(data, function(dat){dat[!far_points]})
                    if(!is.null(closest_points)){
                      removed_data <-  lapply(data, function(dat){dat[closest_points]})
                      data <- lapply(data, function(dat){dat[!closest_points]})
                    }
                    if(any(far_points)){
                      if(!suppress_warnings){
                          warning(paste("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance. The total number of points removed due do being far is",sum(far_points)))
                      }
                    }
                    if(include_distance_to_graph){
                      data[[".distance_to_graph"]] <- norm_XY
                    }


                } else if(tolower(duplicated_strategy) == "jitter"){
                    norm_XY <- NULL
                    fact <- process_factor_unit(private$vertex_unit, private$length_unit)
                    PtE_new <- NULL
                    far_points <- NULL

                    grp_dat <- data[[".group"]]
                    if(length(grp_dat) == 0){
                      grp_dat <- rep(1, nrow(PtE))
                    }

                    for(grp in unique(grp_dat)){
                      idx_grp <- which(grp_dat == grp)
                      PtE_grp <- PtE[idx_grp,, drop=FALSE]
                      dup_points_grp <- duplicated(PtE_grp)
                      while(sum(dup_points_grp)>0){
                        cond_0 <- PtE_grp[dup_points_grp,2] == 0
                        cond_1 <- PtE_grp[dup_points_grp,2] == 1
                        idx_pte0 <- which(cond_0)
                        idx_pte1 <-  which(cond_1)
                        idx_other_pte <- which(!cond_0 & !cond_1)
                        d_points <- which(dup_points_grp)
                        if(length(idx_pte0)>0){
                          PtE_grp[d_points[idx_pte0],2] <- PtE_grp[d_points[idx_pte0],2] + 1e-2 * runif(length(idx_pte0))
                        }
                        if(length(idx_pte1)>0){
                          PtE_grp[d_points[idx_pte1],2] <- PtE_grp[d_points[idx_pte1],2] - 1e-2 * runif(length(idx_pte1))
                        }
                        if(length(idx_other_pte)>0){
                          delta_dif <- min(1e-2, 1-max(PtE_grp[d_points[idx_other_pte],2]), min(PtE_grp[d_points[idx_other_pte],2]))
                          PtE_grp[d_points[idx_other_pte],2] <- PtE_grp[d_points[idx_other_pte],2] + delta_dif * runif(length(idx_other_pte))
                        }
                        dup_points_grp <- duplicated(PtE_grp)
                      }
                      XY_new_grp <- self$coordinates(PtE = PtE_grp, normalized = TRUE)
                      point_coords_grp <- point_coords[idx_grp,, drop=FALSE]
                      norm_XY_grp <- compute_aux_distances(lines = point_coords_grp, points = XY_new_grp, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit, transform = private$transform)
                      far_points_grp <- (norm_XY_grp > tolerance)
                      PtE_grp <- PtE_grp[!far_points_grp,,drop=FALSE]
                      norm_XY_grp <- norm_XY_grp[!far_points_grp]
                      far_points <- c(far_points, far_points_grp)
                      PtE_new <- rbind(PtE_new, PtE_grp)
                      norm_XY <- c(norm_XY, norm_XY_grp)
                    }
                    PtE <- PtE_new
                    data <- lapply(data, function(dat){dat[!far_points]})
                    if(any(far_points)){
                      if(!suppress_warnings){
                              warning(paste("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance. The total number of points removed due do being far is",sum(far_points)))
                      }
                    }
                    # PtE <- PtE[!far_points,,drop=FALSE]
                    if(include_distance_to_graph){
                      data[[".distance_to_graph"]] <- norm_XY
                    }
            } else{
                  stop(paste(duplicated_strategy, "is not a valid duplicated strategy!"))
                }

                rm(far_points)
                rm(norm_XY)

            } else{
                stop("The options for 'data_coords' are 'PtE' and 'spatial'.")
            }
        })

    if(only_return_removed){
      if(!is.null(removed_data)){
        return(as.data.frame(removed_data))
      }
    }

      if(verbose == 2) {
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


    length_data <- unique(unlist(lapply(data, length)))

    if(length(length_data)>1){
      stop("There was a problem when processing the data. The number of observations and observation locations (after projecting on the metric graph) are not matching.")
    }

    # Process the data (find all the different coordinates
    # across the different replicates, and also merge the new data to the old data)
    data <- process_data_add_obs(PtE, new_data = data, old_data = NULL,
                                        group_vector, suppress_warnings = suppress_warnings)

    data <- standardize_df_positions(data, self)
    ## convert to Spoints and add
    group_1 <- data[[".group"]]
    group_1 <- which(group_1 == group_1[1])
    PtE <- cbind(data[[".edge_number"]][group_1],
                 data[[".distance_on_edge"]][group_1])
    spatial_points <- self$coordinates(PtE = PtE, normalized = TRUE)
    data[[".coord_x"]] <- rep(spatial_points[,1], times = n_group)
    data[[".coord_y"]] <- rep(spatial_points[,2], times = n_group)

    if(format == "tibble"){
      data <- tidyr::as_tibble(data)
    }

    if(format == "sf"){
      data <- as.data.frame(data)
      data_geometries <- lapply(1:nrow(data), function(i) sf::st_point(as.numeric(data[i, c('.coord_x', '.coord_y')])))
      data <- sf::st_sf(data, geometry = sf::st_sfc(data_geometries), crs = if(!is.null(private$crs)) private$crs else NULL)      
    } 

    class(data) <- c("metric_graph_data", class(data))
        if(!is.null(group)){
      attr(data, "group_variables") <- group
    } else{
      attr(data, "group_variables") <- ".none"
    }
    })

    if(format == "sp"){
      data <- as.data.frame(data)
      sp::coordinates(data) <- ~ .coord_x + .coord_y
    }

          if(verbose == 2) {
      message(sprintf("time: %.3f s", t[["elapsed"]]))
          }
          return(data)

                              },

  #' @description Add observations to the metric graph.
  #' @param data A `data.frame` or named list containing the observations. In
  #' case of groups, the data.frames for the groups should be stacked vertically,
  #' with a column indicating the index of the group. `data` can also be an `sf` object, a
  #' `SpatialPointsDataFrame` object or an `SSN` object.
  #' in which case `data_coords` will automatically be spatial, and there is no need to specify the `coord_x` or `coord_y` arguments.
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
  #' @param data_coords It decides which
  #' coordinate system to use. If `PtE`, the user must provide `edge_number` and
  #' `distance_on_edge`, otherwise if `spatial`, the user must provide
  #' `coord_x` and `coord_y`. The option `euclidean` is `r lifecycle::badge("deprecated")`. Use `spatial` instead.
  #' @param group Vector. If the data is grouped (for example measured at different time
  #' points), this argument specifies the columns (or entries on the list) in
  #' which the group variables are stored. It will be stored as a single column `.group` with the combined entries.
  #' @param group_sep separator character for creating the new group variable when grouping two or more variables.
  #' @param normalized if TRUE, then the distances in `distance_on_edge` are
  #' assumed to be normalized to (0,1). Default FALSE.
  #' @param clear_obs Should the existing observations be removed before adding the data?
  #' @param tibble Should the data be returned as a `tidyr::tibble`?
  #' @param duplicated_strategy Which strategy to handle observations on the same location on the metric graph (that is, if there are two or more observations projected at the same location).
  #' The options are 'closest' and 'jitter'. If 'closest', only the closest observation will be used. If 'jitter', a small perturbation will be performed on the projected observation location. The default is 'closest'.
  #' @param include_distance_to_graph When `data_coord` is 'spatial', should the distance of the observations to the graph be included as a column?
  #' @param tolerance Parameter to control a warning when adding observations.
  #' If the distance of some location and the closest point on the graph is
  #' greater than the tolerance, the function will display a warning.
  #' This helps detecting mistakes on the input locations when adding new data.
  #' @param tolerance_merge tolerance (in edge_length units) for merging points that are very close and are on a common edge. By default, this tolerance is zero, meaning no merges will be performed.
  #' @param merge_strategy The strategies to handle observations that are within the tolerance. The options are `remove`, `merge`, `average`. The default is `merge`, in which one of the observations will be chosen, and the remaining will be used to try to fill all columns with non-NA values. The second strategy is `remove`, meaning that if two observations are within the tolerance one of them will be removed. Finally, `average` will take the average over the close observations for numerical variables, and will choose one non-NA for non-numerical variables.
  #' @param return_removed Should the removed data (if it exists) due to being projected to the same place when using 'closest' `duplicated_strategy`, or due to some merge strategy, be returned?
  #' @param verbose Print progress of the steps when adding observations. There are 3 levels of verbose, level 0, 1 and 2. In level 0, no messages are printed. In level 1, only messages regarding important steps are printed. Finally, in level 2, messages detailing all the steps are printed. The default is 1.
  #' @param suppress_warnings Suppress warnings related to duplicated observations?
  #' @param Spoints `r lifecycle::badge("deprecated")` Use `data` instead.
  #' @return No return value. Called for its side effects. The observations are
  #' stored in the `data` element of the `metric_graph` object.
  add_observations = function(data = NULL,
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
                              Spoints = lifecycle::deprecated()) {

    merge_strategy <- match.arg(merge_strategy, c("remove", "merge", "average"))
    duplicated_strategy <- match.arg(duplicated_strategy, c("closest", "jitter"))
    data_coords <- match.arg(data_coords, c("PtE", "spatial"))

    if(clear_obs){
      df_temp <- data
      self$clear_observations()
      data <- df_temp
    }

    if(inherits(data, "SSN")){
      data <- data$obs
    }

    if(lifecycle::is_present(Spoints)){
           lifecycle::deprecate_warn("1.2.9000", "add_observations(Spoints)", "add_observations(data)",
             details = c("`Spoints` is deprecated, use `data` instead.")
           )
           data <- Spoints
    }

    removed_data <- NULL
    far_data <- NULL
    strc_data <- FALSE

    if(data_coords == "spatial" && 
       !inherits(data, "sf") && 
       !inherits(data, "metric_graph_data") && 
       !"SpatialPointsDataFrame" %in% is(data) && 
       !is.null(private$crs)) {
      
      if(!is.null(data[[coord_x]]) && !is.null(data[[coord_y]])) {
        warning("The crs of the data is not set, the data will be assumed to be in the same crs as the graph.")
        
        # Create sf object from the data
        data_df <- as.data.frame(data)
        data <- sf::st_as_sf(data_df, 
                               coords = c(coord_x, coord_y), 
                               crs = private$crs)
        
      } else{
        stop("`coord_x` and `coord_y` must be columns of the data when `data_coords` is 'spatial'.")
      }
    }

    if(inherits(data, "sf")){
      if(!inherits(data, "data.frame")){
        stop("No data was found in 'data'")
      }
      data_coords <- "spatial"
      coord_x = ".coord_x"
      coord_y = ".coord_y"

      if(!is.null(private$crs)){
        if(!is.na((sf::st_crs(data)))){
          data <- sf::st_transform(data, crs = private$crs)
        } else{
          warning("The crs of the data is not set, the data will be assumed to be in the same crs as the graph.")
          data <- sf::st_set_crs(data, private$crs)
        }
      }
      coord_tmp <- sf::st_coordinates(data,geometry)
      data <- sf::st_drop_geometry(data)
      data[[".coord_x"]] <- coord_tmp[,1]
      data[[".coord_y"]] <- coord_tmp[,2]
      strc_data <- TRUE
    }

    if("SpatialPointsDataFrame"%in%is(data)){
      data_coords <- "spatial"
      coord_x = ".coord_x"
      coord_y = ".coord_y"
      if(!is.null(private$proj4string)){
        if(!is.na(sp::proj4string(data))){
          data <- sp::spTransform(data,sp::CRS(private$proj4string))
        } else {
          warning("The crs of the data is not set, the data will be assumed to be in the same crs as the graph.")
          sp::proj4string(data) <- private$proj4string
        }
      }
      coord_tmp <- data@coords
      data <- data@data
      data[[".coord_x"]] <- coord_tmp[,1]
      data[[".coord_y"]] <- coord_tmp[,2]
      strc_data <- TRUE
    }

    if(length(tolerance)>1){
      tolerance <- tolerance[[1]]
      warning("'tolerance' had more than one element, only the first one will be used.")
    }

    if(verbose>0){
      message("Adding observations...")
      if(data_coords == "PtE"){
        if(normalized){
          message("Assuming the observations are normalized by the length of the edge.")
        } else{
          message("Assuming the observations are NOT normalized by the length of the edge.")
        }
      }
      if(private$longlat){
        message(paste("The unit for edge lengths is", private$length_unit))
        message(paste0("The current tolerance for removing distant observations is (in ",private$length_unit,"): ", tolerance))
      }
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
            stop("No data provided!")
        }

        if(!is.null(data)){
          if(!is.list(data) && !is.data.frame(data)){
            stop("'data' must be either a list or a data.frame!")
          }
        }


        data <- as.list(data)

        if(!strc_data){
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
        if(verbose > 0){
          if(data_coords == "spatial"){
          message("Converting data to PtE")
          if(private$longlat){
            message("This step may take long. If this step is taking too long consider pruning the vertices to possibly obtain some speed up.")
          }
          }
        }



        ## Check data for repeated observations
        if(data_coords == "spatial"){
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
          if(!suppress_warnings){
            warning("There is at least one 'column' of the data with repeated (possibly different) values at the same location for the same group variable. Only one of these values will be used. Consider using the group variable to differentiate between these values or provide different names for such variables.")
          }
        }

        t <- system.time({
            if(data_coords == "PtE"){
                PtE <- cbind(data[[edge_number]], data[[distance_on_edge]])
                if(!normalized){
                  PtE[, 2] <- PtE[,2] / self$edge_lengths[PtE[, 1]]
                  if(any(PtE[,2] > 1)){
                   stop("There were invalid distances on edges. If your data is normalized, please set the 'normalized' argument to TRUE.") 
                  }
                }
              } else if(data_coords == "spatial"){
                point_coords <- cbind(data[[coord_x]], data[[coord_y]])
                PtE <- self$coordinates(XY = point_coords)

                if(tolower(duplicated_strategy) == "closest"){
                  norm_XY <- NULL
                  fact <- process_factor_unit(private$vertex_unit, private$length_unit)
                  PtE_new <- NULL
                  far_points <- NULL
                  dup_points <- NULL
                  closest_points <- NULL
                  grp_dat <- data[[".group"]]
                  if(length(grp_dat) == 0){
                    grp_dat <- rep(1, nrow(PtE))
                  }

                  for(grp in unique(grp_dat)){
                      idx_grp <- which(grp_dat == grp)
                      PtE_grp <- PtE[idx_grp,, drop=FALSE]
                      XY_new_grp <- self$coordinates(PtE = PtE_grp, normalized = TRUE)
                      point_coords_grp <- point_coords[idx_grp,, drop=FALSE]
                    # norm_XY <- max(sqrt(rowSums( (point_coords-XY_new)^2 )))
                      norm_XY_grp <- compute_aux_distances(lines = point_coords_grp, points = XY_new_grp, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit, transform = private$transform)
                    # norm_XY <- max(norm_XY)
                    # if(norm_XY > tolerance){
                    #   warning("There was at least one point whose location is far from the graph,
                    #     please consider checking the input.")
                    #   }
                      far_points_grp <- (norm_XY_grp > tolerance)
                      PtE_grp <- PtE_grp[!far_points_grp,,drop=FALSE]
                      XY_new_grp <- XY_new_grp[!far_points_grp,,drop=FALSE]
                      point_coords_grp <- point_coords_grp[!far_points_grp,,drop=FALSE]
                      dup_points_grp <- duplicated(XY_new_grp) | duplicated(XY_new_grp, fromLast=TRUE)
                      norm_XY_grp <- norm_XY_grp[!far_points_grp]
                      if(any(dup_points_grp)){
                        old_new_coords <- cbind(point_coords_grp[dup_points_grp,], XY_new_grp[dup_points_grp,], norm_XY_grp[dup_points_grp], which(dup_points_grp))
                        old_new_coords <- as.data.frame(old_new_coords)
                        colnames(old_new_coords) <- c("coordx", "coordy", "pcoordx", "pcoordy", "dist", "idx")
                        old_new_coords <- dplyr::as_tibble(old_new_coords)
                        old_new_coords <- old_new_coords %>% dplyr::group_by(pcoordx, pcoordy) %>% dplyr::mutate(min_dist = min(dist)) %>% dplyr::mutate(min_idx = dist == min_dist, min_idx = get_only_first(min_idx)) %>% dplyr::ungroup()
                        min_dist_idx <- old_new_coords[["idx"]][!old_new_coords[["min_idx"]]]
                        closest_points_grp <- rep(FALSE, length(dup_points_grp))
                        closest_points_grp[min_dist_idx] <- TRUE
                        norm_XY_grp <- norm_XY_grp[!closest_points_grp]
                        PtE_grp <- PtE_grp[!closest_points_grp,,drop=FALSE]
                        closest_points <- c(closest_points, closest_points_grp)
                      } else{
                        closest_points_grp <- rep(FALSE, length(norm_XY_grp))
                        closest_points <- c(closest_points, closest_points_grp)
                      }

                      dup_points <- c(dup_points, dup_points_grp)
                      far_points <- c(far_points, far_points_grp)
                      PtE_new <- rbind(PtE_new, PtE_grp)
                      norm_XY <- c(norm_XY, norm_XY_grp)
                    }
                    PtE <- PtE_new

                    if(sum(dup_points)>0){
                        if(!suppress_warnings){
                          warning("There were points projected at the same location. Only the closest point was kept. To keep all the observations change 'duplicated_strategy' to 'jitter'.")
                        }
                    }
                    if(any(far_points)){
                        if(!suppress_warnings){
                          warning(paste("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance. The total number of points removed due do being far is",sum(far_points)))
                        }
                        far_data <- lapply(data, function(dat){dat[far_points]})
                    }                    
                    data <- lapply(data, function(dat){dat[!far_points]})
                    if(!is.null(closest_points)){
                      removed_data <-  lapply(data, function(dat){dat[closest_points]})
                      data <- lapply(data, function(dat){dat[!closest_points]})
                    }                    
                    if(include_distance_to_graph){
                      data[[".distance_to_graph"]] <- norm_XY
                    }


                } else if(tolower(duplicated_strategy) == "jitter"){
                    norm_XY <- NULL
                    fact <- process_factor_unit(private$vertex_unit, private$length_unit)
                    PtE_new <- NULL
                    far_points <- NULL

                    grp_dat <- data[[".group"]]
                    if(length(grp_dat) == 0){
                      grp_dat <- rep(1, nrow(PtE))
                    }

                    for(grp in unique(grp_dat)){
                      idx_grp <- which(grp_dat == grp)
                      PtE_grp <- PtE[idx_grp,, drop=FALSE]
                      dup_points_grp <- duplicated(PtE_grp)
                      while(sum(dup_points_grp)>0){
                        cond_0 <- PtE_grp[dup_points_grp,2] == 0
                        cond_1 <- PtE_grp[dup_points_grp,2] == 1
                        idx_pte0 <- which(cond_0)
                        idx_pte1 <-  which(cond_1)
                        idx_other_pte <- which(!cond_0 & !cond_1)
                        d_points <- which(dup_points_grp)
                        if(length(idx_pte0)>0){
                          PtE_grp[d_points[idx_pte0],2] <- PtE_grp[d_points[idx_pte0],2] + 1e-2 * runif(length(idx_pte0))
                        }
                        if(length(idx_pte1)>0){
                          PtE_grp[d_points[idx_pte1],2] <- PtE_grp[d_points[idx_pte1],2] - 1e-2 * runif(length(idx_pte1))
                        }
                        if(length(idx_other_pte)>0){
                          delta_dif <- min(1e-2, 1-max(PtE_grp[d_points[idx_other_pte],2]), min(PtE_grp[d_points[idx_other_pte],2]))
                          PtE_grp[d_points[idx_other_pte],2] <- PtE_grp[d_points[idx_other_pte],2] + delta_dif * runif(length(idx_other_pte))
                        }
                        dup_points_grp <- duplicated(PtE_grp)
                      }
                      XY_new_grp <- self$coordinates(PtE = PtE_grp, normalized = TRUE)
                      point_coords_grp <- point_coords[idx_grp,, drop=FALSE]
                      norm_XY_grp <- compute_aux_distances(lines = point_coords_grp, points = XY_new_grp, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit, transform = private$transform)
                      far_points_grp <- (norm_XY_grp > tolerance)
                      PtE_grp <- PtE_grp[!far_points_grp,,drop=FALSE]
                      norm_XY_grp <- norm_XY_grp[!far_points_grp]
                      far_points <- c(far_points, far_points_grp)
                      PtE_new <- rbind(PtE_new, PtE_grp)
                      norm_XY <- c(norm_XY, norm_XY_grp)
                    }
                    PtE <- PtE_new
                    if(any(far_points)){
                      if(!suppress_warnings){
                            warning(paste("There were points that were farther than the tolerance. These points were removed. If you want them projected into the graph, please increase the tolerance. The total number of points removed due do being far is",sum(far_points)))
                      }
                      far_data <- lapply(data, function(dat){dat[far_points]})
                    }
                    data <- lapply(data, function(dat){dat[!far_points]})
                    # PtE <- PtE[!far_points,,drop=FALSE]
                    if(include_distance_to_graph){
                      data[[".distance_to_graph"]] <- norm_XY
                    }
                } else{
                  stop(paste(duplicated_strategy, "is not a valid duplicated strategy!"))
                }

                rm(far_points)
                rm(norm_XY)

            } else{
                stop("The options for 'data_coords' are 'PtE' and 'spatial'.")
            }
        })

      if(verbose == 2) {
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

    ## Filtering data that are very close and are on a common edge
    if(length(tolerance_merge)>1){
      warning("tolerance_merge is not of length 1. Only the first element will be used.")
      tolerance_merge <- tolerance_merge[[1]]
    }
    
    if(is.null(tolerance_merge)){
      warning("tolerance_merge is NULL, so it was set to 0.")
      tolerance_merge <- 0
    }

    if(tolerance_merge < 0){
      stop("tolerance_merge cannot be negative.")
    }

    removed_merge <- NULL

    if(tolerance_merge > 0){
        # data, group_vector and PtE
        if(!is.null(group_vector)){
          ord_idx <- order(PtE[,1], PtE[,2], group_vector)
          group_vector <- group_vector[ord_idx]
          PtE <- PtE[ord_idx,,drop=FALSE]      
        } else{
          ord_idx <- order(PtE[,1], PtE[,2])
          PtE <- PtE[ord_idx,,drop=FALSE]
        }

        data <- lapply(data, function(dat){dat[ord_idx]})

        aux_length <- self$edge_lengths[PtE[,1]] * PtE[,2]

        merge_idx <- get_idx_within_merge_tolerance(PtE, group_vector, aux_length, tolerance_merge)

        removed_merge_idx <- setdiff(1:length(ord_idx), merge_idx)
        
        removed_merge <- lapply(data, function(dat){dat[removed_merge_idx]})
        removed_merge[[".edge_number"]] <- PtE[removed_merge_idx, 1]
        removed_merge[[".distance_on_edge"]] <- PtE[removed_merge_idx, 2]

        data <- lapply(data, function(dat){dat[merge_idx]})
        group_vector <- group_vector[merge_idx]
        PtE <- PtE[merge_idx,,drop=FALSE]

        if(merge_strategy %in% c("average", "merge")){
          merge_idx_map <- setNames(seq_along(merge_idx), as.character(merge_idx))
          ref_idx_merges <- find_merged_indices_for_unselected(merge_idx, length(ord_idx))
          data <- apply_merge_strategy(data, removed_merge, merge_idx_map, ref_idx_merges, merge_strategy)
        }
    }

    # Process the data (find all the different coordinates
    # across the different replicates, and also merge the new data to the old data)
    length_data <- unique(unlist(lapply(data, length)))

    if(length(length_data)>1){
      stop("There was a problem when processing the data. The number of observations and observation locations (after projecting on the metric graph) are not matching.")
    }

    private$data <- process_data_add_obs(PtE, new_data = data, private$data,
                                        group_vector, suppress_warnings = suppress_warnings)

    private$data <- standardize_df_positions(private$data, self, edge_number = ".edge_number", distance_on_edge = ".distance_on_edge")

    ## convert to Spoints and add
    PtE <- self$get_PtE()

    spatial_points <- self$coordinates(PtE = PtE, normalized = TRUE)

    private$data[[".coord_x"]] <- rep(spatial_points[,1], times = n_group)
    private$data[[".coord_y"]] <- rep(spatial_points[,2], times = n_group)
    if(tibble){
      private$data <- tidyr::as_tibble(private$data)
    }
    private$group_col <- group
    # distance_graph_tmp <- private$data[[".distance_to_graph"]]

    class(private$data) <- c("metric_graph_data", class(private$data))
    if(!is.null(group)){
      attr(private$data, "group_variables") <- group
    } else{
      attr(private$data, "group_variables") <- ".none"
    }
    # attr(private$data, "distance_to_graph") <- distance_graph_tmp
    })
          if(verbose == 2) {
      message(sprintf("time: %.3f s", t[["elapsed"]]))
          }

    if(return_removed){
      ret_list <- list()
      if(!is.null(removed_data)){
        ret_list[["removed"]] <- as.data.frame(removed_data)
      }
      if(!is.null(far_data)){
        ret_list[["far_data"]] <- as.data.frame(far_data)
      }
      if(!is.null(removed_merge)){
        ret_list[["removed_merge"]] <- as.data.frame(removed_merge)
      }
      if(length(ret_list) > 0){
        return(ret_list)
      } else {
        return(invisible(NULL))
      }
    }
  },

#' @description Use `dplyr::mutate` function on the internal edge weights object.
#' @param ... Arguments to be passed to `dplyr::mutate()`.
#' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
#' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
#' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble".
#' @details A wrapper to use `dplyr::mutate()` on the internal edge weights object and return the result in the requested format.
#' @return A `tidyr::tibble`, `sf` or `sp` object containing the resulting data list after the mutate.
mutate_weights = function(..., .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
    if (!inherits(private$edge_weights, "tbl_df")) {
        edge_weights_res <- tidyr::as_tibble(private$edge_weights)
    } else {
        edge_weights_res <- private$edge_weights
    }

    if (.drop_all_na) {
        is_tbl <- inherits(edge_weights_res, "tbl_df")
        idx_temp <- idx_not_all_NA(edge_weights_res)
        edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        if (is_tbl) {
            edge_weights_res <- tidyr::as_tibble(edge_weights_res)
        }
    }

    if (.drop_na) {
        if (!inherits(edge_weights_res, "tbl_df")) {
            idx_temp <- idx_not_any_NA(edge_weights_res)
            edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        } else {
            edge_weights_res <- tidyr::drop_na(edge_weights_res)
        }
    }

    edge_weights_res <- dplyr::mutate(.data = edge_weights_res, ...)

    if (!inherits(edge_weights_res, "metric_graph_weights")) {
        class(edge_weights_res) <- c("metric_graph_weights", class(edge_weights_res))
    }

    return(private$format_weights(edge_weights_res, format))
},


#' @description Use `dplyr::select` function on the internal edge weights object.
#' @param ... Arguments to be passed to `dplyr::select()`.
#' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
#' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
#' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble".
#' @details A wrapper to use `dplyr::select()` on the internal edge weights object and return the result in the requested format.
#' @return A `tidyr::tibble`, `sf` or `sp` object containing the resulting data list after the select.
select_weights = function(..., .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
    if (!inherits(private$edge_weights, "tbl_df")) {
        edge_weights_res <- tidyr::as_tibble(private$edge_weights)
    } else {
        edge_weights_res <- private$edge_weights
    }

    edge_weights_res <- dplyr::select(.data = edge_weights_res, ...)

    if (.drop_all_na) {
        is_tbl <- inherits(edge_weights_res, "tbl_df")
        idx_temp <- idx_not_all_NA(edge_weights_res)
        edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        if (is_tbl) {
            edge_weights_res <- tidyr::as_tibble(edge_weights_res)
        }
    }

    if (.drop_na) {
        if (!inherits(edge_weights_res, "tbl_df")) {
            idx_temp <- idx_not_any_NA(edge_weights_res)
            edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        } else {
            edge_weights_res <- tidyr::drop_na(edge_weights_res)
        }
    }

    if (!inherits(edge_weights_res, "metric_graph_weights")) {
        class(edge_weights_res) <- c("metric_graph_weights", class(edge_weights_res))
    }

    return(private$format_weights(edge_weights_res, format))
},


#' @description Use `dplyr::filter` function on the internal edge weights object.
#' @param ... Arguments to be passed to `dplyr::filter()`.
#' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
#' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
#' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble".
#' @details A wrapper to use `dplyr::filter()` on the internal edge weights object and return the result in the requested format.
#' @return A `tidyr::tibble`, `sf` or `sp` object containing the resulting data list after the filter.
filter_weights = function(..., .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
    if (!inherits(private$edge_weights, "tbl_df")) {
        edge_weights_res <- tidyr::as_tibble(private$edge_weights)
    } else {
        edge_weights_res <- private$edge_weights
    }

    if (.drop_all_na) {
        is_tbl <- inherits(edge_weights_res, "tbl_df")
        idx_temp <- idx_not_all_NA(edge_weights_res)
        edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        if (is_tbl) {
            edge_weights_res <- tidyr::as_tibble(edge_weights_res)
        }
    }

    if (.drop_na) {
        if (!inherits(edge_weights_res, "tbl_df")) {
            idx_temp <- idx_not_any_NA(edge_weights_res)
            edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        } else {
            edge_weights_res <- tidyr::drop_na(edge_weights_res)
        }
    }

    edge_weights_res <- dplyr::filter(.data = edge_weights_res, ...)

    if (!inherits(edge_weights_res, "metric_graph_weights")) {
        class(edge_weights_res) <- c("metric_graph_weights", class(edge_weights_res))
    }

    return(private$format_weights(edge_weights_res, format))
},


#' @description Use `dplyr::summarise` function on the internal edge weights object grouped by the edge numbers.
#' @param ... Arguments to be passed to `dplyr::summarise()`.
#' @param .groups A vector of strings containing the names of the columns to be grouped, when computing the summaries. The default is `NULL`.
#' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
#' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
#' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble".
#' @details A wrapper to use `dplyr::summarise()` on the internal edge weights object and return the result in the requested format.
#' @return A `tidyr::tibble`, `sf` or `sp` object containing the resulting data list after the summarise.
summarise_weights = function(..., .groups = NULL, .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
    if (!inherits(private$edge_weights, "tbl_df")) {
        edge_weights_res <- tidyr::as_tibble(private$edge_weights)
    } else {
        edge_weights_res <- private$edge_weights
    }

    if (.drop_all_na) {
        is_tbl <- inherits(edge_weights_res, "tbl_df")
        idx_temp <- idx_not_all_NA(edge_weights_res)
        edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        if (is_tbl) {
            edge_weights_res <- tidyr::as_tibble(edge_weights_res)
        }
    }

    if (.drop_na) {
        if (!inherits(edge_weights_res, "tbl_df")) {
            idx_temp <- idx_not_any_NA(edge_weights_res)
            edge_weights_res <- lapply(edge_weights_res, function(dat) { dat[idx_temp] })
        } else {
            edge_weights_res <- tidyr::drop_na(edge_weights_res)
        }
    }

    edge_weights_res <- dplyr::group_by_at(.tbl = edge_weights_res, .vars = .groups)
    edge_weights_res <- dplyr::summarise(.data = edge_weights_res, ...)
    edge_weights_res <- dplyr::ungroup(edge_weights_res)

    if (!inherits(edge_weights_res, "metric_graph_weights")) {
        class(edge_weights_res) <- c("metric_graph_weights", class(edge_weights_res))
    }

    return(private$format_weights(edge_weights_res, format))
},

#' @description Use `tidyr::drop_na()` function on the internal edge weights object.
#' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble".
#' @param ... Arguments to be passed to `tidyr::drop_na()`.
#' @details A wrapper to use `tidyr::drop_na()` within the internal edge weights object.
#' @return A `tidyr::tibble`, `sf`, or `sp` object containing the resulting data list after the drop_na.
drop_na_weights = function(...,format = "tibble") {
    if (!inherits(private$edge_weights, "tbl_df")) {
        edge_weights_res <- tidyr::as_tibble(private$edge_weights)
    } else {
        edge_weights_res <- private$edge_weights
    }

    edge_weights_res <- tidyr::drop_na(data = edge_weights_res, ...)

    if (!inherits(edge_weights_res, "metric_graph_weights")) {
        class(edge_weights_res) <- c("metric_graph_weights", class(edge_weights_res))
    }

    return(private$format_weights(edge_weights_res, format))
},


 #' @description Use `dplyr::mutate` function on the internal metric graph data object.
#' @param ... Arguments to be passed to `dplyr::mutate()`.
#' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
#' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
#' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble".
#' @details A wrapper to use `dplyr::mutate()` within the internal metric graph data object and return the result in the requested format.
#' @return A `tidyr::tibble`, `sf`, or `sp` object containing the resulting data list after the mutate.
mutate = function(..., .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
    if (!inherits(private$data, "tbl_df")) {
        data_res <- tidyr::as_tibble(private$data)
    } else {
        data_res <- private$data
    }

    if (.drop_all_na) {
        is_tbl <- inherits(data_res, "tbl_df")
        idx_temp <- idx_not_all_NA(data_res)
        data_res <- lapply(data_res, function(dat) { dat[idx_temp] })
        if (is_tbl) {
            data_res <- tidyr::as_tibble(data_res)
        }
    }

    if (.drop_na) {
        if (!inherits(data_res, "tbl_df")) {
            idx_temp <- idx_not_any_NA(data_res)
            data_res <- lapply(data_res, function(dat) { dat[idx_temp] })
        } else {
            data_res <- tidyr::drop_na(data_res)
        }
    }

    data_res <- dplyr::mutate(.data = data_res, ...)

    if (!inherits(data_res, "metric_graph_data")) {
        class(data_res) <- c("metric_graph_data", class(data_res))
    }

    return(private$format_data(data_res, format))
},

  #' @description Use `tidyr::drop_na()` function on the internal metric graph data object.
  #' @param ... Arguments to be passed to `tidyr::drop_na()`.
  #' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble".
  #' @details A wrapper to use `dplyr::drop_na()` within the internal metric graph data object.
  #' @return A `tidyr::tibble` object containing the resulting data list after the drop_na.
  drop_na = function(..., format = "tibble") {
    if(!inherits(private$data, "tbl_df")){
      data_res <- tidyr::as_tibble(private$data)
    } else{
      data_res <- private$data
    }

    data_res <- tidyr::drop_na(data = data_res, ...)

    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }

    return(private$format_data(data_res, format))
  },



  #' @description Use `dplyr::select` function on the internal metric graph data object.
  #' @param ... Arguments to be passed to `dplyr::select()`.
  #' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
  #' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
  #' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble". 
  #' @details A wrapper to use `dplyr::select()` within the internal metric graph data object. Observe that it is a bit different from directly using `dplyr::select()` since it does not allow to remove the internal positions that are needed for the metric_graph methods to work.
  #' @return A `tidyr::tibble` object containing the resulting data list after the selection.
  select = function(..., .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
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
    return(private$format_data(data_res, format))
  },

    #' @description Use `dplyr::filter` function on the internal metric graph data object.
  #' @param ... Arguments to be passed to `dplyr::filter()`.
   #' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
  #' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble". 
  #' @details A wrapper to use `dplyr::filter()` within the internal metric graph data object.
  #' @return A `tidyr::tibble` object containing the resulting data list after the filter.
  filter = function(..., .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
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

    return(private$format_data(data_res, format))
  },


  #' @description Use `dplyr::summarise` function on the internal metric graph data object grouped by the spatial locations and the internal group variable.
  #' @param ... Arguments to be passed to `dplyr::summarise()`.
  #' @param .include_graph_groups Should the internal graph groups be included in the grouping variables? The default is `FALSE`. This means that, when summarising, the data will be grouped by the internal group variable together with the spatial locations.
  #' @param .groups A vector of strings containing the names of the columns to be additionally grouped, when computing the summaries. The default is `NULL`.
   #' @param .drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param .drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
   #' @param format The format of the output: "tibble", "sf", or "sp". Default is "tibble". 
  #' @details A wrapper to use `dplyr::summarise()` within the internal metric graph data object grouped by manually inserted groups (optional), the internal group variable (optional) and the spatial locations. Observe that if the integral group variable was not used as a grouping variable for the summarise, a new column, called `.group`, will be added, with the same value 1 for all rows.
  #' @return A `tidyr::tibble` object containing the resulting data list after the summarise.
  summarise = function(..., .include_graph_groups = FALSE, .groups = NULL, .drop_na = FALSE, .drop_all_na = TRUE, format = "tibble") {
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

    return(private$format_data(data_res, format))

  },

 #' @description Return the internal data with the option to filter by groups.
 #' @param group A vector contaning which groups should be returned? The default is `NULL`, which gives the result for the all groups.
 #' @param format Which format should the data be returned? The options are `tibble` for `tidyr::tibble`, `sf` for `POINT`, `sp` for `SpatialPointsDataFrame` and `list` for the internal list format.
 #' @param drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
 #' @param drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
 #' @param tibble `r lifecycle::badge("deprecated")` Use `format` instead.

  get_data = function(group = NULL, format = c("tibble", "sf", "sp", "list"), drop_na = FALSE, drop_all_na = TRUE, tibble = deprecated()){
    if(is.null(private$data)){
      stop("The graph does not contain data.")
    }

  if (lifecycle::is_present(tibble)) {
      lifecycle::deprecate_warn("1.3.0.9000", "get_edge_weights(tibble)", "get_edge_weights(format)",
        details = c("The argument `tibble` was deprecated in favor of the argument `format`.")
      )
    if(tibble){
      format <- "tibble"
    }
  }       

    format <- format[[1]]

    format <- tolower(format)

    if(!(format %in% c("tibble", "sf", "sp", "list"))){
      stop("The possible formats are 'tibble', 'sf', 'sp' and 'list'.")
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


    if(format == "tibble"){
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

    if(format == "sf"){
      data_temp <- as.data.frame(data_temp)
      data_geometries <- lapply(1:nrow(data_temp), function(i) sf::st_point(as.numeric(data_temp[i, c('.coord_x', '.coord_y')])))
      data_temp <- sf::st_sf(data_temp, geometry = sf::st_sfc(data_geometries), crs = if(!is.null(private$crs)) private$crs else NULL)      
    } 

    if(!inherits(data_temp, "metric_graph_data")){
      class(data_temp) <- c("metric_graph_data", class(data_temp))
    }

    if(format == "sp"){
      data_temp <- as.data.frame(data_temp)
      sp::coordinates(data_temp) <- ~ .coord_x + .coord_y
    }

    return(data_temp)
  },
 #' @description Define the columns to be used for creating the directional vertex
 #' weights. Also possible to supply user defined functions for input and output
 #' to create ones own weights.
 #' @param f_in functions for the input vertex (default `w/sum(w)`) uses the columns of name_column
 #' @param f_out functions for the output vertex (deafult `rep(-1,length(w))`) uses the columns of name_column
 #' @details For more details see paper (that does not exists yet).
 #' @return No return value.
 setDirectionalWeightFunction = function(f_in = NULL,
                                         f_out = NULL){


  if(!is.null(self$C)){
    warning('The constraint matrix has been deleted')
  }
   self$C     = NULL
   self$CoB   = NULL
   self$CoB$T = NULL

   if(!is.null(f_in)){
     self$DirectionalWeightFunction_in = f_in
   }else{
     self$DirectionalWeightFunction_in  = function(w){ w/sum(w)}
   }
   if(!is.null(f_out)){
     self$DirectionalWeightFunction_out = f_out
   }else{
     self$DirectionalWeightFunction_out = function(w){ rep(-1,length(w))}
   }
  },

 #' @description Build directional ODE constraint matrix from edges.
 #' @param alpha how many derivatives the processes has
 #' @param weight weighting for each vertex used in the constraint (E x 2)
 #' @details Currently not implemented for circles (edges that start and end
 #' in the same vertex)
 #' @return No return value. Called for its side effects.
  buildDirectionalConstraints = function(alpha = 1){
    if(alpha %% 1 != 0){
      stop("alpha should be an integer")
    }

    weight <- self$get_edge_weights()
    weight <- as.vector(weight[[private$directional_weights]])
    V_indegree = self$get_degrees("indegree")
    V_outdegree = self$get_degrees("outdegree")
    # index_outdegree <- V_outdegree > 0 & V_indegree >0
    # index_in0      <- V_indegree == 0
    # nC = (sum(V_outdegree[index_outdegree] *(1 + V_indegree[index_outdegree])) + sum(V_outdegree[index_in0]-1)) * alpha
    # i_  =  rep(0, nC)
    # j_  =  rep(0, nC)
    # x_  =  rep(0, nC)
    # Vs <- which(index_outdegree)
    # count_constraint <- 0
    # count <- 0
    # for (v in Vs) {
    #   out_edges   <- which(self$E[, 1] %in% v)
    #   in_edges    <- which(self$E[, 2] %in% v)
    #   #for each out edge
    #   n_in <- length(in_edges)
    #   for(i in 1:length(out_edges)){
    #     for(der in 1:alpha){
    #       i_[count + 1:(n_in+1)] <- count_constraint + 1
    #       j_[count + 1:(n_in+1)] <- c(2 * alpha * (out_edges[i]-1) + der,
    #                                   2 * alpha * (in_edges-1)  + alpha + der)


    #       x_[count + 1:(n_in+1)] <- c(as.matrix(self$DirectionalWeightFunction_out(weight[out_edges[i]])),
    #                                   as.matrix(self$DirectionalWeightFunction_in(weight[in_edges])))

    #       count <- count + (n_in+1)
    #       count_constraint <- count_constraint + 1
    #     }
    #   }
    # }
    # Vs0 <- which(index_in0)
    # for (v in Vs0) {
    #   out_edges   <- which(self$E[, 1] %in% v)
    #   #for each out edge
    #   if(length(out_edges)>1){
    #     for(i in 2:length(out_edges)){
    #       for(der in 1:alpha){
    #         i_[count + 1:2] <- count_constraint + 1
    #         j_[count + 1:2] <- c(2 * alpha * (out_edges[i]-1) + der,
    #                              2 * alpha * (out_edges[i-1]-1)   + der)

    #         x_[count + 1:2] <- c(1,
    #                              -1)
    #         count <- count + 2
    #         count_constraint <- count_constraint + 1
    #       }
    #     }
    #   }
    # }
    # C <- Matrix::sparseMatrix(i = i_[1:count],
    #                           j = j_[1:count],
    #                           x = x_[1:count],
    #                           dims = c(count_constraint, 2*alpha*self$nE))
    # self$C = C
    temp_E <- apply(self$E,2,as.integer)
    self$C <-construct_directional_constraint_matrix(E = temp_E, nV = as.integer(self$nV), nE = as.integer(self$nE), alpha = as.integer(alpha),
    V_indegree = as.integer(V_indegree), V_outdegree = as.integer(V_outdegree), weight = weight, 
    DirectionalWeightFunction_out = self$DirectionalWeightFunction_out, 
    DirectionalWeightFunction_in = self$DirectionalWeightFunction_in)

    self$CoB <- c_basis2(self$C)
    self$CoB$T <- t(self$CoB$T)
    self$CoB$alpha <- 1
  },


  #' @description Build Kirchoff constraint matrix from edges.
  #' @param alpha the type of constraint (currently only supports 2)
  #' @param edge_constraint if TRUE, add constraints on vertices of degree 1
  #' @details Currently not implemented for circles (edges that start and end
  #' in the same vertex)
  #' @return No return value. Called for its side effects.
  buildC = function(alpha = 2, edge_constraint = FALSE) {

    if(alpha==2){
      temp_E <- apply(self$E,2,as.integer)

      self$C <- construct_constraint_matrix(temp_E, as.integer(self$nV), as.integer(edge_constraint))
      self$CoB <- c_basis2(self$C)
      self$CoB$T <- t(self$CoB$T)
      self$CoB$alpha <- 2
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

build_mesh = function(h = NULL, n = NULL, continuous = TRUE,
                         continuous.outs = FALSE, continuous.deg2 = FALSE) {
    if (is.null(h) && is.null(n)) {
      stop("You should specify either h or n!")
    }

    if (!is.null(h)) {
      if (length(h) > 1 || !is.numeric(h)) stop("h should be a single number")
      if (h <= 0) stop("h must be positive!")
    }

    if (!is.null(n)) {
      if (length(n) > 1 || !is.numeric(n)) stop("n should be a single number")
      if (n <= 0) stop("n must be positive!")
      if (n %% 1 != 0) {
        warning("A noninteger n was given, we are rounding it to an integer.")
        n <- round(n)
      }
    }

    mesh <- list(PtE = NULL, V = NULL, E = NULL, n_e = integer(length(self$edges)), h_e = NULL, ind = NULL, VtE = NULL)
    mesh$n_e <- integer(length(self$edges))
    attr(mesh, "continuous") <- continuous

    edge_lengths <- self$edge_lengths
    n_edges <- length(self$edges)

    if (continuous) {
      mesh$V <- self$V
      mesh$ind <- seq_len(self$nV)

      if (is.null(n)) {
        mesh$n_e <- ceiling(edge_lengths / h) + 1 - 2
      } else {
        mesh$n_e <- rep(n, n_edges)
      }

      # Call the Rcpp function
      mesh_data <- generate_mesh(n_edges, edge_lengths, mesh$n_e, self$E, mesh$ind, continuous)

      mesh$PtE <- cbind(mesh_data$PtE_edge, mesh_data$PtE_pos)
      mesh$h_e <- mesh_data$h_e
      mesh$E <- cbind(mesh_data$E_start, mesh_data$E_end)

      mesh$VtE <- rbind(self$VtEfirst(), mesh$PtE)
      if (!is.null(mesh$PtE) && nrow(mesh$PtE) > 0) {
        mesh$V <- rbind(mesh$V, self$coordinates(PtE = mesh$PtE))
      }
      self$mesh <- mesh
    } else {
      mesh$ind <- 0

      if (is.null(n)) {
        mesh$n_e <- ceiling(edge_lengths / h) + 1
      } else {
        mesh$n_e <- rep(n + 2, n_edges)
      }

      # Call the Rcpp function for the non-continuous case
      mesh_data <- generate_mesh(n_edges, edge_lengths, mesh$n_e, self$E, mesh$ind, continuous = FALSE)

      mesh$PtE <- cbind(mesh_data$PtE_edge, mesh_data$PtE_pos)
      mesh$h_e <- mesh_data$h_e
      mesh$E <- cbind(mesh_data$E_start, mesh_data$E_end)

      mesh$VtE <- mesh$PtE
      if(!is.null(mesh$PtE) && nrow(mesh$PtE) > 0) {
        mesh$V <- self$coordinates(PtE = mesh$PtE)
      }
      self$mesh <- mesh
      if (continuous.outs) private$mesh_merge_outs()
      private$move_V_first()
      if (continuous.deg2) private$mesh_merge_deg2()
    }
  }, 

  #' @description Get the version of MetricGraph package used to build the graph
  #' @return A character string with the version number
  get_version = function() {
    return(private$version)
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

    # if(!all(private$get_edge_weights_internal()==1)){
    if(!is.null(private$kirchhoff_weights)){
      for(i in 1:self$nV) {
        if(attr(self$vertices[[i]],"degree") > 1) {
          edges.i <- which(rowSums(self$E==i)>0)
          edges.mesh <- which(rowSums(self$mesh$E==i)>0)
          w <- rep(0,length(edges.mesh))
          h <- rep(0,length(edges.mesh))
          for(j in 1:length(edges.mesh)) {
            V.e <- self$mesh$E[edges.mesh[j],which(self$mesh$E[edges.mesh[j],]!=i)]
            if(length(V.e) == 0){
              V.e <- self$mesh$E[edges.mesh[j],1]
            }
            E.e <- self$mesh$VtE[V.e,1] #the edge the mesh node is on
            # w[j] <- attr(self$edges[[E.e]],"weight")
            kw <- attr(self$edges[[E.e]], "kirchhoff_weight")
            w_tmp <- attr(self$edges[[E.e]],"weight")
            w[j] <- w_tmp[[kw]]
            h[j] <- self$mesh$h_e[edges.mesh[j]]
          }
          for(j in 2:attr(self$vertices[[i]],"degree")){
            self$mesh$K[i,i] <- self$mesh$K[i,i] +  (w[j]/w[1] - 1)/h[j]
          }
        }
      }
    }
    # }

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
  #' number and `newdata` is `NULL`, it will be the index of the group as stored internally and if `newdata` is provided, it will be the index of the group stored in `newdata`. If `group`
  #' is a character, then the group will be chosen by its name.
  #' @param type The type of plot to be returned. The options are `ggplot` (the default), that uses `ggplot2`; `plotly` that uses `plot_ly` for 3D plots, which requires the `plotly` package, and `mapview` that uses the `mapview` function, to build interactive plots, which requires the `mapview` package.
  #' @param interactive Only works for 2d plots. If `TRUE`, an interactive plot will be displayed. Unfortunately, `interactive` is not compatible with `edge_weight` if `add_new_scale_weights` is TRUE.
  #' @param vertex_size Size of the vertices.
  #' @param vertex_color Color of vertices.
  #' @param edge_width Line width for edges. If `edge_width_weight` is not `NULL`, this determines the maximum edge width.
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
  #' @param direction Show the direction of the edges? For `type == "mapview"` the arrows are not shown, only the color of the vertices indicating whether they are problematic or not.
  #' @param arrow_size The size of the arrows if direction is TRUE.
  #' @param edge_weight Which column from edge weights to determine the colors of the edges? If `NULL` edge weights are not plotted. To plot the edge weights when the metric graph `edge_weights` is a vector instead of a `data.frame`, simply set to 1.
  #' `edge_weight` is only available for 2d plots. For 3d plots with edge weights, please use the `plot_function()` method.
    #' @param edge_width_weight Which column from edge weights to determine the edges widths? If `NULL` edge width will be determined from `edge_width`. Currently it is not supported for `type = "mapview"`.
    #' @param scale_color_main Color scale for the data to be plotted.
    #' @param scale_color_weights Color scale for the edge weights. Will only be used if `add_new_scale_weights` is TRUE.
    #' @param scale_color_main_discrete Color scale for the data to be plotted, for discrete data.
    #' @param scale_color_weights_discrete Color scale for discrete edge weights. Will only be used if `add_new_scale_weights` is TRUE.
    #' @param scale_color_degree Color scale for the degrees.
    #' @param add_new_scale_weights Should a new color scale for the edge weights be created?
    #' @param scale_color_mapview Color scale to be applied for data when `type = "mapview"`.
    #' @param scale_color_weights_mapview Color scale to be applied for edge weights when `type = "mapview"`.
    #' @param scale_color_weights_discrete_mapview Color scale to be applied for degrees when `type = "mapview"`. If `NULL` `RColorBrewer::brewer.pal(n = n_weights, "Set1")` will be used where `n_weights` is the number of different degrees.
    #' @param scale_color_degree_mapview Color scale to be applied for degrees when `type = "mapview"`. If `NULL` `RColorBrewer::brewer.pal(n = n_degrees, "Set1")` will be used where `n_degrees` is the number of different degrees.
    #' @param plotly  `r lifecycle::badge("deprecated")` Use `type` instead.
##  # ' @param mutate A string containing the commands to be passed to `dplyr::mutate` function in order to obtain new variables as functions of the existing variables.
##  # ' @param filter A string containing the commands to be passed to `dplyr::filter` function in order to obtain new filtered data frame.
##  # ' @param summarise A string containing the commands to be passed to `dplyr::summarise` function in order to obtain new  data frame containing the summarised variable.
##  # ' @param summarise_group_by A vector of strings containing the names of the columns to be additionally grouped, when computing the summaries. The default is `NULL`.
##  # ' @param summarise_by_graph_group Should the internal graph groups be included in the grouping variables? The default is `FALSE`. This means that, when summarising, the data will be grouped by the internal group variable together with the spatial locations.
  #' @param ... Additional arguments to pass to `ggplot()` or `plot_ly()`
  #' @return A `plot_ly` (if `type = "plotly"`) or `ggplot` object.
  plot = function(data = NULL,
                  newdata = NULL,
                  group = 1,
                  type = c("ggplot", "plotly", "mapview"),
                  interactive = FALSE,
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
                  ...) {

 if (lifecycle::is_present(plotly)) {
      lifecycle::deprecate_warn("1.3.0.9000", "plot(plotly)", "plot(type)",
        details = c("The argument `plotly` was deprecated in favor of the argument `type`.")
      )
    if(plotly){
      type <- "plotly"
    }
  }                 

    if(!is.null(data) && is.null(private$data) && is.null(newdata)) {
      stop("The graph does not contain data.")
    }
    if(is.numeric(group) && !is.null(data) && is.null(newdata)) {
      unique_group <- unique(private$data[[".group"]])
      group <- unique_group[group]
    } else if(!is.null(newdata) && is.numeric(group)){
      unique_group <- unique(newdata[[".group"]])
      group <- unique_group[group]
    }
    if(!is.null(newdata)){
      if(!inherits(newdata, "metric_graph_data")){
        stop("'newdata' must be of class 'metric_graph_data'!")
      }
    }

    type <- type[[1]]

    if(!(type %in% c("ggplot", "plotly", "mapview"))){
      stop("type must be one of the following: 'ggplot', 'plotly' or 'mapview'.")
    }

    if(type == "ggplot") {
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
                           edge_weight = edge_weight,
                           edge_width_weight = edge_width_weight,
                           scale_color_main = scale_color_main,
                           scale_color_weights = scale_color_weights,
                           scale_color_degree = scale_color_degree,
                           add_new_scale_weights = add_new_scale_weights,
                           arrow_size = arrow_size,
                           scale_color_main_discrete = scale_color_main_discrete,
                           scale_color_weights_discrete = scale_color_weights_discrete,                           
                           ...)
      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degree" && !private$transform){
          p <- p + labs(x = "Longitude",  y = "Latitude")
        } else{
          p <- p + labs(x = paste0("x (in ",private$vertex_unit, ")"),  y = paste0("y (in ",private$vertex_unit, ")"))
        }
      }
    } else if(type == "plotly") {
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
                           edge_width_weight = edge_width_weight,
                           ...)
      p <-  plotly::layout(p, scene = list(xaxis = list(autorange = "reversed")))                             
      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degree" && !private$transform){
          p <- plotly::layout(p, scene = list(xaxis = list(title = "Longitude"), yaxis = list(title = "Latitude")))
        } else{
          p <- plotly::layout(p, scene = list(xaxis = list(title = paste0("x (in ",private$vertex_unit, ")")), yaxis = list(title = paste0("y (in ",private$vertex_unit, ")"))))
        }
      }
    } else if(type == "mapview"){
requireNamespace("mapview")
edge_width <- 2.5 * edge_width
edges_sf <- self$get_edges(format = "sf")
if(is.null(newdata) && !is.null(data)){
  data_sf <- self$get_data(format = "sf")
  idx_grp <- (data_sf[[".group"]] == group)
  data_sf <- data_sf[idx_grp, , drop=FALSE]
  class(data_sf) <- setdiff(class(data_sf), "metric_graph_data")
} else if (!is.null(data)){
  data_sf <- as.data.frame(newdata)
  data_geometries <- lapply(1:nrow(data_sf), function(i) sf::st_point(as.numeric(data_sf[i, c('.coord_x', '.coord_y')])))
  data_sf <- sf::st_sf(data_sf, geometry = sf::st_sfc(data_geometries), crs = if(!is.null(private$crs)) private$crs else NULL)     
  idx_grp <- (data_sf[[".group"]] == group)
  data_sf <- data_sf[idx_grp, , drop=FALSE]   
  class(data_sf) <- setdiff(class(data_sf), "metric_graph_data")        
}

mapview_output <- p

if(is.null(p)){
  if (!is.null(edge_weight)) {
    if (!(edge_weight %in% colnames(edges_sf))) {
      stop(paste(edge_weight, "is not a valid column in edges_sf"))
    }

    if(is.character(edges_sf[[edge_weight]]) || is.factor(edges_sf[[edge_weight]])){
      if(is.null(scale_color_degree_mapview)){
        scale_color_weights_discrete_mapview <- RColorBrewer::brewer.pal(n = length(levels(edges_sf[[edge_weight]])), "Set1")
      }         
      scale_weights <- scale_color_weights_discrete_mapview
    } else{
      scale_weights <- scale_color_weights_mapview
    }
    mapview_output <- mapview::mapview(
      x = edges_sf,
      zcol = edge_weight,  
      color = scale_weights,  
      lwd = edge_width,  
      layer.name = "Edges",  
      col.regions = scale_weights,  
      ...
    )
  } else {
    mapview_output <- mapview::mapview(
      x = edges_sf,
      lwd = edge_width,  
      color = edge_color, 
      layer.name = "Edges",
      ...
    )
  }
} else{
    if (!is.null(edge_weight)) {
    if (!(edge_weight %in% colnames(edges_sf))) {
      stop(paste(edge_weight, "is not a valid column in edges_sf"))
    }
    if(is.character(edges_sf[[edge_weight]]) || is.factor(edges_sf[[edge_weight]])){
      if(is.null(scale_color_degree_mapview)){
        scale_color_weights_discrete_mapview <- RColorBrewer::brewer.pal(n = length(levels(edges_sf[[edge_weight]])), "Set1")
      }      
      scale_weights <- scale_color_weights_discrete_mapview
    } else{
      scale_weights <- scale_color_weights_mapview
    }    
    mapview_output <- mapview_output + mapview::mapview(
      x = edges_sf,
      zcol = edge_weight,  
      color = scale_weights,  
      lwd = edge_width,  
      layer.name = "Edges",  
      col.regions = scale_weights,  
      ...
    )
  } else {
    mapview_output <- mapview_output + mapview::mapview(
      x = edges_sf,
      lwd = edge_width,  
      color = edge_color, 
      layer.name = "Edges",
      ...
    )
  }
}



if (degree) {
  vertices_sf <- self$get_vertices(format = "sf")
  vertices_sf$degree <- as.factor(vertices_sf$degree)
  if(is.null(scale_color_degree_mapview)){
    scale_color_degree_mapview <- RColorBrewer::brewer.pal(n = length(levels(vertices_sf$degree)), "Set1")
  }

  mapview_output <- mapview_output + mapview::mapview(
    x = vertices_sf,
    zcol = "degree",
    cex = vertex_size,
    col.regions = scale_color_degree_mapview,
    color = scale_color_degree_mapview,
    layer.name = "Vertices (Degree)",
    ...
  )
} else if (direction) {
  vertices_sf <- self$get_vertices(format = "sf")
  problematic_vertices <- vertices_sf[vertices_sf$problematic == TRUE, ]
  non_problematic_vertices <- vertices_sf[vertices_sf$problematic == FALSE, ]

  mapview_output <- mapview_output + mapview::mapview(
    x = problematic_vertices,
    cex = vertex_size,
    col.regions = "red",
    color = "red",
    layer.name = "Problematic Vertices",
    ...
  )

  mapview_output <- mapview_output + mapview::mapview(
    x = non_problematic_vertices,
    cex = vertex_size,
    col.regions = "green",
    color = "green",
    layer.name = "Non-problematic Vertices",
    ...
  )
} else if(vertex_size > 0){
  vertices_sf <- self$get_vertices(format = "sf")
  mapview_output <- mapview_output + mapview::mapview(
    x = vertices_sf,
    cex = vertex_size,
    col.regions = vertex_color,
    color = vertex_color,
    layer.name = "Vertices",
    ...
  )
}

if (!is.null(data)) {
  if (!(data %in% names(data_sf))) {
    stop(paste(data, "is not an existing column name in the dataset."))
  }
  data_size <- 2 * data_size

  mapview_output <- mapview_output + mapview::mapview(
    x = data_sf,
    zcol = data,
    cex = data_size,
    col.regions = scale_color_mapview,
    color = scale_color_mapview,
    layer.name = data,
    ...
  )
}

if(vertex_size > 0){
  if (mesh) {
    if (is.null(self$mesh)) {
      stop("The metric graph does not contain a mesh.")
    }

    mesh_df <- as.data.frame(self$mesh$V)
    colnames(mesh_df) <- c("X", "Y")
    mesh_sf <- sf::st_as_sf(mesh_df, coords = c("X", "Y"))

    if (!is.null(private$crs)) {
      sf::st_crs(mesh_sf) <- private$crs
    } else {
      sf::st_crs(mesh_sf) <- sf::NA_crs_
    }

    mapview_output <- mapview_output + mapview::mapview(
      x = mesh_sf,
      cex = vertex_size * 0.5,
      col.regions = "gray",
      color = "black",
      layer.name = "Mesh",
      ...
    )
  } else{
    warning("The mesh was not shown since vertex_size is zero.")
  }
} 

if (!is.null(X)) {
  if (is.null(X_loc)) {
    stop("X supplied but not X_loc")
  }
  if (length(X) != nrow(X_loc)) {
    stop("The number of observations does not match the number of locations!")
  }
  
  points_xy <- self$coordinates(PtE = X_loc)
  x_loc_sf <- sf::st_as_sf(data.frame(x = points_xy[, 1], y = points_xy[, 2], val = as.vector(X)),
                           coords = c("x", "y"))

  mapview_output <- mapview_output + mapview::mapview(
    x = x_loc_sf,
    zcol = "val",
    cex = data_size,
    col.regions = scale_color_mapview,
    layer.name = "X Points",
    ...
  )
}

return(mapview_output)

    }
    if(interactive && (type == "ggplot")){
      print(plotly::ggplotly(p))
      return(invisible(p))
    }
    return(p)
  },

  #' @description Plots the connections in the graph
  #' @return No return value. Called for its side effects.
  plot_connections = function(){
        g <- make_graph(edges = c(t(self$E)), directed = FALSE)
        plot(g)
  },

  #' @description Checks if the graph is a tree (without considering directions)
  #' @return TRUE if the graph is a tree and FALSE otherwise.
  is_tree = function(){
        g <- make_graph(edges = c(t(self$E)), directed = FALSE)
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
  #' @param type The type of plot to be returned. The options are `ggplot` (the default), that uses `ggplot2`; `plotly` that uses `plot_ly` for 3D plots, which requires the `plotly` package, and `mapview` that uses the `mapview` function, to build interactive plots, which requires the `mapview` package.
  #' @param continuous Should continuity be assumed when the plot uses `newdata`?
  #' @param interpolate_plot Should the values to be plotted be interpolated?
  #' @param vertex_size Size of the vertices.
  #' @param vertex_color Color of vertices.
  #' @param edge_width Width for edges.
  #' @param edge_weight Which column from edge weights to plot? If `NULL` edge weights are not plotted. To plot the edge weights when the metric graph `edge_weights` is a vector instead of a `data.frame`, simply set to 1.
  #' @param edge_color For 3D plot, color of edges.
  #' @param line_width For 3D plot, line width of the function curve.
  #' @param line_color Color of the function curve.
  #' @param scale_color Color scale to be used for data and weights.
  #' @param scale_color_mapview Color scale to be applied for data when `type = "mapview"`.
  #' @param support_width For 3D plot, width of support lines.
  #' @param support_color For 3D plot, color of support lines.
  #' @param mapview_caption Caption for the function if `type = "mapview"`.
  #' @param p Previous plot to which the new plot should be added.
  #' @param plotly  `r lifecycle::badge("deprecated")` Use `type` instead.
  #' @param improve_plot  `r lifecycle::badge("deprecated")` Use `interpolate` instead. There is no need to use it to improve the edges.
  #' @param ... Additional arguments for `ggplot()` or `plot_ly()`
  #' @return Either a `ggplot` (if `plotly = FALSE`) or a `plot_ly` object.
  plot_function = function(data = NULL,
                           newdata = NULL,
                           group = 1,
                           X = NULL,
                           type = c("ggplot", "plotly", "mapview"),
                           continuous = TRUE,
                           interpolate_plot = TRUE,
                           edge_weight = NULL,
                           vertex_size = 5,
                           vertex_color = "black",
                           edge_width = 1,
                           edge_color = 'black',
                           line_width = NULL,
                           line_color = 'rgb(0,0,200)',
                           scale_color = ggplot2::scale_color_viridis_c(option = "D"),
                           scale_color_mapview = viridis::viridis(100, option = "D"),
                           support_width = 0.5,
                           support_color = "gray",
                           mapview_caption = "Function",
                           p = NULL,
                           plotly = deprecated(),
                           improve_plot = deprecated(),
                           ...){
    if (is.null(line_width)) {
      line_width = edge_width
    }


     if (lifecycle::is_present(plotly)) {
      lifecycle::deprecate_warn("1.3.0.9000", "plot(plotly)", "plot(type)",
        details = c("The argument `plotly` was deprecated in favor of the argument `type`.")
      )
    if(plotly){
      type <- "plotly"
    }
  }                 

    type <- type[[1]]

    if(!(type %in% c("ggplot", "plotly", "mapview"))){
      stop("type must be one of the following: 'ggplot', 'plotly' or 'mapview'.")
    }

    mesh <- FALSE

    if(!is.null(edge_weight)){
      edge_weight <- edge_weight[[1]]
      continuous <- FALSE
      newdata <- do.call(rbind, self$edges)
      newdata <- self$edgeweight_to_data(loc = newdata, weight_col = edge_weight,
                                        data_coords = "spatial",
                                        add = FALSE,
                                        return = TRUE,
                                        verbose = 0,
                                        suppress_warnings = TRUE)
      data <- edge_weight
    }

    if(is.null(data) && is.null(X) && is.null(edge_weight)){
      stop("You should provide either 'data', 'X' or 'edge_weight'.")
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

      if(attr(self$mesh, "continuous")){
        n.v <- dim(self$V)[1]
        XV <- X[1:n.v]
      } 
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

    PtE_edges <- lapply(self$edges, function(edge){attr(edge, "PtE")})

    x.loc <- y.loc <- z.loc <- i.loc <- NULL
    kk = 1
    for (i in 1:self$nE) {
      Vs <- self$E[i, 1]
      Ve <- self$E[i, 2]
      if (mesh) {
        ind <- self$mesh$PtE[, 1] == i

        if (sum(ind)==0) {
          if(attr(self$mesh,"continuous")){
            vals <- rbind(c(0, XV[Vs]),
                        c(1, XV[Ve]))
          } else{
            vals <- NULL
          }

        } else {
          if(attr(self$mesh,"continuous")) {
            vals <- rbind(c(0, XV[Vs]),
                          cbind(self$mesh$PtE[ind, 2], X[n.v + which(ind)]),
                          c(1, XV[Ve]))
          } else {
              vals <- cbind(self$mesh$PtE[ind, 2], X[which(ind)])
          }


        }

        if(interpolate_plot){
          PtE_tmp <- PtE_edges[[i]]
          PtE_tmp <- setdiff(PtE_tmp, vals[,1])
          if(length(PtE_tmp)>0){
                PtE_tmp <- cbind(PtE_tmp, NA)
                vals <- rbind(vals,PtE_tmp)
          }

          # if(nrow(vals)>0){
          #   ord_idx <- order(vals[,1])
          #   vals <- vals[ord_idx,]
          #   if(vals[1,1] > 0){
          #     vals <- rbind(c(0,NA), vals)
          #   }
          #   if(vals[nrow(vals),1] < 1){
          #     vals <- rbind(vals, c(1,NA))
          #   }
          #   max_val <- max(vals[,2], na.rm=TRUE)
          #   min_val <- min(vals[,2], na.rm=TRUE)
          #   vals[,2] <- na.const(pmax(pmin(object = zoo::na.approx(object = vals[,2],
          #                                         x = vals[,1],
          #                                             na.rm=FALSE, ties = "mean"),
          #                                      max_val), min_val))
          #   vals <- vals[(vals[,1] >= 0) & (vals[,1]<=1),]
          # }
          if(nrow(vals) > 0) {
            # Sort by first column
            ord_idx <- order(vals[,1])
            vals <- vals[ord_idx,]

            # Add boundary points if needed
            if(vals[1,1] > 0) {
              vals <- rbind(c(0,NA), vals)
            }
            if(vals[nrow(vals),1] < 1) {
              vals <- rbind(vals, c(1,NA))
            }

            # Only proceed with interpolation if there are any non-NA values
            if(any(!is.na(vals[,2]))) {
              max_val <- max(vals[,2], na.rm=TRUE)
              min_val <- min(vals[,2], na.rm=TRUE)

              # Perform interpolation
              interpolated <- try(
                zoo::na.approx(object = vals[,2],
                               x = vals[,1],
                               na.rm=FALSE, 
                               ties = "mean"),
                silent = TRUE
              )

              if(!inherits(interpolated, "try-error")) {
                vals[,2] <- na.const(pmax(pmin(interpolated, max_val), min_val))
              }
            }
            # Keep only values in [0,1] range
            vals <- vals[(vals[,1] >= 0) & (vals[,1] <= 1),]
          }
        }


      } else {
        X <- as.data.frame(X)
        vals <- X[X[, 1]==i, 2:3, drop = FALSE]

    if(continuous){
      if(nrow(vals)>0){

        if(!interpolate_plot){
            if (max(vals[, 1]) < 1) {
              # Check if we can add end value from another edge
              start_Ei <- which(self$E[, 1] == Ve)  # Edges that start at Ve
              end_Ei <- which(self$E[, 2] == Ve)    # Edges that end at Ve

              min.val <- NA
              max.val <- NA

              # Process edges starting at Ve
              if (length(start_Ei) > 0) {
                valid_X_start <- X[X[, 1] %in% start_Ei, , drop = FALSE]
                ind_start <- which(valid_X_start[, 2] == 0)
                if (length(ind_start) > 0) {
                  ind_min <- which.min(valid_X_start[ind_start, 3])
                  min.val <- valid_X_start[ind_start, 3][ind_min]
                }
              }

              # Process edges ending at Ve
              if (length(end_Ei) > 0) {
                valid_X_end <- X[X[, 1] %in% end_Ei, , drop = FALSE]
                ind_end <- which(valid_X_end[, 2] == 1)
                if (length(ind_end) > 0) {
                  ind_max <- which.max(valid_X_end[ind_end, 3])
                  max.val <- valid_X_end[ind_end, 3][ind_max]
                }
              }

              # Add the closest value to the current vertex
              if (!is.na(min.val) && !is.na(max.val)) {
                if (1 - max.val < min.val) {
                  vals <- rbind(vals, c(1, max.val))
                } else {
                  vals <- rbind(vals, c(1, min.val))
                }
              } else if (!is.na(min.val)) {
                vals <- rbind(vals, c(1, min.val))
              } else if (!is.na(max.val)) {
                vals <- rbind(vals, c(1, max.val))
              }
          } 

          if   (min(vals[, 1]) > 0) {
              # Check if we can add start value from another edge
              start_Ei <- which(self$E[, 1] == Vs)  # Edges that start at Vs
              end_Ei <- which(self$E[, 2] == Vs)    # Edges that end at Vs

              min.val <- NA
              max.val <- NA

              # Process edges starting at Vs
              if (length(start_Ei) > 0) {
                valid_X_start <- X[X[, 1] %in% start_Ei, , drop = FALSE]
                ind_start <- which(valid_X_start[, 2] == 0)
                if (length(ind_start) > 0) {
                  ind_min <- which.min(valid_X_start[ind_start, 3])
                  min.val <- valid_X_start[ind_start, 3][ind_min]
                }
              }

              # Process edges ending at Vs
              if (length(end_Ei) > 0) {
                valid_X_end <- X[X[, 1] %in% end_Ei, , drop = FALSE]
                ind_end <- which(valid_X_end[, 2] == 1)
                if (length(ind_end) > 0) {
                  ind_max <- which.max(valid_X_end[ind_end, 3])
                  max.val <- valid_X_end[ind_end, 3][ind_max]
                }
              }

              # Add the closest value to the current vertex
              if (!is.na(min.val) && !is.na(max.val)) {
                if (1 - max.val < min.val) {
                  vals <- rbind(c(0, max.val), vals)
                } else {
                  vals <- rbind(c(0, min.val), vals)
                }
              } else if (!is.na(min.val)) {
                vals <- rbind(c(0, min.val), vals)
              } else if (!is.na(max.val)) {
                vals <- rbind(c(0, max.val), vals)
              } else if (nrow(vals) > 0) {
                # Add a placeholder value based on the closest existing value in vals
                idx_tmp <- which.min(vals[, 1])
                vals <- rbind(c(0, vals[idx_tmp, 2]), vals)
              }
          } 
        } else {        
          PtE_tmp <- PtE_edges[[i]]

          # Check for edges ending at the starting point of the current edge
          if (any(self$E[, 2] == self$E[i, 1])) {
            edge_new <- which(self$E[, 2] == self$E[i, 1])[1]
            idx_new <- which(X[, 1] == edge_new)
            new_val <- X[idx_new, 2:3, drop = FALSE]
            if (nrow(new_val) > 0) {
              sub_fact <- ifelse(any(vals[, 1] == 0), 1 + 1e-6, max(new_val[, 1]))
              new_val[, 1] <- new_val[, 1] - sub_fact
              vals <- rbind(vals, new_val)
            }
          } else if (any(self$E[-i, 1] == self$E[i, 1])) {
            edge_new <- which(self$E[-i, 1] == self$E[i, 1])[1]
            idx_new <- which(X[, 1] == edge_new)
            new_val <- X[idx_new, 2:3, drop = FALSE]
            if (nrow(new_val) > 0) {
              sum_fact <- ifelse(any(vals[, 1] == 1), -1e-6, min(new_val[, 1]))
              new_val[, 1] <- -new_val[, 1] + sum_fact
              vals <- rbind(vals, new_val)
            }
          }

          # Check for edges starting at the ending point of the current edge
          if (any(self$E[, 1] == self$E[i, 2])) {
            edge_new <- which(self$E[, 1] == self$E[i, 2])[1]
            idx_new <- which(X[, 1] == edge_new)
            new_val <- X[idx_new, 2:3, drop = FALSE]
            if (nrow(new_val) > 0) {
              sum_fact <- ifelse(any(vals[, 1] == 1), 1 + 1e-6, 1 - min(new_val[, 1]))
              new_val[, 1] <- new_val[, 1] + sum_fact
              vals <- rbind(vals, new_val)
            }
          } else if (any(self$E[, 2] == self$E[i, 2])) {
            edge_new <- which(self$E[, 2] == self$E[i, 2])[1]
            idx_new <- which(X[, 1] == edge_new)
            new_val <- X[idx_new, 2:3, drop = FALSE]
            if (nrow(new_val) > 0) {
              sub_fact <- ifelse(any(vals[, 1] == 0), 1 + 1e-6, -max(new_val[, 1]) - 1)
              new_val[, 1] <- -new_val[, 1] - sub_fact
              vals <- rbind(vals, new_val)
            }
          }

          # Add remaining values from PtE_tmp not in vals
          PtE_tmp <- setdiff(PtE_tmp, vals[, 1])
          if (length(PtE_tmp) > 0) {
            PtE_tmp <- cbind(PtE_tmp, NA)
            PtE_tmp <- as.data.frame(PtE_tmp)
            colnames(PtE_tmp) <- c(".distance_on_edge", data)
            vals <- rbind(vals, PtE_tmp)
          }

          # Sort values by the first column
          ord_idx <- order(vals[, 1])
          vals <- vals[ord_idx, ]

          # Interpolate missing values within bounds
          # max_val <- max(vals[, 2], na.rm = TRUE)
          # min_val <- min(vals[, 2], na.rm = TRUE)
          # vals[, 2] <- na.const(
          #   pmax(
          #     pmin(
          #       zoo::na.approx(object = vals[, 2], x = vals[, 1], na.rm = FALSE, ties = "mean"),
          #       max_val
          #     ),
          #     min_val
          #   )
          # )

          # Only proceed if there are any non-NA values
          if(any(!is.na(vals[,2]))) {
            max_val <- max(vals[,2], na.rm = TRUE)
            min_val <- min(vals[,2], na.rm = TRUE)

            # Try interpolation with error handling
            interpolated <- try(
              zoo::na.approx(
                object = vals[,2], 
                x = vals[,1], 
                na.rm = FALSE, 
                ties = "mean"
              ),
              silent = TRUE
            )

            if(!inherits(interpolated, "try-error")) {
              vals[,2] <- na.const(
                pmax(
                  pmin(
                    interpolated,
                    max_val
                  ),
                  min_val
                )
              )
            }
          } # If all values are NA, vals[,2] remains unchanged


          # Filter values to lie within [0, 1]
          vals <- vals[(vals[, 1] >= 0) & (vals[, 1] <= 1), ]

          # Add a starting value if vals is non-empty and no value at 0
          if (nrow(vals) > 0 && !any(vals[, 1] == 0)) {
            vals <- rbind(c(0, vals[1, 2, drop = TRUE]), vals)
          }            
        }
      } else{
        if(interpolate_plot){
          vals <- NULL
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
                vals <- rbind(vals, c(1, min.val))
                  # if(length(min.val)>0){
                  #   vals <- rbind(vals, c(1, min.val[[1]]))
                  # } else{
                  #   vals <- rbind(vals, c(1, X[ind, 3,drop=TRUE]))
                  # }
                  vals <- rbind(vals, c(1, X[ind, 3,drop=TRUE]))
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
                  # if(length(max.val)>0){
                  #   vals <- rbind(vals, c(1, max.val[[1]]))
                  # } else{
                  #   vals <- rbind(vals, c(1, X[ind, 3, drop=TRUE]))
                  # }
                  vals <- rbind(vals, c(1, X[ind, 3, drop=TRUE]))
                }
              }

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
                } else{
                  vals <- rbind(c(0, vals[1, 2, drop=TRUE]), vals)
                }
              }
              if(ncol(vals)<2){
                vals <- NULL
              }
        }
      }
      } else if(interpolate_plot){

          PtE_tmp <- PtE_edges[[i]]
          PtE_tmp <- setdiff(PtE_tmp, vals[,1])
          if(length(PtE_tmp)>0){
                PtE_tmp <- cbind(PtE_tmp, NA)
                PtE_tmp <- as.data.frame(PtE_tmp)
                colnames(PtE_tmp) <- c(".distance_on_edge", data)
                vals <- rbind(vals,PtE_tmp)
          }
            # if(nrow(vals)>0){
            #   ord_idx <- order(vals[,1])
            #   vals <- vals[ord_idx,]
            #   if(vals[1,1] > 0){
            #     vals <- rbind(c(0,NA), vals)
            #   }
            #   if(vals[nrow(vals),1] < 1){
            #     vals <- rbind(vals, c(1,NA))
            #   }
            #   max_val <- max(vals[,2], na.rm=TRUE)
            #   min_val <- min(vals[,2], na.rm=TRUE)
            #   vals[,2] <- na.const(pmax(pmin(object = zoo::na.approx(object = vals[,2],
            #                                         x = vals[,1],
            #                                             na.rm=FALSE, ties = "mean"),
            #                                      max_val), min_val))
            #   vals <- vals[(vals[,1] >= 0) & (vals[,1]<=1),]
            # }
            
            if(nrow(vals) > 0) {
              # Sort by first column
              ord_idx <- order(vals[,1])
              vals <- vals[ord_idx,]

              # Add boundary points if needed
              if(vals[1,1] > 0) {
                vals <- rbind(c(0,NA), vals)
              }
              if(vals[nrow(vals),1] < 1) {
                vals <- rbind(vals, c(1,NA))
              }

              # Only proceed with interpolation if there are any non-NA values
              if(any(!is.na(vals[,2]))) {
                max_val <- max(vals[,2], na.rm=TRUE)
                min_val <- min(vals[,2], na.rm=TRUE)

                # Perform interpolation
                interpolated <- try(
                  zoo::na.approx(object = vals[,2],
                                 x = vals[,1],
                                 na.rm=FALSE, 
                                 ties = "mean"),
                  silent = TRUE
                )

                if(!inherits(interpolated, "try-error")) {
                  vals[,2] <- na.const(pmax(pmin(interpolated, max_val), min_val))
                }
              }
              # Keep only values in [0,1] range
              vals <- vals[(vals[,1] >= 0) & (vals[,1] <= 1),]
            }

      }
      }

        data.to.plot <- vals
        if(!is.null(data.to.plot)){
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



    }

    data <- data.frame(x = x.loc, y = y.loc, i = i.loc, z = z.loc)

    if(type == "plotly"){
      requireNamespace("plotly")
      if(is.null(p)){
        p <- self$plot(type = "plotly",
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

      p <-  plotly::layout(p, scene = list(xaxis = list(autorange = "reversed")))

      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degree" && !private$transform){
          p <- plotly::layout(p, scene = list(xaxis = list(title = "Longitude"), yaxis = list(title = "Latitude")))
        } else{
          p <- plotly::layout(p, scene = list(xaxis = list(title = paste0("x (in ",private$vertex_unit, ")")), yaxis = list(title = paste0("y (in ",private$vertex_unit, ")"))))
        }
      }

    } else if (type == "ggplot"){
      if(is.null(p)) {
          p <- ggplot(data = data) +
          geom_path( mapping = aes(x = x, y = y,
                                     group = i,
                                     colour = z), linewidth = line_width) + labs(colour = "") + scale_color # + scale_color_viridis() +
      } else {
        p <- p + geom_path(data = data, mapping =
                           aes(x = x, y = y,
                               group = i, colour = z),
                           linewidth = line_width) + labs(colour = "") + scale_color # + scale_color_viridis()
      }
          p <- self$plot(edge_width = 0, vertex_size = vertex_size,
                     vertex_color = vertex_color, p = p)
      if(!is.null(private$vertex_unit)){
        if(private$vertex_unit == "degree" && !private$transform){
          p <- p + labs(x = "Longitude",  y = "Latitude")
        } else{
          p <- p + labs(x = paste0("x (in ",private$vertex_unit, ")"),  y = paste0("y (in ",private$vertex_unit, ")"))
        }
      }
    } else if(type == "mapview"){
    requireNamespace("mapview")
    data <- na.omit(data[, c("x", "y", "z", "i")])
    data_split <- split(data, data$i)
    
    linestrings <- do.call(c, lapply(data_split, function(group_data) {
      if (nrow(group_data) < 2) {
          return(NULL)  
      }
      lapply(1:(nrow(group_data) - 1), function(j) {
        sf::st_linestring(as.matrix(group_data[j:(j + 1), c("x", "y")]))
      })
    }))
    
    data_segment <- do.call(rbind, lapply(data_split, function(group_data) {
      data.frame(
        z = group_data$z[-nrow(group_data)], 
        i = group_data$i[-nrow(group_data)]  
      )
    }))
    
    data_sf <- sf::st_as_sf(
      data_segment,
      geometry = sf::st_sfc(linestrings),
      crs = if (!is.null(private$crs)) private$crs else sf::NA_crs_
    )

    data_sf <- sf::st_as_sf(
      data_segment, 
      geometry = sf::st_sfc(linestrings),
      crs = if (!is.null(private$crs)) private$crs else sf::NA_crs_
    )

    edges_sf <- self$get_edges(format = "sf")

    if(is.null(p)){
      mapview_output <- mapview::mapview(
        x = edges_sf,
        lwd = edge_width,  
        color = edge_color, 
        layer.name = "Edges",
        ...
      )
    } else{
      mapview_output <- mapview::mapview(
        x = edges_sf,
        lwd = 1.75 * edge_width,  
        color = edge_color, 
        layer.name = "Edges",
        ...
      )
    }
    edge_width <- 2.5 * edge_width
    mapview_output <- mapview_output + mapview::mapview(
        x = data_sf,
        zcol = "z",
        color = scale_color_mapview, 
        lwd = edge_width, 
        cex = vertex_size, 
        col.regions = scale_color_mapview, 
        layer.name = mapview_caption
      )


    return(mapview_output)
    }
    return(p)
  },

  #' @description Plots a movie of a continuous function evolving on the graph.
  #' @param X A m x T matrix where the ith column represents the function at the
  #' ith time, evaluated at the mesh locations.
  #' @param type Type of plot. Either `"plotly"` or `"ggplot"`.
  #' @param vertex_size Size of the vertices.
  #' @param vertex_color Color of vertices.
  #' @param edge_width Width for edges.
  #' @param edge_color For 3D plot, color of edges.
  #' @param line_width For 3D plot, line width of the function curve.
  #' @param line_color Color of the function curve.
  #' @param ... Additional arguments for ggplot or plot_ly.
  #' @return Either a `ggplot` (if `plotly=FALSE`) or a `plot_ly` object.
  plot_movie = function(X,
                        type = "plotly",
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

    if (type == "plotly") {
      requireNamespace("plotly")
      plotly <- TRUE
    } else if (type == "ggplot") {
      requireNamespace("ggplot2")
      plotly <- FALSE
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
        if(private$vertex_unit == "degree" && !private$transform){
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
     if(verbose == 2) {
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
          if(verbose == 2) {
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
          if(verbose == 2) {
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

    if(verbose == 2) {
      message("Part 1/2")
      bar_line_vertex <- msg_progress_bar(length(self$edges))
    }


    lines <- matrix(nrow = 2*length(self$edges), ncol = 3)
    for(i in 1:length(self$edges)){
      if(verbose == 2) {
        bar_line_vertex$increment()
      }
      points <- self$edges[[i]]
      n <- dim(points)[1]
      lines[2*i-1,] <- c(i, points[1,])
      lines[2*i, ] <- c(i, points[n,])
    }

    #save all vertices that are more than tolerance distance apart
    if(verbose == 2) {
      message("Computing auxiliary distances")
    }

    if(!project_data && longlat){
    if(!project || !longlat){
        fact <- process_factor_unit(vertex_unit, length_unit)
          dists <- compute_aux_distances(lines = lines[,2:3,drop=FALSE], crs = crs, longlat = longlat, proj4string = proj4string, fact = fact, which_longlat = which_longlat, length_unit = private$length_unit, transform = private$transform)
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


    if(verbose == 2) {
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

    if(verbose == 2) {
      message("Part 2/2")
      bar_line_vertex <- msg_progress_bar(max(lines[, 1]))
    }

    lvl <- matrix(0, nrow = max(lines[,1]), 4)
    k=1
    lines_keep_id <- NULL

    for (i in 1:max(lines[, 1])) {

      if(verbose == 2) {
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
      ll <- compute_line_lengths(self$edges[[i]], longlat = longlat, unit = length_unit, crs = crs, proj4string, which_longlat, vertex_unit, project_data, private$transform)
      if(ll > tolerance) {
        lvl[k,] <- c(i, ind1, ind2, ll)
        k=k+1
        lines_keep_id <- c(lines_keep_id, i)
      }
    }

    lvl <- lvl[1:(k-1),,drop = FALSE]
    self$edges <- self$edges[lines_keep_id]

    if(is.vector(private$edge_weights)){
        private$edge_weights <- private$edge_weights[lines_keep_id]
    } else{
        private$edge_weights <- private$edge_weights[lines_keep_id,,drop=FALSE]
    }

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
      if(v1 != v2){
        # if they are different, only consider edges with different vertices as endpoints
        valid_ind <- (rowSums(self$mesh$E == v1) > 0) & (rowSums(self$mesh$E == v2) > 0)
        e <- which(valid_ind & (rowSums(self$mesh$E == v1) + rowSums(self$mesh$E == v2) == 2))
      } else{
        e <- which(rowSums((self$mesh$E == v1)) == 2)
      }
      # Handle the case of multiple edges
      # In the case the edge lengths are different, we can identify
      # In the case they are equal, we cannot identify, but it does not matter, as there is no difference then.
      e_bkp <- e
      if(length(e)>1){
       ind <- which(self$edge_lengths[ei] == self$mesh$h_e[e])
       e <- e[ind]
       e <- e[1]
      }
      # a loop might have been split into smaller parts, in this case, we simply get the first
      if(self$edge_lengths[ei] == sum(self$mesh$h_e[e_bkp])){
        e <- e_bkp[which(self$mesh$E[e_bkp,1] == v1)]
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
                     edge_weight = NULL,
                     edge_width_weight = NULL,
                     scale_color_main = ggplot2::scale_color_viridis_c(option = "D"),
                     scale_color_weights = ggplot2::scale_color_viridis_c(option = "A"),
                     scale_color_degree = ggplot2::scale_color_viridis_d(option = "D"),
                     add_new_scale_weights = TRUE,
                     arrow_size,
                     scale_color_main_discrete,
                     scale_color_weights_discrete,
                     ...){
    xyl <- c()

    nc <- do.call(rbind,lapply(self$edges, function(x) dim(x)[1]))
    xyl <- cbind(do.call(rbind,self$edges), rep(1:length(nc), times = nc))
    df_plot <- data.frame(x = xyl[, 1], y = xyl[, 2], grp = xyl[, 3])

    if(!is.null(edge_weight)){
      edge_weight <- edge_weight[[1]]
      e_weights <- private$get_edge_weights_internal(data.frame = TRUE)
      e_weights <- e_weights[,edge_weight, drop = FALSE]
      colnames(e_weights) <- "weights"
      e_weights["grp"] <- 1:self$nE
      df_plot <- merge(df_plot, e_weights)
    } else{
      df_plot[["weights"]] <- rep(edge_color, nrow(df_plot))
    }
    if(!is.null(edge_width_weight)){
      edge_width_weight <- edge_width_weight[[1]]
      e_weights <- private$get_edge_weights_internal(data.frame = TRUE)
      e_weights <- e_weights[,edge_width_weight, drop = FALSE]
      e_weights[,1] <- e_weights[,1] * line_width / max(e_weights[,1])
      e_weights[,1] <- e_weights[,1]
      colnames(e_weights) <- "widths"
      e_weights["grp"] <- 1:self$nE
      df_plot <- merge(df_plot, e_weights)
    } else{
      df_plot[["widths"]] <- rep(line_width, nrow(df_plot))
    }

    if(is.null(p)){
      if(!is.null(edge_weight)){
        if(is.factor(df_plot[["weights"]]) || is.character(df_plot[["weights"]])){
          scale_weights <- scale_color_weights_discrete
        } else{
          scale_weights <- scale_color_weights
        }
        p <- ggplot() + geom_path(data = df_plot,
                                  mapping = aes(x = x, y = y, group = grp,
                                  colour = weights, linewidth = widths),
                                  ...) + ggplot2::scale_linewidth_identity() + scale_weights + labs(colour = edge_weight)
          if(add_new_scale_weights){
            p <- p + new_scale_color()
          }
      } else{
        p <- ggplot() + geom_path(data = df_plot,
                                  mapping = aes(x = x, y = y, group = grp, linewidth = widths),
                                  color = edge_color,
                                  # linewidth = line_width,
                                  ...) + ggplot2::scale_linewidth_identity()
      }
    } else {
      if(!is.null(edge_weight)){
        if(is.factor(df_plot[["weights"]]) || is.character(df_plot[["weights"]])){
          scale_weights <- scale_color_weights_discrete
        } else{
          scale_weights <- scale_color_weights
        }        
        p <- p + geom_path(data = df_plot,
                           mapping = aes(x = x, y = y, group = grp, colour = weights, linewidth =widths),
                           ...) + ggplot2::scale_linewidth_identity() + scale_weights + labs(colour = edge_weight)
          if(add_new_scale_weights){
            p <- p + new_scale_color()
          }
      } else{
        p <- p + geom_path(data = df_plot,
                           mapping = aes(x = x, y = y, group = grp,  linewidth = widths), color = edge_color, ...) + ggplot2::scale_linewidth_identity()
      }
    }

    if(direction) {
      mid.l <- self$coordinates(PtE = cbind(1:self$nE, rep(0.49,self$nE)))
      mid.u <- self$coordinates(PtE = cbind(1:self$nE, rep(0.5,self$nE)))
      p <- p + geom_path(data = data.frame(x = c(mid.l[,1], mid.u[,1]),
                                           y = c(mid.l[, 2], mid.u[,2]),
                                           edge = c(1:self$nE,1:self$nE)),
                          mapping = aes(x = x, y = y, group = edge),
                         arrow = ggplot2::arrow(length=arrow_size),
                          linewidth = marker_size/2, ...)
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
                            size= marker_size, ...) + scale_color_degree # +
    # scale_color_viridis(discrete = TRUE, guide_legend(title = ""))
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

      if(is.factor(as.vector(y_plot[!is.na(as.vector(y_plot))])) || is.character(as.vector(y_plot[!is.na(as.vector(y_plot))]))){
        scale_main <- scale_color_main_discrete
      } else{
        scale_main <- scale_color_main
      }

      p <- p + geom_point(data = data.frame(x = x[!is.na(as.vector(y_plot))],
                                            y = y[!is.na(as.vector(y_plot))],
                                            val = as.vector(y_plot[!is.na(as.vector(y_plot))])),
                          mapping = aes(x, y, color = val),
                          size = data_size, ...) + labs(color = data) + scale_main #+
        # scale_colour_gradientn(colours = viridis(100), guide_legend(title = ""))

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
      if(is.factor(as.vector(X)) || is.character(as.vector(X))){
        scale_main <- scale_color_main_discrete
      } else{
        scale_main <- scale_color_main
      }
      p <- p + geom_point(data = data.frame(x = x, y = y,
                                            val = as.vector(X)),
                          mapping = aes(x, y, color = val),
                          size = data_size) + labs(colour = "") + scale_main #+
        # scale_color_viridis()
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
                    edge_width_weight = NULL,
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

    if(!is.null(edge_width_weight)){
        edge_width_weight <- edge_width_weight[[1]]
        e_weights <- private$get_edge_weights_internal(data.frame = TRUE)
        e_weights <- e_weights[,edge_width_weight, drop = FALSE]
        e_weights[,1] <- e_weights[,1] * line_width / max(e_weights[,1])
        e_weights[,1] <- e_weights[,1]
        colnames(e_weights) <- "widths"
        e_weights["i"] <- 1:self$nE
        data.plot <- merge(data.plot, e_weights)
    } else{
      data.plot[["widths"]] <- rep(line_width, nrow(data.plot))
    }

    if(is.null(p)) {
      p <- plotly::plot_ly(data=data.plot, x = ~y, y = ~x, z = ~z,...)
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = ~widths,
                                       color = edge_color),
                           split = ~i, showlegend = FALSE)
    } else {
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = ~widths,
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

      if(verbose == 2) {
        message("Computing auxiliary distances")
      }

      within_dist <- t(as.matrix(sf::st_is_within_distance(points_sf, lines_sf, dist = tolerance)))

      if(verbose == 2) {
        message("Done!")
      }

      if(verbose == 2) {
        message("Snapping vertices")
        bar_multiple_snaps <- msg_progress_bar(length(self$edges))
      }

      for(i in 1:length(self$edges)){
        if(verbose == 2) {
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
        if(verbose == 2) {
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
      dists <- compute_aux_distances(lines = self$V, crs = private$crs, longlat = private$longlat, proj4string = private$proj4string, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit, transform = private$transform)
      v.merge <- NULL
      k <- 0

      if(self$nV>1){
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
    if(verbose == 2) {
      message("Removing small circles")
    }

    if(is.null(private$degrees)){
      degrees <- private$compute_degrees(verbose=verbose)$degrees
    } else{
      degrees <- private$degrees$degrees
    }

    if(threshold > 0) {
      loop.ind <- which(self$E[,1] == self$E[,2])
      if(length(loop.ind)>0) {
        loop.size <- self$edge_lengths[loop.ind]
        ind <- loop.ind[loop.size < threshold]
        if(length(ind)>0) {
          v.loop <- self$E[ind,1]
          v.degrees <- degrees[v.loop]

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
          if(!is.null(self$edge_lengths)){
            self$edge_lengths <- self$edge_lengths[-ind]
          }
          self$E <- self$E[-ind,]
          self$nE <- self$nE - length(ind)
          if(is.vector(private$edge_weights)){
            private$edge_weights <- private$edge_weights[-ind]
          } else{
            private$edge_weights <- private$edge_weights[-ind,,drop=FALSE]
          }

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


  remove.first.deg2 = function(res, check_circles) {
    if(!check_circles){
      ind <- which(res$degrees==2 & !res$problematic)
    } else{
      ind <- which(res$degrees==2 & !res$problematic & !res$problematic_circles)
    }
    res.out <- res
    if(length(ind)>0) {
      ind <- ind[1]
      e1 <- which(self$E[,2]==ind)
      e2 <- which(self$E[,1]==ind)

      if(!check_circles || (self$E[e1,1] != self$E[e2,2])){
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
             PtEedge1 <- attr(coords, "PtE")
             PtEedge2 <- attr(tmp, "PtE")

            if(which_line_starts == 1){
               coords <- rbind(coords, tmp[-1,,drop=FALSE])
               E_new <- matrix(c(v1,v2),1,2)
               PtEedge1 <- PtEedge1 * self$edge_lengths[e_rem[1]]
               PtEedge2 <- self$edge_lengths[e_rem[1]] + PtEedge2[-1] * self$edge_lengths[e_rem[2]]
               attr(coords, "PtE") <- c(PtEedge1,PtEedge2)/(self$edge_lengths[e_rem[1]] + self$edge_lengths[e_rem[2]])
            } else{
               coords <- rbind(tmp,coords[-1,,drop=FALSE])
               E_new <- matrix(c(v2,v1),1,2)
               PtEedge2 <- PtEedge2 * self$edge_lengths[e_rem[2]]               
               PtEedge1 <- self$edge_lengths[e_rem[2]] + PtEedge1[-1] * self$edge_lengths[e_rem[1]]
               attr(coords, "PtE") <- c(PtEedge2,PtEedge1)/(self$edge_lengths[e_rem[1]] + self$edge_lengths[e_rem[2]])
            }

          # Updating the merged graph

          res.out$degrees <- res$degrees[-ind]
          res.out$problematic <- res$problematic[-ind]
          res.out$problematic_circles <- res.out$problematic_circles[-ind]

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
            # if(private$edge_weights[e_rem[2]] != private$edge_weights[e_rem[1]]){
            if(compare_with_na(private$edge_weights[e_rem[2]], private$edge_weights[e_rem[1]])){
                private$prune_warning <- TRUE
            }
            private$edge_weights <- private$edge_weights[-e_rem[2]]
          } else{
            # if(any(private$edge_weights[e_rem[2],,drop=FALSE] != private$edge_weights[e_rem[1],,drop=FALSE])){
            if(compare_with_na(private$edge_weights[e_rem[2],,drop=FALSE], 
                              private$edge_weights[e_rem[1],,drop=FALSE],
                              is_matrix = TRUE)){
                private$prune_warning <- TRUE
            }
            private$edge_weights <- private$edge_weights[-e_rem[2],,drop=FALSE]
          }
          }

      } else{
          res.out$problematic_circles[ind] <- TRUE
          res.out$circles_avoided <- TRUE 
      }
    }
    class(self$edges) <- "metric_graph_edges"
    return(res.out)
  },

  # Compute lengths

  compute_lengths = function(longlat, unit, crs, proj4string, which_longlat, vertex_unit, project_data, transform){
          ll <- sapply(self$edges,
          function(edge){compute_line_lengths(edge, longlat = longlat, unit = unit, crs = crs, proj4string, which_longlat, vertex_unit, project_data, transform)})
          return(ll)
  },

  approx_coordinates = function(edge){
    # Calculate the Euclidean distances between consecutive rows (vertices)
    dists <- sqrt(rowSums((edge[-1, ,drop=FALSE] - edge[-nrow(edge), ,drop=FALSE])^2))

    # Get cumulative sum to obtain positions along the edge
    cum_dists <- c(0, cumsum(dists))  # Start from 0 for the first vertex

    # Normalize the cumulative distances to get relative positions between 0 and 1
    rel_positions <- cum_dists / cum_dists[length(cum_dists)]

    return(rel_positions)
  },

  exact_PtE_coordinates = function(edge) {
  # Get distances using compute_aux_distances between consecutive rows
  lines <- edge[-1, , drop = FALSE]      # All rows except the first
  points <- edge[-nrow(edge), , drop = FALSE] # All rows except the last
  
  # Compute distances between each consecutive point
  dists <- compute_aux_distances(lines = lines, crs = private$crs, longlat = private$longlat, 
                                 proj4string = private$proj4string, points = points, 
                                 fact = private$fact, which_longlat = private$which_longlat, 
                                 length_unit = private$length_unit, transform = private$transform)
  
  # Calculate the relative positions using cumulative sums
  total_length <- sum(dists)
  relative_positions <- c(0, cumsum(dists) / total_length)
  
  return(relative_positions)
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


  #'  Gets the edge weights
  #'  data.frame If the edge weights are given as vectors, should the result be returned as a data.frame?
  #'  A vector or `data.frame` containing the edge weights.

  get_edge_weights_internal = function(data.frame = FALSE){
    tmp <- private$edge_weights
    row.names(tmp) <- NULL
    if(!is.data.frame(tmp) && data.frame){
      tmp <- data.frame(.weights = tmp)
    }
    return(tmp)
  },

format_data = function(data_res, format) {

    if(!(format%in%c("tibble", "sf", "sp"))){
      stop("The possible formats are 'tibble', 'sf', 'sp'.")
    }

    if (format == "sf") {
        # Ensure the coordinates exist in the data
        if (!all(c(".coord_x", ".coord_y") %in% colnames(data_res))) {
            stop("Cannot convert to sf: .coord_x and .coord_y columns are missing.")
        }
        
        # Create sf object from coordinates
        data_res_sf <- sf::st_as_sf(data_res, coords = c(".coord_x", ".coord_y"), crs = sf::NA_crs_)
        class(data_res_sf) <- c("metric_graph_data", class(data_res_sf))
        return(data_res_sf)

    } else if (format == "sp") {
        # Ensure the coordinates exist in the data
        if (!all(c(".coord_x", ".coord_y") %in% colnames(data_res))) {
            stop("Cannot convert to sp: .coord_x and .coord_y columns are missing.")
        }

        # Create SpatialPointsDataFrame from coordinates
        coords <- cbind(data_res$.coord_x, data_res$.coord_y)
        sp_data <- sp::SpatialPointsDataFrame(
            coords = coords,
            data = data_res,
            proj4string = sp::CRS(NA_character_)
        )
        return(sp_data)

    } else {
        data_res <- tidyr::as_tibble(data_res)
        class(data_res) <- c("metric_graph_data", class(data_res))
        return(data_res)
    }
},


  format_weights = function(edge_weights_res, format) {

    if(!(format%in%c("tibble", "sf", "sp"))){
      stop("The possible formats are 'tibble', 'sf', 'sp'.")
    }

  if (format == "sf") {
    edges_geometries <- lapply(self$edges, sf::st_linestring)

    if (is.vector(private$edge_weights)) {
      ew_tmp <- data.frame(.weights = private$edge_weights)
    } else {
      ew_tmp <- as.data.frame(private$edge_weights)
    }
    ew_tmp[[".edge_lengths"]] <- self$edge_lengths

    edges_sf <- sf::st_sf(ew_tmp, geometry = sf::st_sfc(edges_geometries), crs = if (!is.null(private$crs)) private$crs else sf::NA_crs_)
    class(edges_sf) <- c("metric_graph_weights", class(edges_sf))
    return(edges_sf)

  } else if (format == "sp") {
    edges_list <- lapply(1:length(self$edges), function(i) {
      sp::Line(coords = matrix(self$edges[[i]], nrow = dim(self$edges[[i]])[1], ncol = dim(self$edges[[i]])[2]))
    })
    sp_edges <- sp::SpatialLines(
      lapply(1:length(edges_list), function(i) sp::Lines(list(edges_list[[i]]), ID = as.character(i))),
      proj4string = if (!is.null(private$crs)) private$proj4string else sp::CRS(NA_character_)
    )
    
    if (is.vector(private$edge_weights)) {
      ew_tmp <- data.frame(.weights = private$edge_weights)
    } else {
      ew_tmp <- as.data.frame(private$edge_weights)
    }
    ew_tmp[[".ID"]] <- as.character(1:length(sp_edges))
    ew_tmp[[".edge_lengths"]] <- self$edge_lengths
    edges_sldf <- sp::SpatialLinesDataFrame(sp_edges, data = ew_tmp, match.ID = ".ID")
    return(edges_sldf)

  } else {
    edge_weights_res <- tidyr::as_tibble(edge_weights_res)
    class(edge_weights_res) <- c("metric_graph_weights", class(edge_weights_res))
    return(edge_weights_res)
  }
},

  #' data List containing data on the metric graph.

  data = NULL,

  # Initial edges added

  initial_edges_added = NULL,

  # Initial graph

  initial_graph = NULL,

  # Version of MetricGraph package used to build the graph
  
  version = as.character(utils::packageVersion("MetricGraph")),

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

  # Kichhoff weights

  kirchhoff_weights = NULL,

  # manual edge lengths

  manual_edge_lengths = FALSE,

  # perform merges?

  perform_merges = NULL,

  #grouping variables when adding data

  group_variables = NULL,

  project_data = NULL,

  degrees = NULL,

  # Warning if edges with different weights have been merged

  prune_warning = FALSE,

  clear_initial_info = function(){
    private$addinfo = FALSE
    private$initial_edges_added = NULL
  },

  # Observe that by construction of get_PtE, all PtEs are sorted.
  split_edge = function(Ei, t_values, indices = NULL) {

    # Calculate distances and determine if new vertices need to be added
    edge <- self$edges[[Ei]]
    edge_length <- self$edge_lengths[Ei]
    PtE_edge <- attr(edge, "PtE")

    # Data frame creation and filling NA values
    val_results <- interpolate2(edge, pos = t_values, normalized = TRUE, get_idx = TRUE)
    idx_positions <- val_results[["idx"]]
    val_lines <- val_results[["coords"]]

    # Initialize new_vertices
    new_vertices <- numeric(length(t_values))

    # Loop through each value in t_values, updating self$V conditionally
    for (i in seq_along(t_values)) { 
        val_line <- matrix(val_lines[i, , drop = FALSE], nrow = 1)
        newV <- self$nV + 1
        self$V <- rbind(self$V, val_line)
        self$nV <- self$nV + 1
        new_vertices[i] <- newV
    }

    # Edge updates
    edge_updates <- self$nE + seq_len(length(t_values))
    
    if (!is.null(private$data) && !is.null(indices)) {
      private$temp_PtE[indices, 1] <- edge_updates
      private$temp_PtE[indices, 2] <- (private$temp_PtE[indices, 2] - t_values) / (1 - t_values)
    }
    edge_updates <- c(Ei, edge_updates)

    coords_list1 <- edge[1:idx_positions[1], , drop = FALSE]
    coords_list1 <- rbind(coords_list1, val_lines[1, , drop = FALSE])
    tmp_vec <- c(PtE_edge[1:idx_positions[1]], t_values[1])

    pos_edge_diff <- tmp_vec - tmp_vec[1]
    norm_factor <- tmp_vec[length(tmp_vec)] - tmp_vec[1]
    tmp_PtE <- pos_edge_diff / norm_factor

    # Add tmp_PtE as an attribute
    attr(coords_list1, "PtE") <- tmp_PtE

    n_t_values <- length(t_values)

    coords_list2 <- vector("list", n_t_values)

    for (i in seq_along(t_values)) {
      val_line_start <- matrix(val_lines[i, , drop = FALSE], nrow = 1)

      if (i < n_t_values) {
        val_line_end <- matrix(val_lines[i + 1, , drop = FALSE], nrow = 1)
        if (idx_positions[i] != idx_positions[i + 1]) {
          coords_list2[[i]] <- rbind(
            val_line_start,
            edge[(idx_positions[i] + 1):idx_positions[i + 1], , drop = FALSE],
            val_line_end
          )
          tmp_vec <- c(t_values[i], PtE_edge[(idx_positions[i] + 1):idx_positions[i + 1]], t_values[i + 1])
        } else {
          coords_list2[[i]] <- rbind(
            val_line_start,
            val_line_end
          )
          tmp_vec <- c(t_values[i], t_values[i + 1])
        }
      } else {
        coords_list2[[i]] <- rbind(
          val_line_start,
          edge[(idx_positions[i] + 1):nrow(edge), , drop = FALSE]
        )
        tmp_vec <- c(t_values[i], PtE_edge[(idx_positions[i] + 1):nrow(edge)])
      }

      pos_edge_diff <- tmp_vec - tmp_vec[1]
      norm_factor <- tmp_vec[length(tmp_vec)] - tmp_vec[1]
      tmp_PtE <- pos_edge_diff / norm_factor

      # Add tmp_PtE as an attribute
      attr(coords_list2[[i]], "PtE") <- tmp_PtE
    }

    # Construct aux_matrix for self$E
    if (length(t_values) == 1) {
      aux_matrix <- matrix(
        c(self$E[Ei, 1], new_vertices[1], new_vertices[1], self$E[Ei, 2]),
        nrow = 2, byrow = TRUE
      )
    } else {
      aux_matrix <- rbind(
        c(self$E[Ei, 1], new_vertices[1]),
        cbind(new_vertices[-length(new_vertices)], new_vertices[-1]),
        c(new_vertices[length(new_vertices)], self$E[Ei, 2])
      )
    }

    # Update self$E, self$edges, and related attributes
    self$E[Ei, ] <- aux_matrix[1, ]
    self$E <- rbind(self$E, aux_matrix[-1, ])
    self$edges[[Ei]] <- coords_list1
    self$edges <- c(self$edges, coords_list2)
    self$nE <- self$nE + length(t_values)

    # Update edge lengths and weights
    segment_lengths <- c(t_values[1], diff(t_values), 1 - t_values[length(t_values)]) * self$edge_lengths[Ei]
    self$edge_lengths <- c(self$edge_lengths, segment_lengths[-1])
    self$edge_lengths[Ei] <- segment_lengths[1]

    if (is.vector(private$edge_weights)) {
      private$edge_weights <- c(private$edge_weights, rep(private$edge_weights[Ei], length(t_values)))
    } else {
      private$edge_weights <- rbind(private$edge_weights, do.call(rbind, replicate(length(t_values), private$edge_weights[Ei, , drop = FALSE], simplify = FALSE)))
    }

    return(new_vertices)
  },

  compute_laplacian_PtE = function(PtE, normalized = TRUE, verbose = verbose) {
      if(verbose == 2){
        message("Processing the graph locations...")
      }
    t <- system.time({
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

      graph.temp$add_observations(data = df_temp, normalized = normalized, verbose = 0,
                  suppress_warnings = TRUE)
    })
      if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
      if(verbose == 2){
        message("Turning observations of the auxiliary graph to vertices...")
      }
    t <- system.time(
      graph.temp$observation_to_vertex(mesh_warning = FALSE)
    )
    if(verbose > 0){
      message("Computing the Laplacian...")
    }
    t <- system.time({
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
    })
    if(verbose == 2){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
       }

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

  if(verbose == 2) {
    bar_line_line <- msg_progress_bar(length(self$edges)-1)
  }
  for(i in 1:(length(self$edges)-1)) {
    if(verbose == 2) {
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
              if(min(compute_aux_distances(lines = self$V, crs=private$crs, longlat=private$longlat, proj4string = private$proj4string, points = p, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit, transform = private$transform))>tol) {
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
            if(min(compute_aux_distances(lines = self$V, crs=private$crs, longlat=private$longlat, proj4string = private$proj4string, points = p, fact = fact, which_longlat = private$which_longlat, length_unit = private$length_unit, transform = private$transform))>tol) {
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
  if(verbose == 2) {
    bar_eu <- msg_progress_bar(length(e.u))
  }
  for (i in 1:length(e.u)) {
    if(verbose == 2) {
      bar_eu$increment()
    }
    dists <- sort(PtE[which(PtE[,1]==e.u[i]),2])
    if(length(dists) > 0){
    if(dists[1] > tolerance){
      private$split_edge(e.u[i], dists[1])
    }
    }
    if(length(dists)>1) {
      dists_up <- dists
      for(j in 2:length(dists)){
        dists_up[j] <- (dists[j] - dists[j-1])/(1 - dists[j-1])
        if(dists_up[j] > tolerance){
          private$split_edge(self$nE, dists_up[j])
        } 
      }
    }
  }
  return(self)
},

 # Function to compute the degrees of the vertices,

  compute_degrees = function(verbose = 0, add = FALSE) {
    if (verbose == 2) {
      message("Computing degrees...")
    }

    # Vectorized computation of in-degrees and out-degrees
    degrees_out <- tabulate(self$E[, 1], nbins = self$nV)
    degrees_in <- tabulate(self$E[, 2], nbins = self$nV)

    # Compute total degrees
    degrees <- degrees_in + degrees_out

    if (add) {
      private$degrees <- list(degrees = degrees, indegrees = degrees_in, outdegrees = degrees_out)
    }

    return(list(degrees = degrees, indegrees = degrees_in, outdegrees = degrees_out))
  },

  # Reference edges for the vertices

  ref_edges = NULL,

  #  Creates/updates the vertices element of the metric graph list

   create_update_vertices = function(verbose = 0) {
    degrees <- private$compute_degrees(verbose = verbose, add = TRUE)

    # Check if the number of rows in self$V matches the length of degrees$degrees
    n_vertices <- nrow(self$V)
    if (length(degrees$degrees) != n_vertices ||
        length(degrees$indegrees) != n_vertices ||
        length(degrees$outdegrees) != n_vertices) {
      stop("Mismatch in the number of vertices and degree information. Please check your graph data.")
    }

    if (verbose == 2) {
      if (is.null(self$vertices)) {
        message("Creating vertices object")
      } else {
        message("Updating vertices object")
      }
      bar_update_attr_edges <- msg_progress_bar(n_vertices)
    }

    # Set column names for self$V
    colnames(self$V) <- c("X", "Y")

    # Compute problematic vertices
    problematic <- (degrees$degrees > 1) & ((degrees$indegrees == 0) | (degrees$outdegrees == 0))
    ids <- seq_len(n_vertices)

    # Prepare longlat (which should always be present)
    longlat <- rep(private$longlat, n_vertices)

    if (is.null(private$crs$input)) {
      # Case when crs is NULL
      vertices_df <- data.frame(
        X = self$V[, 1],
        Y = self$V[, 2],
        degree = degrees$degrees,
        indegree = degrees$indegrees,
        outdegree = degrees$outdegrees,
        problematic = problematic,
        longlat = longlat,
        id = ids
      )

      # Convert each row to a "metric_graph_vertex" object without crs
      self$vertices <- lapply(1:n_vertices, function(i) {
        vert <- as.numeric(vertices_df[i, c("X", "Y")])  # Convert to numeric vector
        attr(vert, "degree") <- vertices_df$degree[i]
        attr(vert, "indegree") <- vertices_df$indegree[i]
        attr(vert, "outdegree") <- vertices_df$outdegree[i]
        attr(vert, "problematic") <- vertices_df$problematic[i]
        attr(vert, "longlat") <- vertices_df$longlat[i]
        attr(vert, "id") <- vertices_df$id[i]
        class(vert) <- "metric_graph_vertex"
        if (verbose == 2) {
          bar_update_attr_edges$increment()
        }
        return(vert)
      })
    } else {
      # Case when crs is not NULL
      crs <- rep(private$crs$input, n_vertices)
      vertices_df <- data.frame(
        X = self$V[, 1],
        Y = self$V[, 2],
        degree = degrees$degrees,
        indegree = degrees$indegrees,
        outdegree = degrees$outdegrees,
        problematic = problematic,
        longlat = longlat,
        crs = crs,
        id = ids
      )

      # Convert each row to a "metric_graph_vertex" object with crs
      self$vertices <- lapply(1:n_vertices, function(i) {
        vert <- as.numeric(vertices_df[i, c("X", "Y")])  # Convert to numeric vector
        attr(vert, "degree") <- vertices_df$degree[i]
        attr(vert, "indegree") <- vertices_df$indegree[i]
        attr(vert, "outdegree") <- vertices_df$outdegree[i]
        attr(vert, "problematic") <- vertices_df$problematic[i]
        attr(vert, "longlat") <- vertices_df$longlat[i]
        attr(vert, "crs") <- vertices_df$crs[i]
        attr(vert, "id") <- vertices_df$id[i]
        class(vert) <- "metric_graph_vertex"
        if (verbose == 2) {
          bar_update_attr_edges$increment()
        }
        return(vert)
      })
    }

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
      if(!is.na(ind) && length(ind)>0) {
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

  # Bounding box

  bounding_box = NULL,

  compute_bounding_box = function(){
    all_coords <- do.call(rbind, self$edges) 

    # Calculate the bounding box
    min_x <- min(all_coords[, 1])
    max_x <- max(all_coords[, 1])
    min_y <- min(all_coords[, 2])
    max_y <- max(all_coords[, 2])

    private$bounding_box <- list(min_x = min_x, max_x = max_x, min_y = min_y, max_y = max_y)
  },

  # Temp PtE

  temp_PtE = NULL,

  # # longlat
  # longlat = NULL,

  # tolerance

  tolerance = NULL,

  # transform to long lat?

  transform = FALSE,

  # edge_weights

  edge_weights = NULL,

  #
  directional_weights = NULL,

    set_first_weights = function(weights = rep(1, self$nE)){
      if(is.null(weights)){
        weights <- rep(1, self$nE)
      }
    if(!is.vector(weights) && !is.data.frame(weights)){
      stop("'weights' must be either a vector or a data.frame!")
    }
    
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

  }

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
  #' @param vertex_unit The unit in which the vertices are specified. The options are 'degree' (the great circle distance in km), 'km', 'm' and 'miles'. The default is `NULL`, which means no unit. However, if you set `length_unit`, you need to set `vertex_unit`.
  #' @param length_unit The unit in which the lengths will be computed. The options are 'km', 'm' and 'miles'. The default is `vertex_unit`. Observe that if `vertex_unit` is `NULL`, `length_unit` can only be `NULL`.
  #' If `vertex_unit` is 'degree', then the default value for `length_unit` is 'km'.
  #' @param longlat If TRUE, then it is assumed that the coordinates are given.
  #' in Longitude/Latitude and that distances should be computed in meters. It takes precedence over
  #' `vertex_unit` and `length_unit`, and is equivalent to `vertex_unit = 'degree'` and `length_unit = 'm'`.
   #' @param tolerance Vertices that are closer than this number are merged when
   #' constructing the graph (default = 1e-10). If `longlat = TRUE`, the
   #' tolerance is given in km.
   #' @param by_length Sort the components by total edge length? If `FALSE`,
   #' the components are sorted by the number of vertices.
    #' @param edge_weights Either a number, a numerical vector with length given by the number of edges, providing the edge weights, or a `data.frame` with the number of rows being equal to the number of edges, where
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
                               check_connected = FALSE, edge_weights = edge_weights,...)
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
          data_tmp_graph[[".edge_number"]] <- match(data_tmp_graph[[".edge_number"]], edge_keep)
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
