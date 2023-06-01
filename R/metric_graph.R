#' @title Metric graph
#' @description Class representing a general metric graph.
#' @details A graph object created from vertex and edge matrices, or from an
#' `sp::SpatialLines` object where each line is representing and edge. For more details,
#'  see the vignette:
#' \code{vignette("metric_graph", package = "MetricGraph")}
#' @return Object of \code{\link[R6]{R6Class}} for creating metric graphs.
#' @examples
#' library(sp)
#' line1 <- Line(rbind(c(0, 0), c(2, 0)))
#' line2 <- Line(rbind(c(2, 0), c(1, 1)))
#' line3 <- Line(rbind(c(1, 1), c(0, 0)))

#' lines <-  SpatialLines(list(Lines(list(line1), ID = "1"),
#'                            Lines(list(line2), ID = "2"),
#'                            Lines(list(line3), ID = "3")))
#' graph <- metric_graph$new(lines)
#' graph$plot()
#'
#' @export
metric_graph <-  R6::R6Class("metric_graph",
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

  #' @field EID Vector with the IDs of the edges in the graph.
  EID = NULL,


  #' @field LtE Matrix with edge positions on the lines.
  LtE = NULL,

  #' @field ELend Vector with the locations of the end points of the edges on
  #' the lines in the graph. The locations are normalized on the line.
  ELend = NULL,

  #' @field ELstart Vector with the locations of the starting points of the
  #' edges on the lines in the graph. The locations are normalized on the line.
  ELstart = NULL,

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

  #' @field lines The lines in the graph.
  lines = NULL,

  #' @field geo_dist Geodesic distances between the vertices in the graph.
  geo_dist = NULL,

  #' @field res_dist Resistance distances between the observation locations.
  res_dist = NULL,

  #' @field Laplacian The weighted graph Laplacian of the vertices in the
  #' graph. The weights are given by the edge lengths.
  Laplacian = NULL,

  #' @description Create a new `metric_graph` object.
  #' @param lines Object of type `SpatialLinesDataFrame` or `SpatialLines`.
  #' @param V n x 2 matrix with Euclidean coordinates of the n vertices.
  #' @param E m x 2 matrix where each row represents one of the m edges.
  #' @param longlat If TRUE, then it is assumed that the coordinates are given.
  #' in Longitude/Latitude and that distances should be computed in km.
  #' @param tolerance List that provides tolerances during the construction of
  #' the graph:
  #' - `vertex_vertex` Vertices that are closer than this number are merged
  #' (default = 1e-7).
  #' - `vertex_line` If a vertex at the end of one line is closer than this
  #' number to another line, this vertex is connected to that line
  #' (default = 1e-7).
  #' - `line_line` If two lines at some point are closer than this number, a new
  #' vertex is added at that point and the two lines are connected (default = 0).
  #'
  #' If `longlat = TRUE`, the tolerances are given in km.
  #' @param check_connected If `TRUE`, it is checked whether the graph is
  #' connected and a warning is given if this is not the case.
  #' @param adjust_lines Set to `TRUE` to adjust the lines object to match the graph
  #' connections. This can take some time for large graphs, so by default it is `TRUE`
  #' for graphs with at most 100 lines, and `FALSE` for larger graphs.
  #' @param remove_deg2 Set to `TRUE` to remove all vertices of degree 2 in the
  #' initialization. Default is `FALSE`.
  #' @param remove_circles All circlular edges with a length smaller than this number
  #' are removed. The default is the `vertex_vertex` tolerance.
  #' @param verbose Print progress of graph creation.
  #' @details A graph object can be initialized in two ways. The first method
  #' is to specify V and E. In this case, all edges are assumed to be straight
  #' lines. The second option is to specify the graph via the `lines` input.
  #' In this case, the vertices are set by the end points of the lines.
  #' Thus, if two lines are intersecting somewhere else, this will not be
  #' viewed as a vertex.
  #' @return A `metric_graph` object.
  initialize = function(lines = NULL,
                        V = NULL,
                        E = NULL,
                        longlat = FALSE,
                        tolerance = list(vertex_vertex = 1e-7,
                                         vertex_line = 1e-7,
                                         line_line = 0),
                        check_connected = TRUE,
                        adjust_lines = NULL,
                        remove_deg2 = FALSE,
                        remove_circles = TRUE,
                        verbose = FALSE) {

    private$longlat <- longlat

    tolerance_default = list(vertex_vertex = 1e-7,
                             vertex_line = 1e-7,
                             line_line = 0)



    for(i in 1:length(tolerance_default)){
      if(!(names(tolerance_default)[i] %in% names(tolerance))){
        tolerance[names(tolerance_default)[i]] <- tolerance_default[i]
      }
    }

    if(is.null(tolerance$buffer_line_line)){
      tolerance$buffer_line_line <- max(tolerance$line_line/2 - 1e-10,0)
    }
    max_tol <- max(c(tolerance$vertex_vertex,
                     tolerance$vertex_line,
                     tolerance$line_line))

    private$tolerance <- tolerance

    PtE_tmp_line_line <- NULL
    PtE_tmp_line_vertex <- NULL

    if(is.null(lines) && is.null(V) && is.null(E)) {
      lines <- logo_lines()
    }
    if(!is.null(lines)){
      if(!is.null(V) || !is.null(E)){
        warning("object initialized from lines, then E and V are ignored")
      }
      if (inherits(lines,"SpatialLinesDataFrame")) {
        self$lines = SpatialLines(lines@lines)
      } else if (inherits(lines,"SpatialLines")) {
        self$lines = lines
      } else {
        stop("lines should be of class SpatialLines or SpatialLinesDataFrame")
      }

    } else {
      if(is.null(V) || is.null(E)){
        stop("You must supply lines or V and E")
      }
      if(ncol(V)!=2 || ncol(E)!=2){
        stop("V and E must have two columns!")
      }
      lines <- list()
      for(i in 1:dim(E)[1]) {
        id <- sprintf("%d", i)
        lines[[i]] <- Lines(list(Line(rbind(V[E[i,1], ],
                                            V[E[i,2], ]))),
                            ID = id)
      }
      self$lines <- SpatialLines(lines)
    }

    if(verbose){
      message("Setup lines and merge close vertices")
    }

    t <- system.time(
      private$line_to_vertex(tolerance = tolerance$vertex_vertex,
                           longlat = longlat)
      )
    if(verbose){
      message(sprintf("time: %.3f s", t[["elapsed"]]))
    }

    if(is.null(adjust_lines)) {
      if(length(self$lines) < 100) {
        adjust_lines = TRUE
      } else {
        adjust_lines = FALSE
      }
    }

    if (tolerance$line_line > 0) {
    private$addinfo <- TRUE

    if(verbose){
      message("Find line-line intersections")
    }

    t <- system.time(
      points_add <- private$find_line_line_points(tol = tolerance$line_line)
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
      private$add_vertices(PtE, tolerance = tolerance$line_line)
      )

      if(verbose){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    }

    if(!is.null(private$initial_added_vertex) && adjust_lines){
      if(verbose){
        message(sprintf("Split %d lines", length(private$initial_added_vertex)))
      }

      t <- system.time(
      for(i in 1:length(private$initial_added_vertex)){

          private$split_line_at_added_vertex(private$initial_line_added[i],
                                            private$initial_added_vertex[i],
                                            private$initial_edges_added[i,])
      })

      if(verbose){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    } else if (!adjust_lines && verbose) {
      message("The lines object is not updated, so plots might not be accurate")
    }
    private$clear_initial_info()
    }

    if(tolerance$vertex_line > 0){
      private$addinfo <- TRUE
      if(verbose){
        message("Snap vertices to close lines")
      }

      t <- system.time(
        PtE_tmp <- private$coordinates_multiple_snaps(XY = self$V,
                                              tolerance = tolerance$vertex_line)
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
          private$add_vertices(PtE_tmp, tolerance = tolerance$vertex_line)
          )

        if(verbose){
          message(sprintf("time: %.3f s", t[["elapsed"]]))
        }
      }

      if(!is.null(private$initial_added_vertex) && adjust_lines){

        if(verbose){
          message(sprintf("Split %d lines: ", length(private$initial_added_vertex)))
        }

        t <- system.time(
          for(i in 1:length(private$initial_added_vertex)){
            private$split_line_at_added_vertex(private$initial_line_added[i],
                                               private$initial_added_vertex[i],
                                               private$initial_edges_added[i,])
          })

        if(verbose){
          message(sprintf("time: %.3f s", t[["elapsed"]]))
        }
      } else if (!adjust_lines && verbose) {
        message("The lines object is not updated, so plots might not be accurate")
      }
      private$clear_initial_info()
    }
    private$merge_close_vertices(tolerance$vertex_vertex, longlat = longlat)
    if(is.logical(remove_circles)){
      private$remove_circles(tolerance$vertex_vertex)
    } else {
      private$remove_circles(remove_circles)
    }


    if (remove_deg2) {
      if (verbose) {
        message("Remove degree 2 vertices")
      }
      t <- system.time(
      private$merge.all.deg2()
      )
      if(verbose){
        message(sprintf("time: %.3f s", t[["elapsed"]]))
      }
    }

    private$initial_graph <- self$clone()

    # Cloning again to add the initial graph to the initial graph
    private$initial_graph <- self$clone()

    # Checking if graph is connected
    if (check_connected) {
      g <- graph(edges = c(t(self$E)), directed = FALSE)
      components <- igraph::clusters(g, mode="weak")
      nc <- components$no
      if(nc>1){
        message("The graph is disconnected. You can use the function 'graph_components' to obtain the different connected components.")
      }
    }

    # Checking if there is some edge with infinite length
    if(any(!is.finite(self$edge_lengths))){
      warning("There is at least one edge of infinite length. Please, consider redefining the graph.")
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
      graph.temp$observation_to_vertex()
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

        graph.temp$observation_to_vertex()
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
  prune_vertices = function(){
   while(sum(self$get_degrees()==2)>0) {
     private$remove.first.deg2()
   }
   if(!is.null(self$data)){
      x_coord <- self$data[["__coord_x"]]
      y_coord <- self$data[["__coord_y"]]
      new_PtE <- self$coordinates(XY = cbind(x_coord, y_coord))
      group_vec <- self$data[["__group"]]
      self$data[["__edge_number"]] <- new_PtE[,1]
      self$data[["__distance_on_edge"]] <- new_PtE[,2]
      order_idx <- order(group_vec, new_PtE[,1], new_PtE[,2])
      self$data <- lapply(self$data, function(dat){dat[order_idx]})
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

  #' @description Gets the spatial points from the data.
  #' @return A `SpatialPoints` object of the observation locations.
  get_Spoints = function(){
     if(is.null(self$data)){
      stop("There is no data!")
    }
    group <- self$data[["__group"]]
    group <- which(group == group[1])
    Spoints <- SpatialPoints(cbind(self$data[["__coord_x"]][group],
                                   self$data[["__coord_y"]][group]))
    return(Spoints)
  },

  #' @description Adds observation locations as vertices in the graph.
  #' @param tolerance Observations locations are merged to a single vertex if
  #' they are closer than this number (given in relative edge distance between
  #' 0 and 1). The default is `1e-15`.
  #' @return No return value. Called for its side effects.
  observation_to_vertex = function(tolerance = 1e-15) {
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
  #' `distance_on_edge`, otherwise if `euclidean`, the user must provide
  #' `coord_x` and `coord_y`.
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
  #' @return No return value. Called for its side effects. The observations are
  #' stored in the `data` element of the `metric_graph` object.
  add_observations = function(Spoints = NULL,
                              data = NULL,
                              edge_number = "edge_number",
                              distance_on_edge = "distance_on_edge",
                              coord_x = "coord_x",
                              coord_y = "coord_y",
                              data_coords = c("PtE", "euclidean"),
                              group = NULL,
                              normalized = FALSE,
                              tolerance = max(self$edge_lengths)/2) {
    data_coords <- data_coords[[1]]
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
          } else if(data_coords == "euclidean"){
            point_coords <- cbind(data[[coord_x]], data[[coord_y]])
            PtE <- self$coordinates(XY = point_coords)
            XY_new <- self$coordinates(PtE = PtE, normalized = TRUE)
            norm_XY <- max(sqrt(rowSums( (point_coords-XY_new)^2 )))
            if(norm_XY > tolerance){
              warning("There was at least one point whose location is far from the graph,
                please consider checking the input.")
              }
        } else{
            stop("The options for 'data_coords' are 'PtE' and 'euclidean'.")
        }
      }
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

    for (i in 1:(dim(self$LtE)[2])) {
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
  #' elements \eqn{B_{ij} = <d\phi_i, \phi_j>}, and the vector with weights
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
    A <- Matrix(0, nrow = n, ncol = dim(self$mesh$V)[1])
    for (i in 1:n) {
      A[i, self$mesh$E[x[i, 1], 1]] = 1 - x[i, 2]
      A[i, self$mesh$E[x[i, 1], 2]] = x[i, 2]
    }
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
    }
    return(p)
  },

  #' @description Plots the connections in the graph
  #' @return No return value. Called for its side effects.
  plot_connections = function(){
        g <- graph(edges = c(t(self$E)), directed = FALSE)
        plot(g)
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
                    PtE_dim," ", dim(self$mesh$V)[1], " and ",
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
      index <- (self$LtE@p[i] + 1) : (self$LtE@p[i + 1])
      LinesPos <- cbind(self$LtE@i[index] + 1,
                        self$LtE@x[index])
      LinesPos <- LinesPos[order(LinesPos[, 2]), , drop = FALSE]
      for (j in 1:length(index)) {
        if (j==1) {
          index_j <- vals[, 1] <= LinesPos[j, 2]
        } else {
          index_j <- (vals[, 1] <= LinesPos[j, 2]) &
            (vals[, 1] > LinesPos[j - 1, 2])
        }
        if (sum(index_j) == 0)
          next
        rel.pos = vals[index_j, 1]
        if (j == 1) {
          rel.pos <- rel.pos / LinesPos[j, 2]
        } else {
          rel.pos <- (rel.pos - LinesPos[j - 1, 2]) /
            (LinesPos[j, 2] - LinesPos[j - 1, 2])
        }
        if (j == dim(LinesPos)[1])
          rel.pos = self$ELend[i] * rel.pos
        if (j==1)
          rel.pos = rel.pos + self$ELstart[i]

        data.to.plot <- cbind(rel.pos,vals[index_j, 2])
        Line_edge <- SpatialLines(list(self$lines@lines[[LinesPos[j, 1]]]))

        data.to.plot.order <- data.to.plot[order(data.to.plot[, 1]), ,
                                           drop = FALSE]

        coords <- interpolate2(Line_edge,
                               pos = data.to.plot.order[, 1, drop = TRUE],
                               normalized = TRUE)

        x.loc <- c(x.loc, coords[, 1])
        y.loc <- c(y.loc, coords[, 2])
        z.loc <- c(z.loc, data.to.plot.order[, 2])
        i.loc <- c(i.loc, rep(kk, length(coords[, 1])))
        kk = kk+1
      }
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
      index <- (self$LtE@p[i] + 1) : (self$LtE@p[i + 1])
      LinesPos <- cbind(self$LtE@i[index] + 1,
                       self$LtE@x[index])
      LinesPos <- LinesPos[order(LinesPos[, 2]), , drop = FALSE]
      for (j in 1:length(index)) {
        if (j==1) {
          index_j <- vals[, 1] <= LinesPos[j, 2]
        } else {
          index_j <- (vals[, 1] <= LinesPos[j, 2]) &
            (vals[, 1] > LinesPos[j - 1, 2])
        }
        if (sum(index_j) == 0)
          next
        rel.pos = vals[index_j, 1]
        if (j == 1) {
          rel.pos <- rel.pos / LinesPos[j, 2]
        } else {
          rel.pos <- (rel.pos - LinesPos[j - 1, 2]) /
            (LinesPos[j, 2] - LinesPos[j - 1, 2])
        }
        if (j == dim(LinesPos)[1])
          rel.pos = self$ELend[i] * rel.pos
        if (j==1)
          rel.pos = rel.pos + self$ELstart[i]

        data.to.plot <- cbind(rel.pos,vals[index_j, 2:(frames+1)])
        Line_edge <- SpatialLines(list(self$lines@lines[[LinesPos[j, 1]]]))

        data.to.plot.order <- data.to.plot[order(data.to.plot[, 1]), ,
                                           drop = FALSE]

        coords <- interpolate2(Line_edge,
                               pos = data.to.plot.order[, 1, drop = TRUE],
                               normalized = TRUE)
        x.loc <- c(x.loc, rep(coords[, 1], frames))
        y.loc <- c(y.loc, rep(coords[, 2], frames))
        z.loc <- c(z.loc, c(data.to.plot.order[, 2:(frames+1)]))
        i.loc <- c(i.loc, rep(rep(kk, length(coords[, 1])), frames))
        f.loc <- c(f.loc, rep(1:frames, each = length(coords[, 1])))
        kk = kk+1
      }
    }
    data <- data.frame(x = x.loc, y = y.loc, z = z.loc, i = i.loc, f = f.loc)

    if(plotly){
      requireNamespace("plotly")
      x <- y <- ei <- NULL
      for (i in 1:self$nE) {
        xi <- self$lines@lines[[i]]@Lines[[1]]@coords[, 1]
        yi <- self$lines@lines[[i]]@Lines[[1]]@coords[, 2]
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
      self$observation_to_vertex()
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

      LT <- private$edge_pos_to_line_pos2(PtE[, 1], PtE[, 2])
      for (i in 1:dim(PtE)[1]) {
        Points[i,] <- interpolate2(self$lines[LT[i, 1], ] ,
                                   pos = LT[i, 2], normalized = TRUE)
      }
      return(Points)
    } else {
      Spoints <- SpatialPoints(XY)
      SP <- snapPointsToLines(Spoints, self$lines)
      coords.old <- as.data.frame(Spoints@coords)
      colnames(coords.old) <- paste(colnames(coords.old), '_old', sep="")
      Spoints@coords = SP@coords
      Spoints@bbox   = SP@bbox
      LtE = cbind(match(SP@data[,1], self$EID), 0)

      for (ind in unique(LtE[, 1])) {
        index.p <- LtE[, 1] == ind
        LtE[index.p,2]=projectVecLine2(self$lines[ind,], Spoints[index.p,],
                                       normalized=TRUE)

      }
      PtE <- LtE
      for (ind in unique(LtE[, 1])) {
        Es_ind <- which(self$LtE[ind, ] > 0)
        index.p <- which(LtE[, 1] == ind)
        for (j in index.p) {
          E_ind <- which.min(replace(self$ELend[Es_ind],
                                     self$ELend[Es_ind] < LtE[j,2], NA))
          PtE[j, 1] <- Es_ind[E_ind]
          PtE[j, 2] <- (LtE[j, 2] - self$ELstart[PtE[j, 1]]) /
            (self$ELend[PtE[j, 1]] - self$ELstart[PtE[j, 1]])
        }
      }
      if (!normalized) {
        PtE <- cbind(PtE[, 1], PtE[, 2] * self$edge_lengths[PtE[, 1]])
      }
      return(PtE)
    }
  }
    ),

  private = list(
  #function for creating Vertex and Edges from self$lines
  line_to_vertex = function(tolerance = 0, longlat = FALSE) {
    lines <- c()
    for(i in 1:length(self$lines)){
      points <- self$lines@lines[[i]]@Lines[[1]]@coords
      n <- dim(points)[1]
      lines <- rbind(lines, c(i, points[1,]), c(i, points[n,]))
    }

    #save all vertices that are more than tolerance distance apart
    dists <- spDists(lines[, 2:3, drop = FALSE], longlat = longlat)
    vertex <- lines[1, , drop = FALSE]
    for (i in 2:dim(lines)[1]) {
      i.min <- which.min(dists[i, 1:(i-1)])
      if (dists[i, i.min] > tolerance) {
        vertex <- rbind(vertex, lines[i, ])
      }
    }

    lvl <- matrix(0, nrow = max(lines[,1]), 4)
    k=1
    lines_keep_id <- NULL
    for (i in 1:max(lines[, 1])) {
      which.line <- sort(which(lines[, 1] == i))
      line <- lines[which.line, ]
      #index of vertex corresponding to the start of the line
      ind1 <- which.min((vertex[, 2] - line[1, 2])^2 +
                          (vertex[, 3] - line[1, 3])^2)
      #index of vertex corresponding to the end of the line
      ind2 <- which.min((vertex[, 2] - line[2, 2])^2 +
                          (vertex[, 3] - line[2, 3])^2)

      self$lines@lines[[i]]@Lines[[1]]@coords[1,] <- vertex[ind1, 2:3]
      i.e <- dim(self$lines@lines[[i]]@Lines[[1]]@coords)[1]
      self$lines@lines[[i]]@Lines[[1]]@coords[i.e,] <- vertex[ind2, 2:3]
      ll <- LineLength(self$lines@lines[[i]]@Lines[[1]], longlat = longlat)
      if(ll > tolerance) {
        lvl[k,] <- c(i, ind1, ind2, ll)
        k=k+1
        lines_keep_id <- c(lines_keep_id, i)
      }
    }

    lvl <- lvl[1:(k-1),,drop = FALSE]
    self$lines <- self$lines[lines_keep_id]
    self$V <- vertex[, 2:3]
    self$E <- lvl[, 2:3, drop = FALSE]
    self$edge_lengths <- lvl[,4]
    self$nV <- dim(self$V)[1]
    self$nE <- dim(self$E)[1]
    self$LtE <- Matrix::sparseMatrix(j = 1:dim(self$E)[1],
                                     i = c(1:length(self$lines)),
                                     x = rep(1,dim(self$E)[1]),
                                     dims = c(dim(self$E)[1],
                                              length(self$lines)))
    self$ELend <- rep(1, dim(self$E)[1])
    self$ELstart <- rep(0, dim(self$E)[1])
    self$EID <- as.vector(sapply(slot(self$lines,"lines"),
                                 function(x) slot(x, "ID")))
  },

  #Compute PtE for mesh given PtE for graph
  PtE_to_mesh = function(PtE){
    VtE <- rbind(self$VtEfirst(), self$mesh$PtE)
    PtE_update <- matrix(0, dim(PtE)[1], 2)
    for (i in 1:dim(PtE)[1]) {
      ei <- PtE[i, 1]
      ind <- which(VtE[,1] == ei)
      if (PtE[i,2]<min(VtE[ind, 2])) { #node on first edge on line
        v1 <- self$E[ei, 1]
        ind2 <- which.min(VtE[ind, 2])
        v2 <- ind[ind2]
        d1 <- 0
        d2 <- VtE[v2, 2]
      } else if (PtE[i,2] > max(VtE[ind,2])) { #node on last edge on line
        ind2 <- which.max(VtE[ind, 2])
        v1 <- ind[ind2]
        v2 <- self$E[ei, 2]
        d1 <- VtE[v1, 2]
        d2 <- 1
      } else {
        ind2 <- sort(sort(abs(VtE[ind, 2] - PtE[i, 2]),
                          index.return = TRUE)$ix[1:2])
        v1 <- ind[ind2[1]]
        v2 <- ind[ind2[2]]
        d1 <- VtE[v1, 2]
        d2 <- VtE[v2, 2]
      }
      #edge on mesh
      e <- which(rowSums((self$mesh$E == v1) + (self$mesh$E == v2)) == 2)
      #distance on edge
      if (self$mesh$E[e, 1] == v1) {
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

    coords <- lapply(coordinates(self$lines), function(x) x[[1]])
    nc <- do.call(rbind,lapply(coords, function(x) dim(x)[1]))
    xyl <- cbind(do.call(rbind,coords), rep(1:length(nc), times = nc))

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
        xi <- self$lines@lines[[i]]@Lines[[1]]@coords[, 1]
        yi <- self$lines@lines[[i]]@Lines[[1]]@coords[, 2]
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

    coordinates_multiple_snaps = function(XY, tolerance){
      Spoints <- SpatialPoints(XY)
      coords_line <- c()
      coords_tmp <- c()
      Spoints_sf <- sf::st_as_sf(Spoints)
      lines_sf <- sf::st_as_sf(self$lines)

      # within_dist <- gWithinDistance(Spoints, self$lines, dist = tolerance, byid = TRUE)
      within_dist <- t(as.matrix(sf::st_is_within_distance(Spoints_sf, lines_sf, dist = tolerance)))


      for(i in 1:length(self$lines)){
        select_points <- matrix(XY[within_dist[i,],], ncol=2)
        if(nrow(select_points) > 0){
          SP_tmp <- SpatialPoints(select_points)
          SP <- snapPointsToLines(SP_tmp, self$lines[i])
          idx_tol <- (SP@data[["snap_dist"]] <= tolerance)
          coords_line <- c(coords_line, SP@data[["nearest_line_id"]][idx_tol])
          coords_tmp <- rbind(coords_tmp, SP@coords[idx_tol,])
        }
      }

      Spoints@coords = coords_tmp

      LtE = cbind(match(coords_line, self$EID), 0)

      for (ind in unique(LtE[, 1])) {
        index.p <- LtE[, 1] == ind
        # LtE[index.p,2]=rgeos::gProject(self$lines[ind,], Spoints[index.p,],
        #                                normalized=TRUE)

        LtE[index.p,2]=projectVecLine2(self$lines[ind,], Spoints[index.p,],
                                       normalized=TRUE)
      }
      PtE <- LtE

      for (ind in unique(LtE[, 1])) {

        Es_ind <- which(self$LtE[ind, ] > 0)

        index.p <- which(LtE[, 1] == ind)

        for (j in index.p) {
          E_ind <- which.min(replace(self$ELend[Es_ind],
                                     self$ELend[Es_ind] < LtE[j,2], NA))
          PtE[j, 1] <- Es_ind[E_ind]
          PtE[j, 2] <- (LtE[j, 2] - self$ELstart[PtE[j, 1]]) /
            (self$ELend[PtE[j, 1]] - self$ELstart[PtE[j, 1]])
        }
      }
      return(PtE)
      },

  # Version 2 edge_pos_to_line_pos
  # Gets relative position on the line
  edge_pos_to_line_pos2 = function(E_i, t_i){
    stopifnot(length(E_i) == length(t_i))
    tmp_matrix <- matrix(NA, nrow=length(E_i), ncol=2)

    for(i in 1:length(E_i)){
      tmp_matrix[i,1] <- which(self$LtE[,E_i[i]] == 1)
    }
      start_points <- self$ELstart[E_i]
      tmp_matrix[,2] <- t_i * (self$ELend[E_i] - self$ELstart[E_i]) + self$ELstart[E_i]
    return(tmp_matrix)
    # line_E_i <- which(self$LtE[,E_i] == 1)
    # start_point <- self$ELstart[E_i]
    # exact_point <- t_i * (self$ELend[E_i] - self$ELstart[E_i]) + self$ELstart[E_i]
    # return(cbind(line_E_i, exact_point))
  },

  #utility function to merge close vertices
  merge_close_vertices = function(tolerance, longlat) {
    if(tolerance > 0) {
      dists <- spDists(self$V, longlat = longlat)
      v.merge <- NULL
      k <- 0
      for (i in 2:self$nV) {
        i.min <- which.min(dists[i, 1:(i-1)])
        if (dists[i, i.min] < tolerance) {
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
          self$V <- self$V[-v.rem,]
          self$nV <- self$nV - 1
          self$E[self$E == v.rem] <- v.keep
          self$E[self$E > v.rem] <- self$E[self$E > v.rem] - 1
        }
      }
    }
  },

  # utility function to remove small circles
  remove_circles = function(threshold) {
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

          self$lines <- self$lines[-ind]
          self$E <- self$E[-ind,]
          self$EID <- self$EID[-ind]
          self$edge_lengths <- self$edge_lengths[-ind]
          self$ELend <- self$ELend[-ind]
          self$ELstart <- self$ELstart[-ind]
          self$nE <- self$nE - length(ind)
          self$LtE <- self$LtE[-ind,-ind]
        }
      }
    }
  },

  # utility function to merge lines connected by degree 2 vertices
  merge.all.deg2 = function() {
   while(sum(self$get_degrees()==2)>0) {
     private$remove.first.deg2()
   }
  },

  remove.first.deg2 = function() {
    ind <- which(self$get_degrees()==2)
    if(length(ind)>0) {
      ind <- ind[1]
      e1 <- which(self$E[,2]==ind)
      e2 <- which(self$E[,1]==ind)
      e_rem <- sort(c(e1,e2))
      v1 <- setdiff(self$E[e_rem[1],],ind)
      v2 <- setdiff(self$E[e_rem[2],],ind)
      if(v1 > ind) {
        v1 <- v1-1
      }
      if(v2 > ind) {
        v2 <- v2 - 1
      }

      e_remidx <- which(self$LtE[,e_rem[2]] == 1)

      if(e_remidx == e_rem[2]){
                line_keep1 <- line_keep2 <- NULL
                if(e_rem[1]>1) {
                  line_keep1 <- self$lines[1:(e_rem[1]-1)]
                }

                line_keep2 <- self$lines[setdiff((e_rem[1]+1):length(self$lines),e_rem[2])]

                line_merge <- list()
                coords <- self$lines@lines[[e_rem[1]]]@Lines[[1]]@coords
                tmp <- self$lines@lines[[e_rem[2]]]@Lines[[1]]@coords
                diff_ss <- norm(as.matrix(coords[1,] - tmp[1,]))
                diff_se <- norm(as.matrix(coords[1,] - tmp[dim(tmp)[1],]))
                diff_es <- norm(as.matrix(coords[dim(coords)[1],] - tmp[1,]))
                diff_ee <- norm(as.matrix(coords[dim(coords)[1],] - tmp[dim(tmp)[1],]))
                diffs <- c(diff_ss, diff_se, diff_es, diff_ee)
                if(which.min(diffs) == 1) {
                  coords <- rbind(coords[rev(1:dim(coords)[1]),], tmp)
                  E_new <- matrix(c(v2,v1),1,2)
                } else if(which.min(diffs)==2){
                  coords <- rbind(tmp,coords)
                  E_new <- matrix(c(v2,v1),1,2)
                } else if(which.min(diffs)==3) {
                  coords <- rbind(coords, tmp)
                  E_new <- matrix(c(v1,v2),1,2)
                } else {
                  coords <- rbind(coords, tmp[rev(1:dim(tmp)[1]),])
                  E_new <- matrix(c(v1,v2),1,2)
                }
                line_merge <-  Lines(list(Line(coords)), ID = sprintf("new%d",1))

                if(!is.null(line_keep1) && !is.null(line_keep2)) {
                  line_new <- SpatialLines(c(line_keep1@lines, line_merge, line_keep2@lines))
                } else if (is.null(line_keep1)) {
                  line_new <- SpatialLines(c(line_merge, line_keep2@lines))
                } else if (is.null(line_keep2)) {
                  line_new <- SpatialLines(c(line_keep1@lines, line_merge))
                } else {
                  line_new <- SpatialLines(c(line_merge))
                }
                for(i in 1:length(line_new)) {
                  slot(line_new@lines[[i]],"ID") <- sprintf("%d",i)
                }

                #update lines
                self$lines <- line_new
                self$EID <- as.vector(sapply(slot(self$lines,"lines"), function(x) slot(x, "ID")))
                self$LtE <- self$LtE[-e_rem[2],-e_rem[2]]
      } else{
                E_new1 <- self$E[e1,1]
                E_new2 <- self$E[e2,2]
                if(E_new1 > ind){
                  E_new1 <- E_new1 - 1
                }
                if(E_new2 > ind){
                  E_new2 <- E_new2 -1
                }
                E_new <- matrix(c(E_new1, E_new2),1,2)
                self$LtE <- self$LtE[,-e_rem[2]]
      }

      #update vertices
      self$V <- self$V[-ind,]
      self$nV <- self$nV - 1

      #update edges
      self$E[self$E >= ind] <- self$E[self$E >= ind] - 1
      self$E <- self$E[-e_rem[2],,drop=FALSE]
      self$E[e_rem[1],] <- E_new
      self$EID <- self$EID[-ind]
      self$edge_lengths[e_rem[1]] <- self$edge_lengths[e_rem[1]] + self$edge_lengths[e_rem[2]]
      self$edge_lengths <- self$edge_lengths[-e_rem[2]]
      self$ELend <- self$ELend[-e_rem[2]]
      self$ELstart <- self$ELstart[-e_rem[2]]
      self$nE <- self$nE - 1
    }
  },

  # Vertex added in the initial processing

  initial_added_vertex = NULL,

  # Initial lines added

  initial_line_added = NULL,

  # Initial edges added

  initial_edges_added = NULL,

  # Initial graph

  initial_graph = NULL,

  # Conjugate_line

  conjugate_initial_line = NULL,

  # pruned

  pruned = FALSE,

  # Split lines

  split_lines_ids_times = NULL,

  # should the information be saved when splitting edges?

  addinfo = FALSE,

  clear_initial_info = function(){
    private$initial_added_vertex = NULL
    private$initial_line_added = NULL
    private$addinfo = FALSE
    private$initial_edges_added = NULL
    private$conjugate_initial_line = NULL
  },

  get_list_coords = function(pos, line1, line2){
          new_lines <- c()
          for(i in 1:length(self$lines)){
            if(i!=pos){
              line_coords <- coordinates(self$lines[i])[[1]][[1]]
              id_line <- self$lines[i]@lines[[1]]@ID
              new_lines <- c(new_lines, Lines(list(Line(line_coords)), ID = id_line))
            } else{
              new_lines <- c(new_lines, line1)
            }
          }
          if(!is.null(line2)){
            new_lines <- c(new_lines, line2)
          }
          return(new_lines)
  },

  split_line_at_added_vertex = function(integer_id_line,
                                        added_vertex_id,
                                        added_edges_id){

      id_conjugate_line <- which(private$conjugate_initial_line[, 1] == integer_id_line)
      id_lines_of_interest <- c(integer_id_line,
                                private$conjugate_initial_line[id_conjugate_line, 2])


      distances_lines <- c()
      idx_min <- c()
      for(id_ in id_lines_of_interest){
          line <- self$lines@lines[id_]
          # tmp_dist <- rgeos::gDistance(SpatialPoints(coords = matrix(self$V[added_vertex_id,],ncol=2)), SpatialLines(line))
          tmp_dist <- distance2(SpatialPoints(coords = matrix(self$V[added_vertex_id,],ncol=2)), SpatialLines(line))
          distances_lines <- c(distances_lines, tmp_dist)
        }

        min_idx <- which.min(distances_lines)

        private$conjugate_initial_line <- rbind(private$conjugate_initial_line,
                                                      cbind(integer_id_line, nrow(self$LtE)+1))

        integer_id_line <- id_lines_of_interest[min_idx][[1]]
        line_coords <- self$lines@lines[integer_id_line][[1]]@Lines[[1]]@coords

        tmp_dist <- sapply(1:nrow(line_coords), function(j){
              (self$V[added_vertex_id,1]-line_coords[j,1])^2 + (self$V[added_vertex_id,2] - line_coords[j,2])^2
          })
        idx_tmp <-  which.min(tmp_dist)

        closest_coord <- idx_tmp

        if(closest_coord == nrow(line_coords)){
          closest_coord <- closest_coord - 1
        }


      slot_id_lines <- sapply(slot(self$lines,"lines"), function(x) slot(x, "ID"))
      id_line <- slot_id_lines[integer_id_line]

      times_split <- sum(which(private$split_lines_ids_times[,1] == integer_id_line))
      private$split_lines_ids_times <- rbind(private$split_lines_ids_times,
                                             cbind(integer_id_line, 1))

        coords1 <- rbind(matrix(line_coords[1:closest_coord,],ncol=2),
                         matrix(self$V[added_vertex_id,],ncol=2))
                    line1 <- Lines(list(Line(coords1)), ID = id_line)

        coords2 <- rbind(matrix(self$V[added_vertex_id,], ncol=2),
                         matrix(line_coords[(closest_coord+1):nrow(line_coords),],
                                ncol=2))
                line2 <- Lines(list(Line(coords2)),
                               ID = paste0(id_line, "__",
                                           as.character(times_split+1)))

        self$lines <- SpatialLines(private$get_list_coords(integer_id_line,
                                                           line1,
                                                           line2))

        self$LtE <- rbind(self$LtE, rep(0, self$nE))
        self$LtE[nrow(self$LtE), added_edges_id[2]] <- 1
        if(length(id_conjugate_line) == 0){
          id_original_line <- integer_id_line
        } else{
          id_original_line <- private$conjugate_initial_line[id_conjugate_line, 1]
        }

        self$LtE[id_original_line, added_edges_id[2]] <- 0
        self$ELend[added_edges_id[1]] <- 1
        self$ELstart[added_edges_id[2]] <- 0
        self$ELend[added_edges_id[2]] <- 1
        self$edge_lengths[added_edges_id[1]] <- LineLength(line1@Lines[[1]])
        self$edge_lengths[added_edges_id[2]] <- LineLength(line2@Lines[[1]])
        self$EID <- c(self$EID, paste0(id_line, "__", as.character(times_split+1)))
  },

  # @description function for splitting lines in the graph
  # @param Ei index of line to split
  # @param t  position on line to split (normalized)
  # @param tolerance tolerance for merging overlapping vertices
  split_edge = function(Ei, t, tolerance = 0) {

    index <- (self$LtE@p[Ei] + 1):(self$LtE@p[Ei + 1])
    LinesPos <- cbind(self$LtE@i[index] + 1, self$LtE@x[index])
    LinesPos <- LinesPos[order(LinesPos[, 2]), , drop = FALSE]
    j <- min(which(t <= LinesPos[, 2]))
    t_mod <- t
    if (j == 1) {
      t_mod <- t_mod/LinesPos[j, 2]
    } else {
      t_mod <- (t_mod - LinesPos[j - 1, 2]) /
        (LinesPos[j, 2] - LinesPos[j - 1, 2])
    }
    mult_ <- 1
    if(j== dim(LinesPos)[1] ){
      mult_ <- self$ELend[Ei]
    }
    if(j==1)
      mult_ = mult_ - self$ELstart[Ei]
    t_mod = mult_ * t_mod

    if(j==1){
      t_mod = t_mod + self$ELstart[Ei]
    }

    Line <- self$lines[LinesPos[j,1], ]
    # val_line <- rgeos::gInterpolate(Line, t_mod, normalized = TRUE)@coords
     val_line <- interpolate2(Line, pos = t_mod, normalized = TRUE)

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
        #change LtE
        self$ELend <- c(self$ELend, self$ELend[Ei])
        self$ELend[Ei] <- t_mod
        self$ELstart <- c(self$ELstart, t_mod)
        LtE.i <- self$LtE@i
        LtE.p <- self$LtE@p
        LtE.x <- self$LtE@x
        LtE.dim <- self$LtE@Dim
        LtE.dim[2] <- LtE.dim[2] + 1
        # add the new column
        LtE.i_new <- LinesPos[j,1] - 1
        LtE.x_new <- LinesPos[j,2]
        large <- LtE.x[index] > t

        if(private$addinfo){
          private$initial_added_vertex <- c(private$initial_added_vertex, newV)
          private$initial_line_added <- c(private$initial_line_added, LtE.i_new+1)
        }


        n.p <- length(self$LtE@p)
        if(sum(large)>1){
          index.rem = index[large]
          index.rem <- index.rem[-length(index.rem)]
          LtE.i_new = c(LtE.i_new,LtE.i[index.rem])
          LtE.x_new = c(LtE.x_new,LtE.x[index.rem])
          LtE.i <- LtE.i[-index.rem]
          LtE.x <- LtE.x[-index.rem]
          self$LtE@p[(Ei+1):n.p] <- self$LtE@p[(Ei+1):n.p] - sum(large) - 1
        }
        LtE.p_new <- self$LtE@p[n.p] + sum(large)
        self$LtE <- Matrix::sparseMatrix(i = c(LtE.i, LtE.i_new)+1,
                                         p = c(LtE.p, LtE.p_new),
                                         x = c(LtE.x, LtE.x_new),
                                         dims = LtE.dim)

        if(add_V){
          self$V <- rbind(self$V, c(val_line))
        }
        l_e <- self$edge_lengths[Ei]
        self$edge_lengths[Ei] <- t * l_e
        self$edge_lengths <- c(self$edge_lengths, (1 - t) * l_e)
        self$nE <- self$nE + 1
        self$E <- rbind(self$E, c(newV, self$E[Ei, 2]))
        self$E[Ei, 2] <- newV
        self$nV <- dim(self$V)[1]

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
    graph.temp$observation_to_vertex()
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

  find_line_line_points = function(tol) {
  lines_sf <- sf::st_as_sf(self$lines)
  # dists <- gWithinDistance(self$lines, dist = tol, byid = TRUE)
  dists <- t(as.matrix(sf::st_is_within_distance(lines_sf, dist = tol)))
  points_add <- NULL
  points_add_PtE <- NULL
  for(i in 1:(length(self$lines)-1)) {
    #lines within tol of line i
    inds <- i+which(as.vector(dists[i, (i+1):length(self$lines)]))
    if(length(inds)>0) {
      for(j in inds) {
        #first check if there are intersections
        # intersect_tmp <- rgeos::gIntersection(self$lines[i], self$lines[j])
        intersect_tmp <- intersection2(self$lines[i], self$lines[j])

        if( "GEOMETRYCOLLECTION" %in% sf::st_geometry_type(intersect_tmp)){
          intersect_tmp <- sf::st_collection_extract(intersect_tmp, type = "LINESTRING")
          tmp_inter <- 1
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
            intersect_tmp <- sf::as_Spatial(intersect_tmp$geom)
            # coord_tmp <-gInterpolate(intersect_tmp, d=0.5, normalized = TRUE)
            coord_tmp <- interpolate2(intersect_tmp, pos=0.5, normalized = TRUE)
            coord_tmp <- matrix(coordinates(coord_tmp),1,2)
          }

            # for(k in 1:length(intersect_tmp)) {
            coord_tmp <- matrix(coord_tmp, ncol=2)
            for(k in 1:nrow(coord_tmp)){
              p <- matrix(coord_tmp[k,],1,2)
              #add points if they are not close to V or previous points
              if(min(spDists(self$V, p))>tol) {
                p_cur <- rbind(p_cur,p)
                p2 <- snapPointsToLines(SpatialPoints(p),self$lines[i])
                points_add <- rbind(points_add, p, coordinates(p2))
                # points_add_PtE <- rbind(points_add_PtE,
                #                         c(i,gProject(self$lines[i],
                #                                      SpatialPoints(p))),
                #                         c(j,gProject(self$lines[j],SpatialPoints(p))))
                points_add_PtE <- rbind(points_add_PtE,
                                        c(i,projectVecLine2(self$lines[i],
                                                     SpatialPoints(p))),
                                        c(j,projectVecLine2(self$lines[j],SpatialPoints(p))))

              }
            }
        }
        #now check if there are intersections with buffer
        # tmp_line <- gBuffer(self$lines[i], width = tol)
        lines_tmp_sf <- sf::st_as_sf(self$lines[i])
        lines2_tmp_sf <- sf::st_as_sf(self$lines[j])
        tmp_line <- sf::st_buffer(lines_tmp_sf, dist = tol)
        # intersect_tmp <- rgeos::gIntersection(tmp_line, self$lines[j])
        # intersect_tmp <- intersection2(tmp_line, self$lines[j])
        intersect_tmp <- intersection3(tmp_line, lines2_tmp_sf)
        if( "GEOMETRYCOLLECTION" %in% sf::st_geometry_type(intersect_tmp)){
          intersect_tmp <- sf::st_collection_extract(intersect_tmp, type = "LINESTRING")
        }
        # if(!is.null(intersect_tmp)) {
        if(nrow(sf::st_coordinates(intersect_tmp))>0){
            # for(k in 1:length(intersect_tmp)) {
            intersect_tmp <- sf::as_Spatial(intersect_tmp$geom)
            # for(k in 1:nrow(coord_tmp)){
            for(k in 1:length(intersect_tmp)) {
            if(inherits(intersect_tmp, "SpatialLines")) {
              # coord_tmp <-gInterpolate(intersect_tmp[k], d=0.5, normalized = TRUE)
              coord_tmp <- interpolate2(intersect_tmp[k], pos=0.5, normalized = TRUE)
              p <- matrix(coordinates(coord_tmp),1,2)
            } else {
              p <- matrix(coordinates(intersect_tmp[k]),1,2)
            }
            #add points if they are not close to V or previous points
            if(min(spDists(self$V, p))>tol) {
              # if(is.null(p_cur) || gDistance(SpatialPoints(p_cur), intersect_tmp[k])>tol) {
              if(is.null(p_cur) || distance2(SpatialPoints(p_cur), intersect_tmp[k])>tol) {
                p2 <- snapPointsToLines(SpatialPoints(p),self$lines[i])
                points_add <- rbind(points_add, p, coordinates(p2))
                # points_add_PtE <- rbind(points_add_PtE,
                #                         c(i,gProject(self$lines[i],SpatialPoints(p))),
                #                         c(j,gProject(self$lines[j],SpatialPoints(p))))
                points_add_PtE <- rbind(points_add_PtE,
                                        c(i,projectVecLine2(self$lines[i],SpatialPoints(p))),
                                        c(j,projectVecLine2(self$lines[j],SpatialPoints(p))))

              }
            }
          }
        }
      }
    }
  }
  p_add <- points_add
  p_pte <- points_add_PtE
  # ret_tmp <- list(points = points_add, PtE = points_add_PtE)
  return(list(points = points_add, PtE = points_add_PtE))
},

add_vertices = function(PtE, tolerance = 1e-10) {

  e.u <- unique(PtE[,1])
  for (i in 1:length(e.u)) {
    dists <- sort(PtE[which(PtE[,1]==e.u[i]),2])
    private$split_edge(e.u[i], dists[1], tolerance)
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

  # longlat
  longlat = NULL,

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
#' line1 <- Line(rbind(c(0, 0), c(1, 0)))
#' line2 <- Line(rbind(c(1, 0), c(2, 0)))
#' line3 <- Line(rbind(c(1, 1), c(2, 1)))

#' Lines <-  SpatialLines(list(Lines(list(line1), ID = "1"),
#'                            Lines(list(line2), ID = "2"),
#'                            Lines(list(line3), ID = "3")))
#' graphs <- graph_components$new(Lines)
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
   #' @param lines Object of type `SpatialLinesDataFrame` or `SpatialLines`.
   #' @param V n x 2 matrix with Euclidean coordinates of the n vertices.
   #' @param E m x 2 matrix where each row represents an edge.
   #' @param longlat If `TRUE`, then it is assumed that the coordinates are given
   #' in Longitude/Latitude and that distances should be computed in km.
   #' @param tolerance Vertices that are closer than this number are merged when
   #' constructing the graph (default = 1e-10). If `longlat = TRUE`, the
   #' tolerance is given in km.
   #' @param by_length Sort the components by total edge length? If `FALSE`,
   #' the components are sorted by the number of vertices.
   #' @param ... Additional arguments used when specifying the graphs
   #' @return A `graph_components` object.
   initialize = function(lines = NULL,
                         V = NULL,
                         E = NULL,
                         by_length = TRUE,
                         ...) {
     graph <- metric_graph$new(lines = lines, V = V, E = E,
                               check_connected = FALSE, ...)

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
         ind_keep <- graph$LtE%*%ind_keep
         self$graphs[[k]] = metric_graph$new(lines = graph$lines[which(ind_keep!=0)],
                                             check_connected = FALSE, ...)
       }
       self$sizes <- components$csize
       self$lengths <- unlist(lapply(1:self$n,
                                     function(x) sum(self$graphs[[x]]$edge_lengths)))
       if(by_length) {
         reo <- sort(self$lengths, decreasing = TRUE, index.return = TRUE)$ix
       } else {
         reo <- sort(self$sizes, decreasing = TRUE, index.return = TRUE)$ix
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
