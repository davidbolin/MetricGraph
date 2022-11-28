#' @title Metric graph object for specification of Gaussian processes
#' @description Class representing a general metric graphs.
#' @details A graph object created from vertex and edge matrices, or from an
#' sp::Lines object where each line is representing and edge. For more details,
#'  see the help vignette:
#' \code{vignette("metric_graph", package = "MetricGraph")}
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
  #' @field V matrix with positions in Euclidean space of the vertices of the
  #' graph
  V = NULL,

  #' @field nV the number of vertices
  nV = 0,

  #' @field E matrix with the edges of the graph, where `E[i,1]` is the vertex
  #' at the start of the ith edge and `E[i,2]` is the vertex at the end of the
  #' edge
  E = NULL,

  #' @field nE the number of edges
  nE= 0,

  #' @field edge_lengths vector with the lengths of the edges in the graph
  edge_lengths = NULL,

  #' @field EID vector with the IDs of the edges in the graph
  EID = NULL,


  #' @field LtE matrix with edge positions on the lines
  LtE = NULL,

  #' @field ELend vector with the locations of the end points of the edges on
  #' the lines in the graph. The locations are normalized on the line
  ELend = NULL,

  #' @field ELstart vector with the locations of the starting points of the
  #' edges on the lines in the graph. The locations are normalized on the line
  ELstart = NULL,

  #' @field C constraint matrix used to set Kirchhoff constraints
  C = NULL,

  #' @field CoB change-of-basis object used for Kirchhoff constraints
  CoB = NULL,

  #' @field points the observations in a SpatialPointsDataFrame
  points  = NULL,

  #' @field data a list containing the data on the graph
  data = NULL,

  #' @field PtE matrix specifying the locations of the observation points on
  #' the edges, where `PtE[i,1]` is the edge index for the ith observation
  #' and  `PtE[,2]` is the normalized distance on the edge
  PtE = NULL,

  #' @field PtV vector with the indices of the vertices which are observation
  #' locations
  PtV  = NULL,

  #' @field mesh mesh object used for plotting
  mesh = NULL,

  #' @field lines the lines in the graph
  lines = NULL,

  #' @field geo_dist geodesic distances between the vertices in the graph
  geo_dist = NULL,

  #' @field res_dist resistance distances between the observation locations
  res_dist = NULL,

  #' @field Laplacian the weighted graph Laplacian of the vertices in the
  #' graph. The weights are given by the edge lengths
  Laplacian = NULL,

  #' @description Create a new `metric_graph` object
  #' @param lines object of type `SpatialLinesDataFrame` or `SpatialLines`
  #' @param V n x 2 matrix with Euclidean coordinates of the n vertices
  #' @param E m x 2 matrix where each row represents an edge
  #' @param longlat If TRUE, then it is assumed that the coordinates are given
  #' in Longitude/Latitude and that distances should be computed in km.
  #' @param tolerance vertices that are closer than this number are merged when
  #' constructing the graph (default = 1e-10). If `longlat = TRUE`, the
  #' tolerance is given in km.
  #' @details A graph object can be initialized in two ways. The first method
  #' is to specify V and E. In this case, all edges are assumed to be straight
  #' lines. The second option is to specify the graph via the `lines` input.
  #' In this case, the vertices are set by the end points of the lines.
  #' Thus, if two lines are intersecting somewhere else, this will not be
  #' viewed as a vertex.
  #' @return A metric_graph object
  initialize = function(lines = NULL,
                        V = NULL,
                        E = NULL,
                        longlat = FALSE,
                        tolerance = 1e-10) {

    if(!is.null(lines)){
      if(!is.null(V) || !is.null(E)){
        warning("object initialized from lines, then E and V are ignored")
      }
      self$lines = lines
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
        lines[[i]] <- Lines(list(Line(rbind(V[E[i,1], ], V[E[i,2], ]))), ID = id)
      }
      self$lines <- SpatialLines(lines)
    }
    self$EID = sapply(slot(self$lines,"lines"), function(x) slot(x, "ID"))
    private$line_to_vertex(tolerance = tolerance, longlat = longlat)
    private$initial_graph <- self$clone()

    # Checking if graph is connected
    g <- graph(edges = c(t(self$E)), directed = FALSE)
    components <- igraph::clusters(g, mode="weak")
    nc <- components$no
    if(nc>1){
      warning("The graph is disconnected. You can use the function 'graph_components' to obtain the different connected components.")
    }

    # Checking if there is some edge with infinite length
    if(any(!is.finite(self$edge_lengths))){
      warning("There is at least one edge of infinite length. Please, consider redefining the graph.")
    }
  },

  #' @description Computes shortest path distances between the vertices in the
  #' graph
  #' @param full Should the geodesic distances be computed for all
  #' the available locations. If `FALSE`, it will be computed
  #' separately for the locations of each replicate.
  #' @param repl vector or list containing which replicates
  #' to compute the distance. If `NULL`, it will be computed
  #' for all replicates.
  compute_geodist = function(full = FALSE, repl = NULL) {
    self$geo_dist <- list()
    if(full){
      g <- graph(edges = c(t(self$E)), directed = FALSE)
      E(g)$weight <- self$edge_lengths
      self$geo_dist[["__complete"]] <- distances(g)
    } else{
      if(is.null(repl)){
        lapply()
      } else{

      }
    }
  },

  #' @description Computes shortest path distances between the vertices in the
  #' mesh
  compute_geodist_mesh = function(full = FALSE, repl=NULL) {
    g <- graph(edges = c(t(self$mesh$E)), directed = FALSE)
    E(g)$weight <- self$mesh$h_e
    self$mesh$geo_dist <- distances(g)
  },

  #' @description Computes the resistance distance between the observation
  #' locations
  #' @param PtE points to compute the metric for, if not provided, the metric
  #' is computed and stored for the observations in the graph
  #' @param normalized are the locations in PtE in normalized distance?
  compute_resdist = function(PtE = NULL, normalized = FALSE) {
    if (is.null(PtE)) {
      graph.temp <- self$clone()
      if(is.null(graph.temp$PtV)) {
        graph.temp$observation_to_vertex()
      }
      if(is.null(graph.temp$geo_dist)){
        graph.temp$compute_geodist()
      }
      L <- Matrix(0, graph.temp$nV, graph.temp$nV)
      for (i in 1:graph.temp$nE) {
        tmp <- -1 / graph.temp$geo_dist[graph.temp$E[i, 1], graph.temp$E[i, 2]]
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

      self$res_dist <- R[graph.temp$PtV, graph.temp$PtV]
      reo <- order(self$PtE[,1],self$PtE[,2])
      self$res_dist[reo, reo] <- as.matrix(self$res_dist)
      self$res_dist <- as.matrix(self$res_dist)
    } else {
      graph.temp <- self$clone()
      graph.temp$clear_observations()
      df_temp <- data.frame(y = rep(0, dim(PtE)[1]),
                          edge_number = PtE[,1],
                          distance_on_edge = PtE[,2])
      graph.temp$add_PtE_observations(data_frame = df_temp,
                                   normalized = normalized)
      graph.temp$compute_resdist()
      return(graph.temp$res_dist)
    }
  },

  #' @description Computes the resistance metric between the vertices in the
  #' mesh
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

  #' @description Computes the weigthed graph Laplacian for the graph
  compute_laplacian = function() {
    Wmat <- Matrix(0,self$nV,self$nV)
    for (i in 1:self$nE) {
      Wmat[self$E[i, 1], self$E[i, 2]] <- 1/self$edge_lengths[i]
      Wmat[self$E[i, 2], self$E[i, 1]] <- 1/self$edge_lengths[i]
    }
    self$Laplacian <- Matrix::Diagonal(self$nV,as.vector(Matrix::rowSums(Wmat))) - Wmat
  },

  #' @description Gets PtE from the data
  
  get_PtE = function() {
    if(is.null(self$data)){
      stop("There is no data!")
    }
    repl <- self$data[["__repl"]]
    repl <- repl[repl == repl[1]]

    PtE <- cbind(self$data[["__edge_number"]][repl], 
                self$data[["__distance_on_edge"]][repl])
    return(PtE)
  },

  get_Spoints <- function(){
     if(is.null(self$data)){
      stop("There is no data!")
    }
    repl <- self$data[["__repl"]]
    repl <- repl[repl == repl[1]]
    Spoints <- SpatialPoints(cbind(self$data[["__coord_x"]][repl], 
                                        self$data[["__coord_y"]][repl]))
    return(Spoints)
  },

  #' @description Adds observation locations as vertices in the graph
  observation_to_vertex = function() {
    # Reordering
    PtE <- get_PtE()
    order_idx <- order(PtE[, 1], PtE[, 2])

    PtE <- PtE[order_idx, ]

    if (length(order_idx) == 1) {
      PtE <- matrix(PtE, ncol = 2)
    }
    l <- length(PtE[, 1])
    self$PtV <- rep(0, l)
    for (i in 1:l) {
        e <- as.vector(PtE[i, 1])
        t <- as.vector(PtE[i, 2])
        l_e <- self$edge_lengths[e]
        if (abs(t) < 10^-10) {
          PtE[i, 2] <- 0
          self$PtV[i] <- self$E[e, 1]
        } else if (t > 1 - 10^-10) {
          PtE[i, 2] <- 1
          self$PtV[i] <- self$E[e, 2]
        } else {
          self$split_edge(e, t)
          self$PtV[i] <- dim(self$V)[1]
        }
    }
    
    # Updates the columns `__edge_number` and `__distance_on_edge`
    # and reorders the data. 

    self$data <- lapply(self$data, function(data){return(dat[order_idx,])})
    self$data[["__edge_number"]] <- PtE[ ,1]
    self$data[["__distance_on_edge"]] <- PtE[ ,2]

    self$mesh$PtE <- self$coordinates(XY = self$mesh$V[(nrow(self$VtEfirst()) + 1):nrow(self$mesh$V), ])

    if (!is.null(self$geo_dist)) {
      self$compute_geodist()
    }
    if (!is.null(self$res_dist)) {
      self$compute_resdist()
    }
    if (!is.null(self$CoB)) {
      self$buildC(2)
    }

    # Now we compute an on the method

  },

  #' @description Clear all observations from the object
  clear_observations = function() {
   self$data <- NULL
  },

  #' @description Add observations to the graph
  #' @param Spoints SpatialPoints or SpatialPointsDataFrame of the observations,
  #' which may include the coordinates only, or the coordinates as well as the
  #' observations.
  #' @param data A data.frame or named list containing the observations. In case of replicates, the data.frames for the replicates should stacked vertically, with a column
  #' indicating the index of the replicate. If `data` is not `NULL`, it takes priority over any eventual data in `Spoints`.
  #' @param edge_number Column (or entry on the list) of the `data` that contains the edge numbers. If not supplied,
  #' the column with name "edge_number" will be chosen.   Will not be used if `Spoints` is not `NULL`.
  #' @param distance_on_edge Column (or entry on the list) of the `data` that contains the edge numbers. If not supplied,
  #' the column with name "distance_on_edge" will be chosen.  Will not be used if `Spoints` is not `NULL`.
  #' @param coord_x Column (or entry on the list) of the `data` that contains the x coordinate. If not supplied,
  #' the column with name "coord_x" will be chosen.  Will not be used if `Spoints` is not `NULL` or if `data_coords` is `PtE`.
   #' @param coord_y Column (or entry on the list) of the `data` that contains the y coordinate. If not supplied,
  #' the column with name "coord_x" will be chosen.  Will not be used if `Spoints` is not `NULL` or if `data_coords` is `PtE`.
  #' @param data_coords To be used only if `Spoints` is `NULL`. Which coordinate system to use? If `PtE`, the user must provide
  #' `edge_number` and `distance_on_edge`, otherwise if `euclidean`, the user must provide `coord_x` and `coord_y`.
  #' @param replicates If the data contains replicates, one must provide the column (or entry on the list) in which the replicate indices are stored.
  #' @param normalized if TRUE, then the distances in `distance_on_edge` are assumed to be
  #' normalized to (0,1). Default FALSE. Will not be used if `Spoints` is not `NULL`.
  add_observations = function(Spoints = NULL, data = NULL, edge_number = "edge_number",
                                          distance_on_edge = "distance_on_edge",
                                          coord_x = "coord_x",
                                          coord_y = "coord_y",
                                          data_coords = c("PtE", "euclidean"),
                                          replicates = NULL, normalized = FALSE) {
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

      lapply(data, function(dat){if(nrow(matrix(Spoints@coords, ncol=2)) != length(dat)){
        stop(paste(dat,"has a different number of elements than the number of coordinates!"))
      }})

      if(!is.null(Spoints)){
        PtE <- self$coordinates(Spoints@coords)
      } else{
        if(data_coords == "PtE"){
          if(normalized){
            PtE <- cbind(data[[edge_number]], data[[distance_on_edge]])
          } else{
            PtE <- cbind(data[[edge_number]], data[[distance_on_edge]] / self$edge_lengths[PtE[, 1]])
          }
        } else if(data_coords == "euclidean"){
            point_coords <- cbind(data[[coord_x]], data[[coord_y]])
            PtE <- self$coordinates(point_coords)
        } else{
          stop("The options for 'data_coords' are 'PtE' and 'euclidean'.")
        }
      }
     if(!is.null(replicates)){
      replicate_vector <- data[[replicates]]
     } else{
      replicate_vector <- NULL
     }

     data[[edge_number]] <- NULL
     data[[distance_on_edge]] <- NULL
     data[[coord_x]] <- NULL
     data[[coord_y]] <- NULL 
     data[[replicates]] <- NULL

      
      ## convert everything to PtE

      self$data <- process_data_add_obs(PtE, new_data = data, self$data, replicate_vector)

      ## convert to Spoints and add

      PtE <- get_PtE()

      Spoints <- self$coordinates(PtE = PtE)

      self$data[["__coord_x"]] <- Spoints@coords[,1]
      self$data[["__coord_y"]] <- Spoints@coords[,2]
  
    },


  #' @description build Kirchoff constraint matrix from edges, currently not
  #' implemented for circles (edges that start and end in the same vertex)
  #' @param alpha the type of constraint (currently only supports 2)
  #' @param edge_constraint if TRUE, add constraints on vertices of degree 1
  buildC = function(alpha = 2, edge_constraint = FALSE) {

    if(alpha==2){
      i_  =  rep(0, 2*self$nE)
      j_  =  rep(0, 2*self$nE)
      x_  =  rep(0, 2*self$nE)

      count_constraint <- 0
      count <- 0
      for (v in 1:self$nV) {
        lower.edges  <- which(self$E[, 1] %in% v)
        upper.edges  <- which(self$E[, 2] %in% v)
        n_e <- length(lower.edges) + length(upper.edges)

        #derivative constraint
        if ((edge_constraint & n_e ==1) | n_e > 1) {
          i_[count + 1:n_e] <- count_constraint + 1
          j_[count + 1:n_e] <- c(4 * (lower.edges-1) + 2, 4 * (upper.edges-1) + 4)
          x_[count + 1:n_e] <- c(rep(1,length(lower.edges)),
                                 rep(-1,length(upper.edges)))
          count <- count + n_e
          count_constraint <- count_constraint + 1
        }
        if (n_e > 1) {
          if (length(upper.edges) == 0) {
            edges <- cbind(lower.edges, 1)
          } else if(length(lower.edges) == 0){
            edges <- cbind(upper.edges, 3)
          }else{
            edges <- rbind(cbind(lower.edges, 1),
                           cbind(upper.edges, 3))
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

  #' @description build mesh object for graph
  #' @param h maximum distance between mesh nodes (should be provided if n is not provided)
  #' @param n maximum number of nodes per edge (should be provided if h is not provided)
  #' @details The mesh is a list with the objects
  #' - PtE which contains the mesh locations excluding the original vertices
  #' - V the verties of the mesh
  #' - E the edges of the mesh
  #' - n_e the number of vertices in the mesh per original edge in the graph
  #' - h_e the mesh width per edge in the graph
  #' - ind the indices of the vertices in the mesh
  #' - VtE all mesh locations including the original vertices
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
    for (i in 1:dim(self$LtE)[1]) {
      if (is.null(n)) {
        #remove boundary points
        self$mesh$n_e[i] <- ceiling(self$edge_lengths[i] / h) + 1 - 2
      } else {
        self$mesh$n_e[i] <- n
      }
      if (self$mesh$n_e[i] > 0) {
        d.e <- seq(from = 0, to = 1, length.out = self$mesh$n_e[i] + 2)
        d.e <- d.e[2:(1+self$mesh$n_e[i])]

        self$mesh$h_e <- c(self$mesh$h_e,
                           rep(self$edge_lengths[i] * d.e[1],
                               self$mesh$n_e[i] + 1))

          Line <- self$lines[i,]
          val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized=TRUE)
          Points <- gInterpolate(Line, d.e, normalized=TRUE)
          self$mesh$V <- rbind(self$mesh$V, Points@coords)

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

    ### Update mesh PtE 
    self$mesh$PtE <- self$coordinates(XY = Points@coords)
  },

  #' @description build mass and stiffness matrices for given mesh object
  compute_fem = function() {
    if (is.null(self$mesh)) {
      stop("no mesh provided")
    }
    nV <- dim(self$mesh$V)[1]
    self$mesh$G <- self$mesh$C <- Matrix(0,nrow=nV,ncol=nV)

    for (e in 1:dim(self$mesh$E)[1]) {
      v1 <- self$mesh$E[e, 1]
      v2 <- self$mesh$E[e, 2]

      self$mesh$C[v1,v1] <- self$mesh$C[v1, v1] + self$mesh$h_e[e]/3
      self$mesh$C[v2,v2] <- self$mesh$C[v2, v2] + self$mesh$h_e[e]/3
      self$mesh$C[v1,v2] <- self$mesh$C[v1, v2] + self$mesh$h_e[e]/6
      self$mesh$C[v2,v1] <- self$mesh$C[v2, v1] + self$mesh$h_e[e]/6

      self$mesh$G[v1,v1] <- self$mesh$G[v1, v1] + 1 / self$mesh$h_e[e]
      self$mesh$G[v2,v2] <- self$mesh$G[v2, v2] + 1 / self$mesh$h_e[e]
      self$mesh$G[v1,v2] <- self$mesh$G[v1, v2] - 1 / self$mesh$h_e[e]
      self$mesh$G[v2,v1] <- self$mesh$G[v2, v1] - 1 / self$mesh$h_e[e]
    }
  },

  #' @description Computes observation matrix for mesh
  #' @param PtE locations given as (edge number in graph, normalized location on edge)
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

  #' @description Find one edge corresponding to each vertex
  #' @return VtE matrix where `VtE[i,1]` is the edge number, `VtE[i,2] = 0`
  #' if the vertex is at the start of the edge and `VtE[i,1] = 1` if the vertex
  #' is at the end of the edge
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

  #' @description plot a metric graph
  #' @param data Which column of the data to plot? If `NULL`, no data will be plotted.
  #' @param repl If there are replicates, which replicate to plot? 
  #' @param plotly use plot_ly for 3D plot (default FALSE). This option requires the 'plotly' package.
  #' @param vertex_size size of the vertices
  #' @param vertex_color color of vertices
  #' @param edge_width line width for edges
  #' @param edge_color color of edges
  #' @param data_size size of markers for data
  #' @param mesh Plot the mesh locations?
  #' @param X Additional values to plot
  #' @param X_loc locations of the additional values in the format
  #' (edge, normalized distance on edge)
  #' @param p existing ggplot or plot_ly object to add the graph to
  #' @param ... additional arguments for ggplot or plot_ly
  #' @return a plot_ly or or ggplot object
  #' @examples
  #' library(sp)
  #' line1 <- Line(rbind(c(0,0),c(1,0)))
  #' line2 <- Line(rbind(c(0,0),c(0,1)))
  #' line3 <- Line(rbind(c(0,1),c(-1,1)))
  #' theta <- seq(from=pi,to=3*pi/2,length.out = 20)
  #' line4 <- Line(cbind(sin(theta),1+ cos(theta)))
  #' Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
  #'                               Lines(list(line2),ID="2"),
  #'                               Lines(list(line3),ID="3"),
  #'                               Lines(list(line4),ID="4")))
  #' graph <- metric_graph$new(lines = Lines)
  #' graph$plot()
  plot = function(data = NULL,
                  repl = 1,
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
                  ...) {
    if(!plotly){
      p <- private$plot_2d(line_width = edge_width,
                           marker_size = vertex_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           data_size = data_size,
                           repl = repl,
                           mesh = mesh,
                           X = X,
                           X_loc = X_loc,
                           p = p,
                           ...)
    } else {
      requireNamespace("plotly")
      p <- private$plot_3d(line_width = edge_width,
                           marker_size = vertex_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           data_size = data_size,
                           repl = repl,
                           mesh = mesh,
                           X = X,
                           X_loc = X_loc,
                           p = p,
                           ...)
    }
    return(p)
  },

  #' @description plot continuous function on the graph
  #' @param X Either an m x 3 matrix with (edge number, position on
  #' curve (in length), value) or a vector with values for the function
  #' evaluated at the mesh in the graph
  #' @param plotly if TRUE, then plot is shown in 3D. This option requires the package 'plotly'.
  #' @param vertex_size (for both 2d and 3d plots) size of the vertices
  #' @param vertex_color color of vertices
  #' @param edge_width width for edges
  #' @param edge_color for 3D plot, color of edges
  #' @param line_width for 3D plot, line width of the function curve.
  #' @param line_color color of the function curve
  #' @param support_width for 3D plot, width of support lines
  #' @param support_color for 3D plot, color of support lines
  #' @param p previous plot in which the new plot should be added.
  #' @param ... additional arguments for ggplot or plot_ly
  #' @return either a ggplot or a plot_ly object
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

      if (length(X) != dim(self$V)[1] + PtE_dim) {
        stop("X does not have the correct size")
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

        index <- (self$LtE@p[i] + 1) :
          (self$LtE@p[i + 1])
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
          p2 <- rgeos::gInterpolate(Line_edge, data.to.plot.order[, 1,
                                                                  drop = FALSE],
                                    normalized = TRUE)
          coords <-p2@coords
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
                                i = rep(1:length(z.loc),2))
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
          geom_path(linewidth = line_width) + scale_color_viridis() + labs(colour = "")
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

  #' @description function for splitting lines in the graph
  #' @param Ei index of line to split
  #' @param t  position on line to split (normalized)
  split_edge = function(Ei, t) {

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
      val_line <- gInterpolate(Line, t_mod, normalized = TRUE)@coords

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

    newV <- self$nV + 1
    self$V <- rbind(self$V, c(val_line))
    l_e <- self$edge_lengths[Ei]
    self$edge_lengths[Ei] <- t * l_e
    self$edge_lengths <- c(self$edge_lengths, (1 - t) * l_e)
    self$nE <- self$nE + 1
    self$E <- rbind(self$E, c(newV, self$E[Ei, 2]))
    self$E[Ei, 2] <- newV

    ind <- which(self$data[["__edge_number"]] %in% Ei)
    for (i in ind) {
      if (self$data[["__distance_on_edge"]][i] >= t - 1e-10) {
        self$data[["__edge_number"]][i] <- self$nE
        self$data[["__distance_on_edge"]][i] <- 
        abs(self$data[["__distance_on_edge"]][i] - t) / (1 - t)
      }
    }

    self$nV <- dim(self$V)[1]
  },

  #' @description Add observations on mesh to the object
  #' @param data A data.frame or named list containing the observations. In case of replicates, the data.frames for the replicates should stacked vertically, with a column
  #' indicating the index of the replicate. If `data_frame` is not `NULL`, it takes priority over any eventual data in `Spoints`.
  #' @param replicates If the data_frame contains replicates, one must provide the column in which the replicate indices are stored.
  add_mesh_observations = function(data = NULL, replicates = NULL) {
    if(is.null(self$mesh)){
      stop("You should have a mesh!")
    }
    Spoints <- self$mesh$V[(nrow(self$VtEfirst()) + 1):nrow(self$mesh$V), ]
    rownames(Spoints) <- 1:nrow(Spoints)
    Spoints <- SpatialPoints(coords = Spoints)
    self$add_observations(Spoints = Spoints, data = data, replicates = replicates)
  },

  #' @description Get a copy of the initial graph
  get_initial_graph = function() {
    return(private$initial_graph$clone())
  },

  #' @description Get the observation/prediction matrix A
  #' @param order Which order should be considered? The options are 'internal' and 'original'.
  #'  The order 'internal' is the order of `graph$y``, for a metric_graph `graph`. The order 'original'
  #' is the order of the user's input.
  #' @param obs_to_vert Should the observations be turned into vertices?

  A = function(order = "internal", obs_to_vert = FALSE){
    if(is.null(private$internal_A) && !obs_to_vert){
        stop("The A matrix was not computed. If you want to compute rerun this method with 'obs_to_vertex=TRUE', in which the observations will be turned to vertices and the A matrix will then be computed")
    } else if(is.null(private$internal_A)){
      self$observation_to_vertex()
    }

    if(order == "internal"){
      return(private$internal_A)
    } else if(order == "original"){
      orig_A <- private$internal_A
      for (i in length(private$reorder_idx):1) {
      idx <- private$reorder_idx[[i]]
      A_tmp <- orig_A[1:length(idx), ]
      if (length(idx)==1) {
        A_tmp <- matrix(A_tmp, ncol = 2)
      }
      orig_A[idx, ] <- A_tmp
    }
    return(orig_A)
    } else{
      stop("The order must be either 'internal' or 'original'!")
    }
  },

  #' @description Convert between locations on the graph and Euclidean coordinates
  #' @param PtE matrix with locations on the graph (edge number and normalized
  #' position on the edge).
  #' @param XY matrix with locations in Euclidean space
  #' @param normalized If TRUE, it is assumed that the positions in `PtE` are
  #' normalized to (0,1), and the object returned if `XY` is specified contains
  #' normalized locations.
  #' @return If `PtE` is specified, then a matrix with Euclidean coordinates of
  #' the locations is returned. If `XY` is provided, then a matrix with the
  #' closest locations on the graph is returned
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
      for (i in 1:dim(PtE)[1]) {
        LT <- private$edge_pos_to_line_pos(PtE[i, 1], PtE[i, 2])
        Line <- self$lines[LT[1, 1], ]
        val_line <- gProject(Line, as(Line, "SpatialPoints"),
                             normalized = TRUE)
        Point <- gInterpolate(Line,LT[1, 2], normalized = TRUE)
        x <- c(x, Point@coords[1])
        y <- c(y, Point@coords[2])
      }
      return(cbind(x,y))
    } else {
      Spoints <- SpatialPoints(XY)
      SP <- snapPointsToLines(Spoints, self$lines)
      coords.old <- as.data.frame(Spoints@coords)
      colnames(coords.old) <- paste(colnames(coords.old) ,'_old',sep="")
      Spoints@coords = SP@coords
      Spoints@bbox   = SP@bbox
      LtE = cbind(match(SP@data[,1], self$EID),0)

      for (ind in unique(LtE[, 1])) {
        index.p <- LtE[, 1] == ind
        LtE[index.p,2]=rgeos::gProject(self$lines[ind,], Spoints[index.p,],
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
    #computes which line and which position t_E on Ei belongs to
    # Ei  (int)   edge index
    # t_e (n x 1) number of positions on Ei
    edge_pos_to_line_pos = function(Ei, t_E) {

      LT <- matrix(0, nrow= length(t_E),2)
      L_index <- (self$LtE@p[Ei]+1):(self$LtE@p[Ei+1])

      #LinPos line number and end relative end of the line on the edge
      LinesPos <- cbind(self$LtE@i[L_index] + 1, self$LtE@x[L_index])
      LinesPos <- LinesPos[order(LinesPos[,2]),,drop=F]
      for(j in 1:length(L_index)){
        if(j==1){
          index_j <-  t_E <= LinesPos[j,2]
        }else{
          index_j <- (t_E <= LinesPos[j,2]) &  (t_E > LinesPos[j-1,2])
        }
        if(sum(index_j) == 0)
          next

        LT[index_j,1] = LinesPos[j,1]
        rel.pos = t_E[index_j]
        if(j == 1){
          rel.pos <- rel.pos/LinesPos[j,2]
        }else{
          rel.pos <- (rel.pos-LinesPos[j-1,2])/(LinesPos[j,2]-LinesPos[j-1,2])
        }

        if(j== dim(LinesPos)[1] )
          rel.pos = self$ELend[Ei]*rel.pos
        if(j==1)
          rel.pos = rel.pos + self$ELstart[Ei]

        LT[index_j,2] = rel.pos

      }
      return(LT)
    },

  #function for creating Vertex and Edges from self$lines
  line_to_vertex = function(tolerance = 0, longlat = FALSE) {
    lines <- c()
    for(i in 1:length(self$lines)){
      tmp <- self$lines@lines[[i]]@Lines[[1]]@coords
      self$lines@lines[[i]]@Lines[[1]]@coords <- tmp
      points <- self$lines@lines[[i]]@Lines[[1]]@coords
      n <- dim(points)[1]
      #lines contain [line index, start point, line length
      #               line index, end point, line length]
      lines <- rbind(lines, c(i, points[1,]), c(i, points[n,]))
    }

    #save all vertices that are more than tolerance distance apart, also
    #add vertices that are closer than that if they are on the same lines
    dists <- spDists(lines[, 2:3, drop = FALSE], longlat = longlat)
    vertex <- lines[1, , drop = FALSE]
    for (i in 2:dim(lines)[1]) {
      i.min <- which.min(dists[i, 1:(i-1)])
      if (dists[i, i.min] > tolerance) {
        vertex <- rbind(vertex, lines[i, ])
      } else if (lines[i, 1] == lines[i.min, 1]) {
        vertex <- rbind(vertex, lines[i, ])
      }
    }
    #lvl = c(line index, vertex number of start, vertex number of end, length])
    lvl <- matrix(0, nrow = max(lines[,1]), 4)
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
      lvl[i,] <- c(i, ind1, ind2, ll)
    }
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
  },

  #Compute PtE for mesh given PtE for graph
  PtE_to_mesh = function(PtE){
    VtE <- rbind(self$VtEfirst(),self$mesh$PtE)
    PtE_update <- matrix(0,dim(PtE)[1],2)
    for(i in 1:dim(PtE)[1]) {
      ei <- PtE[i,1]
      ind <- which(VtE[,1]==ei)
      if(PtE[i,2]<min(VtE[ind,2])){ #node in first edge on line
        v1 <- self$E[ei,1]
        ind2 <- which.min(VtE[ind,2])
        v2 <- ind[ind2]
        d1 <- 0
        d2 <- VtE[v2, 2]
      } else if (PtE[i,2] > max(VtE[ind,2])) { #node in last edge on line
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
      #edge in mesh
      e <- which(rowSums((self$mesh$E==v1) + (self$mesh$E==v2))==2)
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
                     repl = 1,
                     mesh = FALSE,
                     X = NULL,
                     X_loc = NULL,
                     p = NULL,
                     ...){
    xyl <- c()

      coords <- lapply(coordinates(self$lines), function(x) x[[1]])
      nc <- do.call(rbind,lapply(coords, function(x) dim(x)[1]))
      xyl <- cbind(do.call(rbind,coords), rep(1:length(nc), times = nc))

    if(is.null(p)){
      p <- ggplot()+ geom_path(data = data.frame(x = xyl[, 1],
                                                 y = xyl[,2],
                                                 group = xyl[,3]),
                               mapping = aes(x = x, y = y, group = group),
                               linewidth = line_width,
                               colour = edge_color, ...)
    } else {
      p <- p + geom_path(data = data.frame(x = xyl[, 1],
                                                 y = xyl[,2],
                                                 group = xyl[,3]),
                         mapping = aes(x = x, y = y, group = group),
                         linewidth = line_width,
                         colour = edge_color, ...)
    }
    if (marker_size > 0) {
      p <- p + geom_point(data = data.frame(x = self$V[, 1],
                                            y = self$V[, 2]),
                          mapping = aes(x, y),
                          colour = vertex_color,
                          size= marker_size, ...)
    }
    if (!is.null(data)) {
      x <- y <- NULL
      y_plot <- self$data[[repl]][, data]
      for (i in 1:length(y_plot)) {

          LT <- private$edge_pos_to_line_pos(self$PtE[i, 1], self$PtE[i, 2])
          Line <- self$lines[LT[1, 1], ]
          val_line <- gProject(Line, as(Line, "SpatialPoints"),
                               normalized = TRUE)
          Point <- gInterpolate(Line,LT[1, 2], normalized = TRUE)
          x <- c(x, Point@coords[1])
          y <- c(y, Point@coords[2])

      }
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
      for (i in 1:length(as.vector(X))) {
        LT <- private$edge_pos_to_line_pos(X_loc[i, 1], X_loc[i, 2])
        Line <- self$lines[LT[1, 1], ]
        val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized = TRUE)
        Point <- gInterpolate(Line,LT[1, 2], normalized = TRUE)
        x <- c(x, Point@coords[1])
        y <- c(y, Point@coords[2])
      }
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
                     repl = 1,
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
      p <- plotly::plot_ly(data=data.plot, x = ~y, y = ~x, z = ~z)
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = edge_color, ...),
                           split = ~i, showlegend = FALSE)
    } else {
      p <- plotly::add_trace(p, data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = edge_color, ...),
                           split = ~i, showlegend = FALSE)
    }


    if(marker_size > 0) {
      data.plot2 <- data.frame(x = self$V[, 1], y = self$V[, 2],
                               z = rep(0, self$nV))
      p <- plotly::add_trace(p, data = data.plot2, x = ~y, y = ~x, z = ~z,
                           type = "scatter3d", mode = "markers",
                           marker = list(size = marker_size,
                                         color = vertex_color, ...))
    }

    if (!is.null(data)) {
      x <- y <- NULL
      y_plot <- self$data[[repl]][, data]
      for (i in 1:nrow(y_plot)) {
        Line <- self$lines[PtE[i, 1], ]
        val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized = TRUE)
        Point <- gInterpolate(Line, PtE[i, 2], normalized = TRUE)
        x <- c(x, Point@coords[1])
        y <- c(y, Point@coords[2])
      }
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
    xr <- 2*(diff(range(self$V[,1])) + diff(range(self$V[,2])))
    return(p)
  },

  # Initial graph

  initial_graph = NULL,


))




