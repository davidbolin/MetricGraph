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

  #' @field A sparse matrix specifying which vertices are observation locations
  A = NULL,

  #' @field CoB change-of-basis object used for Kirchhoff constraints
  CoB = NULL,

  #' @field points the observations in a SpatialPointsDataFrame
  points  = NULL,

  #' @field y vector with data on the graph
  y = NULL,

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
  #' @param tolerance vertices that are closer than this number in Euclidean
  #' distance are merged when constructing the graph (default = 2e-16).
  #' @details A graph object can be initialized in two ways. The first method
  #' is to specify V and E. In this case, all edges are assumed to be straight
  #' lines. The second option is to specify the graph via the `lines` input.
  #' In this case, the vertices are set by the end points of the lines.
  #' Thus, if two lines are intersecting somewhere else, this will not be
  #' viewed as a vertex.
  #' @return A metric_graph object
  initialize = function(lines = NULL, V = NULL, E = NULL,
                        tolerance = 2e-16) {

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
    private$line_to_vertex(tolerance = tolerance)
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
  compute_geodist = function() {
    g <- graph(edges = c(t(self$E)), directed = FALSE)
    E(g)$weight <- self$edge_lengths
    self$geo_dist <- distances(g)
  },

  #' @description Computes shortest path distances between the vertices in the
  #' mesh
  compute_geodist_mesh = function() {
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
      graph.temp$add_PtE_observations(y = rep(0, dim(PtE)[1]), PtE,
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

  #' @description Adds observation locations as vertices in the graph
  observation_to_vertex = function() {
    # Reordering
    order_idx <- order(self$PtE[, 1], self$PtE[, 2])

    private$reorder_idx <- c(private$reorder_idx, list(order_idx))

    self$PtE <- self$PtE[order_idx, ]

    if (length(order_idx) == 1) {
      self$PtE <- matrix(self$PtE, ncol = 2)
    }
    l <- length(self$PtE[, 1])
    self$PtV <- rep(0, l)
    for (i in 1:l) {
        e <- as.vector(self$PtE[i, 1])
        t <- as.vector(self$PtE[i, 2])
        l_e <- self$edge_lengths[e]
        if (abs(t) < 10^-10) {
          self$PtE[i, 2] <- 0
          self$PtV[i] <- self$E[e, 1]
        } else if (t > 1 - 10^-10) {
          self$PtE[i, 2] <- 1
          self$PtV[i] <- self$E[e, 2]
        } else {
          self$split_edge(e, t)
          self$PtV[i] <- dim(self$V)[1]
        }
    }

    if (!is.null(self$geo_dist)) {
      self$compute_geodist()
    }
    if (!is.null(self$res_dist)) {
      self$compute_resdist()
    }
    if (!is.null(self$CoB)) {
      self$buildC(2)
    }
    self$A <- Matrix::Diagonal(self$nV)[self$PtV, ]

    if (length(self$PtV) == 1) {
      self$A <- matrix(self$A, ncol = 2)
    }

    for (i in length(private$reorder_idx):1) {
      idx <- private$reorder_idx[[i]]
      A_tmp <- self$A[1:length(idx), ]
      if (length(idx)==1) {
        A_tmp <- matrix(A_tmp, ncol = 2)
      }
      self$A[idx, ] <- A_tmp
    }

    self$add_responses(private$raw_y)
  },

  #' @description Clear all observations from the object
  clear_observations = function() {
   self$y <- NULL
   self$PtE <- NULL
   self$points <- NULL
   private$raw_y <- c()
   private$reorder_idx <- list()
  },

  #' @description Add observations to the graph
  #' @param Spoints SpatialPoints or SpatialPointsDataFrame of the observations,
  #' which may include the coordinates only, or the coordinates as well as the
  #' observations
  #' @param y the observations. These are used if provided, and otherwise the
  #' observations are assumed to be in Spoints
  #' @param y.index If `y` is not provided, `y.index` gives the column number
  #' for the data to use in `Spoints@data`. If it is not provided, it is assumed
  #' that the data is in the first column
  add_observations = function(Spoints, y = NULL, y.index = NULL) {

    if (!is.null(y) && is.null(y)) {
      stop("if y is provided, then y.index must be provided as well.")
    }

    if("SpatialPointsDataFrame"%in%is(Spoints)){
    if(is.null(y)){
        if(is.null(y.index)) {
          if (dim(Spoints@data)[2] > 1){
            stop("The data field contains multiple columns,
                 please specify which column to use via y.index")
          } else {
            y.index <- 1
          }
        }
        y <- Spoints@data[,y.index]
    }
    }
    if(is.null(y)){
      y <- rep(NA, nrow(Spoints@coords))
    }
    self$y <- c(self$y, y)
    private$raw_y <- c(private$raw_y, y)

    SP <- snapPointsToLines(Spoints, self$lines)
    coords.old <- as.data.frame(Spoints@coords)
    colnames(coords.old) <- paste(colnames(coords.old) ,'_old',sep="")
    Spoints@coords = SP@coords
    Spoints@bbox   = SP@bbox
    LtE = cbind(match(SP@data[,1], self$EID),0)
    if("SpatialPointsDataFrame"%in%is(Spoints)){
      Spoints@data <- cbind(Spoints@data,coords.old)
    }else{
      Spoints <- SpatialPointsDataFrame(Spoints, data = coords.old)
    }
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
    if (is.null(self$points)) {
      self$points <- Spoints
      self$PtE <- PtE
    } else {
      df1 <- self$points@data
      df2 <- Spoints@data
      df1[setdiff(names(df2), names(df1))] <- NA
      df2[setdiff(names(df1), names(df2))] <- NA
      self$points@data <- df1
      Spoints@data <- df2
      self$points <- rbind(self$points, Spoints)
      self$points <- rbind(self$points, Spoints)
      self$PtE <- rbind(self$PtE, PtE)
    }
  },

  #' @description Add observations to the object
  #' @param y vector with the values of the observations
  #' @param PtE matrix where `PtE[i,1]` is the index of the edge for the ith
  #' observation and `PtE[i,2]` is the distance on the edge where the
  #' observation is located
  #' @param normalized if TRUE, then the distances in `PtE` are assumed to be
  #' normalized to (0,1). Default FALSE.
  #' @param Spoints Optional argument of class `SpatialPoints` or
  #' `SpatialPointsDataFrame` specifying the Euclidean coordinates of the
  #' observation locations. If this is not provided, the coordinates are
  #' calculated internally.
  add_PtE_observations = function(y, PtE, Spoints=NULL, normalized = FALSE) {

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
    self$y <- c(self$y, y)
    private$raw_y <- c(private$raw_y, y)

    if(normalized){
      self$PtE = rbind(self$PtE, PtE)
    } else {
      PtE <- cbind(PtE[, 1], PtE[, 2] / self$edge_lengths[PtE[, 1]])
      self$PtE = rbind(self$PtE, PtE)
    }


      if(!is.null(Spoints)){
              if("SpatialPointsDataFrame"%in%is(Spoints)){
        Spoints@data <- cbind(Spoints@data,as.data.frame(PtE))
      }else{
        Spoints <- SpatialPointsDataFrame(Spoints, data = as.data.frame(PtE))
      }
      } else {
        coords <- c()
        for(i in 1:dim(PtE)[1]){

          LT = private$edge_pos_to_line_pos(PtE[i, 1] , PtE[i, 2])
          points <- rgeos::gInterpolate(self$lines[LT[1,1],],
                                        LT[1,2],
                                        normalized = TRUE)
          coords <- rbind(coords, points@coords)
        }
        rownames(coords) <- 1:dim(coords)[1]
        Spoints <- sp::SpatialPoints(coords)
        Spoints <- SpatialPointsDataFrame(Spoints, data = as.data.frame(PtE))
      }
      if(is.null(self$points)){
        self$points <- Spoints
      } else {
        df1 <- self$points@data
        df2 <- Spoints@data
        df1[setdiff(names(df2), names(df1))] <- NA
        df2[setdiff(names(df1), names(df2))] <- NA
        self$points@data <- df1
        Spoints@data <- df2
        self$points = rbind(self$points, Spoints)
      }
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
        self$mesh$PtE <- rbind(self$mesh$PtE, cbind(rep(i, self$mesh$n_e[i]),
                                                    d.e))

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
  #' @param plotly use plot_ly for 3D plot (default FALSE)
  #' @param line_width line width for edges
  #' @param vertex_size size of the vertices
  #' @param vertex_color color of vertices
  #' @param edge_color color of edges
  #' @param data Plot the data?
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
  plot = function(plotly = FALSE,
                  line_width = 0.3,
                  vertex_size = 3,
                  vertex_color = 'black',
                  edge_color = 'black',
                  data = FALSE,
                  data_size = 1,
                  mesh = FALSE,
                  X = NULL,
                  X_loc = NULL,
                  p = NULL,
                  ...) {
    if(!plotly){
      p <- private$plot_2d(line_width = line_width,
                           marker_size = vertex_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           data_size = data_size,
                           mesh = mesh,
                           X = X,
                           X_loc = X_loc,
                           p = p,
                           ...)
    } else {
      p <- private$plot_3d(line_width = line_width,
                           marker_size = vertex_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           data_size = data_size,
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
  #' @param plotly if TRUE, then plot is shown in 3D
  #' @param graph_color for 3D plot, the color of the graph.
  #' @param line_width for 3D plot, the line width of the graph.
  #' @param vertex_size for 3D plot, the vertex size of the vertices
  #' @param color color of curve
  #' @param p previous plot in which the new plot should be added.
  #' @param ... additional arguments for ggplot or plot_ly
  #' @return either a ggplot or a plot_ly object
  plot_function = function(X,
                           plotly = FALSE,
                           graph_color = 'black',
                           line_width = 1,
                           vertex_size = 10,
                           color = 'rgb(0,0,200)',
                           p = NULL,
                           ...){

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
        X <- c(rep(NA, dim(private$initial_graph$V)[1]), X)
      }

      if (length(X) != dim(private$initial_graph$V)[1] + PtE_dim) {
        stop("X does not have the correct size")
      }
    }

    if (mesh) {
      n.v <- dim(private$initial_graph$V)[1]
      XV <- X[1:n.v]
    }

    x.loc <- y.loc <- z.loc <- i.loc <- NULL
    kk = 1
    for (i in 1:private$initial_graph$nE) {
      Vs <- private$initial_graph$E[i, 1]
      Ve <- private$initial_graph$E[i, 2]
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
      if (is.null(private$initial_graph$lines) == TRUE) {
        data.to.plot.order <- vals[order(vals[, 1]), ]
        V <- private$initial_graph$V[private$initial_graph$E[i, ], ]

        alpha <- data.to.plot.order[,1]
        coords <- cbind((1 - alpha) * V[1, 1] + alpha * V[2, 1],
                        (1 - alpha) * V[1, 2] + alpha * V[2, 2])

        x.loc = c(x.loc, coords[, 1])
        y.loc = c(y.loc, coords[, 2])
        z.loc = c(z.loc, data.to.plot.order[, 2])
        i.loc = c(i.loc, rep(kk, length(coords[, 1])))
        kk = kk+1
      } else {
        index <- (private$initial_graph$LtE@p[i] + 1) :
          (private$initial_graph$LtE@p[i + 1])
        LinesPos <- cbind(private$initial_graph$LtE@i[index] + 1,
                          private$initial_graph$LtE@x[index])
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
            rel.pos = private$initial_graph$ELend[i] * rel.pos
          if (j==1)
            rel.pos = rel.pos + private$initial_graph$ELstart[i]

          data.to.plot <- cbind(rel.pos,vals[index_j, 2])
          Line_edge <- SpatialLines(list(private$initial_graph$lines@lines[[LinesPos[j, 1]]]))

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
    }
    data <- data.frame(x = x.loc, y = y.loc, z = z.loc, i = i.loc)

    if(plotly){
      if(is.null(p)){
        p <- self$plot(plotly = TRUE, color = graph_color,
                       line_width = line_width, vertex_size = vertex_size)
      }
      p <- p %>% add_trace(data = data, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = color, ...),
                           split = ~i, showlegend = FALSE)
    } else {
      if(is.null(p)) {
        p <- ggplot(data = data, aes(x = x, y = y,
                                     group = i,
                                     colour = z),
                    linewidth = line_width) +
          geom_path() + scale_color_viridis() + labs(colour = "")
      } else {
        p <- p + geom_path(data = data,
                           aes(x = x, y = y,
                               group = i, colour = z),
                           linewidth = line_width) +
          scale_color_viridis() + labs(colour = "")
      }

    }

    return(p)
  },

  #' @description plot continuous function on the mesh of a graph
  #' @param X Either an m x 3 matrix with (edge number, position on
  #' curve (in length), value) or a vector with values for the function
  #' evaluated at a precomputed mesh.
  #' @param plotly plot in 2D or 3D?
  #' @param graph_color for 3D plot, the color of the graph.
  #' @param line_width for 3D plot, the line width of the curves.
  #' @param vertex_size for 3D plot, the size of the vertices
  #' @param color Color of the curves
  #' @param p previous ggplot or plot_ly object to add the plot to
  #' @param ... additional arguments for ggplot or plot_ly
  #' @return A ggplot or a plot_ly object
  plot_function_mesh = function(X,
                                plotly = FALSE,
                                graph_color = 'black',
                                line_width = 1,
                                vertex_size = 10,
                                color = 'rgb(0,0,200)',
                                p = NULL,
                                ...){
    self$plot_function(X = X, plotly = plotly, graph_color = graph_color,
                       line_width = line_width, vertex_size = vertex_size,
                       color = color, p = p, ...)
  },

  #' @description function for splitting lines in the graph
  #' @param Ei index of line to split
  #' @param t  position on line to split (normalized)
  split_edge = function(Ei, t) {
    if (!is.null(self$lines)) {
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
    } else {
      V1 <- self$V[self$E[Ei, 1], ]
      V2 <- self$V[self$E[Ei, 2], ]
      val_line <- (1 - t) * V1 + t * V2
    }
    newV <- self$nV + 1
    self$V <- rbind(self$V, c(val_line))
    l_e <- self$edge_lengths[Ei]
    self$edge_lengths[Ei] <- t * l_e
    self$edge_lengths <- c(self$edge_lengths, (1 - t) * l_e)
    self$nE <- self$nE + 1
    self$E <- rbind(self$E, c(newV, self$E[Ei, 2]))
    self$E[Ei, 2] <- newV
    ind <- which(self$PtE[, 1] %in% Ei)
    for (i in ind) {
      if (self$PtE[i, 2] >= t - 1e-10) {
        self$PtE[i, 1] <- self$nE
        self$PtE[i, 2] <- abs(self$PtE[i, 2] - t) / (1 - t)
      }
    }
    self$nV <- dim(self$V)[1]
  },

  #' @description function for adding simulated response variables in the
  #' correct order.
  #' @param y A vector of response variables
  add_responses = function(y) {
  stopifnot(length(y) == nrow(self$PtE))
  private$raw_y <- y
  self$y <- y
  idx <- private$reorder_idx[[1]]
  y_tmp <- self$y[idx]
  self$y[1:length(idx)] <- y_tmp
  if (length(private$reorder_idx) > 1) {
    for (i in 2:length(private$reorder_idx)) {
      idx <- private$reorder_idx[[i]]
      y_tmp <- self$y[1:length(idx)]
      self$y <- y_tmp[idx]
    }
  }
},

  #' @description Add observations on mesh to the object
  #' @param y the observations.
  add_mesh_observations = function(y) {
    if(is.null(self$mesh)){
      stop("You should have a mesh!")
    }
    Spoints <- self$mesh$V[(nrow(self$VtEfirst()) + 1):nrow(self$mesh$V), ]
    rownames(Spoints) <- 1:nrow(Spoints)
    Spoints <- SpatialPoints(coords = Spoints)
    self$add_observations(Spoints = Spoints, y = y)
  },

  #' @description Get a copy of the initial graph
  get_initial_graph = function() {
    return(private$initial_graph$clone())
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
  line_to_vertex = function(tolerance = 0) {
    lines <- c()
    for(i in 1:length(self$lines)){
      points <- self$lines@lines[[i]]@Lines[[1]]@coords
      n <- dim(points)[1]
      #lines contain [line index, start point, line length
      #               line index, end point, line length]
      lines <- rbind(lines,
                     c(i, points[1,],
                       sp::LineLength(self$lines@lines[[i]]@Lines[[1]])),
                     c(i, points[n,],
                       sp::LineLength(self$lines@lines[[i]]@Lines[[1]])))
    }

    #save all vertices that are more than tolerance distance apart
    vertex <- lines[1, , drop = FALSE]
    for (i in 1:dim(lines)[1]) {
      tmp <- vertex[, 2:3] - matrix(rep(1,dim(vertex)[1]),
                                   dim(vertex)[1], 1) %*% lines[i,2:3]
      if (min(rowSums(tmp^2)) > tolerance) {
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
      lvl[i,] <- c(i, ind1, ind2, line[1,4])
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

  plot_2d = function(line_width = 0.1,
                     marker_size = 1,
                     vertex_color = 'black',
                     edge_color = 'black',
                     data = FALSE,
                     data_size = 1,
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
                               colour = edge_color)
    } else {
      p <- p + geom_path(data = data.frame(x = xyl[, 1],
                                                 y = xyl[,2],
                                                 group = xyl[,3]),
                         mapping = aes(x = x, y = y, group = group),
                         linewidth = line_width,
                         colour = edge_color)
    }
    if (marker_size > 0) {
      p <- p + geom_point(data = data.frame(x = self$V[, 1],
                                            y = self$V[, 2]),
                          mapping = aes(x, y),
                          colour = vertex_color,
                          size= marker_size)
    }
    if (data) {
      x <- y <- NULL
      for (i in 1:length(self$y)) {

          LT <- private$edge_pos_to_line_pos(self$PtE[i, 1], self$PtE[i, 2])
          Line <- self$lines[LT[1, 1], ]
          val_line <- gProject(Line, as(Line, "SpatialPoints"),
                               normalized = TRUE)
          Point <- gInterpolate(Line,LT[1, 2], normalized = TRUE)
          x <- c(x, Point@coords[1])
          y <- c(y, Point@coords[2])

      }
      p <- p + geom_point(data = data.frame(x = x, y = y,
                                            val = as.vector(self$y)),
                          mapping = aes(x, y, color = val),
                          size = data_size) +
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

  plot_3d = function(line_width = 1,
                     marker_size = 1,
                     vertex_color = 'rgb(0,0,0)',
                     edge_color = 'rgb(0,0,0)',
                     data = FALSE,
                     data_size = 1,
                     mesh = FALSE,
                     p = NULL,
                     ...){
    if (is.null(self$lines)) {
      data.plot <- data.frame(x = c(self$V[E[, 1], 1], self$V[E[, 2], 1]),
                              y = c(self$V[E[, 1], 2], self$V[E[, 2], 2]),
                              z = rep(0, 2 * self$nE),
                              i = c(1:self$nE, 1:self$nE))
    } else {
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
    }
    if(is.null(p)) {
      p <- plot_ly(data=data.plot, x = ~y, y = ~x, z = ~z)
      p <- p %>% add_trace(data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = edge_color, ...),
                           split = ~i, showlegend = FALSE)
    } else {
      p <- p %>% add_trace(data = data.plot, x = ~y, y = ~x, z = ~z,
                           mode = "lines", type = "scatter3d",
                           line = list(width = line_width,
                                       color = edge_color, ...),
                           split = ~i, showlegend = FALSE)
    }


    if(marker_size > 0) {
      data.plot2 <- data.frame(x = self$V[, 1], y = self$V[, 2],
                               z = rep(0, self$nV))
      p <- p %>% add_trace(data = data.plot2, x = ~y, y = ~x, z = ~z,
                           type = "scatter3d", mode = "markers",
                           marker = list(size = marker_size,
                                         color = vertex_color, ...))
    }

    if (data) {
      x <- y <- NULL
      for (i in 1:length(self$y)) {
        Line <- self$lines[PtE[i, 1], ]
        val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized = TRUE)
        Point <- gInterpolate(Line, PtE[i, 2], normalized = TRUE)
        x <- c(x, Point@coords[1])
        y <- c(y, Point@coords[2])
      }
      data.plot <- data.frame(x = x, y = y,
                              z = rep(0,length(x)),
                              val = self$y)
      p <- p %>% add_trace(data = data.plot, x = ~y, y = ~x, z = ~z,
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
      p <- p %>% add_trace(data = data.plot, x = ~y, y = ~x, z = ~z,
                           type = "scatter3d", mode = "markers",
                           marker = list(size = marker_size/2,
                                         color = 'rgb(100,100,100)'),
                           showlegend = FALSE)
    }
    xr <- 2*(diff(range(self$V[,1])) + diff(range(self$V[,2])))
    return(p)
  },

  # Ordering indexes

  reorder_idx = list(),

  # unordered y

  raw_y = c(),

  # Initial graph

  initial_graph = NULL

  # # Initial number of vertices

  # initial_V = NULL,

  # # Initial edges

  # initial_E = NULL,

  # # Initial lines

  # initial_Lines = NULL,

  # # Initial LtE

  # initial_LtE = NULL,

  # # Initial ELstart

  # initial_ELstart = NULL,

  # # Initial ELend

  # initial_ELend = NULL

))




