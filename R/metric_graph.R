#' @title Metric graph object for specification of Gaussian processes
#' @description Class representing general metric graphs.
#' @details A graph object created from vertex and edge matrices, or from an
#' sp::Lines object where each line is representing and edge.
#'
#' @export
metric_graph <-  R6::R6Class("GPGraph::graph",
  public = list(
  #' @field V Position in Euclidean space of the vertices
  V = NULL,

  #' @field nV number of vertices
  nV = 0,

  #' @field E Edges,  E[i,1] is the vertex at the start of the edge and  E[i,2]
  #' is the vertex at the end of the edge
  E = NULL,

  #' @field nE number of edges
  nE= 0,

  #' @field edge_lengths length of edges
  edge_lengths = NULL,

  #' @field EID ID of edges
  EID = NULL,


  #' @field EtL [i,] - edge i position on lines
  LtE = NULL,

  #' @field ELend normalized end point of edge on line (normalized on the line)
  ELend = NULL,
  #' @field Estart normalized start point of edge on line (normalized on the line)
  ELstart = NULL,

  #' @field C constraint matrix used to set Kirchhoff constraints
  C = NULL,

  #' @field A Matrix specifying which vertices are observation locations
  A = NULL,

  #' @field CBobj svd stuct object used for Kirchhoff constraints
  CBobj = NULL,

  #' @field Points Observations in SpatialPointsDataFrame
  Points  = NULL,

  #' @field y the data connected to P
  y = NULL,

  #' @field PtE Points to Line (connected to Points)
  #'  [,1] edge index
  #'  [,2] normalized distance along the line (i.e. distance to initial point)
  PtE = NULL,

  #' @field PtV Vector with the indices of the vertices which are observation
  #' locations
  PtV  = NULL,

  #' @field mesh mesh object used for plotting
  mesh = NULL,

  #' @field Lines List of Lines object for building the graph
  Lines = NULL,

  #' @field geo.dist Geodesic distance matrix
  geo.dist = NULL,

  #' @field res.dist Resistance distance matrix
  res.dist = NULL,

  #' @field Laplacian The weighted graph Laplacian
  Laplacian = NULL,

  #' @description Create a new gpgraph_graph object
  #' @param Lines sp object SpatialLines DataFrame or SpatialLines
  #' @param P n x 2 matrix with Euclidean coordinates of the n vertices
  #' @param E m x 2 matrix where each line represents an edge
  #' @param edge_lengths m x 1 vector with edge lengths
  #' @details A graph object can be initialized in two ways. The first method
  #' is to specify P and E. In this case, if edge_lengths is not specified, all
  #' edges are assumed to be straight lines. Otherwise the edge lengths set in
  #' edge_lengths are used. The second option is to specify the graph based on
  #' Lines. In this case, the vertices are set by the end points of the lines.
  #' Thus, if two lines are intersecting somewhere else, this will not be
  #' viewed as a vertex.
  #' @return A metric_graph object
  initialize = function(Lines = NULL, P = NULL, E = NULL, edge_lengths = NULL) {
    #We have three different ways of initializing:

    #option 1: initialization from lines
    if(!is.null(Lines)){
      if(!is.null(edge_lengths) || !is.null(P) || !is.null(E)){
        warning("object initialized from lines, then E,P,edge_lengths are
                ignored")
      }
      self$nE = length(Lines)
      self$Lines = Lines

      self$EID = sapply(slot(self$Lines,"lines"), function(x) slot(x, "ID"))
      private$line_to_vertex()
    } else {
      if(is.null(P) || is.null(E)){
        stop("You must supply Lines or P and E")
      }
      self$nE <- dim(E)[1]
      self$V   <- P
      self$E <- E
      self$nV <- dim(self$V)[1]
      if(is.null(edge_lengths)){
        self$edge_lengths <- sqrt((self$V[self$E[,2], 1] -
                                     self$V[self$E[, 1], 1])^2 +
                                    (self$V[self$E[,2], 2] -
                                       self$V[self$E[, 1], 2])^2)
      } else {
        self$edge_lengths = edge_lengths
      }
    }
  },

  #' @description Computes shortest path distances between the vertices in the
  #' graph
  compute_geodist = function(){
    g <- graph(edges = c(t(self$E)), directed=FALSE)
    E(g)$weight <- self$edge_lengths
    self$geo.dist <- distances(g)
  },

  #' @description Computes the resistance metric between the vertices in the
  #' graph
  compute_resdist = function(){
    if(is.null(self$geo.dist)){
      self$compute_geodist()
    }
    L <- Matrix(0,self$nV,self$nV)
    for(i in 1:self$nE){
      tmp <- -1/self$geo.dist[self$E[i,1], self$E[i,2]]
      L[self$E[i,2],self$E[i,1]] <- L[self$E[i,1],self$E[i,2]] <- tmp
    }
    for(i in 1:self$nV){
      L[i,i] <- -sum(L[i,-i])
    }
    L[1,1] <- L[1,1] + 1

    Li <- solve(L)
    self$res.dist <- -2*Li + t(diag(Li))%x%rep(1,self$nV) + t(rep(1,self$nV))%x%diag(Li)
  },

  #' @description Compute graph Laplacian for the graph
  compute_laplacian = function(){
    Wmat <- Matrix(0,self$nV,self$nV)
    for(i in 1:self$nE){
      Wmat[self$E[i, 1], self$E[i, 2]] <- 1/self$edge_lengths[i]
      Wmat[self$E[i, 2], self$E[i, 1]] <- 1/self$edge_lengths[i]
    }
    self$Laplacian <- Diagonal(self$nV,as.vector(Matrix::rowSums(Wmat))) - Wmat
  },

  #' @description Add observation locations as vertices in the graph
  #' @export
  observation_to_vertex = function(){

    l <- length(self$PtE[, 1])
    self$PtV <- rep(0, l)
    for(i in 1:l){
        e <- as.vector(self$PtE[i, 1])
        t <- as.vector(self$PtE[i, 2])
        l_e <- self$edge_lengths[e]
        if(abs(t) < 10^-10){
          self$PtE[i, 2] <- 0
          self$PtV[i] <- self$E[e, 1]
        }else if(t > 1 - 10^-10){
          self$PtE[i, 2] <- 1
          self$PtV[i] <- self$E[e, 2]
        }else{
          self$split_edge(e, t)
          self$PtV[i] <- dim(self$V)[1]
        }
    }
    if(!is.null(self$geo.dist)){
      self$compute_geodist()
    }
    if(!is.null(self$res.dist)){
      self$compute_resdist()
    }
    if(!is.null(self$CBobj)) {
      self$buildC(2)
    }
    self$A <- Diagonal(self$nV)[self$PtV, ]
  },

  #' @description Add observations to the object
  #' @param Spoints SpatialPoints or SpatialPointsDataFrame of the observations
  #' @param y        (n x 1) the value of the observations
  #' @param y.index (string, int) column in Spoints where y is located
  add_observations = function(Spoints, y=NULL, y.index=NULL){

    if(is.null(y)){
        y = Spoints@data[,y.index]
    }
    self$y   = y

    SP <- maptools::snapPointsToLines(Spoints, self$Lines)
    coords.old <- Spoints@coords
    colnames(coords.old) <- paste(colnames(coords.old) ,'_old',sep="")
    Spoints@coords = SP@coords
    Spoints@bbox   = SP@bbox
    LtE = cbind(match(SP@data[,1], self$EID),0)
    if("SpatialPointsDataFrame"%in%is(Spoints)){
      Spoints@data <- cbind(Spoints@data,coords.old)
    }else{
      Spoints <- SpatialPointsDataFrame(Spoints, data = coords.old)
    }
    for(ind in unique(LtE[,1])){
        index.p = LtE[,1]==ind
        LtE[index.p,2]=rgeos::gProject(self$Lines[ind,], Spoints[index.p,])
    }
    PtE = LtE
    for(ind in unique(LtE[,1])){
      Es_ind <- which(self$LtE[ind,]>0)
      index.p = which(LtE[,1]==ind)
      for(j in index.p){
        E_ind <- which.min(replace(self$ELend[Es_ind], self$ELend[Es_ind] < LtE[j,2], NA))
        PtE[j,1] <- Es_ind[E_ind]
        PtE[j,2] <- (LtE[j,2] - self$ELstart[PtE[j,1]])/  (self$ELend[PtE[j,1]]- self$ELstart[PtE[j,1]])

      }
    }
    if(is.null(self$Points)){
      self$Points = Spoints
      self$PtE = PtE
    }else{
      self$Points = rbind(self$Points,Spoints)
      self$PtE = rbind(self$PtE,PtE)
    }
  },

  #' @description add observations to the object
  #' @param Spoints SpatialPoints or SpatialPointsDataFrame of the observations
  #' @param y (n x 1) the value of the observations
  #' @param PtE (n x 2) edge index, distance on index
  #' @param normalized if TRUE, then the distances in PtE are assumed to be
  #' normalized to (0,1), default FALSE.
  add_observations2 = function(y, PtE, Spoints=NULL, normalized = FALSE){

    self$y   = y
    if(normalized){
      self$PtE = PtE
    } else {
      self$PtE = cbind(PtE[, 1], PtE[, 2] / self$edge_lengths[PtE[, 1]])
    }

    if(!is.null(self$Lines)){

      if(is.null(Spoints)){
        coords <- c()
        for(i in 1:dim(PtE)[1]){

          LT = private$edge_pos_to_line_pos(self$PtE[i, 1] , self$PtE[i, 2])
          points <- rgeos::gInterpolate(self$Lines[LT[1,1],],
                                        LT[1,2],
                                        normalized = TRUE)
          coords <- rbind(coords, points@coords)
        }
        rownames(coords) <- 1:dim(coords)[1]
        Spoints <- sp::SpatialPoints(coords)
      }
      self$Points = Spoints
    }
  },



  #' @description build Kirchoff constraint matrix from edges, NOT implemented
  #' for circles (i.e. self closed edges)
  #' @param alpha (int) which type of constraint (currently only 2 implemented)
  #' @param edge_constraint (bool) if true add constraint on vertices of
  #' degree 1.
  buildC = function(alpha, edge_constraint=FALSE){

    if(alpha==2){
      i_  =  rep(0, 2*self$nE)
      j_  =  rep(0, 2*self$nE)
      x_  =  rep(0, 2*self$nE)

      count_constraint = 0
      count            = 0
      for(v in 1:self$nV){
        lower.edges  = which(self$E[,1]%in%v)
        upper.edges  = which(self$E[,2]%in%v)
        n_e = length(lower.edges) + length(upper.edges)
        #derivative constraint
        if((edge_constraint & n_e ==1) | n_e > 1) {
          i_[count + 1:n_e] <- count_constraint + 1
          j_[count + 1:n_e] <- c(4 * (lower.edges-1) + 2, 4 * (upper.edges-1) + 4)
          x_[count + 1:n_e]     <- c( rep(1,length(lower.edges)),
                                      rep(-1,length(upper.edges)))
          count <- count + n_e
          count_constraint <- count_constraint + 1
        }
        if(n_e > 1){
          if(length(upper.edges)==0){
            edges <- cbind(lower.edges, 1)

          }else if(length(lower.edges)==0){
            edges <- cbind(upper.edges, 3)

          }else{
            edges <- rbind(cbind(lower.edges, 1),
                           cbind(upper.edges, 3))

          }
          for(i in 2:n_e){
            i_[count + 1:2]   <- count_constraint + 1
            j_[count + 1:2] <- c(4 * (edges[i-1,1] - 1) + edges[i-1, 2],
                                 4 * (edges[i,1]   - 1) + edges[i,   2])
            x_[count + 1:2]     <- c(1,-1)
            count <- count + 2
            count_constraint <- count_constraint + 1
          }
        }
      }
      C <- Matrix::sparseMatrix(i    = i_[1:count],
                                j    = j_[1:count],
                                x    = x_[1:count],
                                dims = c(count_constraint, 4*self$nE) )
      self$C = C

      self$CBobj <- c_basis2(self$C)
      self$CBobj$T <- t(self$CBobj$T)
    }else{
      error("only alpha=2 implimented")
    }
  },

  #' @description build mesh object for graph
  #' @param h maximum distance between mesh nodes
  #' @param n maximum number of nodes per edge
  #' @details The mesh is a list with the objects
  #' - PtE which contains the mesh locations excluding the original vertices
  #' - V the verties of the mesh
  #' - E the edges of the mesh
  #' - n_e the number of vertices in the mesh per original edge in the graph
  #' - h_e the mesh width per edge in the graph
  #' - ind the indices of the vertices in the mesh
  #' - VtE all mesh locations including the original vertices
  build_mesh = function(h,n=NULL){

    self$mesh <- list(PtE = NULL,
                      V = NULL,
                      E = NULL,
                      n_e = rep(0,self$nV),
                      h_e = NULL,
                      ind = 1:self$nV,
                      VtE = NULL)

    self$mesh$V <- self$V
    for (i in 1:dim(self$E)[1]) {
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
        if(is.null(self$Lines)) {
          self$mesh$V <- rbind(self$mesh$V,
                               cbind(self$V[self$E[i,1], 1]*(1 - d.e) +
                                       d.e*self$V[self$E[i, 2],1],
                                     self$V[self$E[i,1], 2]*(1 - d.e) +
                                       d.e*self$V[self$E[i, 2],2]))
        } else {
          Line <- self$Lines[i,]
          val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized=TRUE)
          Points <- gInterpolate(Line, d.e, normalized=TRUE)
          self$mesh$V <- rbind(self$mesh$V, Points@coords)
        }
        V.int <- (max(self$mesh$ind) + 1):(max(self$mesh$ind) + self$mesh$n_e[i])
        self$mesh$ind <- c(self$mesh$ind, V.int)
        self$mesh$E <- rbind(self$mesh$E, cbind(c(self$E[i, 1], V.int),
                                                c(V.int, self$E[i, 2])))
      } else {
        self$mesh$E <- rbind(self$mesh$E, self$E[i, ])
        self$mesh$h_e <- c(self$mesh$h_e,self$edge_lengths[i])
      }
    }
    self$mesh$VtE <- rbind(self$VtEfirst(),self$mesh$PtE)
  },

  #' @description build mass and stiffness matrices for given mesh object
  compute_fem = function(){
    if(is.null(self$mesh)){
      stop("no mesh provided")
    }

    nV <- dim(self$mesh$V)[1]
    self$mesh$G <- self$mesh$C <- Matrix(0,nrow=nV,ncol=nV)

    for(e in 1:dim(self$mesh$E)[1]){
      v1 <- self$mesh$E[e,1]
      v2 <- self$mesh$E[e,2]

      self$mesh$C[v1,v1] <- self$mesh$C[v1,v1] +  self$mesh$h_e[e]/3
      self$mesh$C[v2,v2] <- self$mesh$C[v2,v2] +  self$mesh$h_e[e]/3
      self$mesh$C[v1,v2] <- self$mesh$C[v1,v2] +  self$mesh$h_e[e]/6
      self$mesh$C[v2,v1] <- self$mesh$C[v2,v1] +  self$mesh$h_e[e]/6

      self$mesh$G[v1,v1] <- self$mesh$G[v1,v1] +  1/self$mesh$h_e[e]
      self$mesh$G[v2,v2] <- self$mesh$G[v2,v2] +  1/self$mesh$h_e[e]
      self$mesh$G[v1,v2] <- self$mesh$G[v1,v2] -  1/self$mesh$h_e[e]
      self$mesh$G[v2,v1] <- self$mesh$G[v2,v1] -  1/self$mesh$h_e[e]
    }
  },
  #' @description Computes observation matrix for mesh
  #' @param PtE locations given as (edge number in graph, location on edge)
  mesh_A = function(PtE){
    if(is.null(self$mesh)){
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

  #' @description Find one Edge correspond to each Vertex (warning very
  #' inefficient implementation)
  #' @return VtE (n.v x 2) [1] edge number, [2] 0=  lower edge, 1=  upper edge
  VtEfirst = function(){
    n.V <- dim(self$V)[1]
    VtE <- matrix(0,n.V, 2)

    for(i in 1:n.V){
      Ei <- which(self$E[,1]==i)[1]
      pos <- 0
      if(is.na(Ei)==1){
        pos <- 1
        Ei <- which(self$E[,2]==i)[1]
      }
      VtE[i,] <- c(Ei, pos)
    }
    return(VtE)
  },

  #' @description plot a metric graph
  #' @param plotly use plot_ly for 3D plot (default FALSE)
  #' @param show show the plot?
  #' @param line_width line width for edges
  #' @param marker_size size of markers for vertices
  #' @param vertex_color color of vertices
  #' @param edge_color color of edges
  #' @param data Plot the data?
  #' @param data_size size of markers for data
  #' @param mesh Plot the mesh locations?
  #' @param ... additional arguments for ggplot or plot_ly
  #' @return a plot_ly or or ggplot object
  #' @examples
  #' line1 <- Line(rbind(c(0,0),c(1,0)))
  #' line2 <- Line(rbind(c(0,0),c(0,1)))
  #' line3 <- Line(rbind(c(0,1),c(-1,1)))
  #' theta <- seq(from=pi,to=3*pi/2,length.out = 20)
  #' line4 <- Line(cbind(sin(theta),1+ cos(theta)))
  #' Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
  #'                               Lines(list(line2),ID="2"),
  #'                               Lines(list(line3),ID="3"),
  #'                               Lines(list(line4),ID="4")))
  #' graph <- metric_graph$new(Lines = Lines)
  #' graph$plot()
  plot = function(plotly = FALSE,
                  show = TRUE,
                  line_width = 0.3,
                  marker_size = 3,
                  vertex_color = 'black',
                  edge_color = 'black',
                  data = FALSE,
                  data_size = 1,
                  mesh = FALSE,
                  ...){
    if(!plotly){
      p <- private$plot_2d(show = show,
                           line_width = line_width,
                           marker_size = marker_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           data_size = data_size,
                           mesh = mesh,
                           ...)
    } else {
      p <- private$plot_3d(show = show,
                           line_width = line_width,
                           marker_size = marker_size,
                           vertex_color = vertex_color,
                           edge_color = edge_color,
                           data = data,
                           data_size = data_size,
                           mesh = mesh,
                           ...)
    }
    return(p)
  },

  #' @description plot function X on the graph
  #' @param X Either an m x 3 matrix with (edge number, position on
  #' curve (in length), value) or a vector with values for the function
  #' evaluated at a precomputed mesh.
  #' @param plotly plot in 2D or 3D?
  #' @param show show the plot?
  #' @param graph_color for 3D plot, the color of the graph.
  #' @param graph_width for 3D plot, the line width of the graph.
  #' @param marker_size for 3D plot, the marker size of the vertices
  #' @param color Color of curve
  #' @param ... additional arguments for ggplot or plot_ly
  #' @export
  plot_function = function(X, plotly = TRUE, show = TRUE,
                       graph_color = 'rgb(0,0,0)',
                       graph_width = 1,
                       marker_size = 10,
                       color = 'rgb(0,0,200)',
                       ...){
    if(flat == FALSE){
      p <- self$plot(color = graph_color, line_width = graph_width,
                      marker_size = marker_size)
    } else {
      p <- NULL
    }
    for(i in 1:self$nE){
      if(!is.matrix(X) || is.matrix(X) && dim(X)[1] == 1) {
        ind <- self$mesh$PtE[,1] == i
        vals <- cbind(self$mesh$PtE[ind,2],X[ind])
      } else {
        vals <- X[X[,1]==i,2:3]
      }
      if(isempty(vals))
        next
      if(is.null(self$Lines) == TRUE) {
        p <- private$plot_straight_curve(vals,
                                         self$V[self$E[i,],],
                                         plotly = plotly,
                                         p = p,
                                         color = color,
                                         ...)

      } else {
        if(!is.null(vals)){
          p<-private$plot_edge_line(vals,
                                 Eindex = i,
                                 p = p,
                                 ploty = ploty,
                                 color = color,
                                 ... )
      }

      }
    }
    if(show){
      print(p)
    }
    return(p)
  },
  #' @description plot function mesh X on the graph
  #' @param X (m x 1) a vector with values for the function
  #' evaluated at a precomputed mesh (V,and PtE)
  #' @param plotly use plot_ly for 3D plot?
  #' @param show show the plot?
  #' @param graph_color for 3D plot, the color of the graph.
  #' @param graph_width for 3D plot, the line width of the graph.
  #' @param marker_size for 3D plot, the marker size of the vertices
  #' @param color Color of curve
  #' @param ... additional arguments for ggplot or plot_ly
  #' @export
  plot_function_mesh = function(X, plotly = FALSE, show = TRUE,
                                graph_color = 'black',
                                graph_width = 1,
                                marker_size = 10,
                                color = 'rgb(0,0,200)',
                                ...){
    if(is.null(self$mesh)) {
      stop("no mesh provided")
    }
    if(length(X) != dim(self$V)[1] + dim(self$mesh$PtE)[1]){
      stop("X does not have the correct size")
    }
    if(plotly){
      p <- self$plot(plotly = plotly, color = graph_color,
                     line_width = graph_width, marker_size = marker_size)
    } else {
      p <- NULL
    }
    n.v <- dim(self$V)[1]
    XV <- X[1:n.v]
    for(i in 1:self$nE){
      ind <- self$mesh$PtE[,1] == i
      if(sum(ind)==0){
        vals <- rbind(c(0, XV[self$E[i,1]]),
                      c(1, XV[self$E[i,2]]))

      }else{
        vals <- rbind(c(0, XV[self$E[i,1]]),
                      cbind(self$mesh$PtE[ind,2],X[n.v + which(ind)]),
                      c(1, XV[self$E[i,2]]))

      }



      if(is.null(self$Lines) == TRUE) {
        p <- private$plot_straight_curve(vals,
                                         self$V[self$E[i,],],
                                         plotly = plotly,
                                         p = p,
                                         color = color, ...)

      } else {
        if(!is.null(vals)){
          p <- private$plot_edge_line(vals,
                                 Eindex = i,
                                 p = p,
                                 plotly = plotly,
                                 color = color,
                                 ... )
        }
      }

    }
    if(show){
      print(p)
    }
    return(p)
  },

  #' @description split_edge Function for splitting lines
  #' @param Ei index of line to split
  #' @param t  position on line to split (normalized)
  #' @export
  split_edge = function(Ei, t){
    if (!is.null(self$Lines)) {

      index <- (self$LtE@p[Ei]+1):(self$LtE@p[Ei+1])
      LinesPos <- cbind(self$LtE@i[index] + 1, self$LtE@x[index])
      LinesPos <- LinesPos[order(LinesPos[,2]),,drop=F]
      j <-  min(which(t  <= LinesPos[,2]))
      t_mod <- t
      if(j == 1){
        t_mod <- t_mod/LinesPos[j,2]
      }else{
        t_mod <- (t_mod-LinesPos[j-1,2])/(LinesPos[j,2]-LinesPos[j-1,2])
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

      Line <- self$Lines[LinesPos[j,1], ]
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

      self$LtE <-Matrix::sparseMatrix(i    = c(LtE.i, LtE.i_new)+1,
                                      p    = c(LtE.p, LtE.p_new),
                                      x    = c(LtE.x, LtE.x_new),
                                      dims = LtE.dim )


    }else{
      V1 <- self$V[self$E[Ei, 1], ]
      V2 <- self$V[self$E[Ei, 2], ]
      val_line <- (1-t)*V1 + t*V2
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
    for(i in ind){
      if (self$PtE[i, 2] >= t - 1e-10) {
        self$PtE[i, 1] <- self$nE
        self$PtE[i, 2] <- abs(self$PtE[i, 2] - t)/(1 - t)
      }
    }
    self$nV <- dim(self$V)[1]


  }

  ),

  private = list(
    #' @description computes which line and which position t_E on Ei belongs to
    #' @param Ei  (int)   edge index
    #' @param t_e (n x 1) number of positions on Ei
    edge_pos_to_line_pos = function(Ei, t_E){

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

    #' @description plotting edges when there exists lines
    #' @param vals    (m x 2) position on edge, value
    #' @param Eindex  (int)   which edge
    #' @param p       (ggplot) previous plot object
    #' @param plotly   (?)
    #' @param color   (?)
    #'
    plot_edge_line = function(vals, Eindex, p, plotly, color,  ... ){
      index <- (self$LtE@p[Eindex]+1):(self$LtE@p[Eindex+1])
      LinesPos <- cbind(self$LtE@i[index] + 1, self$LtE@x[index])
      LinesPos <- LinesPos[order(LinesPos[,2]),,drop=F]
      for(j in 1:length(index)){
        if(j==1){
          index_j <- vals[,1] <= LinesPos[j,2]
        }else{
          index_j <- (vals[,1] <= LinesPos[j,2]) &  (vals[,1] > LinesPos[j-1,2])
        }
        if(sum(index_j) > 0)
          next
        rel.pos = vals[index_j,1]
        if(j == 1){
          rel.pos <- rel.pos/LinesPos[j,2]
        }else{
          rel.pos <- (rel.pos-LinesPos[j-1,2])/(LinesPos[j,2]-LinesPos[j-1,2])
        }

        if(j== dim(LinesPos)[1] )
          rel.pos = self$ELend[Eindex]*rel.pos
        if(j==1)
          rel.pos = rel.pos + self$ELstart[Eindex]

        p <- private$plot_curve(cbind(rel.pos,vals[index_j,2]),
                                SpatialLines(list(self$Lines@lines[[LinesPos[j,1]]])),
                                plotly = plotly,
                                p = p,
                                color = color,
                                ...)

      }
      return(p)
    },
  #'
  #'@description function for creating Vertex and Edges from self$lines
  #'
  #'
  line_to_vertex = function(){
    lines <- c()
    for(i in 1:length(self$Lines)){
      points <- self$Lines@lines[[i]]@Lines[[1]]@coords
      n <- dim(points)[1]
      lines <- rbind(lines,
                     c(i, points[1,],
                       sp::LineLength(self$Lines@lines[[i]]@Lines[[1]])),
                     c(i, points[n,],
                       sp::LineLength(self$Lines@lines[[i]]@Lines[[1]])))
    }

    vertex <- lines[1, , drop = FALSE]
    for (i in 1:dim(lines)[1]) {
      tmp <- vertex[, 2:3] - matrix(rep(1,dim(vertex)[1]),
                                   dim(vertex)[1], 1) %*% lines[i,2:3]
      if (min(rowSums(tmp^2)) > 1e-8) {
        vertex <- rbind(vertex, lines[i, ])
      }
    }
    lvl <- matrix(0, nrow = max(lines[,1]), 4)
    for (i in 1:max(lines[, 1])) {
      which.line <- sort(which(lines[, 1] == i))
      line <- lines[which.line, ]
      #ind1 <- (abs( vertex[,2] - line[1,2] )< 1e-10) * (abs( vertex[,3] - line[1,3] )< 1e-10)==1
      #ind2 <-  (abs( vertex[,2] - line[2,2] )< 1e-10) * (abs( vertex[,3] - line[2,3] )< 1e-10)==1
      ind1 <- which.min((vertex[,2] - line[1,2])^2 + ( vertex[,3] - line[1,3] )^2)
      ind2 <- which.min((vertex[,2] - line[2,2])^2 + ( vertex[,3] - line[2,3] )^2)
      lvl[i,1] <- i
      lvl[i,2] <- ind1#which(ind1)
      lvl[i,3] <- ind2#which(ind2)
      lvl[i,4] <- line[1,4]
    }
    self$V <- vertex[, 2:3]
    self$E <- lvl[, 2:3, drop = FALSE]
    self$edge_lengths <- lvl[,4]
    self$nV <- dim(self$V)[1]

    self$LtE <-Matrix::sparseMatrix(j    = 1:dim(self$E)[1],
                                    i    = c(1:length(self$Lines)),
                                    x    = rep(1,dim(self$E)[1]),
                                    dims = c(dim(self$E)[1], length(self$Lines)) )
    self$ELend = rep(1,dim(self$E)[1])
    self$ELstart = rep(0,dim(self$E)[1])

  },
  plot_curve = function(data.to.plot,
                        Line_edge,
                        plotly = TRUE,
                        normalized = TRUE,
                        color = 'rgb(0,0,200)',
                        p = NULL, ...){

    data.to.plot.order <- data.to.plot[order(data.to.plot[, 1]), , drop = FALSE]
    p2 <- rgeos::gInterpolate(Line_edge, data.to.plot.order[, 1, drop = FALSE],
                              normalized = normalized)
    coords <-p2@coords
    data <- data.frame(x = coords[,1], y = coords[,2],
                       z = data.to.plot.order[,2])
    if(!plotly){
      if(is.null(p)){
        p <- ggplot2::ggplot(data = data,
                             ggplot2::aes(x = x, y = y,
                                          colour = z)) +
          ggplot2::geom_path(...)
      } else {
        p <- p + ggplot2::geom_path(data = data)
      }
    } else {
      if(is.null(p)){
        p <- plot_ly(data=data, x = ~y, y = ~x, z = ~z)
        p <- p %>% add_trace(mode="lines", type="scatter3d",
                             line = list(color = color, ...),
                             showlegend = FALSE)
      } else {
        p <- p %>% add_trace(data=data, x = ~y, y = ~x, z = ~z,
                             mode = "lines", type = "scatter3d",
                             line = list(color = color, ...),
                             showlegend = FALSE)
      }
    }
    return(p)
  },
  plot_straight_curve = function(data.to.plot,
                                 V,
                                 plotly = FALSE,
                                 p = NULL,
                                 line_color,
                                 ...){

    data.to.plot.order <- data.to.plot[order(data.to.plot[, 1]), ]
    l <- sqrt(sum((V[1,] - V[2,])^2))
    alpha <- data.to.plot.order[,1]/l
    coords <- cbind((1 - alpha) * V[1, 1] + alpha * V[2, 1],
                    (1 - alpha) * V[1, 2] + alpha * V[2, 2])
    data <- data.frame(x = coords[, 1], y = coords[, 2],
                       z = data.to.plot.order[, 2])
    if (!plotly) {
      if (is.null(p)) {
        p <- ggplot2::ggplot(data = data,
                             ggplot2::aes(x = x, y = y,
                                          colour = z)) +
          ggplot2::geom_path(...)
      } else {
        p <- p + ggplot2::geom_path(data = data, ...)
      }
    } else {
      if (is.null(p)) {
        p <- plot_ly(data=data, x = ~y, y=~x,z=~z)
        p <- p %>% add_trace(mode="lines", type="scatter3d",
                             line = list(...), showlegend = FALSE)
      } else {
        p <- p %>% add_trace(data=data, x = ~y, y = ~x, z = ~z,
                             mode = "lines", type = "scatter3d",
                             line = list(...), showlegend = FALSE)
      }
    }
    return(p)
  },
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

  plot_2d = function(show = TRUE,
                     line_width = 0.1,
                     marker_size = 1,
                     vertex_color = 'black',
                     data = FALSE,
                     data_size = 1,
                     mesh = FALSE,
                     ...){
    xyl <- c()
    if(is.null(self$Lines)){
      xyl <- cbind(c(self$V[E[, 1], 1], self$V[E[, 2], 1]),
                   c(self$V[E[, 1], 2], self$V[E[, 2], 2]),
                   c(1:self$nE, 1:self$nE))
    } else {
      coords <- lapply(coordinates(self$Lines), function(x) x[[1]])
      nc <- do.call(rbind,lapply(coords, function(x) dim(x)[1]))
      xyl <- cbind(do.call(rbind,coords), rep(1:length(nc), times = nc))
    }

    p <- ggplot()+ geom_path(data = data.frame(x = xyl[, 1],
                                              y = xyl[,2],
                                              group = xyl[,3]),
                             mapping = aes(x = x, y = y, group = group),
                             linewidth = line_width)

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
        LT <- private$edge_pos_to_line_pos(self$PtE[i, 1],self$PtE[i, 2])
        Line <- self$Lines[LT[1, 1], ]
        val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized = TRUE)
        Point <- gInterpolate(Line,LT[1, 2], normalized = TRUE)
        x <- c(x, Point@coords[1])
        y <- c(y, Point@coords[2])
      }
      p <- p + geom_point(data = data.frame(x = x, y = y, val = as.vector(self$y)),
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
    p <- p + coord_fixed()
    if(show){
      print(p)
    }
    return(p)
  },
  plot_3d = function(show = TRUE,
                     line_width = 1,
                     marker_size = 1,
                     vertex_color = 'rgb(0,0,0)',
                     edge_color = 'rgb(0,0,0)',
                     data = FALSE,
                     data_size = 1,
                     mesh = FALSE,
                     ...){
    if (is.null(self$Lines)) {
      data.plot <- data.frame(x = c(self$V[E[, 1], 1], self$V[E[, 2], 1]),
                              y = c(self$V[E[, 1], 2], self$V[E[, 2], 2]),
                              z = rep(0, 2 * self$nE),
                              i = c(1:self$nE, 1:self$nE))
    } else {
      x <- y <- ei <- NULL
      for (i in 1:self$nE) {
        xi <- self$Lines@lines[[i]]@Lines[[1]]@coords[, 1]
        yi <- self$Lines@lines[[i]]@Lines[[1]]@coords[, 2]
        ii <- rep(i,length(xi))
        x <- c(x, xi)
        y <- c(y, yi)
        ei <- c(ei, ii)
      }
      data.plot <- data.frame(x = x, y = y, z = rep(0,length(x)), i = ei)
    }
    p <- plot_ly(data=data.plot, x = ~y, y = ~x, z = ~z)
    p <- p %>% add_trace(data = data.plot, x = ~y, y = ~x, z = ~z,
                         mode = "lines", type = "scatter3d",
                         line = list(width = line_width,
                                     color = edge_color, ...),
                         split = ~i, showlegend = FALSE)

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
        Line <- self$Lines[PtE[i, 1], ]
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

    ax <- list(title = '',
               zeroline = FALSE,
               showgrid = FALSE,
               showticklabels=FALSE)
    p <- p %>% layout(title = '',
                      scene = list(xaxis = ax, yaxis = ax, zaxis = ax,
                                   aspectmode='data'))

    if(show){
      print(p)
    }
    return(p)
  }

))




