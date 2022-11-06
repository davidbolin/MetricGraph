#' @title Metric graph object for specification of Gaussian processes
#' @description Class representing general metric graphs.
#' @details A graph object created from vertex and edge matrices, or from an sp::Lines
#' object where each line is representing and edge.
#'
#' @export
metric_graph <-  R6::R6Class("GPGraph::graph", public = list(

  #' @field V Position in Euclidean space of the vertices
  V = NULL,

  #' @field nV number of vertices
  nV = 0,

  #' @field E Edges,  E[i,1] is the vertex at the start of the edge and  E[i,2] is
  #' the vertex at the end of the edge
  E = NULL,

  #' @field nE number of edges
  nE= 0,

  #' @field edge_lengths length of edges
  edge_lengths = NULL,

  #' @field EID ID of edges
  EID = NULL,

  #' @field C constraint matrix used to set Kirchhoff constraints
  C = NULL,

  #' @field A observation matrix specifying which vertices are observation locations
  A = NULL,

  #' @field CBobj svd stuct obj
  CBobj = NULL,

  #' @field Points Observations in SpatialPointsDataFrame
  Points  = NULL,

  #' @field y the data connected to P
  y = NULL,

  #' @field PtE Points to Line (connected to Points),
  #'  [,1] - edge index,
  #'  [,2] - distance along the line (i.e. distance to initial point)
  PtE = NULL,

  #' @field PtV Points to Vertex observations to vertex
  #'  [,1] - vertex index,
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
  #' @details A graph object can be initialized in two ways. The first method is
  #' to specify P and E. In this case, if edge_lengths is not specified, all edges are
  #' assumed to be straight lines. Otherwise the edge lengths set in edge_lengths are used.
  #' The second option is to specify the graph based on Lines. In this case,
  #' the vertices are set by the end points of the lines. Thus, if two lines are intersecting
  #' somewhere else, this will not be viewed as a vertex.
  #' @return A gpgraph_graph object
  initialize = function(Lines = NULL, P = NULL, E = NULL, edge_lengths = NULL) {
    #We have three different ways of initializing:

    #option 1: initialization from lines
    if(!is.null(Lines)){
      if(!is.null(edge_lengths) || !is.null(P) || !is.null(E)){
        warning("object initialized from lines, then E,P,edge_lengths are ignored")
      }
      self$nE = length(Lines)
      self$Lines = Lines
      self$EID = sapply(slot(self$Lines,"lines"), function(x) slot(x, "ID"))
      list_obj <- vertex.to.line(self$Lines)
      self$V   <- list_obj$V
      self$nV <- dim(self$V)[1]
      self$E <- list_obj$EtV
      self$edge_lengths  <- list_obj$edge_lengths
    } else {
      if(is.null(P) || is.null(E)){
        stop("You must supply Lines or P and E")
      }
      self$nE <- dim(E)[1]
      self$V   <- P
      self$E <- E
      self$nV <- dim(self$V)[1]
      self$edge_lengths <- sqrt((self$V[self$E[,2], 1] - self$V[self$E[, 1], 1])^2 +
                            (self$V[self$E[,2], 2] - self$V[self$E[, 1], 2])^2)
    }

  },

  #' @description Split line by point
  #' @param Ei Index of edge to be split
  #' @param t Normalized distance to first edge
  split_line = function(Ei, t){
      Line <- self$Lines[Ei,]

      val_line <- gProject(Line, as(Line, "SpatialPoints"), normalized=TRUE)
      ind <-  (val_line <= t)
      Point <- gInterpolate(Line, t, normalized=TRUE)
      Line1 <- list(as(Line, "SpatialPoints")[ind, ],Point)
      Line2 <- list(Point, as(Line, "SpatialPoints")[ind==F, ])

      if(sum(is(self$Lines)%in%"SpatialLinesDataFrame") > 0){
        self$Lines <-rbind(self$Lines[1:Ei-1,],
                            SpatialLinesDataFrame(as(do.call(rbind,  Line1), "SpatialLines"),
                                                  data=Line@data,match.ID = FALSE),
                            self$Lines[-(1:Ei),],
                            SpatialLinesDataFrame(as(do.call(rbind,  Line2), "SpatialLines"),
                                                  data=Line@data,match.ID = FALSE))
      }else{
        self$Lines <-rbind(self$Lines[1:Ei-1,],
                           as(do.call(rbind,  Line1), "SpatialLines"),
                           self$Lines[-(1:Ei),],
                           as(do.call(rbind,  Line2), "SpatialLines"))

      }


    newV <- self$nV+1#max(self$V[,1])+1
    self$V <- rbind(self$V,c(Point@coords))


    l_e <- self$edge_lengths[Ei]
    self$edge_lengths[Ei] <- t*l_e
    self$edge_lengths <- c(self$edge_lengths, (1-t)*l_e)
    self$nE <- self$nE + 1
    self$E <- rbind(self$E,
                    c(newV,self$E[Ei, 2]))
    self$E[Ei, 2] <- newV
    ind <- which(self$PtE[,1]%in%Ei)
    for(i in ind){
      if( self$PtE[i,2]>= t*l_e-1e-10){
        self$PtE[i,1] <- length(self$edge_lengths) #self$nE #
        self$PtE[i,2] <- abs(self$PtE[i,2]-t*l_e)
      }
    }
    self$nV <- dim(self$V)[1]
  },

  #' @description Computes shortest path distances between the vertices in the graph
  compute_geodist = function(){
    g <- graph(edges = c(t(self$E)), directed=FALSE)
    E(g)$weight <- self$edge_lengths
    self$geo.dist <- distances(g)
  },

  #' @description Computes the resistance metric between the vertices in the graph
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
      Wmat[self$E[i,1],self$E[i,2]] <- Wmat[self$E[i,2],self$E[i,1]] <- 1/self$edge_lengths[i]
    }
    self$Laplacian <- Diagonal(self$nV,as.vector(Matrix::rowSums(Wmat))) - Wmat
  },

  #' @description Add observation locations as vertices in the graph
  observation_to_vertex = function(){

    l <- length(self$PtE[,1])
    self$PtV <- rep(0,l)
    for(i in 1:l){
        e <- as.vector(self$PtE[i,1])
        t <- as.vector(self$PtE[i,2])
        l_e <- self$edge_lengths[e]
        if(abs(t)<10^-10){
          self$PtE[i,2] <- 0
          self$PtV[i] <- self$E[e,1]
        }else if(t>l_e-10^-10){
          self$PtE[i,2] <- l_e
          self$PtV[i] <- self$E[e,2]
        }else{
          self$split_line(e, t/l_e)
          self$PtV[i] <- dim(self$V)[1]
        }
    }
    if(!is.null(self$geo.dist)){
      self$compute_geodist()
    }
    if(!is.null(self$res.dist)){
      self$compute_resdist()
    }
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

    SP <- snapPointsToLines(Spoints, self$Lines)
    coords.old <- Spoints@coords
    colnames(coords.old) <- paste(colnames(coords.old) ,'_old',sep="")
    Spoints@coords = SP@coords
    Spoints@bbox   = SP@bbox
    PtE = cbind(match(SP@data[,1], self$EID),0)
    if("SpatialPointsDataFrame"%in%is(Spoints)){
      Spoints@data <- cbind(Spoints@data,coords.old)
    }else{
      Spoints <- SpatialPointsDataFrame(Spoints, data = coords.old)
    }
    for(ind in unique(PtE[,1])){
        index.p = PtE[,1]==ind
        PtE[index.p,2]=rgeos::gProject(self$Lines[ind,], Spoints[index.p,])
    }

    self$Points = Spoints
    self$PtE = PtE
  },

  #' @description add observations to the object
  #' @param Spoints SpatialPoints or SpatialPointsDataFrame of the observations
  #' @param y        (n x 1) the value of the observations
  #' @param PtE      (n x 2) edge index, distance on index
  add_observations2 = function(y, PtE, Spoints=NULL){

    self$y   = y
    self$PtE = PtE
    if(is.null(Spoints)){
      Edges <- unique(self$PtE[,1])
      coords <- c()
      for(e in Edges){
        ind <- self$PtE[,1] == e
        points <- rgeos::gInterpolate(self$Lines[e,], self$PtE[e, 2], normalized = F)
        coords <- rbind(coords, points@coords)
      }
      Spoints <- sp::SpatialPoints(coords)
    }
    self$Points = Spoints
  },



  #' @description build Kirchoff constraint matrix from edges, NOT implemented for circles (i.e. self closed edges)
  #' @param alpha (int) which type of constraint (currently only 2 implemented)
  #' @param edge_constraint (bool) if true add constraint on vertices of degree 1.
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

  #' @description build mesh object for plotting
  #' @param h maximum distance between mesh nodes
  #' @param n maximum number of nodes per edge
  build_mesh = function(h,n=NULL){

    self$mesh <- list(V=NULL,
                      E=NULL,
                      n_e = NULL,
                      ind = 1:dim(self$V)[1])

    self$mesh$V <- self$V
    for(i in 1:dim(self$E)[1]){
      if(is.null(n)){
        n.e[i] <- max(ceiling(self$edge_lengths[i]/h)+1,3)
      } else {
        n.e[i] = n
      }
      d.e <- seq(from=0,to=1,length.out=n.e[i])[2:(n.e[i]-1)]
      hi[i] <- L[i]*d.e[1]
      n.e[i] <- n.e[i]-2

      self$mesh$V <- rbind(self$mesh$V,
                      cbind(self$V[self$E[i,1],1]*(1-d.e) + d.e*self$V[self$E[i,2],1],
                            self$V[self$E[i,1],2]*(1-d.e) + d.e*self$V[self$E[i,2],2]))
      V.int <- (max(self$mesh$ind)+1):(max(self$mesh$ind)+n.e[i])
      self$mesh$ind <- c(self$mesh$ind,V.int)
      self$mesh$E <- rbind(self$mesh$E, cbind(c(self$E[i,1],V.int), c(V.int,self$E[i,2])))
    }
  }

))

#' Convert line object to graph information
#' @param Lines object
#' @return A list with elements V, EtV, edge_lengths
#' @export
vertex.to.line <- function(Lines){
  lines <- c()
  for(i in 1:length(Lines)){
    points <- Lines@lines[[i]]@Lines[[1]]@coords
    n <- dim(points)[1]
    lines <- rbind(lines,c(i, points[1,], sp::LineLength( Lines@lines[[i]]@Lines[[1]])),
                         c(i, points[n,], sp::LineLength( Lines@lines[[i]]@Lines[[1]])))
  }

  index.dub <- duplicated(lines[,2:3])
  vertex <- cbind( 1:sum(!index.dub), lines[!index.dub,2:3,drop=F])

  lvl <- matrix(0, nrow= max(lines[,1]), 4)
  for(i in 1:max(lines[,1])){
    which.line <- sort(which(lines[,1]==i))
    line <- lines[which.line,]
    ind1 <- (abs( vertex[,2] - line[1,2] )< 1e-10)* (abs( vertex[,3] - line[1,3] )< 1e-10)==1
    ind2 <-  (abs( vertex[,2] - line[2,2] )< 1e-10)* (abs( vertex[,3] - line[2,3] )< 1e-10)==1

    lvl[i,1] <- i
    lvl[i,2] <- which(ind1)
    lvl[i,3] <- which(ind2)
    lvl[i,4] <- line[1,4]
  }

  return(list(V = vertex[,2:3],
               EtV    = lvl[,2:3, drop=FALSE],
               edge_lengths     = lvl[,4]))
}



