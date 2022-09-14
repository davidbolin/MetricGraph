library(sp)
library(Matrix)
library(rgeos)
library(CB)
#'
#' A general graph objects,
#' builds an object of sp::Lines object each Lines is assumes to be an edge.
#' Thus graphs can only be connected by end points
#' @export
#'
graph.obj <-  R6::R6Class("GPGraph::graph", list(

  #' @field El length of edges
  El = NULL,
  #' @field EID ID of edges
  EID = NULL,

  #' @field A constraint matrix, just to setup Kirchoff constraint
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

  #' @field V poisition in the space [,1] - id [,-1] - point
  #'
  V =  NULL,
  #' @field EtV [,1]- index of Lines, [,2] - vertex lower edge, [,3] - vertex upper edge
  EtV = NULL,
  #' @field nE number of edges
  nE= 0,
  #' @field Lines List of Lines object for building the graph
  Lines = NULL, #sp object contaning the lines

  #' @param Lines.in sp object SpatialLinesDataFrame or SpatialLines
  #'
  initialize = function(Lines.in = NULL) {
    if(is.null(Lines.in))
      return()

    self$Lines=  Lines.in
    self$nE = length(Lines.in)
    self$EID = sapply(slot(self$Lines,"lines"), function(x) slot(x, "ID"))
    list_obj <- vertex.to.line(self$Lines)
    self$V   <- list_obj$V
    self$EtV <- list_obj$EtV
    self$El  <- list_obj$El
  },

  #' add observations to the object
  #' @param Spoints SpatialPoints or SpatialPointsDataFrame of the observations
  #' @param y        (n x 1) the value of the observations
  #' @param y.index (string, int) position in Spoints where y is located
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

  #' add observations to the object
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
        points <- rgeos::gInterpolate(graph$Lines[e,], self$PtE[, 2], normalized = F)
        coords <- rbind(coords, points@coords)
      }
      Spoints <- sp::SpatialPoints(coords)
    }
    self$Points = Spoints
  },



  #' build Kirchoff constraint matrix from edges, NOT implemented for circles (i.e. self closed edges)
  #' @param alpha (int) which type of constraint (currently only 2 implemented)
  #' @param edge_constraint (bool) if true add constraint on vertices of degree 1.
  buildA = function(alpha, edge_constraint=FALSE){

    if(alpha==2){
      nE = dim(self$EtV)[1]
      i_  =  rep(0, 2*dim(self$EtV)[1])
      j_  =  rep(0, 2*dim(self$EtV)[1])
      x_  =  rep(0, 2*dim(self$EtV)[1])

      count_constraint = 0
      count            = 0
      for(v in 1:dim(self$V)[1]){
        lower.edges  = which(self$EtV[,2]%in%v)
        upper.edges  = which(self$EtV[,3]%in%v)
        n_e = length(lower.edges) + length(upper.edges)
        #derivative constraint
        if((edge_constraint & n_e ==1) | n_e > 1) {
          i_[count + 1:n_e] <- count_constraint + 1
          j_[count + 1:n_e] <- c(4 * (lower.edges-1) + 2, 4 * (upper.edges-1) + 4)
          x_[count + 1:n_e]     <- c( 1*length(lower.edges),-
                                      1*length(upper.edges))
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
      A <- Matrix::sparseMatrix(i    = i_[1:count],
                                j    = j_[1:count],
                                x    = x_[1:count],
                                dims = c(count_constraint, 4*nE))
      self$A = A

      self$CBobj <- CB::c_basis2_cpp(self$A)
      self$CBobj$T <- t(self$CBobj$T)
    }else{
      error("only alpha=2 implimented")
    }
  }

))

#' Assume that there is only one line each Lines@Lines
#' @param Lines object
#' @return V, EtV, El
#' from Line to vertex
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
  colnames(vertex) <- c("vert","x.coord","y.coord")


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

  return(list(V = vertex,
               EtV    = lvl[,1:3, drop=FALSE],
               El     = lvl[,4]))
}



