library(sp)
library(Matrix)
#'
#' A general graph objects,
#' builds an object of sp::Lines object each Lines is assumes to be an edge.
#' Thus graphs can only be connected by end points
#' @export
#'
graph.obj <-  R6::R6Class("GPGraph::graph", list(

  #' @field El length of edges
  El = NULL,

  #' @field V poisition in the space
  V =  NULL,
  #' @field EtV [,1]- index of Lines, [,2] - vertex lower edge, [,3] - vertex upper edge
  EtV = NULL,
  #' @field nE number of edges
  nE= 0,
  #' @field Lines List of Lines object for building the graph
  Lines = NULL, #sp object contaning the lines

  #' @param Lines sp object SpatialLinesDataFrame or list(Lines)
  #'
  initialize = function(Lines.in = NULL) {
    if(is.null(Lines.in))
      return()

    if(is(Lines.in)=="Line")
      Lines.in = sp::Lines(Lines.in,'1')

    if(is(Lines.in)=="Lines")
      Lines.in = list(Lines.in)

    self$Lines=  Lines.in
    self$nE = length(Lines.in)
    list_obj <- vertex.to.line(self$Lines)
    self$V   <- list_obj$V
    self$EtV <- list_obj$EtV
    self$El  <- list_obj$El
  },
  get_name = function(){return('GPGraph::graph')}))

#' Assume that there is only one line each Lines@Lines
#' @param Lines object
#' @return V, EtV, El
#' from Line to vertex
#' @export
vertex.to.line <- function(Lines){
  lines <- c()
  for(i in 1:length(Lines)){
    points <- Lines[[i]]@Lines[[1]]@coords
    n <- dim(points)[1]
    lines <- rbind(lines,c(i, points[1,], sp::LineLength( Lines[[i]]@Lines[[1]])),
                         c(i, points[n,], sp::LineLength( Lines[[i]]@Lines[[1]])))
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
               El     = lvl[4]))
}


#'
#' Building the precision matrix for the expontial case on each vertex
#' @param theta - kappa
#' @param V vertex position
#' @param EtV [,2-3] index of upper and lower edge
#' @param El length of each vertex
#' @return Q (precision matrix)
#' @export
Q.exp <- function(theta, V,EtV, El){

  i_ <- j_ <- x_ <- rep(0, length(V)*4)
  nE <- dim(EtV)[1]
  for(i in 1:nE){
    l_e <- El[i]
    c1 <- exp(-theta*l_e)
    c2 <- c1^2
    one_m_c2 = 1-c2
    c_1 = 0.5 + c2/one_m_c2
    c_2 = -c1/one_m_c2
    count <- 0
    if(EtV[i,2]!=EtV[i,3]){

      i_[count + 1] <- EtV[i,2]
      j_[count + 1] <- EtV[i,2]
      x_[count + 1] <- c_1

      i_[count + 2] <- EtV[i,3]
      j_[count + 2] <- EtV[i,3]
      x_[count + 2] <- c_1


      i_[count + 3] <- EtV[i,2]
      j_[count + 3] <- EtV[i,3]
      x_[count + 3] <- c_2

      i_[count + 4] <- EtV[i,3]
      j_[count + 4] <- EtV[i,2]
      x_[count + 4] <- c_2
      count <- count + 4
    }else{
      i_[count + 1] <- EtV[i,2]
      j_[count + 1] <- EtV[i,2]
      x_[count + 1] <- tanh(0.5 * theta * l_e)
      count <- count + 1
    }
  }
  if(EtV[nE,2]!=EtV[nE,3]){
    count <- count - 4
  }else{
    count <- count - 1
  }
  n.v <- dim(V)[1]
  Q <- Matrix::sparseMatrix(i=i_[1:count],
                            j=j_[1:count],
                            x=x_[1:count],
                            dims=list(n.v, n.v))
  return(Q)
}
