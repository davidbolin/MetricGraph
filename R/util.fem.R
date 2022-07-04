
#' Compute edge lengths in graph
#'
#' @param P Vertex matrix
#' @param E Edge matrix
#'
#' @return Vector with lengths
#' @export
compute.lengths <- function(P,E){
  L <- NULL
  for(i in 1:dim(E)[1]){
    L[i] <- sqrt( sum( (P[E[i,1],]-P[E[i,2],])^2))
  }
  return(L)
}

#' compute the mesh with a given mesh width h
#'
#' @param P Vertex matrix of graph
#' @param E Edge matrix of graph
#' @param h desired mesh width
#' @param n desired number of vertices per edge in the graph
#'
#' @return list
#' @export
graph.mesh <- function(P,E,h,n=NULL){
  L <- compute.lengths(P,E)
  n.e <- NULL
  hi <- NULL
  P.edge <- P
  E.edge <- NULL
  V <- 1:dim(P)[1]
  for(i in 1:dim(E)[1]){
    if(is.null(n)){
      n.e[i] <- max(ceiling(L[i]/h)+1,3)
    } else {
      n.e[i] = n
    }
    d.e <- seq(from=0,to=1,length.out=n.e[i])[2:(n.e[i]-1)]
    hi[i] <- L[i]*d.e[1]
    n.e[i] <- n.e[i]-2

    P.edge <- rbind(P.edge,
                    cbind(P[E[i,1],1]*(1-d.e) + d.e*P[E[i,2],1],
                          P[E[i,1],2]*(1-d.e) + d.e*P[E[i,2],2]))
    V.int <- (max(V)+1):(max(V)+n.e[i])
    V <- c(V,V.int)
    E.edge <- rbind(E.edge, cbind(c(E[i,1],V.int), c(V.int,E[i,2])))
  }
  return(list(P = P.edge,V=V, E=E.edge,n.e=n.e))
}


#compute the FEM matrices for a given mesh
graph.fem <- function(mesh){
  nV <- dim(mesh$P)[1]
  G <- C <- Matrix(0,nrow=nV,ncol=nV)

  for(e in 1:dim(mesh$E)[1]){
    e1 <- mesh$E[e,1]
    e2 <- mesh$E[e,2]
    h <- norm(as.matrix(mesh$P[e1,]-mesh$P[e2,]),type="2")
    C[e1,e1] <- C[e1,e1] +  h/3
    C[e2,e2] <- C[e2,e2] +  h/3
    C[e1,e2] <- C[e1,e2] +  h/6
    C[e2,e1] <- C[e2,e1] +  h/6
    G[e1,e1] <- G[e1,e1] +  1/h
    G[e2,e2] <- G[e2,e2] +  1/h
    G[e1,e2] <- G[e1,e2] -  1/h
    G[e2,e1] <- G[e2,e1] -  1/h
  }
  return(list(C=C,G=G))
}
