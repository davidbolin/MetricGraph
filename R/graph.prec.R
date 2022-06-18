#' computes the precision matrix for the vertices for an exponential graph model
#' @param P (n x 2) Vertex locations
#' @param E (m x 2) Edges
#' @param kappa  matern parameter
Q.vertex.exp <- function(P,E,kappa,l=NULL,adjust=FALSE){
  if(is.null(l)){
    l <- sqrt((P[E[,2],1]-P[E[,1],1])^2 + (P[E[,2],2]-P[E[,1],2])^2)
  }
  Q <- Matrix(0,nrow=dim(P)[1],ncol=dim(P)[1])
  for(i in 1:dim(E)[1]){
    Q[E[i,1],E[i,1]] = Q[E[i,1],E[i,1]] + 1/2 + exp(-2*kappa*l[i])/(1-exp(-2*kappa*l[i]))
    Q[E[i,2],E[i,2]] = Q[E[i,2],E[i,2]] + 1/2 + exp(-2*kappa*l[i])/(1-exp(-2*kappa*l[i]))
    Q[E[i,1],E[i,2]] = Q[E[i,1],E[i,2]] - exp(-kappa*l[i])/(1-exp(-2*kappa*l[i]))
    Q[E[i,2],E[i,1]] = Q[E[i,2],E[i,1]] - exp(-kappa*l[i])/(1-exp(-2*kappa*l[i]))
  }
  if(adjust){
    Q.ind <- abs(Q)
    index <- which(rowSums(Q.ind>0)<3)
    diag(Q)[index] <- diag(Q)[index] + 0.5
  }
  return(2*kappa*Q)
}
