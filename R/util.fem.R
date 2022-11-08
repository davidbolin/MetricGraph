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
