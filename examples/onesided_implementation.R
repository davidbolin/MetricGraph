# Function to detect starting points in graph
find.mesh.starts <- function(graph){
  if(attr(graph$mesh,"continuous")) {
    V <- graph$mesh$V[1:graph$nV,]
    deg1 <- which(graph$get_degrees() == 1)
    return(deg1[deg1 %in% graph$E[,1]])
  } else {
    return(which(graph$mesh$PtE[,2]==0))
  }
}

num_in_out <- function(graph) {
  outs <- rep(0,graph$nV)
  ins <- rep(0, graph$nV)
  for(i in 1:graph$nV) {
    outs[i] <- sum(graph$E[,1] == i)
    ins[i] <- sum(graph$E[,2] == i)
  }
  return(list(outs = outs, ins = ins))
}
mesh_merge_outs <- function(graph) {
  outs <- num_in_out(graph)$outs
  ind <- which(outs > 1)
  while(length(ind)>0) {
    #find edges going out
    e.ind <- which(graph$E[,1]==ind[1])
    V.keep <- which(graph$mesh$PtE[,1]==e.ind[1])[1]
    V.rem <- which(graph$mesh$PtE[,1]==e.ind[2])[1]
    graph$mesh$PtE <- graph$mesh$PtE[-V.rem,]
    graph$mesh$V <- graph$mesh$V[-V.rem,]
    graph$mesh$E[graph$mesh$E == V.rem] <- V.keep
    graph$mesh$E[graph$mesh$E>V.rem] <- graph$mesh$E[graph$mesh$E>V.rem] - 1
    outs[ind] <- outs[ind] - 1
    ind <- which(outs > 1)
  }
}

mesh_merge_deg2 <- function(graph) {
  outs <- num_in_out(graph)$outs
  ins <- num_in_out(graph)$outs
  ind <- which(outs == 1 & ins == 1)
  cat(ind)
  for(i in 1:length(ind)) {
    V.keep <- graph$mesh$V[ind[i],]
    V.rem <- which(graph$mesh$V[,1] == V.keep[1] & graph$mesh$V[,2] == V.keep[2])[2]
    graph$mesh$PtE <- graph$mesh$PtE[-V.rem,]
    graph$mesh$V <- graph$mesh$V[-V.rem,]
    graph$mesh$E[graph$mesh$E == V.rem] <- ind[i]
    graph$mesh$E[graph$mesh$E>V.rem] <- graph$mesh$E[graph$mesh$E>V.rem] - 1
  }
}

build_mesh = function(self, h=NULL,n=NULL, continuous = TRUE, merge.outs = FALSE) {

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

  if(continuous) {
    self$mesh <- list(PtE = NULL,
                      V = NULL,
                      E = NULL,
                      n_e = rep(0, self$nV),
                      h_e = NULL,
                      ind = 1:self$nV,
                      VtE = NULL)
    attr(self$mesh, 'continuous') <- TRUE
    self$mesh$V <- self$V

    for (i in 1:length(self$edges)) {
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
  } else {
    self$mesh <- list(PtE = NULL,
                      V = NULL,
                      E = NULL,
                      n_e = NULL,
                      h_e = NULL,
                      ind = 0,
                      VtE = NULL)
    attr(self$mesh, 'continuous') <- FALSE
    for (i in 1:length(self$edges)) {
      if (is.null(n)) {
        #remove boundary points
        self$mesh$n_e[i] <- ceiling(self$edge_lengths[i] / h) + 1
      } else {
        self$mesh$n_e[i] <- n + 2
      }
      if (self$mesh$n_e[i] > 0) {
        d.e <- seq(from = 0, to = 1, length.out = self$mesh$n_e[i])

        self$mesh$PtE <- rbind(self$mesh$PtE, cbind(rep(i, self$mesh$n_e[i]),
                                                    d.e))

        self$mesh$h_e <- c(self$mesh$h_e,
                           rep(self$edge_lengths[i] * d.e[2],
                               self$mesh$n_e[i] - 1))

        V.int <- (max(self$mesh$ind) + 1):(max(self$mesh$ind) + self$mesh$n_e[i])
        self$mesh$ind <- c(self$mesh$ind, V.int)

        self$mesh$E <- rbind(self$mesh$E, cbind(V.int[-length(V.int)], V.int[-1]))

      } else {
        self$mesh$E <- rbind(self$mesh$E, self$E[i, ])
        self$mesh$h_e <- c(self$mesh$h_e,self$edge_lengths[i])
      }
    }
    self$mesh$VtE <- self$mesh$PtE
    self$mesh$V <- self$coordinates(PtE = self$mesh$PtE)
   if(merge.outs) {
     mesh_merge_outs(self)
   }
    move.V.first(self)
    mesh_merge_deg2(self)
  }
}


#debug
edge2 <- rbind(c(2,1), c(1,0))
edge3 <- rbind(c(2,-1), c(1,0))
edge4 <- rbind(c(3,0),c(2,-1))
edge7 <- rbind(c(3,0), c(2,1))
edges = list(edge2, edge3, edge4, edge7)


graph <- metric_graph$new(edges = edges)
build_mesh(graph,h=1, continuous = FALSE, merge.outs = TRUE)
graph$plot(mesh = TRUE, direction = TRUE)
v <- rep(0, dim(graph$mesh$V)[1]); v[1] = 1; graph$plot_function(v, plotly = TRUE)

find.mesh.bc <- function(graph) {
  if(attr(graph$mesh,"continuous")) {
    stop("mesh discontinuous")
  }
  starts <- which(graph$mesh$PtE[,2]==0)
  edges <- graph$mesh$PtE[starts,1]
  degrees <- graph$get_degrees()
  connections <-list()
  for(i in 1:length(starts)) {
    vert <- graph$E[edges[i],1] #the vertex of the start point
    deg.start <- degrees[vert]
    if(deg.start == 1) {
      connections[[i]] <- 0 #true starting value
    } else {
      #find edges ending in the vertex
      edge.in <- which(graph$E[,2] == vert)
      #find edges starting in the vertex
      edge.out <- which(graph$E[,1] == vert)

      mesh.in <- NULL
      for(j in 1:length(edge.in)){
        tmp <- which(graph$mesh$PtE[,1] == edge.in[j]) #mesh nodes on the edge
        mesh.in <- c(mesh.in, tmp[length(tmp)])
      }
      connections[[i]] <- mesh.in
    }
  }
  return(connections)
}

set.petrov.matrices <- function(graph) {
  graph$compute_fem(petrov=TRUE)
  if(attr(graph$mesh,"continuous")) {
    starts <- find.mesh.starts(graph)
    if(dim(graph$mesh$Cpet)[1] != dim(graph$mesh$Cpet)[2]+1){
      stop("Graph is not a tree")
    }
    if(length(starts)>1){

      G <- rbind(sparseMatrix(i=1:length(starts),j=starts,x=rep(0,length(starts)),dims=c(length(starts),dim(graph$mesh$Cpet)[1])),
                 t(graph$mesh$Gpet[,-starts[-1]]))
    } else {
      C <- rbind(sparseMatrix(i=1,j=starts,x=1,dims=c(1,dim(graph$mesh$Cpet)[1])),
                 t(graph$mesh$Cpet))
      G <- rbind(sparseMatrix(i=1,j=starts,x=0,dims=c(1,dim(graph$mesh$Cpet)[1])),
                 t(graph$mesh$Gpet))
    }
  } else {
    bc <- find.mesh.bc(graph)
    starts <- which(graph$mesh$PtE[,2]==0)
    C <- t(graph$mesh$Cpet)
    G <- t(graph$mesh$Gpet)

    for(i in length(bc):1) {
      G <- rbind(sparseMatrix(i=1,j=1,x=0,dims=c(1,dim(graph$mesh$Cpet)[1])), G)
      if(bc[[i]][1] == 0) {
        C <- rbind(sparseMatrix(i=1,j=starts[i],x=1,dims=c(1,dim(graph$mesh$Cpet)[1])), C)
      } else {
        C <- rbind(sparseMatrix(i = rep(1,length(bc[[i]])+1),
                                j = c(starts[i], bc[[i]]),
                                x = c(1, rep(-1/length(bc[[i]]), length(bc[[i]]))),
                                dims = c(1, dim(graph$mesh$Cpet)[1])), C)
      }

    }
  }
  return(list(C = C, G = G, n.bc = length(starts), h0 = which(unlist(lapply(bc,length))>1)))
}

#find one mesh node corresponding to each vertex and move it first
move.V.first <- function(graph) {
  nv <- dim(graph$mesh$V)[1]
  for(i in 1:graph$nV) {
    ind <- which(graph$mesh$V[,1] == graph$V[i,1] & graph$mesh$V[,2] == graph$V[i,2])[1]

    if(ind > i && i < nv) {
      if (i == 1) {
        reo <- c(ind, setdiff(i:nv,ind))
      } else {
        reo <- c(1:(i-1), ind, setdiff(i:nv,ind))
      }
      graph$mesh$V <- graph$mesh$V[reo,]
      graph$mesh$PtE <- graph$mesh$PtE[reo,]
      graph$mesh$VtE <- graph$mesh$VtE[reo,]
      Etmp <- graph$mesh$E
      ind1 <- Etmp == ind
      ind2 <- Etmp >= i & Etmp < ind
      graph$mesh$E[ind1] = i
      graph$mesh$E[ind2] = graph$mesh$E[ind2] + 1
    }
  }
}


edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(1,0),c(2,1))
edge3 <- rbind(c(1,0),c(2,-1))
edge4 <- rbind(c(2,-1),c(3,0))
edge5 <- rbind(c(2,-1),c(3,-2))
edges = list(edge1, edge2, edge3, edge4, edge5)
graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=0.05)
graph$plot(mesh=TRUE)

mat <- set.petrov.matrices(graph)
hfull <- c(1/(1+kappa),graph$mesh$h_e)

L <- kappa*mat$C + mat$G
Sigma <- solve(L,diag(hfull)%*%t(solve(L)))
p <- graph$plot(direction = TRUE)
graph$plot_function(diag(Sigma), p = p)




W <- rnorm(n=length(hfull),mean=0,sd = sqrt(hfull))
u <- solve(L,W)
graph$plot_function(u,plotly=TRUE)


#now test a reversed tree
edge1 <- rbind(c(1,0), c(0,0))
edge2 <- rbind(c(2,1), c(1,0))
edge3 <- rbind(c(2,-1), c(1,0))
edge4 <- rbind(c(3,0),c(2,-1))
edge5 <- rbind(c(3,-2),c(2,-1))
edge6 <- rbind(c(2,-1), c(1,-2))
edges = list(edge1, edge2, edge3, edge4, edge5)


graph <- metric_graph$new(edges = edges)
build_mesh(graph,h=0.05, continuous = FALSE, merge.outs = TRUE)
graph$plot(mesh = TRUE, direction = TRUE)

graph$compute_fem(petrov=TRUE)

kappa <- 1
mat <- set.petrov.matrices(graph)
hfull <- c(rep(1/(1+kappa), mat$n.bc),graph$mesh$h_e)
hfull[mat$h0] = 0

L <- kappa*mat$C + mat$G
Sigma <- solve(L,diag(hfull)%*%t(solve(L)))
p <- graph$plot(direction = TRUE)
graph$plot_function(diag(Sigma), p = p)

W <- rnorm(n=length(hfull),mean=0,sd = sqrt(hfull))
u <- solve(L,W)
graph$plot_function(u,plotly=TRUE)



#now test a reversed tree
edge1 <- rbind(c(1,0), c(0,0))
edge2 <- rbind(c(2,1), c(1,0))
edge3 <- rbind(c(2,-1), c(1,0))
edge4 <- rbind(c(3,0),c(2,-1))
edge5 <- rbind(c(3,-2),c(2,-1))
edge6 <- rbind(c(2,-1), c(1,-2))
edge7 <- rbind(c(3,0), c(2,1))
edges = list(edge1, edge2, edge3, edge4, edge5, edge6, edge7)


graph <- metric_graph$new(edges = edges)
build_mesh(graph,h=0.5, continuous = FALSE, merge.outs = TRUE)
graph$plot(mesh = TRUE, direction = TRUE)

v <- rep(0, dim(graph$mesh$V)[1]); v[23] = 1; graph$plot_function(v, plotly = TRUE)

graph$compute_fem(petrov=TRUE)

kappa <- 1
mat <- set.petrov.matrices(graph)
hfull <- c(rep(1/(1+kappa), mat$n.bc),graph$mesh$h_e)
hfull[mat$h0] = 0

L <- kappa*mat$C + mat$G
Sigma <- solve(L,diag(hfull)%*%t(solve(L)))
p <- graph$plot(direction = TRUE)
graph$plot_function(diag(Sigma), p = p)

W <- rnorm(n=length(hfull),mean=0,sd = sqrt(hfull))
u <- solve(L,W)
graph$plot_function(u,plotly=TRUE)


