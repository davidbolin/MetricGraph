
#'
#' The exponential covariance
#' @param D     - vector, or matrix distance
#' @param theta - kappa, sigma
#' @export
r_1 <- function(D,theta){
  return( ( theta[2]^2/(2 * theta[1]))* exp( -theta[1] * abs(D)))
}

#' plot the exponential covariance for parameter set
#' @param theta - sigma_e, kappa, sigma
#' @export
plot_r_1 <-function(theta, t = NULL){
  if(is.null(t)){
    r_0 <- r_1(0, theta[2:3])
    r_p <- -log(0.05)/theta[2]
    t <- seq(0, r_p,length.out = 100)
  }
  r_  <- r_1(t,theta[2:3])
  if(t[1]==0)
    r_[1] = r_[1] + theta[1]^2
  plot(t,r_,type='l')

}

#'
#' the expontial circular covariance
#' @param t     - locations
#' @param l_e   - length circle
#' @param theta - kappa, sigma
#' @export
r_1_circle <- function(t,l_e, theta){

  kappa <- theta[1]
  c = theta[2]^2 / ( 2 * kappa * sinh(kappa *l_e/2))
  r <- matrix(0, length(t), length(t))
  r_0 <- cosh(-kappa * l_e/2)
  if(length(t)==1){
    r[1,1] = c * r_0
    return(r)
  }

  for(i in 1:(length(t)-1)){

    for(ii in (i+1):length(t)){
      r[i, ii] <- cosh(kappa * ( abs( t[i] - t[ii]) - l_e/2 ))
    }
  }
  r <- r + t(r)
  diag(r) <- r_0
  return(c * r)
}


#'
#' Building the precision matrix for all the vertex in
#'the expontial case
#' @param theta - kappa, sigma
#' @param V vertex position
#' @param EtV [,2-3] index of upper and lower edge
#' @param El length of each vertex
#' @param BC boundary condition =0 neumann, 1 = correct boundary to stationary
#' @param build (bool) if true return Q the precision matrix otherwise return list(i,j,x, nv)
#' @return Q (precision matrix)
#' @export
Q.exp <- function(theta, V, EtV, El, BC = 1, build=T){

  kappa <- theta[1]
  sigma <- theta[2]
  i_ <- j_ <- x_ <- rep(0, dim(V)[1]*4)
  nE <- dim(EtV)[1]
  count <- 0
  for(i in 1:nE){
    l_e <- El[i]
    c1 <- exp(-kappa*l_e)
    c2 <- c1^2
    one_m_c2 = 1-c2
    c_1 = 0.5 + c2/one_m_c2
    c_2 = -c1/one_m_c2

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
  n.v <- dim(V)[1]
  if(BC==1){
    #does this work for circle?
    i.table <- table(i_[1:count])
    index = as.integer(names(which(i.table<3)))
    i_ <- c(i_[1:count], index)
    j_ <- c(j_[1:count], index)
    x_ <- c(x_[1:count], rep(0.5, length(index)))
    count <- count + length(index)
  }
  if(build){
    Q <- Matrix::sparseMatrix(i=i_[1:count],
                              j=j_[1:count],
                              x=(2*kappa / sigma^2)*x_[1:count],
                              dims=c(n.v, n.v))


    return(Q)
  }else{
    return(list(i = i_[1:count],
                j = j_[1:count],
                x  = (2*kappa / sigma^2)*x_[1:count],
                dims=c(n.v, n.v)))
  }
}
#'
#' Compute the precision matrix of observations one a line
#' That is the inverse of the exponential covariance at the observation location t
#'
#' @param theta      - kappa, sigma
#' @param t          - (n x 1) relative position on the line start with 0 end with 1
#' @param l_e        - (double) length of the line
#' @param t_sorted   - (bool)
#' @export
Q.exp.line <- function(theta, t,  t_sorted=FALSE){

  l_t = length(t)
  i_ <- j_ <- x_ <- rep(0, 4*(l_t-1) + 2)
  count = 0
  kappa <- theta[1]
  sigma <- theta[2]
  if(t_sorted==F){
    order_t  <- order(t)
    t        <- t[order_t]
    P <- sparseMatrix(seq_along(order_t), order_t, x=1)
  }

  for(i in 2:l_t){

    c1 <- exp(-kappa* (t[i] - t[i-1]))
    c2 <- c1^2
    one_m_c2 = 1-c2
    c_1 = 0.5 + c2/one_m_c2
    c_2 = -c1/one_m_c2

    i_[count + 1] <- i
    j_[count + 1] <- i - 1
    x_[count + 1] <- c_2

    i_[count + 2] <- i - 1
    j_[count + 2] <- i
    x_[count + 2] <- c_2


    i_[count + 3] <- i
    j_[count + 3] <- i
    x_[count + 3] <- c_1

    i_[count + 4] <- i - 1
    j_[count + 4] <- i - 1
    x_[count + 4] <- c_1
    count <- count + 4
  }
  #boundary term correction
  i_[count + 1] <- 1
  j_[count + 1] <- 1
  x_[count + 1] <- 0.5
  i_[count + 2] <- l_t
  j_[count + 2] <- l_t
  x_[count + 2] <- 0.5
  count <- count + 2
  Q <- Matrix::sparseMatrix(i=i_[1:count],
                            j=j_[1:count],
                            x=(2 * kappa / sigma^2)* x_[1:count],
                            dims=c(l_t, l_t))
  if(t_sorted==F)
    Q <- Matrix::t(P)%*%Q%*%P

  return(Q)
}


#' sample a line given end points
#' @param theta  - (3 x 1) sigma_e, sigma, kappa
#' @param u_e    - (2 x 1) the two end points
#' @param Line   - Spatial.datafram
#' @param l_e    - (1 x 1) line length
#' @param  t     - (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt    - (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py    - (n x 1) observation location
#' @param  y     - (n x 1) observations
#' @param  sample - (bool) if true sample else return mean
#' @export
sample.line.expontial<-function(theta, u_e, l_e, t=NULL, Line=NULL, nt=100,  py=NULL, y=NULL, sample=TRUE){

  if(is.null(t)){
    t  = seq(0,1,length.out=nt)
    #t <- rgeos::gProject(Line,rgeos::gInterpolate(Line, t, normalized = T))
    t <- t * l_e
  }
  t_end <- c(0, l_e)
  t <- unique(t)
  t0 <- t
  if(is.null(y)==F){
    ind_remove = which(py %in% t_end)
    if(length(ind_remove)>0){
      y  <- y[-ind_remove]
      py <- py[-ind_remove]
    }

    ind_remove_t <- which(t %in% c(t_end,py))
    if(length(ind_remove_t)>0)
      t <- t[-ind_remove_t]
  }else{
    ind_remove_t <- which(t %in% t_end)
    if(length(ind_remove_t)>0)
      t <- t[-ind_remove_t]
  }

  t <- c(t_end, t)
  if(is.null(py)==F){
    t <- c(py, t)
  }

  Q <- Q.exp.line(theta[2:3], t)
  Q <- Q
  index_E <- length(py) + 1:2
  Q_X <- Q[-index_E,-index_E]
  mu_X <- as.vector(Matrix::solve(Q_X,-Q[-index_E,index_E] %*% u_e))
  if(is.null(py)==F){
    Matrix::diag(Q_X)[1:length(py)] <-  Matrix::diag(Q_X)[1:length(py)] + 1/theta[1]^2
    AtY = rep(0,dim(Q_X)[1])
    AtY[1:length(py)] = (y - mu_X[1:length(py)])/theta[1]^2
    mu_X = mu_X + as.vector(Matrix::solve(Q_X,AtY))
  }

  x <- rep(0, length(t))

  if(sample){
    R_X <- Matrix::Cholesky(Q_X, LDL = FALSE, perm = TRUE)
    z <- rnorm(dim(R_X)[1])
    x[-index_E] <- mu_X + as.vector(Matrix::solve(R_X,Matrix::solve(R_X,z,system = 'Lt'), system='Pt'))
    x[index_E] <- u_e
  }else{
    x[-index_E] <- mu_X
    x[index_E] <- u_e

  }

  x_out <- matrix(0,nrow=length(t0),2)
  x_out[,1] <- t0
  for(i in 1:length(t0))
  {
    ind <- which(t==t0[i])
    x_out[i,2] <- x[ind]
  }
  return(x_out)

}

#'
#' Computes the posterior expectation for each node in the graph
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph.obj - graphical object
#' @param rem.edge  - remove edge
#' @export
posterior.mean.exp <- function(theta, graph.obj, rem.edge=NULL){
  sigma_e <- theta[1]
  #build Q
  Qp.list <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El, build=F)
  #build BSIGMAB
  Qpmu <- rep(0, nrow(graph.obj$V))

  obs.edges <- unique(graph.obj$PtE[,1])
  if(is.null(rem.edge)==F)
    obs.edges <- setdiff(obs.edges, rem.edge)
  i_ <- j_ <- x_ <- rep(0,4*length(obs.edges))
  count <- 0
  for(e in obs.edges){
    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$El[e]
    D_matrix <- as.matrix(dist(c(0,l,graph.obj$PtE[obs.id,2])))
    S <- r_1(D_matrix,theta[2:3])

    #covariance update see Art p.17
    E.ind         <- c(1:2)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])
    Sigma_i       <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    R <- chol(Sigma_i)
    Sigma_iB      <- solve(Sigma_i, t(Bt))
    BtSinvB       <- Bt %*% Sigma_iB

    E <- graph$EtV[e,2:3]
    if(E[1]==E[2]){
      Qpmu[E[1]]    <- Qpmu[E[1]]     +  sum(t(Sigma_iB)%*%y_i)
      Qp[E[1],E[1]] <-  Qp[E[1],E[1]] +  sum(Bt %*% Sigma_iB)
      i_[count+1] <- E[1]
      j_[count+1] <- E[1]
      x_[count+1] <- sum(Bt %*% Sigma_iB)
      count <- count + 1
    }else{
      i_[count+(1:4)] <- c(E[1],E[1],E[2],E[2])
      j_[count+(1:4)] <- c(E[1],E[2],E[1],E[2])
      x_[count+(1:4)] <- c(BtSinvB[1,1], BtSinvB[1,2], BtSinvB[1,2], BtSinvB[2,2] )
      count <- count + 4
      Qpmu[E]    <- Qpmu[E] + t(Sigma_iB)%*%y_i
    }
  }
  i_ <- c(Qp.list$i, i_[1:count])
  j_ <- c(Qp.list$j, j_[1:count])
  x_ <- c(Qp.list$x, x_[1:count])
  Qp <- Matrix::sparseMatrix(i= i_,
                             j= j_,
                             x= x_,
                             dims=Qp.list$dims)


  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Qpmu,system = 'P'), system='L')))
  Qpmu <- as.vector(Matrix::solve(R,Matrix::solve(R, v,system = 'Lt'), system='Pt'))

  return(Qpmu)

}


#'
#' Computes the posterior expectation for each observation in the graph
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#' @param leave.edge.out - compute the expectation of the graph if the observatrions are not on the edge
#' @export
posterior.mean.obs.exp <- function(theta, graph.obj, leave.edge.out = F){

  sigma_e = theta[1]

  Qp <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  if(leave.edge.out==F)
    V.post <- posterior.mean.exp(theta, graph)


  y_hat <- rep(0, length(graph$y))
  Qpmu <- rep(0, nrow(graph.obj$V))
  obs.edges <- unique(graph.obj$PtE[,1])

  for(e in obs.edges){

    if(leave.edge.out==T)
      V.post <- posterior.mean.exp(theta, graph, rem.edge = e)

    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$El[e]
    D_matrix <- as.matrix(dist(c(0,l,graph.obj$PtE[obs.id,2])))
    S <- r_1(D_matrix,theta[2:3])

    #covariance update see Art p.17
    E.ind         <- c(1:2)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])

    E <- graph$EtV[e,2:3]
    y_hat[obs.id] <- t(Bt)%*%V.post[E]
    if(leave.edge.out==F){
      Sigma_i       <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
      Sigma_noise  <- Sigma_i
      diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2

      y_hat[obs.id] <- y_hat[obs.id] +  Sigma_i%*%solve(Sigma_noise,y_i-y_hat[obs.id] )
    }
  }
  return(y_hat)
  #()
}


#' posterior mean calculation not using sparsity
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#'
#' @return The posterior mean
#' @export
posterior.mean.stupid <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  Sigma <- as.matrix(solve(Q))
  SigmaO <- Sigma[graph.obj$PtV,graph.obj$PtV]
  diag(SigmaO) <- diag(SigmaO)  +  sigma_e^2

  return( Sigma[,graph$PtV]%*%solve(SigmaO,graph.obj$y))
}

#' prediction not using sparsity
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#'
#' @return Leave-one-out predictions
#' @export
posterior.leave.stupid <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  Sigma <- as.matrix(solve(Q))
  SigmaO <- Sigma[graph.obj$PtV,graph.obj$PtV]
  diag(SigmaO) <- diag(SigmaO)  +  sigma_e^2
  y_p <- rep(0,length(graph.obj$y))
  for(i in 1:length(graph.obj$y)){
    mu_p <-  Sigma[,graph.obj$PtV[-i]]%*%solve(SigmaO[-i,-i],graph.obj$y[-i])
    y_p[i] <- mu_p[graph.obj$PtV[i]]
  }

  return( y_p)
}


#' Likelihood evaluation not using sparsity
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#' @return The log-likelihood
#' @export
likelihood.exp.graph.covariance <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  Sigma <- as.matrix(solve(Q))[graph.obj$PtV,graph.obj$PtV]
  diag(Sigma) <- diag(Sigma)  +  sigma_e^2
  R <- chol(Sigma)
  return(-sum(log(diag(R))) - 0.5*t(graph.obj$y)%*%solve(Sigma,graph.obj$y))
}

#' Log-likelihood calculation for alpha=1 without integration
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#'
#' @return The log-likelihood
#' @export
likelihood.exp.graph.v2 <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  n.v <- dim(graph.obj$V)[1]
  Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  A <- Diagonal(n.v,rep(1,n.v))[graph.obj$PtV,]
  Q.p <- Q  + t(A)%*%A/sigma_e^2
  mu.p <- solve(Q.p,as.vector(t(A)%*%graph.obj$y/sigma_e^2))

  R <- chol(Q)
  R.p <- chol(Q.p)

  l <- sum(log(diag(R))) - sum(log(diag(R.p))) - length(graph.obj$y)*log(sigma_e)
  v <- graph.obj$y - A%*%mu.p
  l <- l - 0.5*(t(mu.p)%*%Q%*%mu.p + t(v)%*%v/sigma_e^2)
  return(as.double(l))
}

#' Prediction for models assuming observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, sigma, kappa)
#' @param graph.obj Graph object
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL" gives
#' graph Laplacian, and "isoExp" gives isotropic exponential
#' @param ind Indices for cross validation. It should be a vector of the same length as the
#' data, with integer values representing each group in the cross-validation.
#' If NULL, leave-one-out cross validation is performed.
#'
#' @return Vector with all predictions
#' @export
posterior.crossvalidation <- function(theta,
                                      graph.obj,
                                      model = "alpha1",
                                      ind = NULL)
  {
  #setup matrices for prediction
  sigma_e <- theta[1]
  n.v <- dim(graph$V)[1]
  if(model == "isoExp"){
    Sigma <- theta[3]^2*exp(-theta[2]*graph.obj$res.dist)
    Sigma.o <- theta[3]^2*exp(-theta[2]*graph.obj$res.dist[graph.obj$PtV,graph.obj$PtV])
    diag(Sigma.o) <- diag(Sigma.o) + theta[1]^2
  } else if(model == "alpha1" || model == "GL"){
    if(model == "alpha1"){
      Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
    } else if(model == "GL"){
      Q <- theta[2]*(theta[3]*Diagonal(n.v,1) + graph.obj$Laplacian)
    }
    A <- Diagonal(n.v,rep(1,n.v))[graph.obj$PtV,]
    Q.p <- Q  + t(A)%*%A/sigma_e^2
  } else {
    stop("Wrong model choice.")
  }



  if(is.null(ind)){
    ind <- 1:length(graph.obj$y)
  }
  y_p <- rep(0,length(graph.obj$y))
  for(j in 1:length(unique(ind))){
    i <- which(ind == j)
    if(model == "isoExp"){
        y_p[i] <-Sigma[graph.obj$PtV[i],graph.obj$PtV[-i]]%*%solve(Sigma.o[-i,-i],graph.obj$y[-i])
    } else {
      A <- Diagonal(n.v,rep(1,n.v))[graph.obj$PtV[-i],]
      mu.p <- solve(Q  + t(A)%*%A/sigma_e^2, as.vector(t(A)%*%graph.obj$y[-i]/sigma_e^2))
      y_p[i] <- mu.p[graph.obj$PtV[i]]
    }
  }
  return( y_p)
}

#' Log-likelihood calculation for graph Laplacian model
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#'
#' @return The log-likelihood
#' @export
likelihood.graph_laplacian <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  n.v <- dim(graph.obj$V)[1]

  Q <- theta[2]*(theta[3]*Diagonal(n.v,1) + graph.obj$Laplacian)

  A <- Diagonal(n.v,rep(1,n.v))[graph.obj$PtV,]
  Q.p <- Q  + t(A)%*%A/sigma_e^2
  mu.p <- solve(Q.p,as.vector(t(A)%*%graph.obj$y/sigma_e^2))

  R <- chol(Q)
  R.p <- chol(Q.p)

  l <- sum(log(diag(R))) - sum(log(diag(R.p))) - length(graph$y)*log(sigma_e)
  v <- graph$y - A%*%mu.p
  l <- l - 0.5*(t(mu.p)%*%Q%*%mu.p + t(v)%*%v/sigma_e^2)
  return(as.double(l))
}


#'
#' Computes the log likelihood function fo theta for the graph object
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph.obj - graph object
#' @export
likelihood.exp.graph <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  #Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)

  Q.list <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El, build=F)

  Qp <- Matrix::sparseMatrix(i = Q.list$i,
                             j = Q.list$j,
                             x = Q.list$x,
                             dims=Q.list$dims)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  loglik <- 0.5*Matrix::determinant(R)$modulus[1]

  #build BSIGMAB
  Qpmu <- rep(0, nrow(graph.obj$V))
  obs.edges <- unique(graph.obj$PtE[,1])

  i_ <- j_ <- x_ <- rep(0,4*length(obs.edges))
  count <- 0
  for(e in obs.edges){
    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$El[e]
    D_matrix <- as.matrix(dist(c(0,l,graph.obj$PtE[obs.id,2])))
    S <- r_1(D_matrix,theta[2:3])

    #covariance update see Art p.17
    E.ind         <- c(1:2)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind,drop=F],S[E.ind, Obs.ind,drop=F])
    Sigma_i       <- S[Obs.ind,Obs.ind,drop=F] - S[Obs.ind, E.ind,drop=F] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    R <- chol(Sigma_i)
    Sigma_iB      <- solve(Sigma_i, t(Bt))
    BtSinvB       <- Bt %*% Sigma_iB

    E <- graph$EtV[e,2:3]
    if(E[1]==E[2]){
      Qpmu[E[1]]    <- Qpmu[E[1]]     +  sum(t(Sigma_iB)%*%y_i)
      i_[count+1] <- E[1]
      j_[count+1] <- E[1]
      x_[count+1] <- sum(Bt %*% Sigma_iB)
    }else{
      Qpmu[E]    <- Qpmu[E] + t(Sigma_iB)%*%y_i
      i_[count+(1:4)] <- c(E[1],E[1],E[2],E[2])
      j_[count+(1:4)] <- c(E[1],E[2],E[1],E[2])
      x_[count+(1:4)] <- c(BtSinvB[1,1], BtSinvB[1,2], BtSinvB[1,2], BtSinvB[2,2] )
      count <- count + 4
    }
    loglik <- loglik - 0.5  * t(y_i)%*%solve(Sigma_i,y_i)
    loglik <- loglik - sum(log(diag(R)))
  }

  i_ <- c(Q.list$i, i_[1:count])
  j_ <- c(Q.list$j, j_[1:count])
  x_ <- c(Q.list$x, x_[1:count])
  Qp <- Matrix::sparseMatrix(i= i_,
                             j= j_,
                             x= x_,
                             dims=Q.list$dims)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  loglik <- loglik - 0.5*Matrix::determinant(R)$modulus[1]

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Qpmu,system = 'P'), system='L')))
  #Qpmu <- as.vector(solve(R,solve(R, v,system = 'Lt'), system='Pt'))

  loglik <- loglik + 0.5  * t(v)%*%v
  return(loglik[1])
}

