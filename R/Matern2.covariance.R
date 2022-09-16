
#'
#' The Matern alpha=2 covariance
#' @param D      - vector, or matrix distance
#' @param theta  - kappa, sigma
#' @param deriv  - (0,1,2) no derivative, first, or second
#' @export
r_2 <- function(D,theta, deriv = 0){
  kappa <- theta[1]
  sigma <- theta[2]
  aD <- abs(D)
  c = ( sigma^2/(4 * kappa^3))

  R0 =  exp( -kappa * aD)
  if(deriv==0)
    return(  c * (1 + kappa * aD)* R0)


  d1 = -kappa^2 * c * D* R0

  if(deriv==1)
    return( d1)
  if(deriv==2)
    return(kappa^2 * c *( kappa* aD - 1)*R0)
  error("deriv must be either (0,1,2)")

}



#' plot the Matern alpha=2 covariance for parameter set
#' @param theta - sigma_e, kappa, sigma
#' @export
plot_r_2 <-function(theta, t = NULL){
  if(is.null(t)){
    r_p <- 4.743859/theta[2]
    t <- seq(0, r_p,length.out = 100)
  }
  r_  <- r_2(t,theta[2:3])
  if(t[1]==0)
    r_[1] = r_[1] + theta[1]^2
  plot(t,r_,type='l')

}


#' sample a line given end points
#' @param theta  - (3 x 1) sigma_e, sigma, kappa
#' @param u_e    - (4 x 1) X(0),X'(0), X(l_e),X'(l_e) the two end points and there deriv
#' @param Line   - Spatial.datafram
#' @param l_e    - (1 x 1) line length
#' @param  t     - (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt    - (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py    - (n x 1) observation location
#' @param  y     - (n x 1) observations
#' @param  sample - (bool) if true sample else return mean
#' @export
sample.line.matern2 <-function(theta, u_e, l_e, t=NULL, Line=NULL, nt=100,  py=NULL, y=NULL, sample=TRUE){

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
  Sigma <- matrix(0, length(t)+2, length(t)+2)
  d.index <- c(1,2)
  index_E <- 2 + length(py) + 1:2
  D <- outer (t, t, `-`)
  Sigma[-d.index, -d.index] <- r_2(D,         theta[2:3])
  Sigma[ d.index,  d.index] <- -r_2(as.matrix(dist(c(0,l_e))), theta[2:3], 2)
  Sigma[d.index,  -d.index] <- -r_2(D[index_E-2,],             theta[2:3], 1)
  Sigma[-d.index,  d.index] <- t(Sigma[d.index,  -d.index])


  index_boundary <- c(d.index,index_E)
  u_e <- u_e[c(2,4,1,3)] # aragne in X',X
  SinvS <- solve(Sigma[index_boundary,index_boundary],Sigma[index_boundary, -index_boundary] )
  Sigma_X <- Sigma[-index_boundary,-index_boundary] - Sigma[-index_boundary, index_boundary]%*%SinvS
  mu_X <- - t(SinvS)%*%(0-u_e)

  if(is.null(py)==F){
    index_y <- 1:length(py)
    Sigma_Y <- Sigma_X[index_y,index_y,drop=F]
    Matrix::diag(Sigma_Y) <-  Matrix::diag(Sigma_Y) + theta[1]^2

    SinvS <- solve(Sigma_Y,Sigma_X[index_y,,drop=F] )
    Sigma_X <- Sigma_X - Sigma_X[,index_y,drop=F]%*%SinvS
    mu_X = mu_X + t(SinvS)%*%(y- mu_X[index_y])
  }


  x <- rep(0, length(t))

  if(sample){
    R_X <- chol(Sigma_X)
    z <- rnorm(dim(Sigma_X)[1])
    x[-c(1:2)] <- mu_X + t(R_X)%*%z
    x[c(1:2)] <- u_e[c(3,4)]
  }else{
    x[-c(1:2)] <- mu_X
    x[1:2] <- u_e[c(3,4)]

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
#' Building the precision matrix for all the edges in the matern alpha=2 case
#' The ordering is acording to EtV where for each edge there are four random variables
#' processes and derivate for lower and upper edge end point
#' @param theta - kappa, sigma
#' @param V vertex position
#' @param EtV [,2-3] index of upper and lower edge
#' @param El length of each vertex
#' @param BC boundary condition =0 neumann, 1 = correct boundary to stationary
#' @param build (bool) if true return Q the precision matrix otherwise return list(i,j,x, nv)
#' @return Q (precision matrix)
#' @export
Q.matern2 <- function(theta, V, EtV, El, BC = 1, build=T){

  kappa <- theta[1]
  sigma <- theta[2]
  i_ <- j_ <- x_ <- rep(0, dim(EtV)[1]*16)
  nE <- dim(EtV)[1]
  count <- 0

  R_00 <- matrix(c( r_2(0, c(kappa,sigma), 0),
                   -r_2(0, c(kappa,sigma), 1), #0
                   -r_2(0, c(kappa,sigma), 1), #0
                   -r_2(0, c(kappa,sigma), 2)), 2, 2)
  R_node <- rbind(cbind(R_00, matrix(0,2,2)),
                  cbind( matrix(0,2,2),R_00))
  Ajd <- -0.5 * solve(rbind(cbind(R_00, matrix(0,2,2)), cbind(matrix(0,2,2),R_00)))
  for(i in 1:nE){

    l_e <- El[i]
    #lots of redundant caculations
    d_ <- c(0,l_e)
    D <- outer(d_,d_,"-")
    r_0l <-   r_2(l_e, c(kappa,sigma))
    r_11 <- - r_2(l_e, c(kappa,sigma),2)
    # order by node not derivative
    R_01 <- matrix(c(r_0l, r_2(-l_e, c(kappa,sigma),1), r_2(l_e, c(kappa,sigma),1), r_11),2,2)
    #R_node <- rbind(cbind(R_00, R_01), cbind(t(R_01),R_00))
    R_node[1:2,3:4] <- R_01
    R_node[3:4,1:2] <- t(R_01)
    Q_adj = solve(R_node) + Ajd

    if(EtV[i,2]!=EtV[i,3]){

      #lower edge precision u
      i_[count + 1] <- 4*(i-1) + 1
      j_[count + 1] <- 4*(i-1) + 1
      x_[count + 1] <- Q_adj[1, 1]

      #lower edge  u'
      i_[count + 2] <- 4*(i-1) + 2
      j_[count + 2] <- 4*(i-1) + 2
      x_[count + 2] <- Q_adj[2, 2]

      #upper edge  u
      i_[count + 3] <- 4*(i-1) + 3
      j_[count + 3] <- 4*(i-1) + 3
      x_[count + 3] <- Q_adj[3, 3]

      #upper edge  u'
      i_[count + 4] <- 4*(i-1) + 4
      j_[count + 4] <- 4*(i-1) + 4
      x_[count + 4] <- Q_adj[4, 4]

      #lower edge  u, u'
      i_[count + 5] <- 4*(i-1) + 1
      j_[count + 5] <- 4*(i-1) + 2
      x_[count + 5] <- Q_adj[1, 2]
      i_[count + 6] <- 4*(i-1) + 2
      j_[count + 6] <- 4*(i-1) + 1
      x_[count + 6] <- Q_adj[1, 2]

      #upper edge  u, u'
      i_[count + 7] <- 4*(i-1) + 3
      j_[count + 7] <- 4*(i-1) + 4
      x_[count + 7] <- Q_adj[3, 4]
      i_[count + 8] <- 4*(i-1) + 4
      j_[count + 8] <- 4*(i-1) + 3
      x_[count + 8] <- Q_adj[3, 4]

      #lower edge  u, upper edge  u,
      i_[count + 9]  <- 4*(i-1) + 1
      j_[count + 9]  <- 4*(i-1) + 3
      x_[count + 9]  <- Q_adj[1, 3]
      i_[count + 10] <- 4*(i-1) + 3
      j_[count + 10] <- 4*(i-1) + 1
      x_[count + 10] <- Q_adj[1, 3]

      #lower edge  u, upper edge  u',
      i_[count + 11] <- 4*(i-1) + 1
      j_[count + 11] <- 4*(i-1) + 4
      x_[count + 11] <- Q_adj[1, 4]
      i_[count + 12] <- 4*(i-1) + 4
      j_[count + 12] <- 4*(i-1) + 1
      x_[count + 12] <- Q_adj[1, 4]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4*(i-1) + 2
      j_[count + 13] <- 4*(i-1) + 3
      x_[count + 13] <- Q_adj[2, 3]
      i_[count + 14] <- 4*(i-1) + 3
      j_[count + 14] <- 4*(i-1) + 2
      x_[count + 14] <- Q_adj[2, 3]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4*(i-1) + 2
      j_[count + 13] <- 4*(i-1) + 3
      x_[count + 13] <- Q_adj[2, 3]
      i_[count + 14] <- 4*(i-1) + 3
      j_[count + 14] <- 4*(i-1) + 2
      x_[count + 14] <- Q_adj[2, 3]

      #lower edge  u', upper edge  u',
      i_[count + 15] <- 4*(i-1) + 2
      j_[count + 15] <- 4*(i-1) + 4
      x_[count + 15] <- Q_adj[2, 4]
      i_[count + 16] <- 4*(i-1) + 4
      j_[count + 16] <- 4*(i-1) + 2
      x_[count + 16] <- Q_adj[2, 4]



      count <- count + 16
    }else{
        error("circle is implimented")
    }
  }

  if(BC==1){
    #Vertices with of degree 1
    i.table <- table(c(EtV[,2:3]))
    index = as.integer(names(which(i.table==1)))
    #for this vertices locate position
    lower.edges  = which(EtV[,2]%in%index)
    upper.edges  = which(EtV[,3]%in%index)
    for(le in lower.edges){
      ind <- c(4 * (le-1) + 1, 4 * (le-1) + 2)

      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5*c(1/R_00[1,1], 1/R_00[2,2]))
      count <- count + 2
    }
    for(ue in upper.edges){
      ind <- c(4 * (ue-1) + 3, 4 * (ue-1) + 4)
      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5*c(1/R_00[1,1], 1/R_00[2,2]))
      count <- count + 2
    }
  }
  if(build){
    Q <- Matrix::sparseMatrix(i    = i_[1:count],
                              j    = j_[1:count],
                              x    = x_[1:count],
                              dims = c(4*nE, 4*nE))


    return(Q)
  }else{
    return(list(i = i_[1:count],
                j = j_[1:count],
                x  = x_[1:count],
                dims=c(4*n.E, 4*n.E)))
  }
}


#'
#' Computes the log likelihood function fo theta for the graph object
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph.obj - graphical object
#' @export
likelihood.matern2.graph <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q

  n_const <- length(graph$CBobj$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CBobj$T[-ind.const,]
  Q <- Q.matern2(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  R <- Matrix::Cholesky(Tc%*%Q%*%t(Tc), LDL = FALSE, perm = TRUE)
  loglik <- 0.5*Matrix::determinant(R)$modulus[1]

  #build BSIGMAB
  Qpmu      <- rep(0, 4*nrow(graph.obj$EtV))
  obs.edges <- unique(graph.obj$PtE[,1])

  i_ <- j_ <- x_ <- rep(0,16*length(obs.edges))
  count <- 0
  for(e in obs.edges){
    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$El[e]
    t <- c(0,l,graph.obj$PtE[obs.id,2])

    D <- outer (t, t, `-`)
    S <- matrix(0, length(t)+2, length(t)+2 )

    d.index <- c(1,2)
    S[-d.index, -d.index] <- r_2(D,                          theta[2:3])
    S[ d.index,  d.index] <- -r_2(as.matrix(dist(c(0,l))),   theta[2:3], 2)
    S[d.index,  -d.index] <- -r_2(D[1:2,],                   theta[2:3], 1)
    S[-d.index,  d.index] <- t(S[d.index,  -d.index])

    #covariance update see Art p.17
    E.ind         <- c(1:4)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind,drop=F])
    Sigma_i       <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

    R <- chol(Sigma_i, pivot=T)
    if(attr(R,"rank") < dim(R)[1])
      return(-Inf)
    #Sigma_iB      <- solve(Sigma_i, t(Bt))
    Sigma_iB <- t(Bt)
    Sigma_iB[attr(R,"pivot"),] <- forwardsolve(R, backsolve(R, t(Bt[,attr(R,"pivot")]), transpose = TRUE), upper.tri = TRUE)
    BtSinvB       <- Bt %*% Sigma_iB

    E <- graph$EtV[e,2:3]
    if(E[1]==E[2]){
      error("circle not implemented")
    }else{
      BtSinvB <- BtSinvB[c(3,1,4,2), c(3,1,4,2)]
      Qpmu[4*(e-1)+1:4]    <- Qpmu[4*(e-1)+1:4] + (t(Sigma_iB)%*%y_i)[c(3,1,4,2)]

      #lower edge precision u
      i_[count + 1] <- 4*(e-1) + 1
      j_[count + 1] <- 4*(e-1) + 1
      x_[count + 1] <- BtSinvB[1, 1]

      #lower edge  u'
      i_[count + 2] <- 4*(e-1) + 2
      j_[count + 2] <- 4*(e-1) + 2
      x_[count + 2] <- BtSinvB[2, 2]

      #upper edge  u
      i_[count + 3] <- 4*(e-1) + 3
      j_[count + 3] <- 4*(e-1) + 3
      x_[count + 3] <- BtSinvB[3, 3]

      #upper edge  u'
      i_[count + 4] <- 4*(e-1) + 4
      j_[count + 4] <- 4*(e-1) + 4
      x_[count + 4] <- BtSinvB[4, 4]

      #lower edge  u, u'
      i_[count + 5] <- 4*(e-1) + 1
      j_[count + 5] <- 4*(e-1) + 2
      x_[count + 5] <- BtSinvB[1, 2]
      i_[count + 6] <- 4*(e-1) + 2
      j_[count + 6] <- 4*(e-1) + 1
      x_[count + 6] <- BtSinvB[1, 2]

      #upper edge  u, u'
      i_[count + 7] <- 4*(e-1) + 3
      j_[count + 7] <- 4*(e-1) + 4
      x_[count + 7] <- BtSinvB[3, 4]
      i_[count + 8] <- 4*(e-1) + 4
      j_[count + 8] <- 4*(e-1) + 3
      x_[count + 8] <- BtSinvB[3, 4]

      #lower edge  u, upper edge  u,
      i_[count + 9]  <- 4*(e-1) + 1
      j_[count + 9]  <- 4*(e-1) + 3
      x_[count + 9]  <- BtSinvB[1, 3]
      i_[count + 10] <- 4*(e-1) + 3
      j_[count + 10] <- 4*(e-1) + 1
      x_[count + 10] <- BtSinvB[1, 3]

      #lower edge  u, upper edge  u',
      i_[count + 11] <- 4*(e-1) + 1
      j_[count + 11] <- 4*(e-1) + 4
      x_[count + 11] <- BtSinvB[1, 4]
      i_[count + 12] <- 4*(e-1) + 4
      j_[count + 12] <- 4*(e-1) + 1
      x_[count + 12] <- BtSinvB[1, 4]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4*(e-1) + 2
      j_[count + 13] <- 4*(e-1) + 3
      x_[count + 13] <- BtSinvB[2, 3]
      i_[count + 14] <- 4*(e-1) + 3
      j_[count + 14] <- 4*(e-1) + 2
      x_[count + 14] <- BtSinvB[2, 3]

      #lower edge  u', upper edge  u',
      i_[count + 15] <- 4*(e-1) + 2
      j_[count + 15] <- 4*(e-1) + 4
      x_[count + 15] <- BtSinvB[2, 4]
      i_[count + 16] <- 4*(e-1) + 4
      j_[count + 16] <- 4*(e-1) + 2
      x_[count + 16] <- BtSinvB[2, 4]

      count <- count + 16
    }
    loglik <- loglik - 0.5  * t(y_i)%*%solve(Sigma_i,y_i)
    loglik <- loglik - sum(log(diag(R)))
  }
  i_ <- i_[1:count]
  j_ <- j_[1:count]
  x_ <- x_[1:count]
  BtSB <- Matrix::sparseMatrix(i= i_,
                               j= j_,
                               x= x_,
                               dims=dim(Q))
  Qp <- Q + BtSB
  Qp <- Tc%*%Qp%*%t(Tc)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  loglik <- loglik - 0.5*Matrix::determinant(R)$modulus[1]

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Tc%*%Qpmu,system = 'P'), system='L')))
  #Qpmu <- as.vector(solve(R,solve(R, v,system = 'Lt'), system='Pt'))

  loglik <- loglik + 0.5  * t(v)%*%v
  #prior
  loglik <- loglik - (2)*log(theta[2]) - 0.5/theta[2]

  return(loglik[1])
}



#'
#' Computes the posterior mean
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph.obj - graphical object
#' @export
posterior.mean.matern2 <- function(theta, graph.obj, rem.edge=NULL){
  sigma_e <- theta[1]
  #build Q

  n_const <- length(graph$CBobj$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CBobj$T[-ind.const,]
  Q <- Q.matern2(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)

  #build BSIGMAB
  Qpmu      <- rep(0, 4*nrow(graph.obj$EtV))
  obs.edges <- unique(graph.obj$PtE[,1])
  if(is.null(rem.edge)==F)
    obs.edges <- setdiff(obs.edges, rem.edge)

  i_ <- j_ <- x_ <- rep(0,16*length(obs.edges))
  count <- 0
  for(e in obs.edges){
    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$El[e]
    t <- c(0,l,graph.obj$PtE[obs.id,2])

    D <- outer (t, t, `-`)
    S <- matrix(0, length(t)+2, length(t)+2 )

    d.index <- c(1,2)
    S[-d.index, -d.index] <- r_2(D,                          theta[2:3])
    S[ d.index,  d.index] <- -r_2(as.matrix(dist(c(0,l))),   theta[2:3], 2)
    S[d.index,  -d.index] <- -r_2(D[1:2,],                   theta[2:3], 1)
    S[-d.index,  d.index] <- t(S[d.index,  -d.index])

    #covariance update see Art p.17
    E.ind         <- c(1:4)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind,drop=F])
    Sigma_i       <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

    R <- chol(Sigma_i, pivot=T)
    #Sigma_iB      <- solve(Sigma_i, t(Bt))
    Sigma_iB <- t(Bt)
    Sigma_iB[attr(R,"pivot"),] <- forwardsolve(R, backsolve(R, t(Bt[,attr(R,"pivot")]), transpose = TRUE), upper.tri = TRUE)
    BtSinvB       <- Bt %*% Sigma_iB

    E <- graph$EtV[e,2:3]
    if(E[1]==E[2]){
      error("circle not implemented")
    }else{
      BtSinvB <- BtSinvB[c(3,1,4,2), c(3,1,4,2)]
      Qpmu[4*(e-1)+1:4]    <- Qpmu[4*(e-1)+1:4] + (t(Sigma_iB)%*%y_i)[c(3,1,4,2)]

      #lower edge precision u
      i_[count + 1] <- 4*(e-1) + 1
      j_[count + 1] <- 4*(e-1) + 1
      x_[count + 1] <- BtSinvB[1, 1]

      #lower edge  u'
      i_[count + 2] <- 4*(e-1) + 2
      j_[count + 2] <- 4*(e-1) + 2
      x_[count + 2] <- BtSinvB[2, 2]

      #upper edge  u
      i_[count + 3] <- 4*(e-1) + 3
      j_[count + 3] <- 4*(e-1) + 3
      x_[count + 3] <- BtSinvB[3, 3]

      #upper edge  u'
      i_[count + 4] <- 4*(e-1) + 4
      j_[count + 4] <- 4*(e-1) + 4
      x_[count + 4] <- BtSinvB[4, 4]

      #lower edge  u, u'
      i_[count + 5] <- 4*(e-1) + 1
      j_[count + 5] <- 4*(e-1) + 2
      x_[count + 5] <- BtSinvB[1, 2]
      i_[count + 6] <- 4*(e-1) + 2
      j_[count + 6] <- 4*(e-1) + 1
      x_[count + 6] <- BtSinvB[1, 2]

      #upper edge  u, u'
      i_[count + 7] <- 4*(e-1) + 3
      j_[count + 7] <- 4*(e-1) + 4
      x_[count + 7] <- BtSinvB[3, 4]
      i_[count + 8] <- 4*(e-1) + 4
      j_[count + 8] <- 4*(e-1) + 3
      x_[count + 8] <- BtSinvB[3, 4]

      #lower edge  u, upper edge  u,
      i_[count + 9]  <- 4*(e-1) + 1
      j_[count + 9]  <- 4*(e-1) + 3
      x_[count + 9]  <- BtSinvB[1, 3]
      i_[count + 10] <- 4*(e-1) + 3
      j_[count + 10] <- 4*(e-1) + 1
      x_[count + 10] <- BtSinvB[1, 3]

      #lower edge  u, upper edge  u',
      i_[count + 11] <- 4*(e-1) + 1
      j_[count + 11] <- 4*(e-1) + 4
      x_[count + 11] <- BtSinvB[1, 4]
      i_[count + 12] <- 4*(e-1) + 4
      j_[count + 12] <- 4*(e-1) + 1
      x_[count + 12] <- BtSinvB[1, 4]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4*(e-1) + 2
      j_[count + 13] <- 4*(e-1) + 3
      x_[count + 13] <- BtSinvB[2, 3]
      i_[count + 14] <- 4*(e-1) + 3
      j_[count + 14] <- 4*(e-1) + 2
      x_[count + 14] <- BtSinvB[2, 3]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4*(e-1) + 2
      j_[count + 13] <- 4*(e-1) + 3
      x_[count + 13] <- BtSinvB[2, 3]
      i_[count + 14] <- 4*(e-1) + 3
      j_[count + 14] <- 4*(e-1) + 2
      x_[count + 14] <- BtSinvB[2, 3]

      #lower edge  u', upper edge  u',
      i_[count + 15] <- 4*(e-1) + 2
      j_[count + 15] <- 4*(e-1) + 4
      x_[count + 15] <- BtSinvB[2, 4]
      i_[count + 16] <- 4*(e-1) + 4
      j_[count + 16] <- 4*(e-1) + 2
      x_[count + 16] <- BtSinvB[2, 4]

      count <- count + 16
    }
  }
  i_ <- i_[1:count]
  j_ <- j_[1:count]
  x_ <- x_[1:count]
  BtSB <- Matrix::sparseMatrix(i= i_,
                               j= j_,
                               x= x_,
                               dims=dim(Q))
  Qp <- Q + BtSB
  Qp <- Tc%*%Qp%*%t(Tc)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Tc%*%Qpmu,system = 'P'), system='L')))
  Qpmu <- as.vector(Matrix::solve(R,Matrix::solve(R, v,system = 'Lt'), system='Pt'))


  return(t(Tc)%*%Qpmu)
}


#'
#' Computes the posterior expectation for each observation in the graph
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graphical object
#' @param leave.edge.out - compute the expectation of the graph if the observatrions are not on the edge
#' @export
posterior.mean.obs.matern2 <- function(theta, graph.obj, leave.edge.out = F){

  sigma_e = theta[1]

  Qp <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  if(leave.edge.out==F)
    E.post <- posterior.mean.matern2(theta, graph)


  y_hat <- rep(0, length(graph$y))
  Qpmu <- rep(0, nrow(graph.obj$V))
  obs.edges <- unique(graph.obj$PtE[,1])

  for(e in obs.edges){

    if(leave.edge.out==T)
      E.post <- posterior.mean.matern2(theta, graph, rem.edge = e)

    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$El[e]
    t <- c(0,l,graph.obj$PtE[obs.id,2])
    D <- outer (t, t, `-`)
    S <- matrix(0, length(t)+2, length(t)+2 )

    d.index <- c(1,2)
    S[-d.index, -d.index] <- r_2(D,                          theta[2:3])
    S[ d.index,  d.index] <- -r_2(as.matrix(dist(c(0,l))),   theta[2:3], 2)
    S[d.index,  -d.index] <- -r_2(D[1:2,],                   theta[2:3], 1)
    S[-d.index,  d.index] <- t(S[d.index,  -d.index])



    #covariance update see Art p.17
    E.ind         <- c(1:4)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])




    u_e <- E.post[4*(e-1) + c(2,4,1,2)]
    y_hat[obs.id] <- t(Bt)%*%u_e
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


