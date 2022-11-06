#' The Matern covariance with alpha=2
#' @param D vector or matrix with distances
#' @param theta kappa, sigma
#' @param deriv (0,1,2) no derivative, first, or second order
#' @export
r_2 <- function(D, theta, deriv = 0){
  kappa <- theta[1]
  sigma <- theta[2]
  aD <- abs(D)
  c <- ( sigma^2/(4 * kappa^3))

  R0 <-  exp( -kappa * aD)
  if (deriv == 0)
    return( c * (1 + kappa * aD) * R0)


  d1 <- -kappa^2 * c * D * R0

  if (deriv == 1)
    return( d1)
  if (deriv == 2)
    return(kappa^2 * c * ( kappa* aD - 1) * R0)
  stop("deriv must be either (0,1,2)")

}

#' plot the Matern alpha=2 covariance for parameter set
#' @param theta - sigma_e, kappa, sigma
#' @export
plot_r_2 <-function(theta, t = NULL) {
  if (is.null(t)) {
    r_p <- 4.743859 / theta[2]
    t <- seq(0, r_p, length.out = 100)
  }
  r_  <- r_2(t,theta[2:3])
  if (t[1] == 0)
    r_[1] = r_[1] + theta[1]^2
  plot(t, r_, type = "l")

}


#' sample a line given end points
#' @param theta (3 x 1) sigma_e, sigma, kappa
#' @param u_e (4 x 1) X(0),X'(0), X(l_e),X'(l_e) the process and its derivative
#' at the two end points
#' @param Line Spatial data frame
#' @param l_e (1 x 1) line length
#' @param  t (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py (n x 1) observation location
#' @param  y (n x 1) observations
#' @param  sample (bool) if true sample else return mean
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

    if(length(py)==0)
      py <- NULL
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
  u_e <- u_e[c(2,4,1,3)] # arrange in X',X
  SinvS <- solve(Sigma[index_boundary, index_boundary], Sigma[index_boundary, -index_boundary])
  Sigma_X <- Sigma[-index_boundary, -index_boundary] - Sigma[-index_boundary, index_boundary]%*%SinvS
  mu_X <- - t(SinvS) %*% (0-u_e)

  if(is.null(py)==F){
    index_y <- 1:length(py)
    Sigma_Y <- Sigma_X[index_y, index_y, drop=FALSE]
    Matrix::diag(Sigma_Y) <-  Matrix::diag(Sigma_Y) + theta[1]^2

    SinvS <- solve(Sigma_Y, Sigma_X[index_y,,drop = FALSE])
    Sigma_X <- Sigma_X - Sigma_X[,index_y, drop = FALSE] %*% SinvS
    mu_X <- mu_X + t(SinvS) %*% (y- mu_X[index_y])
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





#' Generates samples of the entire graph
#' @param graph   - graphical object
#' @param theta  - (sigma_e, sigma, kappa)
#' @export
graph_posterior_mean_matern2 <- function(graph,  theta, sample=F){

  X <- c()
  V.post.mean <- posterior.mean.matern2(theta, graph)
  for(i in 1:dim(graph$EtV)[1]){
    V.i <-   V.post.mean[4*(i-1) + 1:4]

    ind <- which(graph$PtE[,1] == i)
    if(length(ind)==0){
      X.i   <- sample.line.matern2(theta,
                           V.i,
                           graph$edge_lengths[i],
                           nt = 100,
                           sample=sample)
    }else{
      X.i   <- sample.line.matern2(theta,
                           V.i,
                           graph$edge_lengths[i],
                           nt = 100,
                           y = graph$y[ind],
                           py =graph$PtE[ind,2],
                           sample=sample)

    }
    X.i[,1] <- X.i[,1]/graph$edge_lengths[i]
    X <- rbind(X, cbind(X.i,i))

  }
  return(X)

}

#'
#' Compute covariance of a point to the entire graph (discretized) for
#' Alpha = 2
#'
#' @param  EP    - (2 x 1) [1] -edge number, [2] -normalized location on the edge
#' @param  theta - (3 x 1) (kappa,sigma)
#' @param  graph  - (graph)
#' @param  n.p    - (int) number of points to compute the covariance on on each edge
#' @return C      - (n.p*numer of edges x 3) [1] edge number [2] lenth from lower edge [3] covarians
#' @export
covariance.point.to.graph.alpha2 <- function(EP, theta, graph, n.p = 50){


  kappa <- theta[1]
  sigma <- theta[2]
  #compute covarains of the two edges of EP[1]
  Q <- Q.matern2(theta, graph$V, graph$EtV, graph$edge_lengths)
  if(is.null(graph$CBobj))
    graph$buildA(2, F)

  n_const <- length(graph$CBobj$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CBobj$T[-ind.const,]
  Q_mod <- Tc%*%Q%*%t(Tc)
  R <- Cholesky(Q_mod, LDL = FALSE, perm = TRUE)
  Vs <- graph$EtV[EP[1],2:3]


  #COV[X,Y] = cov[Xtilde+BZ,Y] = B Cov[Z,Y]
  #the four Z is u(0),u'(0), u(T), u'(T) for each edge
  Z <- matrix(0, nrow=4*dim(graph$EtV)[1], ncol=4)
  Z[cbind(4*(EP[1]-1) + 1:4,1:4 )] = 1
  TZ = Tc%*%Z
  V  <- Matrix::solve(R,Matrix::solve(R,TZ,system = 'P'), system='L')
  TCV <- Matrix::solve(R,Matrix::solve(R,V,system = 'Lt'), system='Pt')
  CZ <- t(Tc)%*%TCV

  #modfing  u(0),u'(0),u(T),u'(T)
  CZ <- CZ[,c(2,4,1,3)]

  # compute covarains between two u(0),u'(0),u(T),u'(T) and the point u(p)
  t_norm <- EP[2]
  l <- graph$edge_lengths[EP[1]]
  Sigma <- matrix(0, 5, 5)
  t <- l*c(0,1,t_norm)
  D <- outer (t, t, `-`)
  d.index <- c(1,2)
  Sigma[-d.index, -d.index] <- r_2(D,         theta)
  Sigma[ d.index,  d.index] <- -r_2(as.matrix(dist(c(0,l))),   theta, 2)
  Sigma[d.index,  -d.index] <- -r_2(D[3:4-2,],                 theta, 1)
  Sigma[-d.index,  d.index] <- t(Sigma[d.index,  -d.index])

  B <- Sigma[1:4,5]%*%solve(Sigma[1:4,1:4])#Sigma[-5,5,drop=F]/Sigma[5, 5]
  CV_P <- CZ%*%t(B)
  C <- matrix(0, nrow = n.p * dim(graph$EtV)[1], ncol=3)
  for(i in 1:length(graph$EtV[,1])){
    l <- graph$edge_lengths[i]

    t_s <- seq(0,1,length.out=n.p)
    if(graph$EtV[i,1] == EP[1]){
      Sigma <- matrix(0, length(t_s)+5, length(t_s)+5)
      d.index <- c(1,2)
      index_boundary <- c(d.index,3:4)
      t <- l*c(0, 1, t_norm, t_s)
      D <- outer (t, t, `-`)
      Sigma[-d.index, -d.index] <- r_2(D,         theta)
      Sigma[ d.index,  d.index] <- -r_2(as.matrix(dist(c(0,l))),   theta, 2)
      Sigma[d.index,  -d.index] <- -r_2(D[3:4-2,],                 theta, 1)
      Sigma[-d.index,  d.index] <- t(Sigma[d.index,  -d.index])

      u_e <- CV_P[4*(i-1) + (1:4)]
      u_e <- u_e[c(2,4,1,3)] # aragne in X',X
      SinvS <- solve(Sigma[index_boundary,index_boundary],Sigma[index_boundary, -index_boundary] )
      Sigma_X <- Sigma[-index_boundary,-index_boundary] - Sigma[-index_boundary, index_boundary]%*%SinvS
      C_P <-  t(SinvS)[-1, ]%*%u_e + Sigma_X[1,-1]
    }else{
      Sigma <- matrix(0, length(t_s)+4, length(t_s)+4)
      d.index <- c(1,2)
      index_boundary <- c(d.index,3:4)
      t <- l*c(0,1,t_s)
      D <- outer (t, t, `-`)
      Sigma[-d.index, -d.index] <- r_2(D,         theta)
      Sigma[ d.index,  d.index] <- -r_2(as.matrix(dist(c(0,l))),   theta, 2)
      Sigma[d.index,  -d.index] <- -r_2(D[3:4-2,],                 theta, 1)
      Sigma[-d.index,  d.index] <- t(Sigma[d.index,  -d.index])

      u_e <- CV_P[4*(i-1) + (1:4)]
      u_e <- u_e[c(2,4,1,3)] # aragne in X',X
      SinvS <- solve(Sigma[index_boundary,index_boundary],Sigma[index_boundary, -index_boundary] )
      C_P <-  t(SinvS)%*%u_e

    }
    C[ (i-1) * n.p + (1:n.p),] <- cbind(graph$EtV[i,1],l*t_s, c(C_P) )
  }


  return(C)
}
