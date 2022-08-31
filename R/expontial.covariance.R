
#'
#' The exponential covariance
#' @param D     - vector, or matrix distance
#' @param theta - kappa, sigma
#' @export
r_1 <- function(D,theta){
  return( ( theta[2]^2/(2 * theta[1]))* exp( -theta[1] * abs(D)))
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
#' @return Q (precision matrix)
#' @export
Q.exp <- function(theta, V,EtV, El){

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
  Q <- Matrix::sparseMatrix(i=i_[1:count],
                            j=j_[1:count],
                            x=(2*kappa / sigma^2)*x_[1:count],
                            dims=c(n.v, n.v))
  return(Q)
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
#' @export
sample.line.expontial<-function(theta, u_e, l_e, t=NULL, Line=NULL, nt=100,  py=NULL, y=NULL){

  if(is.null(t)){
    t  = seq(0,1,length.out=nt+2)[c(-1,-(nt+2))]
    t <- rgeos::gProject(Line,rgeos::gInterpolate(Line, t, normalized = T))
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
  R_X <- Matrix::Cholesky(Q_X)
  x <- rep(0, length(t))
  z <- rnorm(dim(R_X)[1])
  x[-index_E] <- mu_X + as.vector(Matrix::solve(R_X,Matrix::solve(R_X,z,system = 'Lt'), system='Pt'))
  x[index_E] <- u_e

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
#' Computes the log likelihood function fo theta for the graph object
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph.obj - graphical object
#' @export
likelihood.exp.graph <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$El)
  R <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)
  loglik <- 0.5*Matrix::determinant(R)$modulus[1]
  Qp <- Q
  #build BSIGMAB
  Qpmu <- rep(0, nrow(graph.obj$V))
  obs.edges <- unique(graph.obj$PtE[,1])

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
    }else{
      Qpmu[E]    <- Qpmu[E] + t(Sigma_iB)%*%y_i
      Qp[E,E] <- Qp[E,E] + BtSinvB + t(BtSinvB)
      Qp[E,E] <- 0.5*(Qp[E,E] + t(Qp[E,E]))
    }
    loglik <- loglik - 0.5  * t(y_i)%*%solve(Sigma_i,y_i)
    loglik <- loglik - sum(log(diag(R)))
  }
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  loglik <- loglik - 0.5*Matrix::determinant(R)$modulus[1]

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Qpmu,system = 'P'), system='L')))
  #Qpmu <- as.vector(solve(R,solve(R, v,system = 'Lt'), system='Pt'))

  loglik <- loglik + 0.5  * t(v)%*%v
  return(loglik[1])
}

