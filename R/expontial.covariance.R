
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
#' @param graph metric_graph object
#' @param BC boundary condition =0 neumann, 1 = correct boundary to stationary
#' @param build (bool) if true return Q the precision matrix otherwise return list(i,j,x, nv)
#' @return Q (precision matrix)
#' @export
Q.exp <- function(theta, graph, BC = 1, build=TRUE){

  kappa <- theta[1]
  sigma <- theta[2]
  i_ <- j_ <- x_ <- rep(0, dim(graph$V)[1]*4)
  count <- 0
  for(i in 1:graph$nE){
    l_e <- graph$edge_lengths[i]
    c1 <- exp(-kappa*l_e)
    c2 <- c1^2
    one_m_c2 = 1-c2
    c_1 = 0.5 + c2/one_m_c2
    c_2 = -c1/one_m_c2

    if(graph$E[i,1]!=graph$E[i,2]){

      i_[count + 1] <- graph$E[i,1]
      j_[count + 1] <- graph$E[i,1]
      x_[count + 1] <- c_1

      i_[count + 2] <- graph$E[i,2]
      j_[count + 2] <- graph$E[i,2]
      x_[count + 2] <- c_1


      i_[count + 3] <- graph$E[i,1]
      j_[count + 3] <- graph$E[i,2]
      x_[count + 3] <- c_2

      i_[count + 4] <- graph$E[i,2]
      j_[count + 4] <- graph$E[i,1]
      x_[count + 4] <- c_2
      count <- count + 4
    }else{
      i_[count + 1] <- graph$E[i,1]
      j_[count + 1] <- graph$E[i,1]
      x_[count + 1] <- tanh(0.5 * theta * l_e)
      count <- count + 1
    }
  }
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
                              dims=c(graph$nV, graph$nV))


    return(Q)
  }else{
    return(list(i = i_[1:count],
                j = j_[1:count],
                x  = (2*kappa / sigma^2)*x_[1:count],
                dims=c(graph$nV, graph$nV)))
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

  index_E <- length(py) + 1:2
  Q_X <- Q[-index_E,-index_E]
  mu_X <- as.vector(Matrix::solve(Q_X,-Q[-index_E,index_E] %*% u_e))
  if(is.null(py)==FALSE){
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
#' Compute covariance of a point to the entire graph (discretized)
#'
#' @param  EP    - (2 x 1) [1] -edge number, [2] -normalized location on the edge
#' @param  theta - (3 x 1) (kappa,sigma)
#' @param  graph  - (graph)
#' @param  n.p    - (int) number of points to compute the covariance on on each edge
#' @return C      - (n.p*numer of edges x 3) [,1] edge number [,2] length from lower edge [,3] covariance
#' @export
covariance.point.to.graph.exp <- function(EP, theta, graph, n.p = 50){

  kappa <- theta[1]
  sigma <- theta[2]
  #compute covarains of the two edges of EP[1]
  Q <- Q.exp(theta, graph$V, graph$E, graph$edge_lengths)
  R <- Cholesky(Q, LDL = FALSE, perm = TRUE)
  Vs <- graph$E[EP[1],]
  Z <- matrix(0, nrow=dim(graph$V)[1], ncol=2)
  Z[Vs[1],1] = 1
  Z[Vs[2],2] = 1
  V  <- solve(R,solve(R,Z,system = 'P'), system='L')
  CV <- solve(R,solve(R,V,system = 'Lt'), system='Pt')

  # compute covarains between two edges and the point
  t_norm <- EP[2]
  l <- graph$edge_lengths[EP[1]]
  Q_line <- as.matrix(Q.exp.line(theta, c(0, l*t_norm,l)))
  Q_AB <- Q_line[2,c(1,3),drop=FALSE]
  Q_AA <- Q_line[2,2]
  B <- -Q_AB/Q_AA

  # What is the covariance of a point to an edge
  #COV[X,Y] = cov[Xtilde+BZ,Y] = B Cov[Z,Y]
  CV_P <- CV%*%t(B)
  C <- matrix(0, nrow = n.p * graph$nE, ncol=3)
  for(i in 1:graph$nE){
    l <- graph$edge_lengths[i]
    t_s <- seq(0,1,length.out=n.p)
    if(i == EP[1]){
      D_matrix <- as.matrix(dist(c(0,l,l*t_norm,l*t_s)))
      S <- r_1(D_matrix, c(kappa,sigma))

      #covariance update see Art p.17
      E.ind         <- c(1:2)
      Obs.ind       <- -E.ind
      Bt            <- solve(S[E.ind, E.ind,drop=FALSE], S[E.ind, Obs.ind,drop=FALSE])
      Sigma_i       <- S[Obs.ind,Obs.ind,drop=FALSE] - S[Obs.ind, E.ind,drop=FALSE] %*% Bt
      C_P <- CV_P[graph$E[i,]]%*%Bt[,-1] + Sigma_i[1,-1]
    }else{

      D_matrix <- as.matrix(dist(c(0,l,l*t_s)))
      S <- r_1(D_matrix,c(kappa,sigma))

      #covariance update see Art p.17
      E.ind         <- c(1:2)
      Obs.ind       <- -E.ind
      Bt            <- solve(S[E.ind, E.ind,drop=FALSE],S[E.ind, Obs.ind,drop=FALSE])
      C_P <- CV_P[graph$E[i,]]%*%Bt

    }
    C[ (i-1) * n.p + (1:n.p),] <- cbind(i,l*t_s, c(C_P) )
  }
  return(C)
}
