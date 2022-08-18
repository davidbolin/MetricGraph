#'
#' Building the precision matrix for the expontial case on each vertex
#' @param theta - kappa
#' @param V vertex position
#' @param EtV [,2-3] index of upper and lower edge
#' @param El length of each vertex
#' @return Q (precision matrix)
#' @export
Q.exp <- function(theta, V,EtV, El){

  i_ <- j_ <- x_ <- rep(0, dim(V)[1]*4)
  nE <- dim(EtV)[1]
  count <- 0
  for(i in 1:nE){
    l_e <- El[i]
    c1 <- exp(-theta*l_e)
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
                            x=x_[1:count],
                            dims=c(n.v, n.v))
  return(Q)
}
#'
#' Compute the precision matrix of observations one a line
#'
#' @param theta      - kappa
#' @param t          - (n x 1) relative position on the line start with 0 end with 1
#' @param l_e        - (double) length of the line
#' @param t_sorted   - (bool)
#' @export
Q.exp.line <- function(theta, t,  t_sorted=FALSE){

  l_t = length(t)
  i_ <- j_ <- x_ <- rep(0, 4*(l_t-1))
  count = 0

  if(t_sorted==F){
    order_t  <- order(t)
    t        <- t[order_t]
    P <- sparseMatrix(seq_along(order_t), order_t, x=1)
  }

  for(i in 2:l_t){

    c1 <- exp(-theta* (t[i] - t[i-1]))
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
  Q <- Matrix::sparseMatrix(i=i_[1:count],
                            j=j_[1:count],
                            x=x_[1:count],
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

  Q <- Q.exp.line(theta[3], t)
  Q <- Q/theta[2]
  index_E <- length(py) + 1:2
  Q_X <- Q[-index_E,-index_E]
  mu_X <- as.vector(Matrix::solve(Q_X,-Q[-index_E,index_E] %*% u_e))
  if(is.null(py)==F){
    Matrix::diag(Q_X)[1:length(py)] <-  Matrix::diag(Q_X)[1:length(py)] + 1/theta[1]^2
    AtY = rep(0,dim(Q_X)[1])
    AtY[1:length(py)] = (y - mu_X[1:length(py)])/theta[1]^2
    mu_X = mu_X + as.vector(Matrix::solve(Q_X,AtY))
  }
  R_X <- Matrix::chol(Q_X)
  x <- rep(0, length(t))
  z <- rnorm(dim(R_X)[1])
  x[-index_E] <- mu_X + as.vector(solve(R_X,solve(R_X,z,system = 'Lt'), system='Pt'))
  x[index_E] <- u_e
  return(cbind(t,x))

}
