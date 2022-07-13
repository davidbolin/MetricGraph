#' computes the precision matrix for the vertices for an exponential graph model
#' @param P (n x 2) Vertex locations
#' @param E (m x 2) Edges
#' @param kappa  matern parameter
#' @export
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


#' computes matrices needed for graph for beta=1
#' @param P vertices
#' @param E edges
#' @param kappa  matern parameter
#' @param sigma  matern parameter
#' @param n  number of nodes to include per edge
#' @return a list with elements:
#'     Sigma : Covariance matrix of 'free' processes and derivatives on the edges
#'     A : Matrix that implements the Kirchoff constraints
#'     rep.ind : the indices of the repeated nodes
#'     p.ind : the indices in Sigma corresponding to the process
#'     d.ind : the indices in Sigma corresponding to the derivative
#'     P : matrix with the locations of the nodes where Sigma is evaluated
#' @export
build.Sigma.beta1 <- function(P,E,kappa,sigma,n){

  #edge lengths
  L <- sqrt((P[E[,2],1]-P[E[,1],1])^2 + (P[E[,2],2]-P[E[,1],2])^2)

  #number of vertices
  n.v <- dim(P)[1]

  #number of edges
  n.e <- dim(E)[1]

  Sigma <- matrix(0,nrow=2*n*length(L),ncol=2*n*length(L))

  p.inds <- vector("list",n.v) #indices of process at vertices
  d.inds <- vector("list",n.v) #indices of derivative at vertices
  d.sign <- vector("list",n.v) #sign of derivative
  proc.index <- NULL #indices for process nodes
  deri.index <- NULL #indices for derivative nodes
  P.sigma <- NULL #coorinates of nodes

  #loop over edges and construct free covariances
  for(i in 1:n.e){
    C <- build.C.beta1(L[i], kappa=kappa, sigma=sigma, nu=3/2)
    x <- seq(from=0,to=L[i],length.out = n)
    P.sigma <- cbind(P.sigma, P[E[i,1],]+outer(P[E[i,2],]-P[E[i,1],],seq(from=0,to=1,length.out=n)))
    S1 <- matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = L[i])
    S2 <- matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = L[i], deriv = c(0,1))
    S3 <- matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = L[i], deriv = c(1,1))
    Sigma[(i-1)*2*n + (1:(2*n)),(i-1)*2*n + (1:(2*n))]  <- rbind(cbind(S1, S2), cbind(t(S2),S3))

    p.inds[[E[i,1]]] <- c(p.inds[[E[i,1]]],(i-1)*2*n + 1)
    p.inds[[E[i,2]]] <- c(p.inds[[E[i,2]]],(i-1)*2*n + n)

    d.inds[[E[i,1]]] <- c(d.inds[[E[i,1]]],(i-1)*2*n + n + 1)
    d.inds[[E[i,2]]] <- c(d.inds[[E[i,2]]],(i-1)*2*n + 2*n)

    d.sign[[E[i,1]]] <- c(d.sign[[E[i,1]]],1)
    d.sign[[E[i,2]]] <- c(d.sign[[E[i,2]]],-1)

    proc.index <- c(proc.index,(i-1)*2*n + (1:n))
    deri.index <- c(deri.index,(i-1)*2*n + n + (1:n))
  }

  #now build A matrix
  i <- j <- x <- NULL
  n.c = 0 #number of constraints
  for(k in 1:n.v){
    n.c <- n.c+1
    d <- length(p.inds[[k]]) #degree of vertex
    i <- c(i,rep(n.c,d))
    j <- c(j,d.inds[[k]])
    x <- c(x,d.sign[[k]])
    if(d>1){
      for(kk in 1:(d-1)){
        n.c <- n.c + 1
        i <- c(i,rep(n.c,2))
        j <- c(j,p.inds[[k]][kk:(kk+1)])
        x <- c(x,c(1,-1))
      }
    }
  }
  A <- Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(n.c,2*n*n.e))

  #fix index vectors with the indices
  rep.ind <- NULL
  for(k in 1:dim(P)[1]){
    d <-length(p.inds[[k]])
    if(d>1){
      rep.ind <- c(rep.ind, p.inds[[k]][2:d],d.inds[[k]][2:d])
      proc.index <- setdiff(proc.index,p.inds[[k]][2:d])
      deri.index <- setdiff(deri.index,p.inds[[k]][2:d])
    }

  }

  return(list(Sigma=Sigma,
              A=A,
              rep.ind = rep.ind,
              p.ind = proc.index,
              d.ind = deri.index,
              P = t(P.sigma[,!duplicated(t(P.sigma))])))
}

#' Title precision for alpha=2 process on interval
#'
#' @param l length of interval
#' @param kappa range parameter
#' @param tau precision parameter
#'
#' @return precision matrix
#' @export
Q.alpha2.base <- function(l,kappa,tau){
  if(1){
    C <- GPGraph:::build.C.beta1(0.5, kappa=kappa, sigma=1, nu=3/2)
    S1 <- GPGraph:::matern.neumann.free2(x, x, C, kappa=kappa, sigma=1, nu=3/2, L = l)
    S2 <- GPGraph:::matern.neumann.free2(x, x, C, kappa=kappa, sigma=1, nu=3/2, L = l, deriv = c(0,1))
    S3 <- GPGraph:::matern.neumann.free2(x, x, C, kappa=kappa, sigma=1, nu=3/2, L = l, deriv = c(1,1))
    Sigma  <- rbind(cbind(S1, S2), cbind(t(S2),S3))
    Q <- solve(Sigma)[c(1,3,2,4),c(1,3,2,4)]
  } else{
    T = kappa*l
    c <- 2*kappa*tau^2/((1-2*T^2)^2 - 2*exp(2*T)*(2*T^2+1)+exp(4*T))
    q1 <- exp(4*T) - (1-2*T^2)^2 + 4*T*exp(2*T)
    q2 <- 4*T*exp(2*T)
    q3 <- 2*exp(T)*(2*T^2*(T-1)-T-exp(2*T)*(T+1)+1)
    q4 <- 2*T*exp(T)*(2*T^2-exp(2*T)-1)
    q5 <- -kappa^2*(1-2*T^2)^2 + exp(2*T)*(2*kappa^2+4*(kappa^2-1)*T^2 - 4*T - 2) - (kappa^2-2)*exp(4*T)
    q6 <- 2*exp(T)*(-2*T^3 - 2*T^2 + T + exp(2*T)*(T-1)+1)
    Q <- c*matrix(c(kappa^2*q1,kappa*q2,kappa^2*q3,kappa*q4,
                    kappa*q2,q5,kappa*q4,q6,
                    kappa^2*q3, kappa*q4, kappa^2*q1, kappa*q2,
                    kappa*q4, q6, kappa*q2, q5),4,4)
  }

  return(Q)
}

#' Build components for Q on line for alpha=2
#'
#' @param loc locations
#' @param kappa range parameter
#' @param tau precision parameter
#'
#' @return List with elements for the precision
#' @importFrom Matrix t
#' @export
Q.alpha2.line <- function(loc,kappa,tau){
  l <- diff(loc)
  n.int <- length(loc)-1
  for(i in 1:n.int){
    A.i <- Matrix(0,nrow=2,ncol=4*n.int)
    if(i == 1){
      Q <- Q.alpha2.base(l[i],kappa,tau)
      A.i[1,4*(i-1)+3] = 1
      A.i[1,4*i+1] = -1
      A.i[2,4*(i-1)+4] = 1
      A.i[2,4*i+2] = 1
      A <- A.i
    } else {
      Q <- bdiag(Q,Q.alpha2.base(l[i],kappa,tau))
      if(i < n.int){
        A.i[1,4*(i-1)+3] = 1
        A.i[1,4*i+1] = -1
        A.i[2,4*(i-1)+4] = 1
        A.i[2,4*i+2] = 1
        A <- rbind(A,A.i)
      }
    }


  }
  ind.proc <- c(seq(from=1,to=4*n.int,by=4),4*n.int-1)
  ind.der <- c(seq(from=2,to=4*n.int,by=4),4*n.int)
  Q.i <- solve(Q)
  Sigma <- Q.i - Q.i%*%t(A)%*%solve(A%*%Q.i%*%t(A))%*%A%*%Q.i
  Sigma <- Sigma[c(ind.proc,ind.der),c(ind.proc,ind.der)]
  return(list(Q=Q,A=A,ind.proc = ind.proc, ind.der = ind.der,
              Sigma = Sigma))
}
