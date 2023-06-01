#Internal function to assemble operator matrix
#' @noRd
make.L <- function(n,kappa,C,G) {
  L0 <- G + kappa^2*C
  Ci <- Diagonal(1/rowSums(C),n=dim(C)[1])
  if(n==0) {
    L <- C
  } else if (n==1) {
    L <- L0
  } else {
    L <- L0
    for(i in 2:n) {
      L <- L%*%Ci%*%L0
    }
  }
  return(L)
}

#Internal function to assemble stiffness matrix
#' @noRd
make.G <- function(n,C,G) {
  Ci <- Diagonal(1/rowSums(C),n=dim(C)[1])
  if(n==0) {
    Gk <- C
  } else if(n==1) {
    Gk <- G
  } else {
    Gk <- G
    for(i in 2:n) {
      Gk <- Gk%*%Ci%*%G
    }
  }
  return(Gk)
}

#' Space-time precision operator discretization
#'
#' The precision matrix for all vertices for space-time field.
#'
#' @param graph A `metric_graph` object.
#' @param t Vector of time points.
#' @param kappa Spatial range parameter.
#' @param rho Drift parameter.
#' @param gamma Temporal range parameter.
#' @param alpha Smoothness parameter (integer) for spatial operator.
#' @param beta Smoothness parameter (integer) for Q-Wiener process.
#' @param sigma Variance parameter.
#' @return Precision matrix.
#' @export
make_Q_spacetime <- function(graph,t,kappa, rho, gamma, alpha, beta, sigma) {

  G <- graph$mesh$G
  C <- graph$mesh$C
  B <- graph$mesh$B
  n <- dim(C)[1]
  Cd <- Diagonal(rowSums(C),n=n)
  Ci <- Diagonal(1/rowSums(C),n=n)
  n <- dim(C)[1]

  nt <- length(t)
  d <- c(Inf, diff(t))
  dm1 <- c(d[2:nt], Inf)
  Gt <- -bandSparse(n = nt, m = nt, k = c(-1, 0, 1),
                    diagonals = cbind(1 / dm1, -(1 / dm1 + 1 / d), 1 / dm1))
  Ct <- bandSparse(n = nt, m = nt, k = c(-1, 0, 1),
                   diagonals = cbind(dm1 / 6, (dm1 + d) / 3, c(d[2:nt],Inf) / 6))
  Ct[1, 1:2] <- c(d[2], d[2] / 2) / 3
  Ct[nt, (nt - 1):nt] <- c(d[nt] / 2, d[nt]) / 3
  Bt <- bandSparse(n = nt, m = nt, k = c(-1, 0, 1),
                   diagonals = cbind(rep(0.5,nt), rep(0,nt), rep(-0.5,nt)))
  Bt[1,1] = -0.5
  Bt[nt,nt] = 0.5
  B0 <- Diagonal(n=nt,0)
  B0[1,1] <- B0[nt,nt] <- 1/2

  if(alpha == 0){
    Q <- kronecker(Gt + gamma^2*Ct, make.L(beta,kappa,C,G))
  } else if (alpha==2) {
    Q <- kronecker(Gt,make.L(beta,kappa,Cd,G)) + 2*gamma*kronecker(B0,make.L(beta+2,kappa,Cd,G))
    Q <- Q + gamma^2*kronecker(Ct, make.L(beta+4,kappa,Cd,G)) + gamma^2*rho^4*kronecker(Ct, G%*%Ci%*%G)
    Q <- Q + 6*gamma^2*rho^2*kronecker(Ct, make.L(beta+2,kappa,Cd,G)%*%Ci%*%G)
    M2 <- make.L(beta+1,kappa,Cd,G)%*%Ci%*%B
    Q <- Q - 2*rho*gamma*(kronecker(t(Bt),M2) + kronecker(Bt,t(M2)))
  } else {
    Q <- kronecker(Gt,make.L(beta,kappa,Cd,G)) + 2*gamma*kronecker(B0,make.L(beta+alpha,kappa,Cd,G))

    for(k in 0:alpha) {
      M1 <- make.L(beta+2*(alpha-k),kappa,Cd,G)%*%Ci%*%make.G(k,Cd,G)
      M2 <- make.L(beta+alpha-k,kappa,Cd,G)%*%Ci%*%make.G(floor(k/2),Cd,G)%*%Ci%*%B
      Q <- Q + gamma^2*choose(alpha,k)*rho^(2*k)*kronecker(Ct,M1)
      Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho^(k)*(kronecker(t(Bt),M2) + kronecker(Bt,t(M2)))
    }
  }
  return(Q/sigma^2)
}

#' Space-time precision operator Euler discretization
#'
#' The precision matrix for all vertices for space-time field
#'
#' @param graph A `metric_graph` object.
#' @param t Vector of time points.
#' @param kappa Spatial range parameter.
#' @param sigma Variance parameter.
#' @param theta Parameter theta for the Euler scheme.
#' @param rho Drift parameter.
#' @param gamma Temporal range parameter.
#' @param alpha Smoothness parameter (integer) for spatial operator.
#' @param beta Smoothness parameter (integer) for Q-Wiener process.
#' @return Precision matrix.
#' @export
make_Q_euler <- function(graph,t,kappa,rho,gamma,alpha,beta,sigma, theta = 1) {

  G <- graph$mesh$G
  C <- graph$mesh$C
  B <- graph$mesh$B

  dt <- t[2] - t[1]
  T <- length(t)
  n <- dim(C)[1]
  if(alpha == 0) {
    Lalpha <- C
  } else if (alpha == 1) {
    Lalpha <- G + kappa^2*C + rho*B
  } else {
    stop("not implemented")
  }

  LL <- Matrix(0,nrow = T*n, ncol = T*n)
  L0 <- 2*gamma*make.L(alpha+beta,kappa,C, G)/sigma^2
  if(beta==0) {
    Ltmp <- C - dt*(1-theta)*gamma*Lalpha
    CC = kronecker(Diagonal(n=T), sigma^2*dt*Ltmp%*%solve(C,t(Ltmp)))
    C.inv <- FALSE
    CC[1:n,1:n] <- solve(L0)
  } else {
    CC = kronecker(Diagonal(n=T), make.L(beta,kappa,G,C)/(dt*sigma^2))
    C.inv <- TRUE
    CC[1:n,1:n] <- L0
  }

  LL[1:n,1:n] <- Diagonal(n=dim(L0)[1])
  for(t in 1:(T-1)) {
    LL[(t*n+1):((t+1)*n),(t*n+1):((t+1)*n)] <- C + dt*theta*gamma*Lalpha
    LL[(t*n+1):((t+1)*n),((t-1)*n+1):(t*n)] <- -C + dt*(1-theta)*gamma*Lalpha
  }
  return(list(L = LL, C = CC, C.inv = C.inv))
}

#' space-time simulation based on implicit Euler discretization in time
#'
#' Simulation with starting value u0
#'
#' @param graph A `metric_graph` object.
#' @param t Vector of time points.
#' @param kappa Spatial range parameter.
#' @param rho Drift parameter.
#' @param gamma Temporal range parameter.
#' @param alpha Smoothness parameter (integer) for spatial operator.
#' @param beta Smoothness parameter (integer) for Q-Wiener process.
#' @param sigma Variance parameter.
#' @param u0 Starting value.
#' @param BC Which boundary condition to use (0,1). Here, 0 is no adjustment on
#' the boundary and 1 results in making the boundary condition stationary.
#' @return Precision matrix.
#' @export
simulate_spacetime <- function(graph, t, kappa, rho, gamma, alpha,
                               beta, sigma, u0, BC = 0) {
  n <- length(u0)
  Cd <- Diagonal(rowSums(graph$mesh$C),n=n)
  Cd[1:graph$nV] <- (graph$get_degrees()==1)*BC
  if(alpha == 1) {
    L <- graph$mesh$G + Cd + kappa^2*graph$mesh$C + rho*graph$mesh$B
  } else if(alpha == 2) {
    Ci <- Diagonal(1/rowSums(graph$mesh$C),n=n)
    L <- graph$mesh$G + Cd + kappa^2*graph$mesh$C + rho*graph$mesh$B
    L <- L%*%Ci%*%L
  } else {
    stop("not implemented")
  }
  if(beta == 0) {
    R <- chol(graph$mesh$C)
  } else {
    stop("not implemented")
  }
  U <- matrix(0,nrow=n,ncol = length(t))
  U[,1] <- u0
  for(i in 1:(length(t)-1)){
    dt <- t[i+1]-t[i]
    U[,i+1] <- as.vector(solve(graph$mesh$C + dt*gamma*L,
                               graph$mesh$C%*%U[,i] + sqrt(dt)*sigma*R%*%rnorm(n)))
  }
  return(U)
}
