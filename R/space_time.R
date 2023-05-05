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
#' The precision matrix for all vertices for space-time field
#' 
#' @param graph metric_graph object
#' @param t vector of time points
#' @param kappa range parameter kappa
#' @param rho drift parameter
#' @param gamma temporal range parameter
#' @param alpha smoothness parameter (integer) for spatial operator
#' @param beta smoothness parameter (integer) for Q-Wiener process
#' @param sigma variance parameter
#' @return Precision matrix
#' @export
make.Q.spacetime <- function(graph,t,kappa, rho, gamma, alpha, beta, sigma) {

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
#' @param graph metric_graph object
#' @param t vector of time points
#' @param kappa range parameter kappa
#' @param rho drift parameter
#' @param gamma temporal range parameter
#' @param alpha smoothness parameter (integer) for spatial operator
#' @param beta smoothness parameter (integer) for Q-Wiener process
#' @return Precision matrix
#' @export


make.Q.euler <- function(graph,t,kappa,rho,gamma,alpha,beta,sigma, theta = 1) {

  G <- graph$mesh$G
  C <- graph$mesh$C
  B <- graph$mesh$B

  dt <- t[2] - t[1]
  T <- length(t)
  n <- dim(L)[1]
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
#' @param graph metric_graph object
#' @param t vector of time points
#' @param kappa range parameter kappa
#' @param rho drift parameter
#' @param gamma temporal range parameter
#' @param alpha smoothness parameter (integer) for spatial operator
#' @param beta smoothness parameter (integer) for Q-Wiener process
#' @param sigma variance parameter
#' @param u0 starting value
#' @return Precision matrix
#' @export
simulate.spacetime <- function(graph, t, kappa, rho, gamma, alpha, beta, sigma, u0, BC = 0) {
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

#' Function for plotting space-time covariances
#' 
#' Plots of marginal spatial and temporal covariances of spacetime model
#' 
#' @param graph metric_graph object
#' @param Q space-time precision matrix
#' @param Qs Spatial precision matrix as reference
#' @param t.ind time indices to plot covariances for
#' @param s.ind space indices to plot covariances for
#' @param t.shift time shifts to plot, i.e., covariances C(u(s(s.ind),t),u(*,t+t.shift))
#' are shown.
#' @param show.temporal Plot the marginal temporal covariances?
#' @param t vector with the time points for which Q is computed
#' @export
plot.spacetime.covariances <- function(graph,
                                       Q,
                                       Qs = NULL,
                                       t.ind,
                                       s.ind,
                                       t.shift=0,
                                       show.temporal = FALSE,
                                       t) {
  if(is.list(Q)) {
    L <- Q$L
    N <- dim(L)[1]
  } else {
    N <- dim(Q)[1]
  }

  n <- N/length(t)

  T <- N/n
  if(length(t.ind)>4)
    stop("max 4 curves allowed")
  if(s.ind > n)
    stop("too large space index")
  if(max(t.ind)>T-1)
    stop("too large time index")

  cols <- c("green", "cyan","blue","red")[(4-length(t.ind)+1):4]

  if(!is.null(Qs)){
    v <- rep(0,dim(Qs)[1]); v[s.ind] <- 1
    c.spatial <- solve(Qs,v)
    p <- graph$plot_function(as.vector(c.spatial), plotly = TRUE, support_width = 1, line_color = "black")
  }

  time.index <- n*(0:(T-1)) + s.ind
  ct <- matrix(0,nrow = length(t.ind),ncol = T)
  for(i in 1:length(t.ind)) {
    v <- rep(0,N)
    v[(t.ind[i]-1)*n+s.ind] <- 1
    if(is.list(Q)){
      if(Q$C.inv) {
        tmp <- solve(Q$L,solve(Q$C,solve(t(Q$L),v)))
      } else {
        tmp <- solve(Q$L,Q$C%*%solve(t(Q$L),v))
      }
    } else {
      tmp <- solve(Q,v)
    }
    ct[i,] <- tmp[time.index]
    for(j in 1:length(t.shift)) {
      ind <- ((t.ind[i]-t.shift[j]-1)*n+1):((t.ind[i]-t.shift[j])*n)
      c <- tmp[ind]
      if(length(t.shift)>1) {
        col <- cols[j]
      } else {
        col <- cols[i]
      }
      if(i == 1 && is.null(Qs)) {
        p <- graph$plot_function(as.vector(c), plotly = TRUE, support_width = 0, line_color = col)
      } else {
        p <- graph$plot_function(as.vector(c), plotly = TRUE, p = p, support_width = 0, line_color = col)
      }

    }
  }
  if(show.temporal){
    df <- data.frame(t=rep(t,length(t.ind)),y=c(t(ct)), i=rep(1:length(t.ind), each=length(t)))
    pt <- plot_ly(df, x = ~t, y = ~y, split = ~i, type = 'scatter', mode = 'lines')
    fig <- subplot(p,pt) %>% layout(title = "Covariances",
                                    scene = list(domain=list(x=c(0,0.5),y=c(0,1))),
                                    scene2 = list(domain=list(x=c(0.5,1),y=c(0,1))))
    fig$x$layout <- fig$x$layout[grep('NA', names(fig$x$layout), invert = TRUE)]
  } else {
    fig <- p
  }

  print(fig)
  return(ct)
}
