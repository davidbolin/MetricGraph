library(sp)
library(Matrix)
library(MetricGraph)
library(htmlwidgets)
library(plotly)

plot.covariances <- function(graph,
                             Q,
                             Qs,
                             t.ind,
                             s.ind,
                             t.shift=0,
                             t) {
  if(is.list(Q)) {
    L <- Q$L
    N <- dim(L)[1]
  } else {
    N <- dim(Q)[1]
  }

  n <- dim(Qs)[1]

  T <- N/n
  if(length(t.ind)>4)
    stop("max 4 curves allowed")
  if(s.ind > n)
    stop("too large space index")
  if(max(t.ind)>T-1)
    stop("too large time index")

  cols <- c("green", "cyan","blue","red")[(4-length(t.ind)+1):4]

  v <- rep(0,dim(Qs)[1]); v[s.ind] <- 1

  c.spatial <- solve(Qs,v)

  p <- graph$plot_function(as.vector(c.spatial), plotly = TRUE, support_width = 1, line_color = "black")
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
      p <- graph$plot_function(as.vector(c), plotly = TRUE, p = p, support_width = 0, line_color = col)
    }
  }
  cat(max(c.spatial)/max(c))
  df <- data.frame(t=rep(t,length(t.ind)),y=c(t(ct)), i=rep(1:length(t.ind), each=length(t)))
  pt <- plot_ly(df, x = ~t, y = ~y, split = ~i, type = 'scatter', mode = 'lines')
  fig <- subplot(p,pt) %>% layout(title = "Covariances",
                                  scene = list(domain=list(x=c(0,0.5),y=c(0,1))),
                                  scene2 = list(domain=list(x=c(0.5,1),y=c(0,1))))
  fig$x$layout <- fig$x$layout[grep('NA', names(fig$x$layout), invert = TRUE)]
  print(fig)
  return(ct)
}




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




line1 <- Line(rbind(c(1,0),c(0,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,0),c(0,-1)))
line4 <- Line(rbind(c(0,1),c(1,0)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))

graph <- metric_graph$new(lines = lines)
#graph$plot(direction = TRUE)
h = 0.01
graph$build_mesh(h = h)
graph$compute_fem()

nt = 200
T = 0.5
t <- seq(from=0, to = T, length.out = nt)


#precision operator discretization
make.Q.spacetime <- function(graph,t,kappa, rho, gamma, alpha, beta) {

  G <- graph$mesh$G
  C <- graph$mesh$C
  B <- graph$mesh$B
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
  return(Q)
}
kappa <- 20
rho <- -10
gamma <- 0.0005
sigma <- 50
alpha <- 1
beta <- 0
Q <- make.Q.spacetime(graph,t,kappa,rho,gamma, alpha, beta)
if(alpha==1){
  Qs <- 2*gamma*make.L(beta+alpha,kappa,graph$mesh$C,graph$mesh$G)/sigma^2
}else if (alpha == 2) {
  C <- graph$mesh$C
  G <- graph$mesh$G
  Cd <- Diagonal(rowSums(C),n=n)
  Ci <- Diagonal(1/rowSums(C),n=n)
  Qs <- 2*gamma*(kappa^4*Cd + (2*kappa^2+rho^2)*G + G%*%Ci%*%G)/sigma^2
}

ct.dir <- plot.covariances(graph,Q = Q/sigma^2,Qs = Qs,
                           t.ind = c(1,nt/2,nt-2), s.ind = 50,t.shift = c(0),t = t)








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
kappa <- 20
rho <- -10
gamma <- 0.01
sigma <- 50
alpha <- 1
beta <- 0
nt = 200
T = 0.25
t <- seq(from=0, to = T, length.out = nt)

Q <- make.Q.euler(graph, t, kappa, rho, gamma, alpha, beta,sigma,theta=1)
ct.eul <- plot.covariances(graph,Q, Qs = 2*gamma*make.L(alpha+beta,kappa,graph$mesh$C,
                                                        graph$mesh$G)/sigma^2,
                 t.ind = c(1, nt/2, nt-2), s.ind = 50,t.shift = 0,t = t)


#precision operator discretization
make.Q.direct <- function(graph,t,kappa, rho, gamma) {
  G <- graph$mesh$G
  C <- graph$mesh$C
  B <- graph$mesh$B
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
  #B0[1,2] <- B0[nt,nt-1] <- -1
  Cd <- Diagonal(rowSums(C),n=n)
  L <- kappa^2*C + G
  Q <- kronecker(Gt, Cd) +  2*gamma*kronecker(B0,L) + gamma^2*kronecker(Ct, L%*%solve(Cd,L))
  Q <- Q + gamma^2*rho^2*kronecker(Ct,G) -gamma*rho*(kronecker(Bt,t(B))+kronecker(t(Bt),B))
  return(Q)
}
Q <- make.Q.direct(graph,t,kappa,rho,gamma)
ct.dir <- plot.covariances(graph,Q = Q/sigma^2,Qs = 2*gamma*L0/sigma^2,
                           t.ind = c(12,nt/2,nt-2), s.ind = 50,t.shift = c(0),t = t)

###########
kappa <- 10
rho <- 0
sigma <- 1.0
dt <- 0.01
n <- dim(graph$mesh$C)[1]
R <- chol(graph$mesh$C)
L <- graph$mesh$G + kappa^2*graph$mesh$C + rho*graph$mesh$B

u0 <- rep(0,n)
u0[1] <- 1
T <- 20
U <- matrix(0,nrow=n,ncol = T)
U[,1] <- u0

for(i in 1:(T-1)){
  U[,i+1] <- as.vector(solve((C + dt*L), C%*%U[,i] + sqrt(dt)*sigma*R%*%rnorm(n)))
}
fig <- graph$plot_movie(U)
fig
if(save.plot){
  htmlwidgets::saveWidget(fig, file = "spacetime.HTML", selfcontained = TRUE)
}




