library(sp)
library(Matrix)
library(MetricGraph)
library(htmlwidgets)
library(plotly)

plot.covariances <- function(graph,Q=NULL,L,C,Ls,Cs, t.ind, s.ind,t.shift=0,t) {
  n <- dim(Cs)[1]
  if(is.null(Q)) {
    N <- dim(L)[1]
  } else {
    N <- dim(Q)[1]
  }
  T <- N/n
  if(length(t.ind)>4)
    stop("max 4 curves allowed")
  if(s.ind > n)
    stop("too large space index")
  if(max(t.ind)>T-1)
    stop("too large time index")

  cols <- c("green", "cyan","blue","red")[(4-length(t.ind)+1):4]

  v <- rep(0,dim(Ls)[1]); v[s.ind] <- 1
  c.spatial <- solve(Ls,Cs%*%solve(t(Ls),v))
  p <- graph$plot_function(as.vector(c.spatial), plotly = TRUE, support_width = 1, line_color = "black")
  time.index <- n*(0:(T-1)) + s.ind
  ct <- matrix(0,nrow = length(t.ind),ncol = T)
  for(i in 1:length(t.ind)) {
    v <- rep(0,N)
    v[(t.ind[i]-1)*n+s.ind] <- 1
    if(is.null(Q)){
      tmp <- solve(L,C%*%solve(t(L),v))
    } else {
      tmp <- solve(Q,v)
    }
    ct[i,] <- tmp[time.index]
    for(j in 1:length(t.shift)) {
      ind <- ((t.ind[i]-t.shift[j]-1)*n+1):((t.ind[i]-t.shift[j])*n)
      c <- tmp[ind]
      p <- graph$plot_function(as.vector(c), plotly = TRUE, p = p, support_width = 0, line_color = cols[i])
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
}



save.plot <- FALSE

line1 <- Line(rbind(c(1,0),c(0,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,0),c(0,-1)))
line4 <- Line(rbind(c(0,1),c(1,0)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))

graph <- metric_graph$new(lines = lines)
graph$plot(direction = TRUE)
h = 0.01
graph$build_mesh(h = h)
graph$compute_fem()


#precision operator discretization
make.Q.direct <- function(graph,t,kappa, rho, kappa.t) {
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
  Cd <- Diagonal(rowSums(C),n=n)
  L <- kappa^2*C + G
  Q <- kappa.t^2*kronecker(Gt, Cd) +  kronecker(Ct, L%*%solve(Cd,L))
  Q <- Q + rho^2*kronecker(Ct,G) -kappa.t*rho*(kronecker(Bt,t(B))+kronecker(t(Bt),B))
  return(Q)
}

nt = 400
T = 2#h^2*nt*20
t <- seq(from=0, to = T, length.out = nt)
kappa <- 10
rho <- -20
kappa.t <- 30
sigma <- 1
n <- dim(graph$mesh$C)[1]
Q <- make.Q.direct(graph,t,kappa,rho,kappa.t)
L0 <- graph$mesh$G + kappa^2*Diagonal(rowSums(graph$mesh$C),n=n)
plot.covariances(graph,Q = Q/sigma^2,Ls = L0,Cs = sigma^2*L0/(2*kappa.t),
                 t.ind = c(nt/2), s.ind = 50,t.shift = 10,t = t)

#symmetric implementation
make.Q.sym <- function(L,C,t,kappa.t) {
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
  Qt <- kappa.t^2*Gt+Ct
  Qt[1,1] <- Qt[1,1] + kappa.t
  Qt[nt,nt] <- Qt[nt,nt] + kappa.t
  C <- Diagonal(rowSums(graph$mesh$C),n=n)
  Q <- kronecker(Qt, L)
}
sigma <- 1
kappa.t <- 0.05
Q <- make.Q.sym(L, graph$mesh$C, t, kappa.t)

plot.covariances(graph,Q = Q/sigma^2,Ls = L,Cs = sigma^2*L/(2*kappa.t),
                 t.ind = c(nt/2,nt-2), s.ind = 50,t.shift = 0,t = t)


make.Q.euler <- function(L0,L,C, t,kappa.t) {
  dt <- t[2] - t[1]
  T <- length(t)
  n <- dim(L)[1]
  Lbar = kappa.t*C + dt*L
  CC = kronecker(Diagonal(n=T), dt*C)
  LL <- Matrix(0,nrow = T*n, ncol = T*n)
  LL[1:n,1:n] <- L0
  CC[1:n,1:n] <- C
  for(t in 1:(T-1)) {
    LL[(t*n+1):((t+1)*n),(t*n+1):((t+1)*n)] <- Lbar
    LL[(t*n+1):((t+1)*n),((t-1)*n+1):(t*n)] <- -kappa.t*C
  }
  return(list(L = LL, C = CC))
}

Q <- make.Q.euler(L0, L, graph$mesh$C, t,kappa.t)
plot.covariances(graph,L = Q$L,C = sigma^2*Q$C, Ls = L,Cs = 4*sigma^2*L/(kappa*sqrt(kappa.t)),
                 t.ind = c(nt/2,nt-2), s.ind = 25,t.shift = 0,t = t)




make.Q.petrov <- function(L0,L,C, t,kappa.t) {
  n <- dim(C)[1]
  h <- diff(t)
  nt <- length(h)
  Gt = as(bandSparse(n=nt,m=nt+1,k=c(0,1),diagonals=cbind(-rep(1,nt), rep(1,nt))),"dgCMatrix")
  Ct <- as(bandSparse(n=nt,m=nt+1,k=c(0,1),diagonals=cbind(0.5*h, 0.5*h)),"dgCMatrix")
  Ct0 = sparseMatrix(i=1:(nt+1),j=1:(nt+1),x=c(h,0) + c(0,h),dims = c(nt+1,nt+1))
  LL <- -kappa.t*kronecker(Gt,C) + kronecker(Ct,L)
  CC <- kronecker(Ct0, C)
  return(list(L = LL, C = CC))
}

Q <- make.Q.petrov(L0, L, Diagonal(rowSums(graph$mesh$C),n=dim(graph$mesh$C)[1]), t,kappa.t)
Q <- Q$L%*%solve(Q$C, t(Q$L))
plot.covariances(graph,Q = Q,Ls = L,Cs = 4*sigma^2*L/(kappa*sqrt(kappa.t)),
                 t.ind = c(nt/2,nt-2), s.ind = 50,t.shift = 0,t = t[-1])







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




