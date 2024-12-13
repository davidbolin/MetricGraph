#' Looking if one can genearlize MTP_2 for multivariate
#'
#' Start by examining a single edge then go to three
#'
#'
#'
#'
library(MetricGraph)
graphics.off()
if(0){
  #create a orhonomal basis
  edge1 <- rbind(c(1,0),c(1,1))
  edge2 <- rbind(c(1,1),c(0,1))
  edge3 <- rbind(c(0,1),c(0,0))
  edge4 <- rbind(c(0,0),c(1,0))
  edges = list(edge1,edge2,edge3,edge4)
  graph <- metric_graph$new(edges = edges)
  print(plot(graph))
  kappa <- 0.5
  tau   <- 1.2
  Q <- spde_precision(kappa = kappa, tau = tau,
                      alpha = 2, graph = graph)
  #create a basis such that it is true
  # (1, 1/kappa)
  # ()
  # gram schmidt

  graph$buildC(2, FALSE)

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
  Q_mod <- Tc %*% Q %*% t(Tc)
  Sigma <- t(Tc) %*% solve(Q_mod)%*% (Tc)
  index <- c(1:4,7:8,11:12)
  Sigma <- Sigma[index,index]
  Q0 <- solve(Sigma)
  Q <- round(Q0,4)
  index0 <- c(1:2,5:6)

  Qs <- list()
  Sigmas <- list()
  Sigmas[[1]] <- Sigma
  Qs[[1]] <- Q[index0, index0]
  Qs[[2]] <- solve(Sigma[-3,-3])[c(1,2,4,5),c(1,2,4,5)]

  #(1/Q0[3,3])*Q[index0, index0]%*%(Sigma[index0,3])%*%t(Sigma[index0,3])%*%Q[index0,index0]
  # A^-1 + A^-1 B (A - BD^-1B)^-1 B
  # Q + Q
  Sigmas[[1]] <- Sigma[-3,-3]
  print('**** diff')
  print(round(Qs[[1]]-Qs[[2]],4)[c(1,2),c(3,4)])
  print('**** cond')
  print(round(solve(Qs[[2]][1:2,1:2],Qs[[2]][c(1:2),c(3:4)]),3))
  index <- c(3,4,7,8)
  for(i in 2:length(index)){
    Qs[[i+1]] <-  solve(Sigma[-index[1:i],-index[1:i]])[1:4,1:4]
    Sigmas[[i+1]] <- Sigma[-index[1:i],-index[1:i]]
    print(paste('i = ',i))
    print('**** diff')
    print(round(Qs[[i+1]]-Qs[[i]],4)[c(1,2),c(3,4)])
    print('**** cond')
    print(round(solve(Qs[[i+1]][1:2,1:2],Qs[[i+1]][c(1:2),c(3:4)]),3))
  }
}

if(1){
  edge1 <- rbind(c(1,0),c(1,1))
  edge2 <- rbind(c(1,1),c(1,2))
  edge3 <- rbind(c(1,2),c(1,3))
  edges = list(edge1,edge2,edge3)
  graph <- metric_graph$new(edges = edges)
  print(plot(graph))
  kappa <- 0.5
  tau   <- 1.2
  Q <- spde_precision(kappa = kappa, tau = tau,
                      alpha = 2, graph = graph)

  graph$buildC(2, FALSE)
  b1 <- c(1, 1/kappa)
  b1 <- b1/sqrt(t(b1)%*%b1)
  b2 <- c(1,1)- (t(c(1,1))%*%b1)*b1
  b2 <- b2/sqrt(t(b2)%*%b2)
  B <- rbind(b1,b2)

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
  Q_mod <- Tc %*% Q %*% t(Tc)
  Sigma <- t(Tc) %*% solve(Q_mod)%*% (Tc)
  index <- c(1:4,7:8,11:12)
  Sigma <- Sigma[index,index]
  Q <- round(solve(Sigma),4)
  print(round(solve(Sigma[-c(3:4),-c(3:4)]),4))
  print(Q[-c(3:4),-c(3:4)])
}
if(0){
  library(MetricGraph)
  edge1 <- rbind(c(1,0),c(0.99,0))
  edge2 <- rbind(c(1+sqrt(0.5),sqrt(0.5)),c(1,0))
  edge3 <- rbind(c(1+sqrt(0.5),-sqrt(0.5)),c(1,0))
  edges = list(edge1,edge2,edge3)
  graph <- metric_graph$new(edges = edges)

  kappa <- 0.1
  tau   <- 1
  P1 <- c(1, 0.5)
  P2 <- c(3, 0.5)
  Q <- spde_precision(kappa = kappa, tau = tau,
                      alpha = 2, graph = graph)

  graph$buildC(2, FALSE)

  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
  Q_mod <- Tc %*% Q %*% t(Tc)
  Sigma <- t(Tc) %*% solve(Q_mod)%*% (Tc)
  index <- c(3,4,5,6,9,10)
  print(solve(Sigma[index,index])[c(5,6),c(5,6)])
}



