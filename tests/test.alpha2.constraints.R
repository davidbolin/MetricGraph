
library(GPGraph)
library(sp)
#just random input params
kappa <- 0.07
sigma   <- 0.1
line.line <- Line(rbind(c(120,80),c(30,80)))
line.line2 <- Line(rbind(c(30,80),c(30,00)))

graph <-  graph.obj$new(sp::SpatialLines(list(Lines(list(line.line),ID="1"),
                                              Lines(list(line.line2),ID="2"))))
Q <- Q.matern2(c(kappa,sigma), graph$V, graph$EtV, graph$El, BC = 1)
graph$buildA(2, F)
Qtilde <- (graph$CBobj$T)%*%Q%*%t(graph$CBobj$T)
Qtilde <- Qtilde[-c(1:2),-c(1:2)]
Sigma.overdetermined  = t(graph$CBobj$T[-c(1:2),])%*%solve(Qtilde)%*%(graph$CBobj$T[-c(1:2),])

t <- c(0,graph$El[1],graph$El[1],graph$El[2]+ graph$El[1])

D <-  rep(1,4)%*%t(t)-t(rep(1,4)%*%t(t))
R00 <- r_2(D,c(kappa,sigma))
R01 <- -r_2(D,c(kappa,sigma),1)
R11 <- -r_2(D,c(kappa,sigma),2)
R <- cbind(rbind(R00,R01),rbind(t(R01),R11))
ind <- rep(0:3,each=2)+rep(c(1,5),times=4) #reoder
R <- R[ind,ind]



test_that("test if the over determined covariances are correct",{
  testthat::expect_equal( c(as.matrix(R)), c(as.matrix(Sigma.overdetermined)), tol=1e-9)
})
