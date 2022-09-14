library(GPGraph)
library(rSPDE)
#library(testthat)
set.seed(1)
kappa <- 0.1*runif(1)+0.1
sigma <- runif(1)+1

c <- 1/(4*kappa^3)

test_that("Check agrement covariance function agrees", {

  x <- seq(0,1,length.out=10)
  testthat::expect_equal( GPGraph::r_2(x, c(kappa,sigma)), rSPDE::matern.covariance(x, kappa, 3/2, sigma)[,1]*c, tol=1e-9)
})


test_that("Check agrement derivative covariance function agrees", {
  x <- seq(-1,1,length.out=20)
  testthat::expect_equal( GPGraph::r_2(x, c(kappa,sigma),1), GPGraph::matern.derivative(x, kappa, 3/2, sigma)[,1]*c, tol=1e-9)
})
test_that("Check agrement derivative covariance function agrees", {
  x <- seq(-1,1,length.out=20)
  testthat::expect_equal( GPGraph::r_2(x, c(kappa,sigma),2), GPGraph::matern.derivative(x, kappa, 3/2, sigma,2)[,1]*c, tol=1e-9)
})

test_that("Check agrement covariance matrix", {
  l_e <- runif(1) + 0.5
  x_ <- c(0,l_e)
  D <- outer(x_,x_,"-")
  r <- rSPDE::matern.covariance(D, kappa=kappa, nu=3/2, sigma=sigma)
  r1 <- -GPGraph::matern.derivative(D,  kappa=kappa, sigma=sigma, nu=3/2, deriv = 1)
  r2 <- -GPGraph::matern.derivative(D,  kappa=kappa, sigma=sigma, nu=3/2, deriv = 2)
  Sigma.0 <- rbind(cbind(r, r1), cbind(t(r1),r2))*c
  r_00 <- GPGraph::r_2(D, c(kappa,sigma))
  r_01 <- - GPGraph::r_2(D, c(kappa,sigma),1)
  r_11 <- - GPGraph::r_2(D, c(kappa,sigma),2)

  Sigma_ <- rbind(cbind(r_00, r_01), cbind(t(r_01),r_11))
  testthat::expect_equal( c(Sigma.0), c(Sigma_), tol=1e-9)
})
l_e <- runif(1) + 0.5
x_ <- c(0,l_e)
D <- outer(x_,x_,"-")
r_00 <- GPGraph::r_2(D, c(kappa,sigma))
r_01 <- - GPGraph::r_2(D, c(kappa,sigma),1)
r_11 <- - GPGraph::r_2(D, c(kappa,sigma),2)
# order by node not derivative
R_00 <- matrix(c(r_00[1], r_01[1,1], r_01[1,1], r_11[1,1]),2,2)
R_01 <- matrix(c(r_00[2], r_01[2,1], r_01[1,2], r_11[2,1]),2,2)
R_node <- rbind(cbind(R_00, R_01), cbind(t(R_01),R_00))
Q_adj = solve(R_node)-0.5 * solve(rbind(cbind(R_00, matrix(0,2,2)), cbind(matrix(0,2,2),R_00)))


build.C.beta1 <- function(L, kappa, sigma){
  C_0 <- matern.neumann.free(c(0,L),c(0,L),kappa,sigma=1, nu=3/2, L=L,deriv=c(1,1))
  return(sigma^2*solve(solve(C_0) -0.5*diag(2)/kappa^2))
}
C <- build.C.beta1(l_e, kappa, sigma)
r_free.2 <- matern.neumann.free2(x_, x_,C, kappa, sigma, nu=3/2, L = l_e)
rd1_free.2 <- matern.neumann.free2(x_, x_,C, kappa, sigma, nu=3/2, L = l_e, deriv = c(0,1))
rd2_free.2 <- matern.neumann.free2(x_, x_,C, kappa, sigma, nu=3/2, L = l_e, deriv = c(1,1))
Sigma.2  <- rbind(cbind(r_free.2, rd1_free.2), cbind(t(rd1_free.2),rd2_free.2))*c
Sigma.2 <- Sigma.2[c(1,3,2,4), c(1,3,2,4)] #order by nodes

test_that("test agrement precision matrix and article method", {
  testthat::expect_equal( c(Q_adj), c(solve(Sigma.2)), tol=1e-9)
})


test_that("test adjusted covariance matrix against article formula",{
  R00R0l <-rbind(cbind(R_00,R_01),cbind(t(R_01),R_00))
  Adj <- solve(rbind(cbind(R_00,-R_01),cbind(-t(R_01),R_00)))
  R_adj <- R_node + R00R0l%*%Adj%*%R00R0l
  testthat::expect_equal( c(Sigma.2), c(R_adj), tol=1e-9)
})

#T_ <- l_e*kappa
#C <- (2*kappa*sigma^(-2))/((1- 2 * T_^2)^2 - 2*exp(2*T_) * (2 * T_^2 + 1) + exp(4* T_))
#q_1 <- exp(4*T_) - (1- 2 * T_)^2 + 4  * T_ * exp(2*T_)
#q_2 = 4* T_ * exp(2*T_)
#Q_11 <- C*q_1*kappa^2
#Q_12 <- C * q_2 * kappa


