#verify claim that one can solve circular edges using same expression as for non circular
# see likelihood inference in article
rm(list=ls())
library(GPGraph)
library(testthat)
if(0){
loc <- c(0.5,0.7,0.9)
sigma <- 1
kappa  <- 2
theta <- c(kappa, sigma)
l_e <- 1
Cov <- r_1_circle(c(0,loc), l_e, theta)
R0  <- r_1(as.matrix(dist(c(0,1,loc))), theta)

Sigma_inv <- solve(R0)[-c(1,2),-c(1,2)]
Sigma_Circle_inv <- solve(Cov)[-1,-1]
B_Circle <- Cov[1,-1]/Cov[1,1]
B = R0[-c(1,2),c(1,2)]%*%solve(R0[c(1,2),c(1,2)])
t(B_Circle)%*%Sigma_Circle%*%B_Circle -t(B%*%c(1,1))%*%Sigma_inv%*%(B%*%c(1,1))

}
#' test agreement between Q and R
test_that("Check agrement beteen covariance and precision matrix formulation", {
  theta <- c(1, 1)
  t <- 0:3
  Q0  <- Q.exp.line(theta, t)
  R0  <- r_1(as.matrix(dist(t)),theta)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))),c(as.matrix(Q0)), tolerance=1e-10)
  theta <- c(1.5, 0.5)
  t <- 0:3
  Q0  <- Q.exp.line(theta, t)
  R0  <- r_1(as.matrix(dist(t)),theta)
  R0_ <- solve(Q0)
  expect_equal(c(as.matrix(solve(R0))),c(as.matrix(Q0)), tolerance=1e-10)
})
