
test_that("Check agrement covariance function agrees", {
  set.seed(1)
  kappa <- 0.1*runif(1)+0.1
  sigma <- runif(1)+1
  c <- 1/(4*kappa^3)
  x <- seq(0,1,length.out=10)
  expect_equal(GPGraph:::r_2(x, c(kappa,sigma)),
               rSPDE::matern.covariance(x, kappa, 3/2, sigma)[,1]*c, tol=1e-9)
})


test_that("Check agrement derivative covariance function agrees", {
  set.seed(1)
  kappa <- 0.1*runif(1)+0.1
  sigma <- runif(1)+1
  c <- 1/(4*kappa^3)
  x <- seq(-1,1,length.out=20)
  expect_equal(GPGraph:::r_2(x, c(kappa,sigma),1),
               matern.derivative(x, kappa, 3/2, sigma)[,1]*c, tol=1e-9)
})
test_that("Check agrement derivative covariance function agrees", {
  set.seed(1)
  kappa <- 0.1*runif(1)+0.1
  sigma <- runif(1)+1
  c <- 1/(4*kappa^3)
  x <- seq(-1,1,length.out=20)
  expect_equal(GPGraph:::r_2(x, c(kappa,sigma),2),
               matern.derivative(x, kappa, 3/2, sigma,2)[,1]*c, tol=1e-9)
})

test_that("Check agrement covariance matrix", {
  set.seed(1)
  kappa <- 0.1*runif(1)+0.1
  sigma <- runif(1)+1
  c <- 1/(4*kappa^3)
  l_e <- runif(1) + 0.5
  x_ <- c(0,l_e)
  D <- outer(x_,x_,"-")
  r <- rSPDE::matern.covariance(D, kappa=kappa, nu=3/2, sigma=sigma)
  r1 <- -matern.derivative(D,  kappa=kappa, sigma=sigma, nu=3/2, deriv = 1)
  r2 <- -matern.derivative(D,  kappa=kappa, sigma=sigma, nu=3/2, deriv = 2)
  Sigma.0 <- rbind(cbind(r, r1), cbind(t(r1),r2))*c
  r_00 <- GPGraph:::r_2(D, c(kappa,sigma))
  r_01 <- - GPGraph:::r_2(D, c(kappa,sigma),1)
  r_11 <- - GPGraph:::r_2(D, c(kappa,sigma),2)

  Sigma_ <- rbind(cbind(r_00, r_01), cbind(t(r_01),r_11))
  testthat::expect_equal( c(Sigma.0), c(Sigma_), tol=1e-9)
})



test_that("test agrement precision matrix and article method", {
  set.seed(1)
  kappa <- 0.1*runif(1)+0.1
  sigma <- runif(1)+1
  c <- 1/(4*kappa^3)
  l_e <- runif(1) + 0.5
  x_ <- c(0,l_e)
  D <- outer(x_,x_,"-")
  r_00 <- GPGraph:::r_2(D, c(kappa,sigma))
  r_01 <- - GPGraph:::r_2(D, c(kappa,sigma),1)
  r_11 <- - GPGraph:::r_2(D, c(kappa,sigma),2)
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

  testthat::expect_equal( c(Q_adj), c(solve(Sigma.2)), tol=1e-7)

  #test adjusted covariance matrix against article formula
  R00R0l <-rbind(cbind(R_00,R_01),cbind(t(R_01),R_00))
  Adj <- solve(rbind(cbind(R_00,-R_01),cbind(-t(R_01),R_00)))
  R_adj <- R_node + R00R0l%*%Adj%*%R00R0l
  testthat::expect_equal( c(Sigma.2), c(R_adj), tol=1e-9)
})

#
# T = kappa*l_e
# tau = 1/sigma
# Q11 = (2*kappa^3 * tau^2*(2*T+sinh (2*T)))/(-2 * T^2+cosh (2*T)-1)
# Q12 = Q21 = (4*kappa^2*tau^2*T^2)/(-2*T^2+cosh(2 *  T)-1)
# Q13 = Q31 = (4 *kappa^3 * tau^2*(sinh (T)+T *cosh (T)))/(2*T^2-cosh (2*T)+1)
# Q14 = Q41 = (4 *kappa^2 * tau ^2 * T *sinh (T))/(-2 * T^2+cosh (2*T)-1)
# Q22 = (2*kappa  *tau^2 *(sinh(2*T)-2 * T))/(-2 * T^2+cosh (2*T)-1)
# Q44 = Q22
# Q23 = Q32 = -(4*kappa^2 * tau^2 * T* sinh(T))/(-2 * T^2+cosh (2 * T)-1)
# Q24 = Q42 = (4 * kappa * tau^2*(sinh (T)-T*cosh (T)))/(2*T^2-cosh (2*T)+1)
# Q33 = (2 * kappa^3*tau^2*(2*T+sinh(2 * T)))/(-2 * T^2+cosh (2*T)-1)
# Q34 = Q43 = (4 * kappa^2*tau^2*T^2)/(2*T^2-cosh(2*T)+1)
#
# Q <- matrix(c(Q11,Q21,Q31,Q41,Q12,Q22,Q32,Q42,
#               Q13,Q23,Q33,Q43,Q14,Q24,Q34,Q44),4,4)
#
# Q <- Q[c(1,3,2,4),c(1,3,2,4)]


test_that("test if computing covariance are equivalent",{
  set.seed(13)
  nt <- 40
  kappa <- 0.3
  sigma_e <- 0.1
  sigma   <- 1
  theta <-  c(sigma_e,kappa,sigma)
  line.line2 <- Line(rbind(c(30,80),c(140,80)))
  line.line <- Line(rbind(c(30,00),c(30,80)))

  graph.temp <-  metric_graph$new(Lines = sp::SpatialLines(list(Lines(list(line.line),ID="1"),
                                                     Lines(list(line.line2),ID="2"))))
  Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 2, graph = graph.temp, BC = 1)
  graph.temp$buildC(2, FALSE)
  Qmod <- (graph.temp$CBobj$T)%*%Q%*%t(graph.temp$CBobj$T)
  Qtilde <- Qmod
  Qtilde <- Qtilde[-c(1:2),-c(1:2)]
  R <- Cholesky(Qtilde,LDL = FALSE, perm = TRUE)
  V0 <- as.vector(solve(R, solve(R,rnorm(6), system = 'Lt')
                        , system = 'Pt'))
  u_e <- t(graph.temp$CBobj$T)%*%c(0,0,V0)
  X <- c()
  for(i in 1:length(graph.temp$edge_lengths)){
    X <- rbind(X,cbind(sample.line.matern2(theta,
                                           u_e[4*(i-1) +1:4],
                                           Line = graph.temp$Lines[i,],
                                           graph.temp$edge_lengths[i],
                                           nt = nt),i))
  }
  X[,2] <- X[,2] + sigma_e*rnorm(nt)

  graph.temp$add_observations2(y = X[,2], PtE = X[,c(3,1)])
  graph.temp$buildC(2, FALSE)
  lik <- likelihood_graph_spde(theta = theta, graph = graph.temp, alpha = 2)
  graph.temp$observation_to_vertex()
  graph.temp$buildC(2, FALSE)
  lik2 <-likelihood_graph_covariance(theta = theta, graph = graph.temp, model = "alpha2")

  expect_equal(as.matrix(lik),as.matrix(lik2), tolerance=1e-10)

  pm <- spde_posterior_mean(theta, alpha = 2, obs = TRUE, graph = graph.temp)# posterior mean
  n.o <- length(graph.temp$y)
  n.v <- dim(graph.temp$V)[1]
  n.c <- 1:length(graph.temp$CBobj$S)
  Q <- spde_precision(kappa = kappa, sigma =sigma, alpha = 2, graph = graph.temp, BC = 1)
  Qtilde <- (graph.temp$CBobj$T)%*%Q%*%t(graph.temp$CBobj$T)
  Qtilde <- Qtilde[-n.c,-n.c]
  Sigma.overdetermined  = t(graph.temp$CBobj$T[-n.c,])%*%solve(Qtilde)%*%(graph.temp$CBobj$T[-n.c,])
  index.obs <-  4*(graph.temp$PtE[,1]-1) + (1 * (abs(graph.temp$PtE[,2])<1e-14)) + (3 * (abs(graph.temp$PtE[,2])>1e-14))
  Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  Sigma.Y <- Sigma
  diag(Sigma.Y) <- diag(Sigma.Y) + theta[1]^2
  pm2 <- Sigma%*%solve(Sigma.Y,graph.temp$y)

  expect_equal(as.matrix(pm),as.matrix(pm2), tolerance=1e-10)
  pm_ <- spde_posterior_mean(theta, alpha = 2, graph = graph.temp, obs = FALSE)
  pm2_ <- Sigma.overdetermined[,index.obs]%*%solve(Sigma.Y,graph.temp$y)
  expect_equal(as.matrix(pm_),as.matrix(pm2_), tolerance=1e-10)
})
