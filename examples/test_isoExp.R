library(sp)
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))
graph <- metric_graph$new(lines = Lines)

graph$build_mesh(h = 0.05)


kappa <- 10
sigma <- 1
sigma_e <- 0.1
theta <-  c(sigma_e, sigma, kappa)

#sample some observation locations
n.obs.per.edge <- 100
PtE <- NULL
for(i in 1:graph$nE){
  if(0){
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge)))
  } else {
    PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge),
                            seq(from=0,to = 1,
                                length.out = n.obs.per.edge+2)[2:(n.obs.per.edge+1)]))
  }

}

D <- graph$compute_resdist(PtE)
Sigma <- sigma^2*exp(-kappa*D)
R <- chol(Sigma)
u <- as.vector(t(R)%*%rnorm(dim(R)[1]))
y <- u + sigma_e*rnorm(n.obs.per.edge * graph$nE)

graph$add_PtE_observations(y,PtE, normalized = TRUE)
graph$plot(data = TRUE)
graph$compute_resdist()
sigma_e_start <- sigma_e
sigma_start <- sigma
kappa_start <- kappa
theta0 <- c(sigma_e_start, sigma_start, kappa_start)
res.exp <- optim(log(theta0),
                 function(x) -likelihood_graph_covariance(exp(x), graph,
                                                          model = "isoExp"))
theta.exp <- exp(res.exp$par)

res.exp <- optim(log(theta.exp),
                 function(x) -likelihood_graph_covariance(exp(x), graph,
                                                          model = "isoExp"))
theta.exp <- exp(res.exp$par)

u_est_exp <- posterior_mean_covariance(theta.exp, graph, model = "isoExp")
graph$plot(X = u_est_exp, X_loc = graph$PtE)

graph$compute_resdist_mesh()

C.exp <- theta.exp[2]^2*exp(-theta.exp[3]*graph$mesh$res_dist[18,])
C.true <- sigma^2*exp(-kappa*graph$mesh$res_dist[18,])

p <- graph$plot_function_mesh(C.true, plotly = TRUE)
graph$plot_function_mesh(C.exp, plotly = TRUE, p = p, color = 'rgb(100,0,0)')
