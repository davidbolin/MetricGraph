library(MetricGraph)
 
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(-1,1))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
edge4 <- cbind(sin(theta),1+ cos(theta))
edges = list(edge1, edge2, edge3, edge4)
graph_bru <- metric_graph$new(edges = edges)
obs_per_edge <- 50
obs_loc <- NULL
for(i in 1:(graph_bru$nE)) {
  obs_loc <- rbind(obs_loc,
                   cbind(rep(i,obs_per_edge), 
                   runif(obs_per_edge)))
}

sigma <- 2
alpha <- 1
nu <- alpha - 0.5
r <- 0.15 # r stands for range
 
u <- sample_spde(range = r, sigma = sigma, alpha = alpha,
                     graph = graph_bru, PtE = obs_loc)

n_obs <- length(u)
sigma.e <- 0.1
 
y <- u + sigma.e * rnorm(n_obs)

 df_graph <- data.frame(y = y, edge_number = obs_loc[,1],
                          distance_on_edge = obs_loc[,2])
# Adding observations and turning them to vertices
graph_bru$add_observations(data = df_graph, normalized=TRUE)

library(INLA)
library(inlabru)
kappa_start <- 15
tau_start <- 1
spde_model_bru <- graph_spde(graph_bru, directional=TRUE, parameterization = "spde", 
stationary_endpoints = "all", start_kappa=kappa_start,start_tau=tau_start, alpha=1)

tmp <- inla.cgeneric.q(spde_model_bru)

Q0 <- tmp$Q

kappa0 <- exp(tmp$theta[2])

tau0 <- exp(-tmp$theta[1])

Q1 <- MetricGraph:::Qalpha1_edges(theta = c(tau0,kappa0), graph = spde_model_bru$graph_spde, w=0,
stationary_points = "all")

Tc <- spde_model_bru$Tc
Qnew <- Tc%*%Q1%*%t(Tc)

sum((Q0-Qnew)^2)




