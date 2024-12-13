library(MetricGraph)
 
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(-1,1))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
edge4 <- cbind(sin(theta),1+ cos(theta))
edges = list(edge1, edge2, edge3, edge4)
graph_bru <- metric_graph$new(edges = edges)
obs_per_edge <- 100
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

sigma_e <- 0.1
n_obs <- nrow(obs_loc)

beta_1 <- obs_loc[,1]
beta_2 <- obs_loc[,2]

y <- 2 - beta_1 + beta_2 + u + rnorm(n_obs, 0, sigma_e)

df_graph <- data.frame(y = as.vector(y), edge_number = obs_loc[,1],
                          distance_on_edge = obs_loc[,2],
                          beta_1 = beta_1,
                          beta_2 = beta_2,
                          dummy = 1:length(beta_1))

# Adding observations 
graph_bru$add_observations(data = df_graph, normalized=TRUE, clear_obs=TRUE)

spde_model_bru <- graph_spde(graph_bru, alpha=alpha, directional=FALSE)

data_spde <- graph_data_spde(graph_spde = spde_model_bru, 
                            loc_name = "loc")

cmp <- y ~ -1 + Intercept(1) + beta_1(beta_1) + beta_2(beta_2) + 
    field(loc, model = spde_model_bru) 

library(inlabru)

spde_bru_fit <-
    bru(cmp, data=data_spde[["data"]])

spde_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model_bru)

summary(spde_result)
summary(spde_bru_fit)


##################################
### alpha = 1 + Directional ######
##################################


sigma <- 2
alpha <- 1
nu <- alpha - 0.5
r <- 0.15 # r stands for range

kappa <- sqrt(8 * nu) / r
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
(4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))

spde_model_bru <- graph_spde(graph_bru, alpha=1, directional=TRUE)

Q_dir <- MetricGraph:::Qalpha1_edges(theta = c(tau,kappa), graph = spde_model_bru$graph_spde, w = 0, stationary_points = "all")

Tc <- spde_model_bru$Tc
Q_dir<- Tc%*%Q_dir%*%t(Tc)

Z <- rnorm(nrow(Q_dir))

LQ <- chol(Q_dir)
u <- solve(LQ, Z)

sigma_e <- 0.1
n_obs <- nrow(obs_loc)

beta_1 <- obs_loc[,1]
beta_2 <- obs_loc[,2]

ord <- spde_model_bru$graph_spde$.__enclos_env__$private$data$dummy

u_field <- spde_model_bru$A %*% u

u_field[ord] <- u_field

y <- 2 - beta_1 + beta_2 + u_field + rnorm(n_obs, 0, sigma_e)

df_graph <- data.frame(y = as.vector(y), edge_number = obs_loc[,1],
                          distance_on_edge = obs_loc[,2],
                          beta_1 = beta_1,
                          beta_2 = beta_2,
                          dummy = 1:length(beta_1))

# Adding observations 
graph_bru$add_observations(data = df_graph, normalized=TRUE, clear_obs=TRUE)

spde_model_bru <- graph_spde(graph_bru, alpha=1, directional=TRUE)

cmp <- y ~ -1 + Intercept(1) + beta_1(beta_1) + beta_2(beta_2) + 
    field(loc, model = spde_model_bru) 

data_spde <- graph_data_spde(graph_spde = spde_model_bru, 
                            loc_name = "loc")

library(inlabru)

spde_bru_fit <-
    bru(cmp, data=data_spde[["data"]])

spde_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model_bru)

summary(spde_result)

res <- graph_lme(y ~ beta_1 + beta_2, model = "WMD1", graph = graph_bru)
res
