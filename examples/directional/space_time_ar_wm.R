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

dim_time <- 15

indep_field <- matrix(nrow = nrow(obs_loc), ncol = dim_time)

for(ii in 1:dim_time){
    indep_field[,ii] <- sample_spde(range = r, sigma = sigma, alpha = alpha,
                     graph = graph_bru, PtE = obs_loc)
}


# correlation AR
rho_ar <- 0.65


## stationary initial condition
field_time <- indep_field
for (ii in 2:dim_time) {
  field_time[, ii] <- rho_ar * field_time[, ii - 1] + sqrt(1 - rho_ar^2) * indep_field[, ii]
}

sigma_e <- 0.1
n_obs <- nrow(field_time)

beta_1 <- obs_loc[,1]
beta_2 <- obs_loc[,2]

y <- 2 - beta_1 + beta_2 + field_time + rnorm(n_obs * dim_time, 0, sigma_e)

df_graph <- data.frame(y = as.vector(y), edge_number = rep(obs_loc[,1], dim_time),
                          distance_on_edge = rep(obs_loc[,2],dim_time),
                          beta_1 = rep(beta_1, dim_time),
                          beta_2 = rep(beta_2, dim_time),
                          time = rep(1:dim_time, each = n_obs),
                          dummy = 1:length(beta_1))

# Adding observations and turning them to vertices
graph_bru$add_observations(data = df_graph, normalized=TRUE, group = "time")

spde_model_bru <- graph_spde(graph_bru, alpha=1, directional=FALSE)

cmp <- y ~ -1 + Intercept(1) + beta_1(beta_1) + beta_2(beta_2) + 
    field(loc, model = spde_model_bru, group = time, control.group = list(model = 'ar1')) 

data_spde <- graph_data_spde(graph_spde = spde_model_bru, 
                            loc_name = "loc", group_col = "time")

library(inlabru)

spde_bru_fit <-
    bru(cmp, data=data_spde[["data"]])

spde_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model_bru)

summary(spde_result)
summary(spde_bru_fit)

############################
### Alpha = 2  #############
############################

sigma <- 2
alpha <- 2
nu <- alpha - 0.5
r <- 0.15 # r stands for range


indep_field <- matrix(nrow = nrow(obs_loc), ncol = dim_time)

for(ii in 1:dim_time){
    indep_field[,ii] <- sample_spde(range = r, sigma = sigma, alpha = alpha,
                     graph = graph_bru, PtE = obs_loc)
}


# correlation AR
rho_ar <- 0.65


## stationary initial condition
field_time <- indep_field
for (ii in 2:dim_time) {
  field_time[, ii] <- rho_ar * field_time[, ii - 1] + sqrt(1 - rho_ar^2) * indep_field[, ii]
}

sigma_e <- 0.1
n_obs <- nrow(field_time)

beta_1 <- obs_loc[,1]
beta_2 <- obs_loc[,2]

y <- 2 - beta_1 + beta_2 + field_time + rnorm(n_obs * dim_time, 0, sigma_e)

df_graph <- data.frame(y = as.vector(y), edge_number = rep(obs_loc[,1], dim_time),
                          distance_on_edge = rep(obs_loc[,2],dim_time),
                          beta_1 = rep(beta_1, dim_time),
                          beta_2 = rep(beta_2, dim_time),
                          time = rep(1:dim_time, each = n_obs),
                          dummy = 1:length(beta_1))

# Adding observations and turning them to vertices
graph_bru$add_observations(data = df_graph, normalized=TRUE, group = "time", clear_obs=TRUE)

spde_model_bru <- graph_spde(graph_bru, alpha=2, directional=FALSE)

cmp <- y ~ -1 + Intercept(1) + beta_1(beta_1) + beta_2(beta_2) + 
    field(loc, model = spde_model_bru, group = time, control.group = list(model = 'ar1')) 

data_spde <- graph_data_spde(graph_spde = spde_model_bru, 
                            loc_name = "loc", group_col = "time")

library(inlabru)

spde_bru_fit <-
    bru(cmp, data=data_spde[["data"]])

summary(spde_bru_fit)

spde_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model_bru)

summary(spde_result)


######################################
## alpha = 1 + directional = TRUE ####
######################################


sigma <- 2
alpha <- 1
nu <- alpha - 0.5
r <- 0.15 # r stands for range

kappa <- sqrt(8 * nu) / r
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
(4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))

dim_time <- 15

spde_model_bru <- graph_spde(graph_bru, alpha=1, directional=TRUE)

Q_dir <- MetricGraph:::Qalpha1_edges(theta = c(tau,kappa), graph = spde_model_bru$graph_spde, w = 0, stationary_points = "all")

Tc <- spde_model_bru$Tc
Q_dir<- Tc%*%Q_dir%*%t(Tc)

indep_field <- matrix(nrow = nrow(obs_loc), ncol = dim_time)
ord <- spde_model_bru$graph_spde$.__enclos_env__$private$data$dummy

LQ <- chol(Q_dir)

for(ii in 1:dim_time){
    Z <- rnorm(nrow(Q_dir)) 
    indep_field[,ii] <- as.vector((spde_model_bru$A)%*%solve(LQ, Z))
    indep_field[ord,ii] <- indep_field[,ii]
}

# correlation AR
rho_ar <- 0.65

## stationary initial condition
field_time <- indep_field
for (ii in 2:dim_time) {
  field_time[, ii] <- rho_ar * field_time[, ii - 1] + sqrt(1 - rho_ar^2) * indep_field[, ii]
}

sigma_e <- 0.1
n_obs <- nrow(obs_loc)

beta_1 <- obs_loc[,1]
beta_2 <- obs_loc[,2]

y <- 2 - beta_1 + beta_2 + field_time + rnorm(n_obs * dim_time, 0, sigma_e)

df_graph <- data.frame(y = as.vector(y), edge_number = rep(obs_loc[,1], dim_time),
                          distance_on_edge = rep(obs_loc[,2],dim_time),
                          beta_1 = rep(beta_1, dim_time),
                          beta_2 = rep(beta_2, dim_time),
                          time = rep(1:dim_time, each = n_obs),
                          dummy = 1:length(beta_1))

# Adding observations and turning them to vertices
graph_bru$add_observations(data = df_graph, normalized=TRUE, group = "time", clear_obs=TRUE)

spde_model_bru <- graph_spde(graph_bru, alpha=1, directional=TRUE)

cmp <- y ~ -1 + Intercept(1) + beta_1(beta_1) + beta_2(beta_2) + 
    field(loc, model = spde_model_bru, group = time, control.group = list(model = 'ar1')) 

data_spde <- graph_data_spde(graph_spde = spde_model_bru, 
                            loc_name = "loc", group_col = "time")

library(inlabru)

spde_bru_fit <-
    bru(cmp, data=data_spde[["data"]])

spde_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model_bru)

summary(spde_result)
summary(spde_bru_fit)
