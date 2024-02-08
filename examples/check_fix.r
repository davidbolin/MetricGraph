pems_graph <- metric_graph$new(edges=pems$edges, longlat=TRUE)
pems_graph$add_observations(data=pems$data, normalized=TRUE)
pems_graph$prune_vertices()
pems_graph$build_mesh(h=0.1)

res <- graph_lme(y ~ -1, graph = pems_graph, model = 'GL1')

summary(res)


res <- graph_lme(y ~ -1, graph = pems_graph, model = 'isoExp')

pred_aug <- augment(res, data.frame(edge_number = pems_graph$mesh$PtE[,1],
                        distance_on_edge = pems_graph$mesh$PtE[,2]), normalized = TRUE,
                        se_fit = TRUE, no_nugget = TRUE, check_euclidean = TRUE)


 pred_aug <- augment(res, normalized = TRUE,
                            se_fit = TRUE, no_nugget = FALSE, check_euclidean = TRUE)

cov_mat <- function(h,theta){
    nu <- 1.5
    sigma <- theta[1]
    kappa <- theta[2]
    (sigma^2/(2^(nu - 1) * gamma(nu))) * ((kappa * abs(h))^nu) * 
            besselK(kappa * abs(h), nu)
        }


    res <- graph_lme(y ~ -1, graph = pems_graph, model = list(type = "isoCov", cov_function = cov_mat), model_options = list(start_par_vec = c(1,1)))


edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(-1,1))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
edge4 <- cbind(sin(theta),1+ cos(theta))
edges = list(edge1, edge2, edge3, edge4)
graph <- metric_graph$new(edges = edges)
graph$plot()

graph$clear_observations()
range <- 0.15
sigma <- 2
sigma_e <- 0.1
theta <-  c(sigma_e, sigma, kappa)

n.obs.per.edge <- 75
PtE <- NULL
for(i in 1:graph$nE){
  #add locations sampled at random to each edge
  PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge)))
}

u <- sample_spde(range = range, sigma = sigma, alpha = 1,
                 graph = graph, PtE = PtE)

beta <- c(2,1)

X_cov <- cbind(1, runif(nrow(PtE)))

y <- X_cov %*% beta +  u + sigma_e*rnorm(n.obs.per.edge * graph$nE)


df_graph <- data.frame(y=y, x1 = X_cov[,2], edge_number = PtE[,1], distance_on_edge=PtE[,2])
graph$add_observations(data=df_graph, normalized = TRUE)
fit_cov <- graph_lme(y ~ x1, graph = graph, model = "WM1")
