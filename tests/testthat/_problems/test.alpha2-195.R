# Extracted from test.alpha2.R:195

# test -------------------------------------------------------------------------
set.seed(13)
nt <- 40
kappa <- 0.3
sigma_e <- 0.1
sigma   <- 1
theta <-  c(sigma_e,sigma,kappa)
edge2 <- rbind(c(30, 80), c(140, 80))
edge1 <- rbind(c(30, 00), c(30, 80))
edges <- list(edge1, edge2)
graph <- metric_graph$new(edges = edges)
Q <- spde_precision(kappa = kappa, tau = 1/sigma,
                      alpha = 2, graph = graph, BC = 1)
graph$buildC(2, FALSE)
Qmod <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
Qtilde <- Qmod
Qtilde <- Qtilde[-c(1:2),-c(1:2)]
R <- Cholesky(Qtilde,LDL = FALSE, perm = TRUE)
V0 <- as.vector(Matrix::solve(R, Matrix::solve(R,rnorm(6), system = 'Lt')
                        , system = 'Pt'))
u_e <- t(graph$CoB$T) %*% c(0, 0, V0)
X <- c()
for(i in 1:length(graph$edge_lengths)){
    X <- rbind(X,cbind(MetricGraph:::sample_alpha2_line(kappa = kappa,
                                          tau = 1/sigma,
                                          sigma_e = sigma_e,
                                          u_e = u_e[4*(i-1) +1:4],
                                          l_e = graph$edge_lengths[i],
                                          nt = nt),i))
  }
X[,2] <- X[,2] + sigma_e*rnorm(2*nt)
df_test <- data.frame(y = X[,2], edge_number = X[,3], distance_on_edge = X[,1])
graph$add_observations(data = df_test, normalized = FALSE)
graph$buildC(2, FALSE)
lik <- -MetricGraph:::likelihood_alpha2(theta = theta, graph = graph, data_name = "y", 
                             X_cov = NULL, repl = NULL, BC = 1, parameterization = "spde")
graph2 <- graph$clone()
graph2$observation_to_vertex()
graph2$buildC(2, FALSE)
lik2 <-MetricGraph:::likelihood_graph_covariance(graph = graph2,
                                     model = "WM2", repl = NULL, y_graph = graph2$get_data()[["y"]],
                                     log_scale = FALSE, X_cov = NULL)
lik2 <- lik2(exp(theta))
