#library(MetricGraph)
set.seed(1)
kappa <- 1
sigma <- 1.3
tau <- 1/sigma
sigma_e <- 0.1
alpha <- 1
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(1,0),c(2,1))
edge3 <- rbind(c(1,0),c(2,-1))
edge4 <- rbind(c(2,-1),c(3,0))
edge5 <- rbind(c(2,-1),c(3,-2))
edges = list(edge1, edge2, edge3, edge4, edge5)
graph <- metric_graph$new(edges = edges)

n.obs.per.edge <- 25
PtE <- NULL
for(i in 1:graph$nE){
  #add locations sampled at random to each edge
  PtE <- rbind(PtE, cbind(rep(i, n.obs.per.edge), runif(n.obs.per.edge)))
}


####
# sampling for alpha 1
###
####
# sampling edge vertices
####
Q_edges <- MetricGraph:::Qalpha1_edges(c(tau,kappa), graph, w = 0,BC=1, build=TRUE)
graph$buildDirectionalConstraints(alpha = 1)
n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tc <- graph$CoB$T[-ind.const, ]
Q_T <- Matrix::forceSymmetric(Tc%*%Q_edges%*%t(Tc))
R <- Matrix::Cholesky(Q_T,LDL = FALSE, perm = TRUE)
Z <- rnorm(dim(R)[1])
u_T <- as.vector(Matrix::solve(R, Matrix::solve(R,Z,
                                                system = 'Lt'), system = 'Pt'))
u_e <- t(Tc)%*%u_T

####
# sampling given edge vertices
####
#simulate
PtY <- rep(0, length(PtE[,1]))
for(e in unique(PtE[,1])){
  l <- graph$edge_lengths[e]
  obs.id = PtE[,1]==e
  t_ <- PtE[obs.id, 2]
  samp <- MetricGraph:::sample_alpha1_line(kappa = kappa, tau = tau,
                                          u_e = u_e[2*(e-1)+(1:2)], t = t_,
                                          l_e = l)
  u_y <- samp[,2]
  PtY[obs.id] <- u_y + sigma_e*rnorm(length(u_y))
}
df_data <- data.frame(y = PtY, edge_number = PtE[,1],
                      distance_on_edge = PtE[,2])
graph$add_observations(data = df_data, normalized = TRUE)
print(graph$plot(data = "y"))
theta_true <- c(tau,kappa,sigma_e)
#MetricGraph:::likelihood_alpha1_directional(log(theta_true),graph, parameterization="spde",data_name='y')
res <- graph_lme(y ~ -1, graph = graph, model = 'wmd1', optim_method='Nelder-Mead')

summary(res)


graph$observation_to_vertex()
graph$buildDirectionalConstraints(alpha = 1)
Data_full <- graph$get_data()
negg.likelihood <- function(theta){
  tau   <- exp(theta[1])
  kappa <- exp(theta[2])
  sigma_e <- exp(theta[3])
  Q_edges <- MetricGraph:::Qalpha1_edges(c(tau,kappa), graph, w = 0,BC=1, build=TRUE)
  n_const <- length(graph$CoB$S)
  ind.const <- c(1:n_const)
  Tc <- graph$CoB$T[-ind.const, ]
  Q_T <- Matrix::forceSymmetric(Tc%*%Q_edges%*%t(Tc))
  Sigma_big <- as.matrix(t(Tc)%*%solve(Q_T)%*%Tc)
  index <- 2*(Data_full$.edge_number - 1) + (Data_full$.distance_on_edge> 1-0.001) + 1
  Sigma_y <- Sigma_big[index,index]
  diag(Sigma_y) <- diag(Sigma_y) + sigma_e^2
  R <- chol(Sigma_y)
  v <- solve(t(R),Data_full$y)
  lik  <- -sum(log(diag(R))) - sum(v^2)/2 - 0.5*length(PtY)*log(2*pi)
  return(-lik)
}

res.lik <- optim(log(c(tau,kappa,sigma_e)), negg.likelihood)
theta_est <- c(res$coeff$random_effects, res$coeff$measurement_error)
negg.likelihood(log(theta_est))
cat('true:(tau, kappa, sigma_e) = ',round(c(tau,kappa,sigma_e),3),'\n')
cat('lm  :(tau, kappa, sigma_e) = ',round((theta_est),3),'\n')
cat('lik :(tau, kappa, sigma_e) = ',round(exp(res.lik$par),3),'\n')
