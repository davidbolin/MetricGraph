#integrating out the effect of edges of degree one
edge1 <- rbind(c(1,0),c(0,0))
edge2 <- rbind(c(1+sqrt(0.5),sqrt(0.5)),c(1,0))
edge3 <- rbind(c(1+sqrt(0.5),-sqrt(0.5)),c(1,0))
edges = list(edge1,edge2,edge3)
graph <- metric_graph$new(edges = edges)
graph$plot(direction = T)
graph$build_mesh(h=0.01)
kappa <- 1
tau   <- 1
alpha <- 1
graph$setDirectionalWeightFunction(f_in = function(x){(x/sum(x))})
Eu <- MetricGraph:::posterior_mean_obs_alpha1(c(0,tau, kappa),
                                              graph,
                                              resp, #resp must be in the graph's internal order
                                              PtE_resp,
                                              graph$mesh$PtE,
                                              type = "PtE",
                                              directional = T)
graph$plot_function(X = Eu, plotly = TRUE)
graph$buildDirectionalConstraints(alpha = 1)
Q <- MetricGraph:::Qalpha1_edges(c(tau,kappa),
                                 graph,
                                 w = 0,
                                 BC=1, build=TRUE)
n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tc <- as.matrix(graph$CoB$T[-ind.const, ])
Q_ <- Matrix::forceSymmetric(Tc%*%Q%*%t(Tc))
Cov <- t(Tc)%*%as.matrix(solve(Q_))%*%Tc

# integrating out the two end points to generate a smaller graf!
Q2 <- Q[c(1:2,4,6),c(1:2,4,6)]
Q2[3,3] <- Q2[4,4] <- Q[4,4] - Q[3,4]^2/Q[3,3] #boundary adjustement
Tc2 <- Tc[-c(4,5),-c(3,5)]
Q_2 <- Matrix::forceSymmetric(Tc2%*%Q2%*%t(Tc2))
Cov2 <- t(Tc2)%*%as.matrix(solve(Q_2))%*%Tc2
