library('sf')
library(GPGraph)

if(0){
  Lines <- as_Spatial(read_sf('examples/data.pems/lines.shp'))
} else if(0){
  line1 <- Line(rbind(c(0,0),c(1,0)))
  line2 <- Line(rbind(c(0,0),c(0,1)))
  line3 <- Line(rbind(c(0,1),c(-1,1)))
  theta <- seq(from=pi,to=3*pi/2,length.out = 20)
  line4 <- Line(cbind(sin(theta),1+ cos(theta)))
  Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                                Lines(list(line2),ID="2"),
                                Lines(list(line3),ID="3"),
                                Lines(list(line4),ID="4")))
} else {
  V <- rbind(c(0,0),c(1,0),c(0,1),c(1,1))
  E = rbind(c(1,2),c(1,3), c(1,4), c(3,4))

}

kappa <- 10
sigma <- 1
sigma_e <- 0.1
theta <-  c(sigma_e, kappa, sigma)

N <- 2
n.obs <- c(100, 500, 1000, 1500)

results <- data.frame(v1 = rep(0,length(n.obs)),
                      v2 = rep(0,length(n.obs)),
                      co = rep(0,length(n.obs)),
                      v1.std = rep(0,length(n.obs)),
                      v2.std = rep(0,length(n.obs)),
                      co.std = rep(0,length(n.obs)))

for(j in 1:length(n.obs)){
  lik <- lik.v2 <- lik.cov <- rep(0,N)
  for(i in 1:N){
    cat(j,i,"\n")
    #graph <-  metric_graph$new(Lines)
    graph <-  metric_graph$new(P = V, E = E)
    p.edge <- graph$edge_lengths/sum(graph$edge_lengths)
    edge.samp <- sample(x = 1:graph$nE, n.obs[j], replace = TRUE, prob = p.edge)
    edge.pos <- runif(n.obs[j])
    PtE <- cbind(edge.samp,graph$edge_lengths[edge.samp]*edge.pos)
    y <- rnorm(n.obs[j])
    graph$add_observations2(y = y, PtE = PtE)

    lik[i] <- system.time(likelihood_graph_spde(theta,graph, alpha = 1, version = 1))[["elapsed"]]
    graph$observation_to_vertex()
    lik.v2[i] <- system.time(likelihood_graph_spde(theta, graph, alpha = 1, version = 2))[["elapsed"]]
    lik.cov[i] <- system.time(likelihood_graph_covariance(theta, graph, model = "alpha1"))[["elapsed"]]
  }
  results$v1[j] = mean(lik)
  results$v2[j] = mean(lik.v2)
  results$co[j] = mean(lik.cov)
  results$v1.std[j] = std(lik)
  results$v2.std[j] = std(lik.v2)
  results$co.std[j] = std(lik.cov)

}

print(results)
plot(n.obs,results$v1,col=1,type="l", ylim = c(0,max(c(max(results$v1),max(results$v2),max(results$co)))))
lines(n.obs,results$v2, col = 2)
lines(n.obs,results$co, col = 3)
