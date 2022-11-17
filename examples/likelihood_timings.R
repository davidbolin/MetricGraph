library('sf')
library(GPGraph)

graph.type = 2
if(graph.type == 1){
  Lines <- as_Spatial(read_sf('examples/data.pems/lines.shp'))
} else if(graph.type == 2){
  line1 <- Line(rbind(c(0,0),c(1,0)))
  line2 <- Line(rbind(c(0,0),c(0,1)))
  line3 <- Line(rbind(c(0,1),c(-1,1)))
  line4 <- Line(rbind(c(-1,1),c(-1,0)))
  line5 <- Line(rbind(c(-1,0),c(0,0)))
  Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                                Lines(list(line2),ID="2"),
                                Lines(list(line3),ID="3"),
                                Lines(list(line4),ID="4"),
                                Lines(list(line5),ID="5")))
} else {
  V <- rbind(c(0,0),c(1,0),c(0,1),c(1,1))
  E = rbind(c(1,2),c(1,3), c(1,4), c(3,4))
}

kappa <- 20
sigma <- 1
sigma_e <- 0.1
alpha <- 2
theta <-  c(sigma_e, kappa, sigma)

N <- 5
if(graph.type == 1){
  n.obs <- c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000)
} else {
  n.obs <- c(40, 80, 120, 160, 200, 240, 280, 320, 360, 400)
}


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
    if(graph.type == 1 || graph.type == 2){
      graph <-  metric_graph$new(Lines)
    } else {
      graph <-  metric_graph$new(P = V, E = E)
    }

    if(graph.type == 2){
      nE <- graph$nE
      edge.pos <- seq(from=0,to=1,length.out = n.obs[j]+2)[2:(n.obs[j]+1)]
      edge.num <- rep(1:graph$nE,each=n.obs[j])
      PtE <- cbind(edge.num,graph$edge_lengths[edge.num]*edge.pos)
      y <- rnorm(n.obs[j]*nE)
    } else if (graph.type == 1) {
      nE <- 1
      p.edge <- graph$edge_lengths/sum(graph$edge_lengths)
      edge.samp <- sample(x = 1:graph$nE, n.obs[j], replace = TRUE, prob = p.edge)
      edge.pos <- rep(0,n.obs[j])
      edge.pos[!duplicated(edge.samp)] <- 0.5
      if(sum(duplicated(edge.samp))>0){
        d.u <- unique(edge.samp[duplicated(edge.samp)])
        for(k in 1:length(d.u)) {
          edge.pos[edge.samp==d.u[k]] = seq(from=0,to=1,length.out = sum(edge.samp==d.u[k])+2)[2:(sum(edge.samp==d.u[k])+1)]
        }
      }
      PtE <- cbind(edge.samp,graph$edge_lengths[edge.samp]*edge.pos)
      y <- rnorm(n.obs[j])
    }


    graph$add_observations2(y = y, PtE = PtE)
    graph$buildC(2)
    lik[i] <- system.time(likelihood_graph_spde(theta,graph, alpha = alpha, version = 1))[["elapsed"]]
    graph$observation_to_vertex()
    if(alpha == 1){
      lik.v2[i] <- system.time(likelihood_graph_spde(theta, graph, alpha = alpha, version = 2))[["elapsed"]]
    } else {
      lik.v2[i] <- system.time(likelihood_graph_spde(theta, graph, alpha = alpha, version = 1))[["elapsed"]]
    }

    if(alpha == 1) {
      lik.cov[i] <- system.time(likelihood_graph_covariance(theta, graph, model = "alpha1"))[["elapsed"]]
    } else {
      lik.cov[i] <- system.time(likelihood_graph_covariance(theta, graph, model = "alpha2"))[["elapsed"]]
    }
  }
  results$v1[j] = mean(lik)
  results$v2[j] = mean(lik.v2)
  results$co[j] = mean(lik.cov)
  results$v1.std[j] = std(lik)
  results$v2.std[j] = std(lik.v2)
  results$co.std[j] = std(lik.cov)

}

c = 0.1
print(results)
plot(nE*n.obs,results$v1,col=1,
     type="l",
     ylim = c(0,max(c(max(results$v1),max(results$v2),max(results$co*c)))),
     xlab = "n",
     ylab = "Time for likelihood evaluation (s)",
     lwd = 2)
points(nE*n.obs,results$v1,col=1, pch=1)
lines(nE*n.obs,results$v2, col = 2,lwd = 2)
points(nE*n.obs,results$v2,col=2, pch=2)
lines(nE*n.obs,results$co*c, col = 3,lwd = 2)
points(nE*n.obs,results$co*c,col=3, pch=5)
legend("topleft", legend = c("Bridge method",
                             "Extended graph",
                             "Covariance-based (times 0.1)"),
       col=1:3, lty=1, pch=c(1,2,5),lwd=2)

