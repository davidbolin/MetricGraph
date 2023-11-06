library(ggplot2)
library(viridis)
graphics.off()
#simple example
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(1,0),c(2,1))
edge3 <- rbind(c(1,0),c(2,-1))
edge4 <- rbind(c(2,-1),c(3,0))
edge5 <- rbind(c(2,-1),c(3,-2))
edges = list(edge1, edge2, edge3, edge4, edge5)
graph <- metric_graph$new(edges = edges)
graph$build_mesh(h=1, continuous = FALSE, continuous.outs = TRUE, continuous.deg2 = TRUE)
print(graph$plot(mesh=TRUE, direction = TRUE))
tau = 0.5
kappa <- 1
Q= MetricGraph:::Qalpha1_v2(c(tau,kappa), graph, w = 0 ,BC = 2, build = TRUE)
