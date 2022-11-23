#build graph Figure covariance between points

library(GPGraph)

theta <- c(2, 5) #kappa, sigma
P <- rbind(c(0,0),
           c(1,0),
           c(2,0),
           c(1,1),
           c(2,1),
           c(0,1),
           c(1,-1))

#specify edges
E <- rbind(c(1,2),
           c(2,3),
           c(2,4),
           c(4,5),
           c(5,3),
           c(1,6),
           c(4,6),
           c(2,7))

graph <-  metric_graph$new(V=V,E=E)
graph$plot()
#if Line is null assume straight line
#If we have a line keep
#compute covarians for a given point
P1 <- c(1,0.1) #edge 1, 0.5 length in

C <- covariance.point.to.graph.exp(P1, theta, graph, n.p = 50)

gg <- plot_X_to_graph(C, graph)
print(gg)

gg <- plot_X_to_graph(C, graph, flat = FALSE, color='rgb(0,0,200)')
print(gg)

