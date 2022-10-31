#build graph Figure covariance between points

library(GPGraph)

theta <- c(0.1, 1) #kappa, sigma
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

graph <-  graph.obj$new()
graph$EtV   <- cbind(1:dim(E)[1],E)
graph$V <- P
graph$El    <- sqrt((graph$V[ graph$EtV[,3], 1] - graph$V[ graph$EtV[, 2], 1])^2 +
                    (graph$V[ graph$EtV[,3], 2] - graph$V[ graph$EtV[, 2], 2])^2)


#if Line is null assume straight line
#If we have a line keep
#compute covarians for a given point
P1 <- c(8,0.1) #edge 1, 0.5 length in

C <- covariance.point.to.graph(P1, theta, graph, n.p = 50)
#g <- plot_straight_curve(C[C[,1]==1,2:3], graph$V[graph$EtV[1,2:3],],)
gg <- plot_X_to_graph(C, graph)
print(gg)
