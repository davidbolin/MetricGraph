test_that("Testing if splitting edges, and adding obsevations casuse error", {

edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(-1,1))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
edge4 <- cbind(sin(theta),1+ cos(theta))
edges = list(edge1, edge2, edge3, edge4)
graph <- metric_graph$new(edges = edges)

graph$.__enclos_env__$private$split_edge(1,0.5)

#graph$LtE
graph$plot()
PtE <- cbind(c(1,2),c(0.2,0.8))
df_graph <- data.frame(y = c(1,2), edge_number = PtE[,1],
                        distance_on_edge=PtE[,2])
graph$add_observations(data = df_graph)
graph$.__enclos_env__$private$split_edge(1,0.5)
graph$plot(data="y")

xc = c(-0.15,-0.51, 0.35,0.45)
yc = c(0.4, 0.25,0,0.4)
Spoints = sp::SpatialPoints(cbind(xc, yc))
Spoints = sp::SpatialPointsDataFrame(Spoints,  data.frame(a=1:4))

df_2 <- data.frame(y = c(1,2,3,4))
graph$add_observations(Spoints, data = df_2)
expect_equal(42,42)
})

