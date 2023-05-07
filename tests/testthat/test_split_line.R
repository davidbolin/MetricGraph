test_that("Testing if splitting edges, and adding obsevations casuse error", {

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line4),ID="3"),
                              Lines(list(line3),ID="4")))
graph <- metric_graph$new(lines = Lines)

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
Spoints = SpatialPoints(cbind(xc, yc))
Spoints = SpatialPointsDataFrame(Spoints,  data.frame(a=1:4))

df_2 <- data.frame(y = c(1,2,3,4))
graph$add_observations(Spoints, data = df_2)
expect_equal(42,42)
})

