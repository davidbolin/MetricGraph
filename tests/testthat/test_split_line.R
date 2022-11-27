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

graph$plot()

graph$split_edge(1,0.5)

#graph$LtE
graph$plot()
PtE <- cbind(c(1,2),c(0.2,0.8))
df_graph <- data.frame(y = c(1,2), edge_number = PtE[,1],
                        distance_on_edge=PtE[,2])
graph$add_PtE_observations(data_frame = df_graph)
graph$split_edge(1,0.5)
graph$plot()

xc = c(-0.5,-0.5, 0.75,0.45)
yc = c(1, 0.25,0,0)
Spoints = SpatialPoints(cbind(xc, yc))
Spoints = SpatialPointsDataFrame(Spoints,  data.frame(a=1:4))

Spoints2 = Spoints = SpatialPoints(cbind(xc+1, yc+1))
df_2 <- data.frame(y = c(1,2,3,4))
graph$add_observations(Spoints, data_frame = df_2)
expect_equal(42,42)
})

