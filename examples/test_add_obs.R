library(sp)
library(GPGraph)
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line4),ID="3"),
                              Lines(list(line3),ID="4")))
graph <- metric_graph$new(Lines = Lines)
obs.per.edge <- 4
obs.loc <- NULL
for(i in 1:(graph$nE)) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), runif(obs.per.edge) *
                           graph$edge_lengths[i]))
}
y <- rep(NA, obs.per.edge * graph$nE)
graph$add_observations2(y,obs.loc)
graph$observation_to_vertex()
graph$plot(line_width = 0.3)
