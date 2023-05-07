# testing if Kirchoff is unique with respect to direction
# Date: 2023-04-18
library(MetricGraph)
library(sp)
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,0),c(-1,0)))
line4 <- Line(rbind(c(0,-1),c(0,0)))
line5 <- Line(rbind(c(0,-1),c(-1,0)))
line6 <- Line(rbind(c(-1,0),c(0,1)))
line7 <- Line(rbind(c(0,1),c(1,0)))
line8 <- Line(rbind(c(1,0),c(0,-1)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4"),
                              Lines(list(line5),ID="5"),
                              Lines(list(line6),ID="6"),
                              Lines(list(line7),ID="7"),
                              Lines(list(line8),ID="8")))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1")))
#Lines = sp::SpatialLines(list(Lines(list(line5),ID="5"),
#                              Lines(list(line6),ID="6"),
#                              Lines(list(line7),ID="7"),
#                              Lines(list(line8),ID="8")))
graph <- metric_graph$new(lines = Lines)

Q0 <- Qalpha2(c(1,1,1), graph,w=0.25,BC=0)
graph$buildC(2,TRUE)

n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tc <- graph$CoB$T[-ind.const, ]
Sigmaw025 <- t(Tc)%*%solve(Tc%*%Q0%*%t(Tc))%*%Tc


Q0_ <- Qalpha2(c(1,1,1), graph,w=0.5,BC=0)
Sigmaw05 <- t(Tc)%*%solve(Tc%*%Q0_%*%t(Tc))%*%Tc
