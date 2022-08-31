
graphics.off()
library(GPGraph)
library(sp)
library(maptools)
dbh = 20
kappa <- 0.5
sp.dbh = spCircle(dbh/2, centerPoint=c(x=30,y=80), spID='tree.1')
line.cirlce <- Line(sp.dbh$spCircle@polygons[[1]]@Polygons$circPlot@coords)
line.line <- Line(rbind(c(30,80),c(40,80)))
graph <- graph.obj$new(sp::SpatialLines(list(Lines(list(line.cirlce),ID="1"))))
graph2 <-  graph.obj$new(sp::SpatialLines(list(Lines(list(line.cirlce),ID="1"),
                                               Lines(list(line.line),ID="2"))))


Q <- Q.exp(kappa, graph$V, graph$EtV, graph$El)
Q2 <- Q.exp(kappa, graph2$V, graph2$EtV, graph2$El)

#https://tmieno2.github.io/R-as-GIS-for-Economists/

#simulate Expontial on a Line given end points
#simulate Expontial on a Line given end points and observations

Q_ <- Q.exp.line(kappa, c(0,30,64))

xc = c(1.2,31)
yc = c(1.5,80)
Spoints = SpatialPoints(cbind(xc, yc))
Spoints = SpatialPointsDataFrame(Spoints,  data.frame(a=1:2))
SP <- snapPointsToLines(Spoints, graph2$Lines)

graph2$add_observations(Spoints, c(1,2))

X <- sample.line.expontial(c(1,1,kappa),
                           c(0,0),
                           Line = graph2$Lines[1,],
                           graph2$El[1])
fig <- plot_curve(X, graph2$Lines[1,])
print(fig)

X2 <- sample.line.expontial(c(0.0001,1,kappa),
                           c(0,0),
                           Line = graph2$Lines[1,],
                           graph2$El[1],
                           py = graph2$PtE[1,2],
                           y  = 10
                           )

fig <- plot_curve(X2, graph2$Lines[1,])
print(fig)
