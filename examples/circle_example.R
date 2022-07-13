library(GPGraph)
library(sp)
dbh = 20
kappa <- 0.01
sp.dbh = spCircle(dbh/2, centerPoint=c(x=30,y=80), spID='tree.1')
line.cirlce <- Line(sp.dbh$spCircle@polygons[[1]]@Polygons$circPlot@coords)
graph <- graph.obj$new(line.cirlce)

Q <- Q.exp(kappa, graph$V, graph$EtV, graph$El)
