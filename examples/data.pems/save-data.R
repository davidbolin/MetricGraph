Lines <- read_sf('examples/data.pems/lines.shp')
lines <- as_Spatial(Lines)

EtV <- read.csv('examples/data.pems/E.csv',header=T, row.names = NULL)
PtE <- read.csv('examples/data.pems/PtE.csv',header=T, row.names = NULL)
PtE[,1] <- PtE[,1] + 1
Y <- read.csv('examples/data.pems/Y.csv',header=T, row.names = NULL)
Y <- colMeans(as.matrix(Y[,-1]))
edge_length_m <- EtV[,4]
PtE[,2] = PtE[,2]/edge_length_m[PtE[,1]]

pems <- list(lines = lines, PtE = PtE, Y = Y)
save(pems, file = "pems.rda")
