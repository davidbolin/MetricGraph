library(GPGraph)
library(sp)
library(INLA)

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

tmp <- gpgraph_spde(graph)

inla.cgeneric.q(tmp)

u <- sample_spde_mesh(kappa = 10, sigma = 2, graph = graph)

A <- diag(c(1,1,1,1))



data_tmp <- data.frame(y = u)

index <- list()

index[["field"]] <- 1:ncol(A)

stk.dat <- inla.stack(data = list(y=u), 
                        A = A, 
                        effects = c(
      index,
      list(Intercept = 1)
    ))

f.s <- y ~ -1 + Intercept + f(field, model = tmp)

inla(f.s, data = inla.stack.data(stk.dat), verbose = TRUE)



line1 <- Line(rbind(c(0,0),c(0,1)))
line2 <- Line(rbind(c(0,1),c(0,2.5)))
Lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2")))
graph <- metric_graph$new(Lines = Lines)



