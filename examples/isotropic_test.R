d1 <- 1
d2 <- 2
d3 <- 0.1
d4 <- 1
sigma = 2
kappa = 2
r1 <- sigma^2*exp(-kappa*d1)
r2 <- sigma^2*exp(-kappa*d2)
r3 <- sigma^2*exp(-kappa*d3)
r4 <- sigma^2*exp(-kappa*d4)
Sigma <- rbind(c(sigma^2, r1, r2, r3, r3*r4/sigma^2),
               c(r1, sigma^2, r1*r2/sigma^2, r1*r3/sigma^2, r1*r3*r4/sigma^4),
               c(r2, r1*r2/sigma^2, sigma^2, r2*r3/sigma^2, r2*r3*r4/sigma^4),
               c(r3, r1*r3/sigma^2, r2*r3/sigma^2, sigma^2, r4),
               c(r3*r4/sigma^2, r1*r3*r4/sigma^4, r2*r3*r4/sigma^4, r4, sigma^2))


A <- 1/(1*3)
Q1 <- exp(-2*kappa*d1)/(1-exp(-2*kappa*d1))
Q2 <- exp(-2*kappa*d2)/(1-exp(-2*kappa*d2))
Q3 <- exp(-2*kappa*d3)/(1-exp(-2*kappa*d3))
Q4 <- exp(-2*kappa*d4)/(1-exp(-2*kappa*d4))
q1 <- -exp(-kappa*d1)/(1-exp(-2*kappa*d1))
q2 <- -exp(-kappa*d2)/(1-exp(-2*kappa*d2))
q3 <- -exp(-kappa*d3)/(1-exp(-2*kappa*d3))
q4 <- -exp(-kappa*d4)/(1-exp(-2*kappa*d4))
Q <- sigma^(-2)*rbind(c(3*A+Q1+Q2+Q3, q1, q2, q3,0),
                      c(q1, 1+Q1, 0, 0, 0),
                      c(q2, 0, 1+Q2, 0, 0),
                      c(q3, 0, 0, 1+Q3+Q4, q4),
                      c(0,  0, 0, q4, 1+Q4))

solve(Sigma)-Q

#test fem implementation with generalized kirchhoff
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,0),c(0,-1)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3")))

line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,0),c(0,-1)))
line4 <- Line(rbind(c(0,1),c(1,0)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4")))

line1 <- Line(rbind(c(0,0),c(1,0), c(0,1), c(0,0)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1")))

graph <- metric_graph$new(lines = lines)
graph$plot()
h = 0.005
graph$build_mesh(h = h)
graph$compute_fem()
G <- graph$mesh$G
C <- graph$mesh$C
B <- graph$mesh$B
Cb<- Diagonal(0,n=dim(G)[1])
diag(Cb)[1:graph$nV] <- 2-graph$get_degrees()
Cb[1,1] <- tanh(kappa*(1+1+sqrt(2))/2)/2
kappa <- 1
Q <- (G + kappa*Cb + kappa^2*C)/(2*kappa)
graph$plot_function(diag(solve(Q)), vertex_size = 0)
