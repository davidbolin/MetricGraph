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
line3 <- Line(rbind(c(1,0),c(1,1)))
line4 <- Line(rbind(c(0,1),c(1,1)))
#line5 <- Line(rbind(c(0,0),c(1,1)))
line6 <- Line(rbind(c(0,1),c(1,0)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4"),
 #                             Lines(list(line5),ID="5"),
                              Lines(list(line6),ID="6")))

graph <- metric_graph$new(lines = lines, tolerance = list(line_line=0.1))
graph$plot()

#tree
line1 <- Line(rbind(c(0,0),c(-1,1)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,0),c(1,1)))
line4 <- Line(rbind(c(-1,1),c(-1.5,2)))
line4 <- Line(rbind(c(-1,1),c(-1,2)))
line5 <- Line(rbind(c(0,1),c(-0.5,2)))
line6 <- Line(rbind(c(0,1),c(0,2)))
line7 <- Line(rbind(c(1,1),c(0.5,2)))
line8 <- Line(rbind(c(1,1),c(1,2)))
lines = sp::SpatialLines(list(Lines(list(line1),ID="1"),
                              Lines(list(line2),ID="2"),
                              Lines(list(line3),ID="3"),
                              Lines(list(line4),ID="4"),
                              Lines(list(line5),ID="5"),
                              Lines(list(line6),ID="6"),
                              Lines(list(line7),ID="7"),
                              Lines(list(line8),ID="8")))

graph <- metric_graph$new(lines = lines)
graph$plot()
h = 0.005
graph$build_mesh(h = h)
graph$compute_fem()
G <- graph$mesh$G
C <- graph$mesh$C
Cb<- Diagonal(0,n=dim(G)[1])
diag(Cb)[1:graph$nV] <- 2-graph$get_degrees()
kappa <- 0.1
d.max <-max(graph$get_degrees())
l.min <- min(graph$edge_lengths)
if(d.max<3){
  beta0 <- Inf
} else if (l.min < sqrt(2)/kappa){
  beta0 <- kappa^2*l.min*d.max/(4*(d.max-2))
} else {
  beta0 <- kappa*d.max/(2*sqrt(2)*(d.max-2))
}
beta <- kappa/3.75
Q <- (G + beta*Cb + kappa^2*C)/(2*kappa)
graph$plot_function(diag(solve(Q)), vertex_size = 0)
cat("kappa = ",kappa, "beta = ",beta, "beta0 = ", beta0, "min eig = ", min(eigen(Q)$values))



#line
len <- 1
h = 0.0005
graph <- metric_graph$new(lines = sp::SpatialLines(list(Lines(list(Line(rbind(c(0,0),
                                                                              c(len,0)))),ID="1"))))

graph$build_mesh(h = h)
graph$compute_fem()
kappa <- 2
cr <- 2
Cb<- Diagonal(0,n=dim(graph$mesh$G)[1])
diag(Cb)[1:graph$nV] <- kappa*c(-cr,1)
Q <- (graph$mesh$G + Cb + cr*kappa^2*graph$mesh$C)/(2*kappa)
min(eigen(Q)$values)

v <- rep(0,dim(Q)[1]);v[1] <- 1
cv1 <- solve(Q,v)
v <- rep(0,dim(Q)[1]);v[dim(Q)[1]] <- 1
cv2 <- solve(Q,v)
v <- rep(0,dim(Q)[1]);v[round(dim(Q)[1]/2)] <- 1
cv3 <- solve(Q,v)
p <- graph$plot_function(cv1,plotly=TRUE,line_color = "red",support_width = 0)
p <- graph$plot_function(cv2,p=p,plotly=TRUE,line_color = "blue",support_width = 0)
graph$plot_function(cv3,p=p,plotly=TRUE,line_color = "black",support_width = 0)

t <- graph$mesh$V[,1]
k <- 0
ek <- sin(k*pi*t/len) - k*pi*cos(k*pi*t/len)/(len*kappa)
ek <- exp(-kappa*t)
p <- graph$plot_function(ef[,dim(graph$mesh$V)[1]-k],plotly = TRUE)
graph$plot_function(ek/sqrt(sum(ek^2)),p=p,plotly = TRUE,line_color = "red")



len <- 2
c = 2
beta=3
x <- seq(from=-2*c*beta^2,to=0,length.out=1000)
plot(x,tanh(sqrt(abs(x)*len)),type="l",ylim=c(0,2))
lines(x,beta*sqrt(abs(x))*(c-1)/(x+c*beta^2),col=2,ylim=c(-2,2))
lines(-c*beta^2*c(1,1),c(-10,10),col=3)
