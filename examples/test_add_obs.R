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
A <- graph$A
sigma <- 1
sigma.e <- 0.1
nu <- 0.5
r <- 0.2
kappa <- sqrt(8 * nu) / r
theta <- c(sigma, kappa)
Qalpha1 <- function(theta, graph, BC = 1, build = TRUE) {

  kappa <- theta[2]
  sigma <- theta[1]
  i_ <- j_ <- x_ <- rep(0, dim(graph$V)[1]*4)
  count <- 0
  for(i in 1:graph$nE){
    l_e <- graph$edge_lengths[i]
    c1 <- exp(-kappa*l_e)
    c2 <- c1^2
    one_m_c2 = 1-c2
    c_1 = 0.5 + c2/one_m_c2
    c_2 = -c1/one_m_c2

    if (graph$E[i, 1] != graph$E[i, 2]) {

      i_[count + 1] <- graph$E[i, 1]
      j_[count + 1] <- graph$E[i, 1]
      x_[count + 1] <- c_1

      i_[count + 2] <- graph$E[i, 2]
      j_[count + 2] <- graph$E[i, 2]
      x_[count + 2] <- c_1


      i_[count + 3] <- graph$E[i, 1]
      j_[count + 3] <- graph$E[i, 2]
      x_[count + 3] <- c_2

      i_[count + 4] <- graph$E[i, 2]
      j_[count + 4] <- graph$E[i, 1]
      x_[count + 4] <- c_2
      count <- count + 4
    }else{
      i_[count + 1] <- graph$E[i, 1]
      j_[count + 1] <- graph$E[i, 1]
      x_[count + 1] <- tanh(0.5 * kappa * l_e)
      count <- count + 1
    }
  }
  if(BC == 1){
    #does this work for circle?
    i.table <- table(i_[1:count])
    index = as.integer(names(which(i.table < 3)))
    i_ <- c(i_[1:count], index)
    j_ <- c(j_[1:count], index)
    x_ <- c(x_[1:count], rep(0.5, length(index)))
    count <- count + length(index)
  }
  if(build){
    Q <- Matrix::sparseMatrix(i = i_[1:count],
                              j = j_[1:count],
                              x = (2 * kappa / sigma^2) * x_[1:count],
                              dims = c(graph$nV, graph$nV))


    return(Q)
  } else {
    return(list(i = i_[1:count],
                j = j_[1:count],
                x = (2 * kappa / sigma^2) * x_[1:count],
                dims = c(graph$nV, graph$nV)))
  }
}
Q <- Qalpha1(theta, graph)
sizeQ <- nrow(Q)
nsim <- 1
Z <- rnorm(sizeQ * nsim)
dim(Z) <- c(sizeQ, nsim)
sizeA <- nrow(A)
eps <- rnorm(sizeA * nsim)
dim(eps) <- c(sizeA, nsim)
LQ <- chol(Q)
u <- solve(LQ, Z)
y <- A%*%u + sigma.e * eps
graph$y <- y
graph$plot(line_width = 0.3, data=TRUE)
