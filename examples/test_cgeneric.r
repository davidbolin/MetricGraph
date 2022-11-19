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

A <- graph$A

sigma <- 1

sigma.e <- 0.1

nu <- 0.5

r <- 0.2

kappa <- sqrt(8 * nu) / r

theta <- c(sigma, kappa)

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

spde_model <- gpgraph_spde(graph)

spde_model_check <- gpgraph_spde(graph, start_kappa = kappa,
                                    start_sigma = sigma,
                                    parameterization = "spde")

Q_chk <- inla.cgeneric.q(spde_model_check)$Q

sum((Q_chk@i - Q@i)^2)
sum((Q_chk@p - Q@p)^2)
sum((Q_chk@x-Q@x)^2)

data_list <- list(y = as.vector(y), obs.loc = obs.loc)

library(inlabru)

cmp <-
    y ~ -1 + Intercept(1) + field(obs.loc, model = spde_model)

spde_bru_fit <-
    bru(cmp,
        data=data_list,
      options=list(
      family = "gaussian",
      inla.mode = "experimental"))



spde.index <- graph_spde_make_index(name="field", graph=graph)

stk.dat <- inla.stack(data = list(y=as.vector(y)), 
                        A = A, 
                        effects = c(
      spde.index,
      list(Intercept = 1)
    ))

f.s <- y ~ -1 + Intercept + f(field, model = spde_model)

spde_fit <- inla(f.s, data = inla.stack.data(stk.dat), verbose = TRUE)

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)

summary(spde_result)
