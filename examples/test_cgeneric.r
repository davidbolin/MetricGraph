library(sp)
library(INLA)
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

obs.per.edge <- 50
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

nsim <- 30

Z <- rnorm(sizeQ * nsim)
dim(Z) <- c(sizeQ, nsim)

sizeA <- nrow(A)

eps <- rnorm(sizeA * nsim)
dim(eps) <- c(sizeA, nsim)

LQ <- chol(Q)
u <- solve(LQ, Z)

y <- A%*%u + sigma.e * eps

graph$y <- as.vector(y)

# graph$plot(line_width = 0.3, data=TRUE)

spde_model <- gpgraph_spde(graph)

spde_model_check <- gpgraph_spde(graph, start_kappa = kappa,
                                    start_sigma = sigma,
                                    parameterization = "spde")

Q_chk <- inla.cgeneric.q(spde_model_check)$Q

sum((Q_chk@i - Q@i)^2)
sum((Q_chk@p - Q@p)^2)
sum((Q_chk@x-Q@x)^2)

obs.loc.rep <- obs.loc
for(i in 2:nsim){
    obs.loc.rep <- rbind(obs.loc.rep, obs.loc)
}
data_list <- list(y = as.vector(y),
                            loc = obs.loc.rep)

# data_list <- list(y = as.vector(y), obs.loc = obs.loc)

library(inlabru)

repl <- rep(1:nsim, each=200)

cmp <-
    y ~ -1 + Intercept(1) + field(loc, model = spde_model,
    replicate = repl)

spde_bru_fit <-
    bru(cmp,
        data=data_list,
      options=list(
      family = "gaussian",
      inla.mode = "experimental"))



spde_bru_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model)

summary(spde_bru_result)



spde.index <- graph_spde_make_index(name="field", graph=graph, n.repl = nsim)
A <- graph_spde_make_A(graph, n.repl = nsim)

stk.dat <- inla.stack(data = list(y=as.vector(y)), 
                        A = A, 
                        effects = c(
      spde.index,
      list(Intercept = 1)
    ))

f.s <- y ~ -1 + Intercept + f(field, model = spde_model, replicate = field.repl)

data_stk <- graph_stack(stk.dat, "field")

spde_fit <- inla(f.s, data = data_stk)

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)

summary(spde_result)


##############

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

obs.per.edge <- 50
obs.loc <- NULL
for(i in 1:(graph$nE)) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), runif(obs.per.edge) * 
                          graph$edge_lengths[i]))
}

y <- rep(NA, obs.per.edge * graph$nE)

n.obs <- obs.per.edge * graph$nE

# We will also add equally spaced nodes where we will do predictions:A

obs.loc2 <- NULL
for(i in 1:(graph$nE)) {
  obs.loc2 <- rbind(obs.loc2,
                   cbind(rep(i,obs.per.edge), 1:obs.per.edge/(obs.per.edge+1) * 
                          graph$edge_lengths[i]))
}

obs.loc <- rbind(obs.loc, obs.loc2)

y <- c(y, rep(NA, obs.per.edge * graph$nE))

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

y <- as.vector(y)

y <- y[1:n.obs]

y <- c(y, rep(NA, n.obs))

# idx <- as.vector(graph$A %*% 1:(2*n.obs+4))

# graph$y[idx-4] <- y

graph$add_responses(y)

graph$plot(line_width = 0.3, data=TRUE)

spde_model <- gpgraph_spde(graph)

spde_model_check <- gpgraph_spde(graph, start_kappa = kappa,
                                    start_sigma = sigma,
                                    parameterization = "spde")

Q_chk <- inla.cgeneric.q(spde_model_check)$Q

sum((Q_chk@i - Q@i)^2)
sum((Q_chk@p - Q@p)^2)
sum((Q_chk@x-Q@x)^2)

obs.loc.rep <- obs.loc
if(nsim>1){
    for(i in 2:nsim){
        obs.loc.rep <- rbind(obs.loc.rep, obs.loc)
    }
}
data_list <- list(y = as.vector(y),
                            loc = obs.loc)

# data_list <- list(y = as.vector(y), obs.loc = obs.loc)

library(inlabru)

cmp <-
    y ~ -1 + Intercept(1) + field(loc, model = spde_model)

spde_bru_fit <-
    bru(cmp,
        data=data_list,
      options=list(
      family = "gaussian",
      inla.mode = "experimental"))



spde_bru_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model)

summary(spde_bru_result)



spde.index <- graph_spde_make_index(name="field", graph=graph, n.repl = nsim)
A <- graph_spde_make_A(graph, n.repl = nsim)

stk.dat <- inla.stack(data = list(y=as.vector(y)), 
                        A = A, 
                        effects = c(
      spde.index,
      list(Intercept = 1)
    ))

f.s <- y ~ -1 + Intercept + f(field, model = spde_model)

data_stk <- graph_stack(stk.dat, "field")

spde_fit <- inla(f.s, data = data_stk)

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)

summary(spde_result)

m.prd <- spde_fit$summary.fitted.values$mean[(n.obs+1):(2*n.obs)]

prd.y <- c(y[1:n.obs], m.prd)

graph$add_responses(prd.y)

graph$plot(data=TRUE)

m.prd.matrix <- cbind(obs.loc2, m.prd)

graph$plot(X = as.vector(m.prd), X_loc = as.matrix(obs.loc2), show=FALSE)

graph$plot_function(X = m.prd.matrix, marker_size = 3, plotly=TRUE)



##############

library(sp)
library(INLA)
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

obs.per.edge <- 50
obs.loc <- NULL
for(i in 1:(graph$nE)) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), obs.per.edge:1/(obs.per.edge+1) * 
                          graph$edge_lengths[i]))

#   obs.loc <- rbind(obs.loc,
#                    cbind(rep(i,obs.per.edge), 1:obs.per.edge/(obs.per.edge+1) * 
#                           graph$edge_lengths[i]))
}

n.obs <- obs.per.edge * graph$nE

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

nsim <- 30

Z <- rnorm(sizeQ * nsim)
dim(Z) <- c(sizeQ, nsim)

sizeA <- nrow(A)

eps <- rnorm(sizeA * nsim)
dim(eps) <- c(sizeA, nsim)

LQ <- chol(Q)
u <- solve(LQ, Z)

y <- A%*%u + sigma.e * eps

graph$y <- as.vector(y)

# graph$plot(line_width = 0.3, data=TRUE)

spde_model <- gpgraph_spde(graph)

spde_model_check <- gpgraph_spde(graph, start_kappa = kappa,
                                    start_sigma = sigma,
                                    parameterization = "spde")

Q_chk <- inla.cgeneric.q(spde_model_check)$Q

sum((Q_chk@i - Q@i)^2)
sum((Q_chk@p - Q@p)^2)
sum((Q_chk@x-Q@x)^2)

obs.loc.rep <- obs.loc
if(nsim>1){
    for(i in 2:nsim){
        obs.loc.rep <- rbind(obs.loc.rep, obs.loc)
    }   
}
data_list <- list(y = as.vector(y),
                            loc = obs.loc.rep)

# data_list <- list(y = as.vector(y), obs.loc = obs.loc)

library(inlabru)

repl <- rep(1:nsim, each=n.obs)

cmp <-
    y ~ -1 + Intercept(1) + field(loc, model = spde_model,
    replicate = repl)

spde_bru_fit <-
    bru(cmp,
        data=data_list,
      options=list(
      family = "gaussian",
      inla.mode = "experimental"))



spde_bru_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model)

summary(spde_bru_result)



spde.index <- graph_spde_make_index(name="field", graph=graph, n.repl = nsim)
A <- graph_spde_make_A(graph, n.repl = nsim)

stk.dat <- inla.stack(data = list(y=as.vector(y)), 
                        A = A, 
                        effects = c(
      spde.index,
      list(Intercept = 1)
    ))

f.s <- y ~ -1 + Intercept + f(field, model = spde_model, replicate = field.repl)

data_stk <- graph_stack(stk.dat, "field", spde.index)

spde_fit <- inla(f.s, data = data_stk)

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)

summary(spde_result)


##### 
#####
#####
#####


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

obs.per.edge <- 2
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

## Get a copy and do prediction on the copy

graph_prd <- graph$clone()

n.obs <- nrow(graph_prd$PtE)

# We will also add equally spaced nodes where we will do predictions:A

obs.loc2 <- NULL

obs.per.edge.prd <-
graph_prd$edge_lengths/min(graph_prd$edge_lengths)

max_subdiv <- 1

obs.per.edge.prd <- sapply(obs.per.edge.prd, 
        function(x){min(max_subdiv,round(x))})

obs.per.edge.prd <- obs.per.edge.prd * (graph_prd$edge_lengths > 1/50)

for(i in 1:(graph_prd$nE)) {
    if(obs.per.edge.prd[i]>0){
        obs.loc2 <- rbind(obs.loc2,
                         cbind(rep(i,obs.per.edge.prd[i]), 1:obs.per.edge.prd[i]/(obs.per.edge.prd[i]+1) * 
                                graph$edge_lengths[i]))
    }
}

y <- rep(NA, nrow(obs.loc2))

graph_prd$add_observations2(y,obs.loc2)

graph_prd$observation_to_vertex()

graph_prd$plot(line_width = 0.3)

A <- graph_prd$A

sigma <- 1

sigma.e <- 0.1

nu <- 0.5

r <- 0.2

kappa <- sqrt(8 * nu) / r

theta <- c(sigma, kappa)

Q <- Qalpha1(theta, graph_prd)

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

y <- as.vector(y)

y <- y[1:n.obs]

graph$add_responses(y)

graph$plot(data=TRUE)

n.full <- nrow(graph_prd$PtE)

y <- c(y, rep(NA, n.full - n.obs))

# idx <- as.vector(graph$A %*% 1:(2*n.obs+4))

# graph$y[idx-4] <- y

graph_prd$add_responses(y)

graph_prd$plot(line_width = 0.3, data=TRUE)

graph$plot(line_width = 0.3, data=TRUE)

spde_model <- gpgraph_spde(graph)

spde_model_check <- gpgraph_spde(graph, start_kappa = kappa,
                                    start_sigma = sigma,
                                    parameterization = "spde")

Q_chk <- inla.cgeneric.q(spde_model_check)$Q

sum((Q_chk@i - Q@i)^2)
sum((Q_chk@p - Q@p)^2)
sum((Q_chk@x-Q@x)^2)

obs.loc.rep <- obs.loc
if(nsim>1){
    for(i in 2:nsim){
        obs.loc.rep <- rbind(obs.loc.rep, obs.loc)
    }
}
data_list <- list(y = as.vector(y),
                            loc = obs.loc)

# data_list <- list(y = as.vector(y), obs.loc = obs.loc)

library(inlabru)

cmp <-
    y ~ -1 + Intercept(1) + field(loc, model = spde_model)

spde_bru_fit <-
    bru(cmp,
        data=data_list,
      options=list(
      family = "gaussian",
      inla.mode = "experimental"))



spde_bru_result <- spde_metric_graph_result(spde_bru_fit, "field", spde_model)

summary(spde_bru_result)



spde.index <- graph_spde_make_index(name="field", graph=graph, n.repl = nsim)
A <- graph_spde_make_A(graph, n.repl = nsim)

stk.dat <- inla.stack(data = list(y=as.vector(y)), 
                        A = A, 
                        effects = c(
      spde.index,
      list(Intercept = 1)
    ))

f.s <- y ~ -1 + Intercept + f(field, model = spde_model)

data_stk <- graph_stack(stk.dat, "field")

spde_fit <- inla(f.s, data = data_stk)

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)

summary(spde_result)

m.prd <- spde_fit$summary.fitted.values$mean[(n.obs+1):(2*n.obs)]

prd.y <- c(y[1:n.obs], m.prd)

graph$add_responses(prd.y)

graph$plot(data=TRUE)

m.prd.matrix <- cbind(obs.loc2, m.prd)

graph$plot(X = as.vector(m.prd), X_loc = as.matrix(obs.loc2), show=FALSE)

graph$plot_function(X = m.prd.matrix, marker_size = 3, plotly=TRUE)


