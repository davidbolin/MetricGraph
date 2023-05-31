test_that("Check cgeneric precision matrices", {
set.seed(1)
testthat::skip_on_cran()
  if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")


line1 <- sp::Line(rbind(c(0,0),c(1,0)))
line2 <- sp::Line(rbind(c(0,0),c(0,1)))
line3 <- sp::Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
line4 <- sp::Line(cbind(sin(theta),1+ cos(theta)))
Lines = sp::SpatialLines(list(sp::Lines(list(line1),ID="1"),
                              sp::Lines(list(line2),ID="2"),
                              sp::Lines(list(line4),ID="3"),
                              sp::Lines(list(line3),ID="4")))
graph <- metric_graph$new(lines = Lines)

obs.per.edge <- 250
obs.loc <- NULL
for(i in 1:(graph$nE)) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), runif(obs.per.edge) * 
                          graph$edge_lengths[i]))
}

n.obs <- nrow(obs.loc)

y <- rep(NA, obs.per.edge * graph$nE)

df_temp <- data.frame(y = y, edge_number = obs.loc[,1], distance_on_edge = obs.loc[,2])

graph$add_observations(data = df_temp, normalized = FALSE)

graph$observation_to_vertex()

A <- graph$A()

tau <- 0.05

sigma.e <- 0.1

nu <- 0.5

r <- 0.3

kappa <- sqrt(8 * nu) / r

theta <- c(sigma, kappa)

Q <- spde_precision(kappa = kappa, tau = tau, alpha = 1, graph = graph)

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

graph$data$y <- as.vector(y)

# graph$plot(line_width = 0.3, data=TRUE)

spde_model <- graph_spde(graph)

spde_model_check <- graph_spde(graph, start_kappa = kappa,
                                    start_sigma = 1/tau,
                                    parameterization = "spde")

Q_chk <- INLA::inla.cgeneric.q(spde_model_check)$Q

expect_equal(sum((Q_chk@i - Q@i)^2), 0)
expect_equal(sum((Q_chk@p - Q@p)^2), 0)
expect_equal(sum((Q_chk@x-Q@x)^2), 0)

  INLA::inla.setOption(num.threads = old_threads)
})

