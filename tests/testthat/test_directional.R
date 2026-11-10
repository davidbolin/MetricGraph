

test_that("simple posterior mean check", {
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(1,0),c(1,1))
edges = list(edge1,edge2)
graph <- metric_graph$new(edges = edges)
theta <- c(1,1,1)
Eu <- MetricGraph:::posterior_mean_alpha1_directional(theta, graph, c(1),
                                                 matrix(c(1,0.5),nrow=1,ncol=2))
Eu_t <- c(rep(r_1(0.5,1,1)/(r_1(0.,1,1)+1),3),r_1(1.5,1,1)/(r_1(0.,1,1)+1))
expect_equal(as.vector(Eu), Eu_t, tolerance = 1e-10)
#no noise case
theta[1] <- 0
Eu <- MetricGraph:::posterior_mean_alpha1_directional(theta, graph, c(1),
                                                      matrix(c(1,0.5),nrow=1,ncol=2))
Eu_t <- c(rep(r_1(0.5,1,1)/r_1(0.,1,1),3),r_1(1.5,1,1)/r_1(0.,1,1))
expect_equal(as.vector(Eu), Eu_t, tolerance = 1e-10)
})
test_that("posterior mean discontinous check", {

  edge1 <- rbind(c(1,0),c(0,0))
  edge2 <- rbind(c(1,1),c(1,0))
  edge3 <- rbind(c(1,-1),c(1,0))
  edges = list(edge1,edge2,edge3)
  graph <- metric_graph$new(edges = edges)
  #graph$plot(direction = T)
  theta <- c(1,1,1)
  theta[1] <- 0
  graph$setDirectionalWeightFunction( f_in = function(x){(x/sum(x))})
  Eu <- MetricGraph:::posterior_mean_alpha1_directional(theta, graph, c(1),
                                                        matrix(c(3,0.5),nrow=1,ncol=2))
  Eu_t <- c(0.5*r_1(0.5,1,1)/r_1(0.,1,1),
            0.5*r_1(1.5,1,1)/r_1(0.,1,1),
            0,
            0,
            r_1(0.5,1,1)/r_1(0.,1,1),
            r_1(0.5,1,1)/r_1(0.,1,1))
  expect_equal(as.vector(Eu), Eu_t, tolerance = 1e-10)
  expect_warning(graph$setDirectionalWeightFunction( f_in = function(x){sqrt(x/sum(x))}))
  Eu <- MetricGraph:::posterior_mean_alpha1_directional(theta, graph, c(1),
                                                        matrix(c(3,0.5),nrow=1,ncol=2))
  Eu_t <- c(sqrt(0.5)*r_1(0.5,1,1)/r_1(0.,1,1),
            sqrt(0.5)*r_1(1.5,1,1)/r_1(0.,1,1),
            0,
            0,
            r_1(0.5,1,1)/r_1(0.,1,1),
            r_1(0.5,1,1)/r_1(0.,1,1))
  expect_equal(as.vector(Eu), Eu_t, tolerance = 1e-10)

  Eu <- MetricGraph:::posterior_mean_alpha1_directional(theta, graph, c(1),
                                                        matrix(c(1,0.5),nrow=1,ncol=2))
  Eu_t <- c(r_1(0.5,1,1)/r_1(0.,1,1),
            r_1(0.5,1,1)/r_1(0.,1,1),
            sqrt(0.5)*r_1(1.5,1,1)/r_1(0.,1,1),
            sqrt(0.5)*r_1(.5,1,1)/r_1(0.,1,1),
            sqrt(0.5)*r_1(1.5,1,1)/r_1(0.,1,1),
            sqrt(0.5)*r_1(.5,1,1)/r_1(0.,1,1))
  expect_equal(as.vector(Eu), Eu_t, tolerance = 1e-10)
})

