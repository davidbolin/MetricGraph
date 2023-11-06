library(rSPDE)
library(MetricGraph)
library(INLA)
graph <- metric_graph$new(tolerance = list(vertex_vertex = 1e-1, vertex_line = 1e-3, line_line = 1e-3))
graph$prune_vertices()
graph$plot()

graph$build_mesh(h = 0.1)
graph$plot(mesh=TRUE)
graph$compute_fem()

sigma <- 0.5
range <- 1.5
alpha <- 2
n.rep <- 10
lgcp_sample <- graph_lgcp(n = n.rep, intercept = 1, sigma = sigma,
                          range = range, alpha = alpha,
                          graph = graph)

graph$plot_function(exp(lgcp_sample[[1]]$u), vertex_size = 0)


df <- data.frame(edge_number = unlist(sapply(lgcp_sample, function(x) x$edge_numbers)),
                 distance_on_edge = unlist(sapply(lgcp_sample, function(x) x$edge_loc)))
df$y <- rep(1,length(df$edge_number))
n.samp <- unlist(sapply(lgcp_sample, function(x) length(x$edge_numbers)))

Atilde <- kronecker(diag(n.rep),graph$fem_basis(graph$mesh$VtE))
atilde <- rep(graph$mesh$weights, n.rep)
ind.tilde <- kronecker(1:n.rep,rep(1,length(graph$mesh$weights)))
ind <- NULL
for(i in 1:n.rep){
  ind <- c(ind, rep(i,n.samp[i]))
}

df$rep <- ind
graph$add_observations(data = df,normalized = TRUE, group = "rep")
graph$plot(vertex_size = 0, data = "y")


graph$add_observations(data = data.frame(y = rep(0,length(atilde)),
                                         e = atilde,
                                         edge_number = rep(graph$mesh$VtE[,1],n.rep),
                                         distance_on_edge = rep(graph$mesh$VtE[,2],n.rep),
                                         rep = ind.tilde),
                       group = "rep",
                       normalized = TRUE)

rspde_model <- rspde.metric_graph(graph, nu = alpha - 1/2)
spde_index <- rspde.make.index(name="field", mesh=graph, nu = alpha - 1/2, n.repl = n.rep)
A <- rspde.make.A(mesh = graph, nu = alpha - 1/2)


stk <- inla.stack(data = graph_data_spde(rspde_model),
                  A = A,
                  effects = c(spde_index, list(Intercept = 1)))

spde_fit <- inla(y ~ -1 + Intercept + f(field, model = rspde_model),
                 family = "poisson", data = inla.stack.data(stk),
                 control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
                 E = inla.stack.data(stk)$e)
spde_result <- rspde.result(spde_fit, "field", rspde_model)

summary(spde_result)
result_df <- data.frame(
  parameter = c("std.dev", "range"),
  true = c(sigma, range),
  mean = c(
    spde_result$summary.std.dev$mean,
    spde_result$summary.range$mean
  ),
  mode = c(
    spde_result$summary.std.dev$mode,
    spde_result$summary.range$mode
  )
)
print(result_df)
posterior_df_fit <- gg_df(spde_result)

library(ggplot2)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) +
  facet_wrap(~parameter, scales = "free") + labs(y = "Density")

n.obs <- length(graph$data$y)
n.field <- dim(graph$mesh$VtE)[1]
u_posterior <- spde_fit$summary.linear.predictor$mean[(n.obs+1):(n.obs+n.field)]
graph$plot_function(u_posterior, vertex_size = 0)

graph$plot_function(lgcp_sample$u, vertex_size = 0)


## References
