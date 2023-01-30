library(MetricGraph)
library(INLA)
library(inlabru)
library(ggplot2)
data(pems)
graph <-  metric_graph$new(lines = pems$lines, longlat = TRUE)
graph$plot()

# Creating the data frame
df_graph <- data.frame(y = pems$Y,
                       edge_number = pems$PtE[,1],
                       distance_on_edge = pems$PtE[,2])

graph$add_observations(data = df_graph, normalized=TRUE)

graph$plot(data="y", vertex_size = 0, data_size = 2)

#create the spde model
spde_model <- graph_spde(graph, parameterization = "spde")

# Fit the model using inlabru
cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_model)

spde_fit <- bru(cmp, data = graph_data_spde(spde_model, loc = "loc"))

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)
summary(spde_result)

posterior_df_fit <- gg_df(spde_result)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) +
  facet_wrap(~parameter, scales = "free") + labs(y = "Density")


#kriging
graph$build_mesh(h = 0.1)

graph$plot(mesh=TRUE)

field_pred <- predict(spde_model,
                      cmp,
                      spde_fit,
                      data = graph$get_mesh_locations(bru = TRUE, loc = "loc"),
                      formula = ~field)

plot(field_pred)


p1 <- ggmap(map2)
p1 <- p1 + xlab(NULL) + ylab(NULL)
p1 <- graph$plot_function(fitted_values$mean, p = p1, edge_width = 0.75)
p1 <- graph$plot(data = TRUE, p = p1, vertex_size = 0, edge_width = 0,
                 data_size = 2)
p1 <- p1 + coord_cartesian(xlim =c(-121.905,-121.88), ylim = c(37.316, 37.328))
p1 <- p1 + xlab(NULL) + ylab(NULL)
p1 <- p1 + scale_colour_continuous(type = "viridis") +
  theme(legend.position = "none")
p1

p2 <- ggmap(map3)
p2 <- p2 + xlab(NULL) + ylab(NULL)
p2 <- graph$plot_function(fitted_values$mean, p = p2, edge_width = 0.5)
p2 <- graph$plot(data = TRUE, p = p2, vertex_size = 0, edge_width = 0,
                 data_size = 2)
p2 <- p2 + coord_cartesian(xlim =c(-121.94,-121.88), ylim = c(37.35, 37.375))
p2 <- p2 + xlab(NULL) + ylab(NULL)
p2 <- p2 + scale_colour_continuous(type = "viridis")#, limits = c(44,53))
p2


#fit fractional model
rm(list=ls())
library(MetricGraph)
library(INLA)
library(inlabru)
library(rSPDE)
data(pems)
graph <-  metric_graph$new(lines = pems$lines, longlat = TRUE)
df_graph <- data.frame(y = pems$Y,
                       edge_number = pems$PtE[,1],
                       distance_on_edge = pems$PtE[,2])

graph$add_observations(data = df_graph, normalized=TRUE)

graph$build_mesh(h = 0.1)
graph$compute_fem()
rspde_model <- rspde.metric_graph(graph)
cmp <- y ~ -1 + Intercept(1) + field(loc, model = rspde_model)
rspde_fit <- bru(cmp, data = graph_data_spde(rspde_model, loc = "loc"))


result_fit <- rspde.result(rspde_fit, "field", rspde_model)
summary(result_fit)


library(ggmap)
map2=get_map(c(left = -121.91, bottom = 37.3,
               right = -121.85, top = 37.35),
             maptype = "toner-lite")
map3=get_map(c(left = -121.96, bottom = 37.34,
               right = -121.85, top = 37.38),
             maptype = "toner-lite")
map1 <- get_map(c(left = -122.2, bottom = 37.05,
                  right = -121.7, top = 37.6),
                maptype = "toner-lite")

