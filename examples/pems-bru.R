library(MetricGraph)
library(INLA)
library(inlabru)
data(pems)
graph <-  metric_graph$new(lines = pems$lines, longlat = TRUE)
graph$plot(vertex_size=1)

# Creating the data frame
df_graph <- data.frame(y = pems$Y,
                       edge_number = pems$PtE[,1],
                       distance_on_edge = pems$PtE[,2])

graph$add_observations(data = df_graph, normalized=TRUE)

graph$plot(data="y", vertex_size = 0, data_size = 2)

#create the spde model
spde_model <- graph_spde(graph)

# Fit the model using inlabru
cmp <- y ~ -1 + Intercept(1) + field(NULL, model = spde_model)

spde_fit <- bru(cmp, data=graph_data_spde(spde_model))

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)
summary(spde_result)

#kriging
graph$build_mesh(h = 0.1)

graph$plot(mesh=TRUE)


fitted_values <- inlabru_predict(spde_model,
                                 spde_fit,
                                 cmp,
                                 PtE = graph$mesh$PtE)
m_prd_bru <- fitted_values$mean
graph_bru$plot_function(X = m_prd_bru, vertex_size = 0)
