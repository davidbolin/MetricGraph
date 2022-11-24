library(GPGraph)
library(ggplot2)
library(maptools)
library(Matrix)
library(igraph)
library(ggplot2)
library(osmdata)
set_overpass_url("https://maps.mail.ru/osm/tools/overpass/api/interpreter")

call <- opq(bbox = c(12.5,55.7,12.6,55.8))
call <- add_osm_feature(call, key = "highway",value=c("motorway",
                                                      "primary",
                                                      "secondary",
                                                      "tertiary"))

data <- osmdata_sp(call)

Lines <- largest_component(SpatialLines(data$osm_lines@lines))

graph <- metric_graph$new(Lines = Lines)

#split circular edges (only needed for alpha = 2)
if(0){
  circ.ind <- which(graph$E[,1]==graph$E[,2])
  PtE <- cbind(circ.ind, rep(0.5,length(circ.ind)))
  graph$add_PtE_observations(y = rep(0, length(circ.ind)), PtE, normalized = TRUE)
  graph$observation_to_vertex()
  graph$clear_observations()
}

p <- graph$plot(marker_size = 0)

p <- p + theme_void() + theme(legend.position="none")
p

graph$build_mesh(h = 0.001)

PtE_tmp <- rbind(graph$VtEfirst(),graph$mesh$PtE)

type = "alpha1"
if(type == "isoExp") { #isotropic exponential
  graph$compute_resdist_mesh()
  Sigma <- exp(-50*graph$mesh$res_dist)
  u <- as.vector(t(chol(forceSymmetric(Sigma)))%*%rnorm(dim(Sigma)[1]))
} else if (type == "alpha1") { #alpha =1 model looks strange
  u <- sample_spde(kappa = 50, sigma = 1, alpha = 1,
                   graph = graph, type="mesh")
} else if (type == "alpha2") { #not working
  u <- sample_spde(kappa = 50, sigma = 1, alpha = 2,
                   graph = graph,
                   type = "mesh")
}

graph$plot_function_mesh(X = u)

graph$plot(X=u, X_loc = PtE_tmp)

n_obs <- length(as.vector(u))
sigma.e <- 0.1
y <- u + sigma.e * rnorm(n_obs)

y_mesh <- y[(graph$nV+1):(graph$nV+nrow(graph$mesh$PtE))]

graph$add_mesh_observations(y_mesh)

graph$observation_to_vertex()

spde_model_bru <- gpgraph_spde(graph, parameterization = "spde")

data_list <- list(y = as.vector(y_mesh), loc = graph$PtE)

cmp <-
    y ~ -1 + Intercept(1) + field(loc,
                            model = spde_model_bru)
library(inlabru)

spde_bru_fit <-
    bru(cmp, data=data_list)

spde_bru_result <- spde_metric_graph_result(spde_bru_fit,
                    "field", spde_model_bru)

summary(spde_bru_result)

# m_prd <- spde_bru_fit$summary.fitted.values$mean

# graph$plot_function_mesh(m_prd[1:n_obs_mesh])

