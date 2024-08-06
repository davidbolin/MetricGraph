library(SSN2)
library(MetricGraph)
#for examples, copy MiddleFork04.ssn directory to R's temporary directory
copy_lsn_to_temp()
mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
                    predpts = c("pred1km.shp", "Knapp"),
                    overwrite = TRUE)

names(mf04p)
#mf04p <- additive.function(mf04p, "h2oAreaKm2", "computed.afv")
ssn_create_distmat(mf04p,
                   predpts = "pred1km", overwrite = TRUE,
                   among_predpts = TRUE
)
ssn_create_distmat(mf04p,
                   predpts = "Knapp", overwrite = TRUE,
                   among_predpts = TRUE, only_predpts = TRUE
)
#distObs <- getStreamDistMat(mf04p)

mf04.glmssn0 <- ssn_glm(Summer_mn ~ ELEV_DEM + SLOPE, mf04p,
                       CorModels = NULL, use.nugget = TRUE, family = "Gaussian")
summary(mf04.glmssn0)

mf04.glmssn1 <- ssn_lm(Summer_mn ~ ELEV_DEM + SLOPE + factor(netID), mf04p,
                       tailup_type = "exponential",
                       additive = "afvArea",
                       estmethod = "ml")
summary(mf04.glmssn1)
res.cross.val <- loocv(mf04.glmssn1)

graph <- metric_graph$new(mf04p$edges$geometry, longlat = FALSE)
data <- data.frame(coordx      = mf04p$obs$NEAR_X,
                   coordy = mf04p$obs$NEAR_Y,
                   y                = mf04p$obs$Summer_mn,
                   elev             = mf04p$obs$ELEV_DEM,
                   slope            = mf04p$obs$SLOPE,
                   netid            = as.factor(mf04p$obs$netID))

graph$add_observations(data = data,
                       data_coords = "spatial",
                       coord_x = "coordx",
                       coord_y = "coordy")
graph$observation_to_vertex()
graph$plot(data = "y", vertex_size = 0.5)

res_lm <- graph_lme(y ~ 1, graph = graph)
res <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'WM1')
res.cross.1 <-posterior_crossvalidation(res)
res2 <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'WM2',improve_hessian=T)
res.cross.2 <-posterior_crossvalidation(res2)

res.dir <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'wmd1', optim_method='Nelder-Mead')
res.cross.dir <-posterior_crossvalidation(res.dir)



graph <- metric_graph$new(mf04p$edges$geometry, longlat = FALSE)
graph$set_edge_weights(weights=data.frame(out.weight= -1,
                                          in.weight = rep(1, graph$nE)))

graph$add_observations(data = data,
                       data_coords = "spatial",
                       coord_x = "coordx",
                       coord_y = "coordy")
graph$observation_to_vertex()
V_indegree = graph$get_degrees("indegree")
V_outdegree = graph$get_degrees("outdegree")
index_outdegree <- V_outdegree > 0 & V_indegree >0
index_in0      <- V_indegree == 0
Vs <- which(index_outdegree)
W <- graph$get_edge_weights()

for (v in Vs) {
  in_edges    <- which(graph$E[, 2] %in% v)
  W[in_edges,2] <- sqrt(W[in_edges,2]/sum(W[in_edges,2]))
}
graph$set_edge_weights(W)
res.dir.root <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'wmd1', optim_method='Nelder-Mead')
res.cross.dir.root <-posterior_crossvalidation(res.dir.root)

reciprocal_tau <- exp(res.dir.root$mle_par_orig[2])
kappa = exp(res.dir.root$mle_par_orig[3])
Q <- MetricGraph:::Qalpha1_edges(c( 1/reciprocal_tau,kappa),
                   graph,
                   w = 0,
                   BC=1, build=T)

graph$buildDirectionalConstraints(alpha = 1)
n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tc <- graph$CoB$T[-ind.const, ]
Sigma <- as.matrix(solve(Tc%*%Q%*%t(Tc)))
Sigma.dir.root <- t(as.matrix(Tc))%*%Sigma%*%(as.matrix(Tc))
###
# creating vertex weights
#
###
graph <- metric_graph$new(mf04p$edges$geometry, longlat = FALSE)
graph$set_edge_weights(weights=data.frame(out.weight= -1,
                                          in.weight = mf04p$edges$h2oAreaKm2))


graph$add_observations(data = data,
                       data_coords = "spatial",
                       coord_x = "coordx",
                       coord_y = "coordy")
graph$observation_to_vertex()
V_indegree = graph$get_degrees("indegree")
V_outdegree = graph$get_degrees("outdegree")
index_outdegree <- V_outdegree > 0 & V_indegree >0
index_in0      <- V_indegree == 0
Vs <- which(index_outdegree)
W <- graph$get_edge_weights()
W_old <- W
for (v in Vs) {
  in_edges    <- which(graph$E[, 2] %in% v)
  W[in_edges,2] <- sqrt(W_old[in_edges,2]/sum(W_old[in_edges,2]))
}


graph$set_edge_weights(W)
res.dir_weightroot <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'wmd1', optim_method='Nelder-Mead')
res.cross.dir_weightroot <-posterior_crossvalidation(res.dir_weightroot)

reciprocal_tau <- exp(res.dir_weightroot$mle_par_orig[2])
kappa = exp(res.dir_weightroot$mle_par_orig[3])
Q <- MetricGraph:::Qalpha1_edges(c( 1/reciprocal_tau,kappa),
                   graph,
                   w = 0,
                   BC=1, build=T)

graph$buildDirectionalConstraints(alpha = 1)
n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tc <- graph$CoB$T[-ind.const, ]
Sigma <- as.matrix(solve(Tc%*%Q%*%t(Tc)))
Sigma.dir.weightroot <- t(as.matrix(Tc))%*%Sigma%*%(as.matrix(Tc))

graph <- metric_graph$new(mf04p$edges$geometry, longlat = FALSE)
graph$set_edge_weights(weights=data.frame(out.weight= -1,
                                          in.weight = mf04p$edges$h2oAreaKm2))


graph$add_observations(data = data,
                       data_coords = "spatial",
                       coord_x = "coordx",
                       coord_y = "coordy")
graph$observation_to_vertex()
for (v in Vs) {
  in_edges    <- which(graph$E[, 2] %in% v)
  W[in_edges,2] <- (W_old[in_edges,2]/sum(W_old[in_edges,2]))
}
graph$set_edge_weights(W)
#graph$buildDirectionalConstraints(alpha = 1)
res.dir_weight <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'wmd1', optim_method='Nelder-Mead')
res.cross.dir_weight <- posterior_crossvalidation(res.dir_weight)


cross.val.res <- rbind(res.cross.1$scores,
                       res.cross.2$scores,
                       res.cross.dir$scores,
                       res.cross.dir.root$scores,
                       res.cross.dir_weight$scores,
                       res.cross.dir_weightroot$scores)
