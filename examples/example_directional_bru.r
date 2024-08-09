library(MetricGraph)
library(SSN2)

copy_lsn_to_temp()
path <- file.path(tempdir(), "MiddleFork04.ssn")

mf04p <- ssn_import(
  path = path,
  predpts = c("pred1km", "CapeHorn"),
  overwrite = TRUE
)

graph <- metric_graph$new(mf04p$edges$geometry, longlat = FALSE)

obs_coords <- sf::st_coordinates(mf04p$obs)

data <- data.frame(coordx      = obs_coords[,1],
                   coordy      = obs_coords[,2],
                   y                = mf04p$obs$Summer_mn,
                   elev             = mf04p$obs$ELEV_DEM,
                   slope            = mf04p$obs$SLOPE,
                   netid            = as.factor(mf04p$obs$netID))
graph$add_observations(data = data,
                       data_coords = "spatial",
                       coord_x = "coordx",
                       coord_y = "coordy")


res.wm1.dir <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'wmd1', optim_method='BFGS')

spde_model_bru <- graph_spde(graph, alpha=1, directional=TRUE)

cmp <-
    y ~ -1 + Intercept(1) + elev(elev) + slope(slope) + netid(netid) + field(loc,
                    model = spde_model_bru)

data_spde_dir <- graph_data_spde(graph_spde = spde_model_bru, 
                            loc_name = "loc")

library(inlabru)

spde_bru_fit_dir <-
    bru(cmp, data=data_spde_dir[["data"]])

for(i in 1:10){
  spde_bru_fit_dir <- bru_rerun(spde_bru_fit_dir)
}

spde_result <- spde_metric_graph_result(spde_bru_fit_dir, "field", spde_model_bru)

summary(spde_result)
summary(spde_bru_fit_dir)

summary(res.wm1.dir)













spde_model <- graph_spde(graph, alpha=1,directional=TRUE)

data_spde <- graph_data_spde(graph_spde = spde_model, name = "field")

f.s <- y ~ -1 + Intercept + elev + slope + netid + f(field, model = spde_model)

stk_dat <- inla.stack(data = data_spde[["data"]], 
                        A = data_spde[["basis"]], 
                        effects = c(
      data_spde[["index"]],
      list(Intercept = 1)
    ))

data_stk <- inla.stack.data(stk_dat)

spde_fit <- inla(f.s, data = data_stk, control.predictor=list(A=inla.stack.A(stk_dat)))

spde_result <- spde_metric_graph_result(spde_fit, "field", spde_model)

summary(spde_result)
