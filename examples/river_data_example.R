library(SSN)
library(MetricGraph)
#for examples, copy MiddleFork04.ssn directory to R's temporary directory
copyLSN2temp()
mf04p <- importSSN(paste0(tempdir(),'/MiddleFork04.ssn'), predpts = "pred1km",
                   o.write = TRUE)
mf04p <- importPredpts(target = mf04p, predpts = "CapeHorn", obj.type = "ssn")
mf04p <- importPredpts(target = mf04p, predpts = "Knapp", obj.type = "ssn")
names(mf04p)
mf04p <- additive.function(mf04p, "h2oAreaKm2", "computed.afv")
createDistMat(mf04p, predpts = "Knapp", o.write = TRUE,
              amongpreds = TRUE)
createDistMat(mf04p, predpts = "CapeHorn", o.write = TRUE,
              amongpreds = TRUE)

distObs <- getStreamDistMat(mf04p)
plot(mf04p, lwdLineCol = "afvArea", lwdLineEx = 10, lineCol = "blue",
     pch = 19, xlab = "x-coordinate (m)", ylab = "y-coordinate (m)",
     asp = 1)

brks <- plot(mf04p, "Summer_mn", lwdLineCol = "afvArea",
             lwdLineEx = 15, lineCol = "black", xlab = "x-coordinate" ,
             ylab = "y-coordinate", asp=1 )

mf04.glmssn0 <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE, mf04p,
                       CorModels = NULL, use.nugget = TRUE)
summary(mf04.glmssn0)

mf04.glmssn1 <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE + netID, mf04p,
                       CorModels = c("Exponential.tailup", "Exponential.taildown",
                                     "Exponential.Euclid"), addfunccol = "afvArea")
summary(mf04.glmssn1)
mf04c.resid2 <- residuals(mf04.glmssn1,
                          cross.validation = TRUE)
abs.res <- mean(abs(getSSNdata.frame(mf04c.resid2)[, "_resid.crossv_"]))
graph <- metric_graph$new(mf04p, longlat = FALSE)

Y.pos <- graph$coordinates(XY = mf04p@obspoints@SSNPoints[[1]]@point.coords)
data <- data.frame(edge_number      = Y.pos[,1],
                   distance_on_edge = Y.pos[,2],
                   y                = mf04p@obspoints@SSNPoints[[1]]@point.data$Summer_mn,
                   elev             = mf04p@obspoints@SSNPoints[[1]]@point.data$ELEV_DEM,
                   slope            = mf04p@obspoints@SSNPoints[[1]]@point.data$SLOPE,
                   netid            = mf04p@obspoints@SSNPoints[[1]]@point.data$netID)

graph$add_observations(data = data, normalized = TRUE)
graph$observation_to_vertex()
graph$plot(data = "y", vertex_size = 0.5)
res_lm <- graph_lme(y ~ 1, graph = graph)
res <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'WM1')
res.cross.1 <-posterior_crossvalidation(res)
res2 <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'WM2',improve_hessian=T)
res.cross.2 <-posterior_crossvalidation(res2)
res.dir <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'wmd1', optim_method='Nelder-Mead')
res.cross.dir <-posterior_crossvalidation(res.dir)
