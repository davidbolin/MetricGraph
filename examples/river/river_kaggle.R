#install.packages('SSNbayes', dependencies = T)

#Or from Github
#remotes::install_github("EdgarSantos-Fernandez/SSNbayes", ref = "HEAD",
#                         upgrade_dependencies = F,
#                         dependencies = F)

#remotes::install_github("EdgarSantos-Fernandez/SSNdata", ref = "HEAD",
#                        upgrade_dependencies = F,
#                        dependencies = F)
library('SSNbayes')
library('SSNdata')
library(RSQLite)
library('SSN2')
library('SSNbayes')
library(ymd)
library(ggplot2)
library(MetricGraph)
## Import SpatialStreamNetwork object
path <- system.file("extdata/clearwater.ssn", package = "SSNdata")
n <- ssn_import(path, predpts = "preds", overwrite = T)

## Import data.frame containing observations and covariates
clear <- readRDS(system.file("extdata//clear_obs.RDS", package = "SSNdata"))

seed <- '202103'
set.seed(seed)

# let's hold 30% of the data for testing the models
ids <- clear[!is.na(clear$temp),]$pid
locs <- sample(ids, round(length(ids) * 0.30), replace = F )
locs <- sort(locs)

clear$y <- clear$temp_backup <- clear$temp # response variable
clear[clear$pid %in% locs,]$y <- NA
clear$dataset <- 'train' # creates a train dataset
clear$dataset[locs] <- 'test' # creates a testing dataset


clear.df <- SSNbayes::collapse(n)
clear.df$addfunccol_cat <- cut(clear.df$afvArea,
                               breaks = seq(min(clear.df$afvArea),
                                            max(clear.df$afvArea),
                                            length.out=5),
                               labels = 1:4,
                               include.lowest = T)
col <- 'gray'

fig <- ggplot(clear.df) +
  geom_path(aes(X1, X2, group = slot, size = addfunccol_cat), lineend = 'round', linejoin = 'round', col = col)+
  geom_point(data = dplyr::filter(clear, clear$date == ymd('2012-08-01')) ,
             aes(x = UTM_Xcoord, y = UTM_Ycoord, col = temp_backup), size = 1.75)+
  geom_text(data = dplyr::filter(clear, clear$date == ymd('2012-08-01')),
            aes(x = UTM_Xcoord, y = UTM_Ycoord+1000, label = locID),size = 2)+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_viridis(option = 'C')+
  scale_shape_manual(values = c(16))+
  ylab("Latitude") +
  xlab("Longitude")+
  coord_fixed()+
  theme_bw()+
  guides(size = 'none')+
  labs(size="",colour = "Temperature(°C)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=13),
        legend.title=element_text(size=13),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background =element_rect(fill='white'))

print(fig)


metric.obj <- metric_graph$new(n$edges$geometry)
clear.dat <- dplyr::filter(clear, clear$date == ymd('2012-08-01'))
clear.dat <- clear.dat[is.na(clear.dat$temp_backup)==F,]
Y.pos <- metric.obj$coordinates(XY = cbind(clear.dat$NEAR_X,clear.dat$NEAR_Y))
data <- data.frame(edge_number      = Y.pos[,1],
                   distance_on_edge = Y.pos[,2],
                   y                = clear.dat$temp_backup,
                   elev             = clear.dat$elev,
                   slope            = clear.dat$SLOPE)

metric.obj$add_observations(data = data, normalized = TRUE)
metric.obj$observation_to_vertex()
metric.obj$plot(data = "y", vertex_size = 0.5)

res <- graph_lme(y ~ elev+slope , graph = metric.obj, model = 'WM1')
res.cross.1 <-posterior_crossvalidation(res)
#res2 <- graph_lme(y ~ elev+slope + netid, graph = graph, model = 'WM2',improve_hessian=T)
#res.cross.2 <-posterior_crossvalidation(res2)
res.dir <- graph_lme(y ~ elev+slope , graph = metric.obj, model = 'wmd1', optim_method='Nelder-Mead')
res.cross.dir <-posterior_crossvalidation(res.dir)
time <- length(unique(clear$date))
rho <- 0.75
Qtime <- Matrix::sparseMatrix(i = c(1:time, 1:(time - 1)),
                              j = c(1:time, 2:time),
                              x = c(rep(1+rho^2,time),rep(-rho,time-1))/(1-rho^2), symmetric = TRUE,
                     repr = "T")
Qtime[1,1] <- 1/(1-rho^2)
Qtime[time,time] <- 1/(1-rho^2)
solve(Qtime)
Q <- Matrix::kronecker(Qtime,Qtime)
metric.obj$buildDirectionalConstraints(alpha = 1)
reciprocal_tau <- 1
kappa <- 1
Q.list <- MetricGraph:::Qalpha1_edges(c( 1/reciprocal_tau,kappa),
                        metric.obj,
                        w = 0,
                        BC=1, build=FALSE)
n_const <- length(metric.obj$CoB$S)
ind.const <- c(1:n_const)
Tc <- metric.obj$CoB$T[-ind.const, ]
Q <- Matrix::sparseMatrix(i = Q.list$i,
                          j = Q.list$j,
                          x = Q.list$x,
                          dims = Q.list$dims)
Q_space <- forceSymmetric(Tc%*%Q%*%t(Tc))
Q_0 <- kronecker(Qtime,Q_space)
R <- Matrix::Cholesky(Q_0,
                      LDL = FALSE, perm = TRUE)
det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]

