
library('tidyverse')
library('abind')
library('SSN2')
library('bayesplot')
library('ymd')
library('viridis')
library('ggrepel')
library('devtools')
library('RColorBrewer')

# install_github("jayverhoef/SSN")

RNGkind(sample.kind = "Rounding")

library('SSNbayes')
library('SSNdata')

path <- system.file("extdata/clearwater.ssn", package = "SSNdata")
n <- ssn_import(path, predpts = "preds", overwrite = T)

## Import data.frame containing observations and covariates
clear <- readRDS(system.file("extdata//clear_obs.RDS", package = "SSNdata"))

seed <- '202103'
set.seed(seed)

ids <- clear[!is.na(clear$temp),]$pid
locs <- sample(ids, round(length(ids) * 0.30), replace = F )
locs <- sort(locs)

clear$y <- clear$temp_backup <- clear$temp # response variable 
clear[clear$pid %in% locs,]$y <- NA
clear$dataset <- 'train' # creates a train dataset
clear$dataset[locs] <- 'test' # creates a testing dataset


library(MetricGraph)

metric.obj <- metric_graph$new(n$edges$geometry)
# clear.dat <- clear[!is.na(clear$temp_backup),]
clear.dat <- clear

Y.pos <- metric.obj$coordinates(XY = cbind(clear.dat$NEAR_X,clear.dat$NEAR_Y))

clear.dat[["edge_number"]] <- Y.pos[,1]
clear.dat[["distance_on_edge"]] <- Y.pos[,2]

clear.dat[["fact_date"]] <- as.factor(clear.dat[,"date"])
levels(clear.dat[["fact_date"]]) <- 1:length(unique(clear.dat[["fact_date"]]))

metric.obj$add_observations(data = clear.dat, normalized = TRUE, group = "fact_date", clear_obs = TRUE)

## Fitting a model first treating time as replicates 

spde_model_bru <- graph_spde(metric.obj, alpha=1, directional=TRUE)

cmp <- y ~ -1 + Intercept(1) + SLOPE(SLOPE) + elev(elev) + h2o_area(h2o_area) + air_temp(air_temp) + sin_cov(sin) + cos_cov(cos) + 
    field(loc, model = spde_model_bru, replicate = fact_date)

data_spde_dir <- graph_data_spde(graph_spde = spde_model_bru, 
                            loc_name = "loc", repl_col = "fact_date")

library(inlabru)

spde_bru_fit_repl <-
    bru(cmp, data=data_spde_dir[["data"]])

spde_result <- spde_metric_graph_result(spde_bru_fit_repl, "field", spde_model_bru)

summary(spde_result)

sigma_start <- summary(spde_result)[,"mean"][1]
range_start <- summary(spde_result)[,"mean"][2]

spde_model_bru_time <- graph_spde(metric.obj, alpha=1, directional=TRUE)

cmp_time <- y ~ -1 + Intercept(1) + SLOPE(SLOPE) + elev(elev) + h2o_area(h2o_area) + air_temp(air_temp) + sin_cov(sin) + cos_cov(cos) + 
    field(loc, model = spde_model_bru_time, group = fact_date, control.group = list(model = 'ar1')) 

data_spde_time <- graph_data_spde(graph_spde = spde_model_bru_time, 
                            loc_name = "loc", group_col = "fact_date")

spde_bru_fit_time <-
    bru(cmp_time, data=data_spde_time[["data"]], options=list(verbose=TRUE))


spde_result <- spde_metric_graph_result(spde_bru_fit_time, "field", spde_model_bru_time)

summary(spde_result)



####

library(INLA)

data_spde <- graph_data_spde(graph_spde = spde_model_bru_time, name = "field", group_col = "fact_date", covariates = c("SLOPE", "elev", "h2o_area","air_temp", "sin","cos"))
f.s <- y ~ f(SLOPE,model="linear") + f(elev,model="linear") + f(h2o_area,model="linear") + 
f(air_temp,model="linear") + f(sin ,model="linear") + f(cos,model="linear") + f(field, model = spde_model_bru_time, group = field.group, control.group = list(model = 'ar1')) 

stk_dat <- inla.stack(data = data_spde[["data"]], 
                        A = data_spde[["basis"]], 
                        effects = 
      data_spde[["index"]]
    )
data_stk <- inla.stack.data(stk_dat)

spde_fit <- inla(f.s, data = data_stk, control.predictor=list(A=inla.stack.A(stk_dat)))

spde_result2 <- spde_metric_graph_result(spde_fit, "field", spde_model_bru_time)

summary(spde_result2)






