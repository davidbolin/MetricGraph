
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

metric.obj <- metric_graph$new(n$edges$geometry, perform_merges=FALSE)
# clear.dat <- clear[!is.na(clear$temp_backup),]
clear.dat <- clear

Y.pos <- metric.obj$coordinates(XY = cbind(clear.dat$NEAR_X,clear.dat$NEAR_Y))

clear.dat[["edge_number"]] <- Y.pos[,1]
clear.dat[["distance_on_edge"]] <- Y.pos[,2]

clear.dat[["fact_date"]] <- as.factor(clear.dat[,"date"])
levels(clear.dat[["fact_date"]]) <- 1:length(unique(clear.dat[["fact_date"]]))

library(fastDummies)
clear.dat[["month"]] <- factor(as.character(as.numeric(clear.dat[["fact_date"]]) %% 12))

clear.dat <- dummy_cols(clear.dat, select_columns = "month")

metric.obj$add_observations(data = clear.dat, normalized = TRUE, group = "fact_date", clear_obs = TRUE)

## Fitting a model first treating time as replicates

res.wm1.dir <- graph_lme(y ~ SLOPE + elev + h2o_area + air_temp + sin + cos,
    model = list(type = "WhittleMatern", fem = FALSE, alpha = 1, version = 1, directional=0), graph = metric.obj)

summary(res.wm1.dir)

spde_model_bru <- graph_spde(metric.obj, alpha=1, directional=TRUE)

data_spde_dir <- graph_data_spde(graph_spde = spde_model_bru,
                            loc_name = "loc", repl_col = "fact_date")

repl <- data_spde_dir[["repl"]]

cmp <- y ~ -1 + Intercept(1) + month_1(month_1) +
    month_2(month_2) + month_5(month_5) + month_8(month_8) +
    month_3(month_3) + month_6(month_6) + month_9(month_9) +
    month_4(month_4) + month_7(month_7) + month_10(month_10) +
    month_11(month_11) +
    SLOPE(SLOPE) + elev(elev) + h2o_area(h2o_area) + air_temp(air_temp) + sin_cov(sin) + cos_cov(cos) +
    field(loc, model = spde_model_bru, replicate = fact_date)

library(inlabru)

spde_bru_fit_repl <-
    bru(cmp, data=data_spde_dir[["data"]], options = list(verbose=TRUE))

spde_result <- spde_metric_graph_result(spde_bru_fit_repl, "field", spde_model_bru)

summary(spde_result)
summary(spde_bru_fit_repl)

## Prediction

graph_data <- metric.obj$get_data()
idx_pred <- which(graph_data$dataset == "test")
data_list <- lapply(graph_data, function(i){i[idx_pred]})
data_list[["loc"]] <- cbind(graph_data[idx_pred, ".edge_number"], graph_data[idx_pred,".distance_on_edge"])


field_pred <- predict(spde_model_bru,
                                cmp,
                                spde_bru_fit_repl,
                                newdata = data_list,
                                repl_col = "fact_date",
                                formula = ~ Intercept + SLOPE + elev + h2o_area + air_temp + sin_cov + cos_cov + field)







## Fitting the space-time model

spde_model_bru_time <- graph_spde(metric.obj, alpha=1, directional=TRUE)

data_spde_time <- graph_data_spde(graph_spde = spde_model_bru_time,
                            loc_name = "loc", group_col = "fact_date")

group <- data_spde_time[["group"]]

# cmp_time <- y ~ -1 + Intercept(1) + month_1(month_1) +
#     month_2(month_2) + month_5(month_5) + month_8(month_8) +
#     month_3(month_3) + month_6(month_6) + month_9(month_9) +
#     month_4(month_4) + month_7(month_7) + month_10(month_10) +
#     month_11(month_11)  + SLOPE(SLOPE) + elev(elev) + h2o_area(h2o_area) + air_temp(air_temp) + sin_cov(sin) + cos_cov(cos) +
#     field(loc, model = spde_model_bru_time, group = fact_date, control.group = list(model = 'ar1'))

cmp_time <- y ~ -1 + Intercept(1) + SLOPE(SLOPE) + elev(elev) + h2o_area(h2o_area) + air_temp(air_temp) + sin_cov(sin) + cos_cov(cos) +
    field(loc, model = spde_model_bru_time, group = fact_date, control.group = list(model = 'ar1'))

spde_bru_fit_time <-
    bru(cmp_time, data=data_spde_time[["data"]], options=list(verbose=TRUE))

spde_result <- spde_metric_graph_result(spde_bru_fit_time, "field", spde_model_bru_time)

summary(spde_result)
summary(spde_bru_fit_time)

### Prediction

data_pred <- clear

data_spde_time <- graph_data_spde(graph_spde = spde_model_bru_time,
                            loc_name = "loc", group_col = "fact_date")

graph_data <- metric.obj$get_data()
idx_pred <- which(graph_data$dataset == "test")
data_list <- lapply(graph_data, function(i){i[idx_pred]})
data_list[["loc"]] <- cbind(graph_data[idx_pred, ".edge_number"], graph_data[idx_pred,".distance_on_edge"])


group <- data_list$fact_date
# field_pred <- predict(spde_model_bru_time,
#                                 cmp_time,
#                                 spde_bru_fit_time,
#                                 newdata = data_list,
#                                 group_col = "fact_date",
#                                 formula = ~ Intercept + month_1 + month_2 + month_3 + month_4 + 
#                                 month_5 + month_6 + month_7 + month_8 + month_9 + month_10 + 
#                                 month_11 + SLOPE + elev + h2o_area + air_temp + sin_cov + cos_cov + field)

field_pred <- predict(spde_model_bru_time,
                                cmp_time,
                                spde_bru_fit_time,
                                newdata = data_list,
                                group_col = "fact_date",
                                formula = ~ Intercept + SLOPE + elev + h2o_area + air_temp + sin_cov + cos_cov + field)



SSN2::ssn_create_distmat(n, "preds" , overwrite=TRUE, among_predpts = TRUE)
net_num <- as.numeric(as.character(distinct(rbind(
distinct(SSN2::ssn_get_data(n), netID),
distinct(SSN2::ssn_get_data(n, name = "preds"), netID))) %>% pull()))

use_saved <- T #NB; set to False if want to run the model 


if(use_saved == T) { # loading the model fit from Github
    url <- "https://github.com/EdgarSantos-Fernandez/SSNdata/raw/main/inst/extdata/fit_ar.rds"
    fit_ar <- readRDS(url(url, method="libcurl"))
    #fit_ar <- readRDS(system.file("extdata//fit_ar.rds", package = "SSNdata"))
}

library('StanHeaders')
library('rstan')

fit <- fit_ar
class(fit) <- c("stanfit") # changing the class, so we can use functions for stanfit objects 

#Creating a summary statistics (man, SD, percentiles, etc.)
stats <- summary(fit)$summary %>%
                       data.frame()
stats <- data.frame(Variable = rownames(stats), stats)



stats[1:15,]

ypred <- data.frame(stats[grep("y\\[", row.names(stats)),])
ypred$ytrue <- clear$temp_backup
ypred$date <- rep(unique(clear$date), each = length(unique(clear$locID)) )
ypred$locID <- rep(1:length(unique(clear$locID)), times = length(unique(clear$date)))
ypred$dataset <- ifelse(ypred$sd == 0, 'obs', 'pred')
ypred$ypred <- ypred$mean
ypred$y <- clear$temp_backup
ypred$temp_pred_2.5 <- stats[grep('y',rownames(stats)),'X2.5.']
ypred$temp_pred_97.5 <- stats[grep('y',rownames(stats)),'X97.5.']
ypred$pid <- 1:nrow(ypred)

idx_pred <- match(data_list$temp_backup, ypred$y)

pred_ssn <- ypred[idx_pred,]$mean
pred_arwdm <- field_pred$pred$mean
pred_replwdm <- field_pred_repl$pred$mean

sum((pred_ssn-data_list$temp_backup)^2)
sum((pred_arwdm-data_list$temp_backup)^2)
sum((pred_replwdm-data_list$temp_backup)^2)
