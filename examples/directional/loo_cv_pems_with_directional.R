## --------------------------------------------------------------------------
## Leave-one-out cross-validation comparison on the pems data, extended
## with the directional Whittle--Matern model (WMD1).
##
## The osm-direction graph  is built by the companion script
## `pems_build_osm_graph.R`. Run that script first.
## --------------------------------------------------------------------------

library(MetricGraph)
set.seed(1)


pems_graph <- metric_graph$new(edges = pems$edges, verbose = 0)
pems_graph$add_observations(data = pems$data, normalized = TRUE,
                            verbose = 0)


cat("Fitting WM alpha=1 (BC=0)...\n")
fit_alpha1 <- graph_lme(y ~ 1, graph = pems_graph, BC = 0,
                        model = list(type = "WhittleMatern", alpha = 1),
                        optim_method = "Nelder-Mead")

cat("Fitting WM alpha=2 (BC=0)...\n")
fit_alpha2 <- graph_lme(y ~ 1, graph = pems_graph, BC = 0,
                        model = list(type = "WhittleMatern", alpha = 2),
                        optim_method = "Nelder-Mead")

cat("Fitting WM alpha=1 (BC=1)...\n")
fit_alpha1_bc <- graph_lme(y ~ 1, graph = pems_graph, BC = 1,
                           model = list(type = "WhittleMatern", alpha = 1),
                           optim_method = "Nelder-Mead")

cat("Fitting WM alpha=2 (BC=1)...\n")
fit_alpha2_bc <- graph_lme(y ~ 1, graph = pems_graph, BC = 1,
                           model = list(type = "WhittleMatern", alpha = 2),
                           optim_method = "Nelder-Mead")

cat("Fitting GL alpha=1...\n")
fit_GL1 <- graph_lme(y ~ -1, graph = pems_graph,
                     model = list(type = "graphLaplacian", alpha = 1),
                     optim_method = "Nelder-Mead")
fit_GL1 <- graph_lme(y ~ 1, graph = pems_graph,
                     model = list(type = "graphLaplacian", alpha = 1),
                     previous_fit = fit_GL1,
                     optim_method = "Nelder-Mead")

cat("Fitting GL alpha=2...\n")
fit_GL2 <- graph_lme(y ~ 1, graph = pems_graph,
                     model = list(type = "graphLaplacian", alpha = 2),
                     previous_fit = fit_GL1,
                     optim_method = "Nelder-Mead")

cat("Fitting isoExp...\n")
fit_isoexp <- graph_lme(y ~ 1, graph = pems_graph,
                        model = list(type = "isoCov"),
                        optim_method = "Nelder-Mead")

cat("Fitting WMD alpha=1 (raw graph, default vertex condition)...\n")
fit_wmd1 <- graph_lme(y ~ 1, graph = pems_graph, BC = 0, model = "WMD1",
                      optim_method = "Nelder-Mead")

cat("Fitting WMD alpha=1 (raw graph, variance-stationary)...\n")
pems_graph_vs <- pems_graph$clone()
pems_graph_vs$setDirectionalWeightFunction(
  f_in = function(x) sqrt(x / sum(x)))
fit_wmd1_vs <- graph_lme(y ~ 1, graph = pems_graph_vs, BC = 0, model = "WMD1",
                         optim_method = "Nelder-Mead")

# Directional WMD1 on the OSM-direction graph

osm_meta_path     <- file.path("examples", "directional", "pems_osm.rds")
osmdir_edges_path <- file.path("examples", "directional",
                               "pems_osmdir_edges.rds")

osm_weight_map <- c(motorway       = 8,
                    trunk          = 6,
                    primary        = 4,
                    secondary      = 2,
                    motorway_link  = 1,
                    trunk_link     = 1,
                    primary_link   = 1,
                    secondary_link = 1,
                    unmatched      = 1)


osm_meta     <- readRDS(osm_meta_path)
osmdir_edges <- readRDS(osmdir_edges_path)
pems_graph_osmdir <- metric_graph$new(edges = osmdir_edges,
                                      longlat = TRUE, verbose = 0)
pems_graph_osmdir$add_observations(data = pems$data, normalized = TRUE,
                                   verbose = 0)

cat("Fitting WMD alpha=1 (osm-direction, default vertex condition)...\n")
fit_wmd1_osmdir <- graph_lme(y ~ 1, graph = pems_graph_osmdir,
                             BC = 0, model = "WMD1",
                             optim_method = "Nelder-Mead")

cat("Fitting WMD alpha=1 (osm-direction, variance-stationary)...\n")
pems_graph_osmdir_vs <- pems_graph_osmdir$clone()
pems_graph_osmdir_vs$setDirectionalWeightFunction(
  f_in = function(x) sqrt(x / sum(x)))
fit_wmd1_osmdir_vs <- graph_lme(y ~ 1, graph = pems_graph_osmdir_vs,
                                BC = 0, model = "WMD1",
                                optim_method = "Nelder-Mead")

cat("Fitting WMD alpha=1 (osm-direction, length-weighted)...\n")
pems_graph_osmdir_len <- pems_graph_osmdir$clone()
pems_graph_osmdir_len$set_edge_weights(
  weights = as.numeric(pems_graph_osmdir_len$edge_lengths),
  verbose = 0)
fit_wmd1_osmdir_len <- graph_lme(y ~ 1, graph = pems_graph_osmdir_len,
                                 BC = 0, model = "WMD1",
                                 optim_method = "Nelder-Mead")

cat("Fitting WMD alpha=1 (osm-direction, road-type-weighted)...\n")
pems_graph_osmdir_type <- pems_graph_osmdir$clone()
edge_len_osmdir   <- as.numeric(pems_graph_osmdir_type$edge_lengths)
type_weight_osmdir <- ifelse(edge_len_osmdir > 0.5, 5, 1)
pems_graph_osmdir_type$set_edge_weights(weights = type_weight_osmdir,
                                        verbose = 0)
fit_wmd1_osmdir_type <- graph_lme(y ~ 1, graph = pems_graph_osmdir_type,
                                  BC = 0, model = "WMD1",
                                  optim_method = "Nelder-Mead")

cat("Fitting WMD alpha=1 (osm-direction, OSM-class-weighted)...\n")
osm_w <- osm_weight_map[osm_meta$highway]; osm_w[is.na(osm_w)] <- 1
pems_graph_osmdir_osmw <- pems_graph_osmdir$clone()
pems_graph_osmdir_osmw$set_edge_weights(weights = as.numeric(osm_w),
                                        verbose = 0)
fit_wmd1_osmdir_osmw <- graph_lme(y ~ 1, graph = pems_graph_osmdir_osmw,
                                  BC = 0, model = "WMD1",
                                  optim_method = "Nelder-Mead")


# Collect models and run LOO CV

fitted_models_list <- list(
  "isoExp"      = fit_isoexp,
  "GL1"         = fit_GL1,
  "GL2"         = fit_GL2,
  "alpha=1"     = fit_alpha1,
  "alpha=1 bc"  = fit_alpha1_bc,
  "alpha=2"     = fit_alpha2,
  "alpha=2 bc"  = fit_alpha2_bc,
  "WMD1"                = fit_wmd1,
  "WMD1 var-stationary" = fit_wmd1_vs)

fitted_models_list[["WMD1 osm-direction"]] <- fit_wmd1_osmdir
fitted_models_list[["WMD1 osm-direction var-stationary"]] <- fit_wmd1_osmdir_vs
fitted_models_list[["WMD1 osm-direction len-weighted"]] <- fit_wmd1_osmdir_len
fitted_models_list[["WMD1 osm-direction type-weighted"]] <- fit_wmd1_osmdir_type
fitted_models_list[["WMD1 osm-direction osm-weighted"]] <- fit_wmd1_osmdir_osmw

cat("\n--- Negative log-likelihoods ---\n")
print(-sapply(fitted_models_list, logLik))

cat("\nRunning posterior_crossvalidation_loo...\n")
res <- posterior_crossvalidation_loo(fitted_models_list, tibble = FALSE)

cat("\n--- LOO CV scores (lower is better) ---\n")
print(round(res$scores, 3))
