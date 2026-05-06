## ---------------------------------------------------------------------------
## Cross-validation for inlabru fits on metric graphs.
##
## Mirrors the API of rSPDE::cross_validation. For models containing an
## inla_metric_graph_spde component (exact, non-FEM SPDE on metric graphs)
## this dispatches to a metric-graph-aware refit/sample path that follows
## the augmented-graph pattern used by predict.inla_metric_graph_spde.
## For other components (rSPDE / FEM / generic), it uses the standard
## inlabru::bru_set_missing + bru_rerun + inlabru::generate path.
## ---------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Helpers shared with rSPDE (copied verbatim to avoid ::: usage).
# -----------------------------------------------------------------------------

#' @noRd
.cv_process_formula <- function(bru_result) {
  form <- bru_result$bru_info$model$formula[3]
  form <- as.character(form)
  form <- strsplit(form, "f\\(")
  form <- form[[1]]
  form <- form[-1]
  form_proc <- sub(",.*", "", strsplit(form, "f\\(")[1])
  if (length(form) > 1) {
    for (i in 2:(length(form))) {
      form_proc <- paste(form_proc, " + ", sub(",.*", "", strsplit(form, "f\\(")[i]))
    }
  }
  form_proc <- paste("~", "linkfuninv(", form_proc, ")")
  return(stats::as.formula(form_proc))
}

#' @noRd
.cv_process_formula_lhoods <- function(bru_result, like_number) {
  form <- inlabru::as_bru_obs_list(bru_result)[[like_number]]$formula[3]
  form <- as.character(form)
  if (form == ".") {
    return(.cv_process_formula(bru_result))
  }
  form_proc <- paste("~", "linkfuninv(", form, ")")
  return(stats::as.formula(form_proc))
}

#' @noRd
.cv_process_link <- function(link_name) {
  switch(link_name,
    "log" = function(x) INLA::inla.link.log(x, inverse = TRUE),
    "invlog" = function(x) INLA::inla.link.invlog(x, inverse = TRUE),
    "logit" = function(x) INLA::inla.link.logit(x, inverse = TRUE),
    "invlogit" = function(x) INLA::inla.link.invlogit(x, inverse = TRUE),
    "probit" = function(x) INLA::inla.link.probit(x, inverse = TRUE),
    "invprobit" = function(x) INLA::inla.link.invprobit(x, inverse = TRUE),
    "cloglog" = function(x) INLA::inla.link.cloglog(x, inverse = TRUE),
    "invcloglog" = function(x) INLA::inla.link.invcloglog(x, inverse = TRUE),
    "tan" = function(x) INLA::inla.link.tan(x, inverse = TRUE),
    "invtan" = function(x) INLA::inla.link.invtan(x, inverse = TRUE),
    "identity" = function(x) INLA::inla.link.identity(x, inverse = TRUE),
    "invidentity" = function(x) INLA::inla.link.invidentity(x, inverse = TRUE)
  )
}

#' @noRd
.cv_default_linkinv <- function(model_family) {
  if (model_family %in% c("gaussian", "t")) {
    function(x) x
  } else if (model_family %in% c("gamma", "poisson")) {
    function(x) exp(x)
  } else if (model_family == "binomial") {
    function(x) exp(x) / (1 + exp(x))
  } else {
    stop(sprintf(
      "The family '%s' is not supported by cross_validation in MetricGraph.",
      model_family
    ))
  }
}

#' @noRd
.cv_map_models_to_strings <- function(bru_fit) {
  families <- bru_fit$.args$family
  mapping <- list(
    "gaussian" = "Precision for the Gaussian observations",
    "gamma" = "Precision-parameter for the Gamma observations",
    "poisson" = ".none",
    "binomial" = ".none",
    "t" = c(
      "precision for the student-t observations",
      "degrees of freedom for student-t"
    )
  )
  result <- vector("list", length(families))
  for (i in seq_along(families)) {
    family <- families[i]
    if (!(family %in% names(mapping))) {
      stop(sprintf(
        "The family '%s' is not supported by cross_validation in MetricGraph.",
        family
      ))
    }
    base_string <- mapping[[family]]
    if (i == 1 || base_string[1] == ".none") {
      result[[i]] <- base_string
    } else {
      result[[i]] <- paste0(base_string, "[", i, "]")
    }
  }
  result
}


# -----------------------------------------------------------------------------
# Detection of a metric_graph SPDE component inside a bru fit.
# -----------------------------------------------------------------------------

#' @noRd
.cv_find_mg_component <- function(bru_fit) {
  effects <- bru_fit$bru_info$model$effects
  if (is.null(effects) || length(effects) == 0) {
    return(NULL)
  }
  candidates <- list(
    function(eff) eff$main$model,
    function(eff) eff$main$mapper$model,
    function(eff) eff$env$model,
    function(eff) eff$mapper$model
  )
  for (eff_name in names(effects)) {
    eff <- effects[[eff_name]]
    model <- NULL
    for (getter in candidates) {
      m <- tryCatch(getter(eff), error = function(e) NULL)
      if (!is.null(m) && inherits(m, "inla_metric_graph_spde")) {
        model <- m
        break
      }
    }
    if (is.null(model)) {
      mapper <- tryCatch(eff$mapper, error = function(e) NULL)
      if (!is.null(mapper) &&
        inherits(mapper, "bru_mapper_inla_metric_graph_spde")) {
        model <- mapper$model
      }
    }
    if (!is.null(model)) {
      input <- tryCatch(eff$main$input$input, error = function(e) NULL)
      location_var <- if (is.null(input)) "loc" else .bru_input_to_name(input)
      return(list(
        effect_name = eff_name,
        location_var = location_var,
        spde_model = model
      ))
    }
  }
  return(NULL)
}

#' @noRd
.cv_find_spde_var_name <- function(bru_fit) {
  fenv <- environment(bru_fit$bru_info$model$formula)
  if (is.null(fenv)) return(NULL)
  candidates <- character(0)
  envs_to_search <- list(fenv)
  parent <- parent.env(fenv)
  while (!identical(parent, emptyenv()) &&
    !identical(parent, globalenv()) &&
    length(envs_to_search) < 20) {
    envs_to_search[[length(envs_to_search) + 1]] <- parent
    parent <- tryCatch(parent.env(parent), error = function(e) emptyenv())
  }
  envs_to_search[[length(envs_to_search) + 1]] <- globalenv()
  for (env in envs_to_search) {
    nms <- ls(env, all.names = TRUE)
    for (nm in nms) {
      val <- tryCatch(get(nm, envir = env, inherits = FALSE),
        error = function(e) NULL
      )
      if (inherits(val, "inla_metric_graph_spde")) {
        candidates <- c(candidates, nm)
      }
    }
    if (length(candidates) > 0) break
  }
  if (length(candidates) == 0) return(NULL)
  candidates[1]
}


# -----------------------------------------------------------------------------
# Refit on training data only.
#  - For metric-graph SPDE components: build a fresh spde object on a clone
#    of the original graph with NA at test indices, rebuild the data list,
#    rebuild cmp by replacing the spde model name, and call inlabru::bru().
#  - Otherwise: use bru_set_missing + bru_rerun.
# -----------------------------------------------------------------------------

#' @noRd
.cv_bru_rerun_generic <- function(bru_fit, idx_data, true_CV, fit_verbose,
                                  model_options_bru) {
  if (is.null(model_options_bru)) model_options_bru <- list()
  options <- model_options_bru
  if (!true_CV) {
    options$control.mode <- list(theta = bru_fit$mode$theta, fixed = TRUE)
  }
  options$verbose <- isTRUE(fit_verbose)
  info <- inlabru::as_bru_info(bru_fit)
  options <- inlabru::bru_options(
    info[["options"]],
    inlabru::as.bru_options(options)
  )
  result <- inlabru::bru_set_missing(bru_fit, keep = idx_data)
  result <- inlabru::bru_rerun(result, options = options)
  result
}

#' @noRd
.cv_bru_rerun_metric_graph <- function(bru_fit, idx_data_per_lik, true_CV,
                                       fit_verbose, model_options_bru,
                                       mg_info) {
  spde_model <- mg_info$spde_model
  loc_name <- mg_info$location_var

  # Find variable name of the spde model in the formula's environment so we
  # can splice a fresh model into the cmp via text substitution. The fallback
  # value matches predict.inla_metric_graph_spde's substitution target.
  spde_var_name <- .cv_find_spde_var_name(bru_fit)

  orig_graph <- spde_model$graph_spde
  graph_tmp <- orig_graph$get_initial_graph()
  graph_tmp$clear_observations()

  original_data <- spde_model$.__enclos_env__$private$data
  if (is.null(original_data)) {
    original_data <- bru_fit$bru_info$lhoods[[1]]$data
  }

  group_variables <- attr(
    spde_model$graph_spde$.__enclos_env__$private$data,
    "group_variable"
  )
  if (is.null(group_variables)) group_variables <- ".none"

  if (group_variables == ".none") {
    graph_tmp$add_observations(
      data = original_data,
      edge_number = ".edge_number",
      distance_on_edge = ".distance_on_edge",
      data_coords = "PtE",
      normalized = TRUE,
      verbose = 0,
      suppress_warnings = TRUE
    )
  } else {
    graph_tmp$add_observations(
      data = original_data,
      edge_number = ".edge_number",
      distance_on_edge = ".distance_on_edge",
      data_coords = "PtE",
      normalized = TRUE,
      group = group_variables,
      verbose = 0,
      suppress_warnings = TRUE
    )
  }

  graph_tmp$observation_to_vertex(mesh_warning = FALSE)

  # Mark held-out responses as NA in the cloned graph data.
  responses <- bru_fit$.args$family
  if (length(idx_data_per_lik) != length(responses)) {
    if (length(responses) == 1L) {
      idx_data_per_lik <- list(idx_data_per_lik[[1]])
    }
  }
  lhoods <- inlabru::as_bru_obs_list(bru_fit)
  for (i_lik in seq_along(lhoods)) {
    response_name <- as.character(lhoods[[i_lik]]$formula[[2]])
    if (length(response_name) > 1) response_name <- response_name[1]
    if (!is.null(graph_tmp$.__enclos_env__$private$data[[response_name]])) {
      keep_idx <- idx_data_per_lik[[i_lik]]
      n_obs_lik <- length(graph_tmp$.__enclos_env__$private$data[[response_name]])
      drop_idx <- setdiff(seq_len(n_obs_lik), keep_idx)
      graph_tmp$.__enclos_env__$private$data[[response_name]][drop_idx] <- NA
    }
  }

  spde____model <- graph_spde(graph_tmp,
    alpha = spde_model$alpha,
    directional = spde_model$directional
  )

  cmp_orig <- bru_fit$bru_info$model$formula
  cmp_c <- as.character(cmp_orig)
  if (!is.null(spde_var_name)) {
    pattern <- paste0("(?<![A-Za-z0-9_.])", spde_var_name, "(?![A-Za-z0-9_.])")
    cmp_c[3] <- sub(pattern, "spde____model", cmp_c[3], perl = TRUE)
  } else {
    cmp_c[3] <- sub(
      "model\\s*=\\s*[A-Za-z_.][A-Za-z0-9_.]*",
      "model = spde____model",
      cmp_c[3]
    )
  }
  cmp_new <- stats::as.formula(paste(cmp_c[2], cmp_c[1], cmp_c[3]))
  environment(cmp_new) <- new.env(parent = environment(cmp_orig))
  assign("spde____model", spde____model, envir = environment(cmp_new))

  data_spde_new <- graph_data_spde(spde____model,
    loc_name = loc_name,
    drop_all_na = FALSE, drop_na = FALSE
  )[["data"]]

  options <- if (is.null(model_options_bru)) list() else model_options_bru
  if (!true_CV) {
    options$control.mode <- list(theta = bru_fit$mode$theta, fixed = TRUE)
  }
  options$verbose <- isTRUE(fit_verbose)
  info <- bru_fit[["bru_info"]]
  bru_opts <- inlabru::bru_options(
    info[["options"]],
    inlabru::as.bru_options(options)
  )

  bru_fit_new <- inlabru::bru(cmp_new,
    data = data_spde_new,
    options = bru_opts,
    allow_combine = FALSE
  )

  attr(bru_fit_new, "mg_info") <- list(
    spde_model_train = spde____model,
    location_var = loc_name,
    graph_train = graph_tmp,
    cmp = cmp_new
  )
  bru_fit_new
}


# -----------------------------------------------------------------------------
# Posterior linear-predictor sampling at test locations.
# -----------------------------------------------------------------------------

#' @noRd
.cv_sample_post_lp_generic <- function(bru_fit, i_lik, test_list,
                                       n_samples, print) {
  link_name <- bru_fit$.args$control.family[[i_lik]]$link
  model_family <- bru_fit$.args$family[[i_lik]]
  if (link_name == "default") {
    linkfuninv <- .cv_default_linkinv(model_family)
  } else {
    linkfuninv <- .cv_process_link(link_name)
  }
  formula_tmp <- .cv_process_formula_lhoods(bru_fit, i_lik)
  env_tmp <- environment(formula_tmp)
  assign("linkfuninv", linkfuninv, envir = env_tmp)
  if (print) cat("Generating samples...\n")
  data <- inlabru::as_bru_obs_list(bru_fit)[[i_lik]]$data
  post <- inlabru::generate(bru_fit,
    newdata = data,
    formula = formula_tmp, n.samples = n_samples
  )
  post[test_list[[i_lik]], , drop = FALSE]
}

#' @noRd
.cv_sample_post_lp_metric_graph <- function(bru_fit_train, i_lik, test_list,
                                            n_samples, print, full_bru_fit) {
  link_name <- full_bru_fit$.args$control.family[[i_lik]]$link
  model_family <- full_bru_fit$.args$family[[i_lik]]
  if (link_name == "default") {
    linkfuninv <- .cv_default_linkinv(model_family)
  } else {
    linkfuninv <- .cv_process_link(link_name)
  }
  formula_tmp <- .cv_process_formula_lhoods(full_bru_fit, i_lik)
  env_tmp <- environment(formula_tmp)
  assign("linkfuninv", linkfuninv, envir = env_tmp)
  if (print) cat("Generating samples...\n")

  data <- inlabru::as_bru_obs_list(bru_fit_train)[[i_lik]]$data
  post <- inlabru::generate(bru_fit_train,
    newdata = data,
    formula = formula_tmp, n.samples = n_samples
  )
  post[test_list[[i_lik]], , drop = FALSE]
}


# -----------------------------------------------------------------------------
# Posterior response samples (linear predictor + likelihood noise).
# -----------------------------------------------------------------------------

#' @noRd
.cv_get_response_samples <- function(post_linear_predictors, fit_for_hyper,
                                     i_lik, n_samples_total, print) {
  model_family <- fit_for_hyper$.args$family[[i_lik]]
  family_mappings <- .cv_map_models_to_strings(fit_for_hyper)

  if (family_mappings[[i_lik]][[1]] != ".none") {
    meas_err_par <- lapply(family_mappings[[i_lik]], function(param) {
      hyper_sample <- INLA::inla.hyperpar.sample(n_samples_total,
        fit_for_hyper,
        improve.marginals = TRUE
      )
      hyper_sample[, param]
    })
  }

  if (model_family == "gaussian") {
    sd_sample <- 1 / sqrt(as.vector(meas_err_par[[1]]))
    Y_sample <- lapply(seq_len(nrow(post_linear_predictors)), function(i) {
      post_linear_predictors[i, ] + sd_sample * stats::rnorm(n_samples_total)
    })
  } else if (model_family == "gamma") {
    phi_sample <- as.vector(meas_err_par[[1]])
    Y_sample <- lapply(seq_len(nrow(post_linear_predictors)), function(i) {
      scale_temp <- post_linear_predictors[i, ] / phi_sample
      stats::rgamma(n_samples_total, shape = phi_sample, scale = scale_temp)
    })
  } else if (model_family == "t") {
    sd_sample <- 1 / sqrt(as.vector(meas_err_par[[1]]))
    deg_sample <- as.vector(meas_err_par[[2]])
    Y_sample <- lapply(seq_len(nrow(post_linear_predictors)), function(i) {
      post_linear_predictors[i, ] +
        sd_sample * stats::rt(n_samples_total, df = deg_sample)
    })
  } else if (model_family == "poisson") {
    Y_sample <- lapply(seq_len(nrow(post_linear_predictors)), function(i) {
      stats::rpois(n_samples_total, post_linear_predictors[i, ])
    })
  } else if (model_family == "binomial") {
    Y_sample <- lapply(seq_len(nrow(post_linear_predictors)), function(i) {
      stats::rbinom(n_samples_total,
        size = 1, prob = post_linear_predictors[i, ]
      )
    })
  } else {
    stop(sprintf(
      "The family '%s' is not supported by cross_validation in MetricGraph.",
      model_family
    ))
  }
  if (print) cat("Samples generated!\n")
  Y_sample
}


# -----------------------------------------------------------------------------
# Main user-facing function.
# -----------------------------------------------------------------------------

#' Perform cross-validation on a list of fitted inlabru models on metric graphs.
#'
#' Mirrors [rSPDE::cross_validation()] for `bru` fits (output from `inlabru::bru()`),
#' with built-in support for the exact, non-FEM SPDE models in
#' \pkg{MetricGraph} (objects of class `inla_metric_graph_spde`). For models
#' fit with such a component, the function rebuilds the SPDE on a graph clone
#' with held-out responses set to NA, refits when `true_CV = TRUE`, and draws
#' posterior samples via the `inlabru::generate` path used by
#' [`predict.inla_metric_graph_spde`]. For models without a metric-graph SPDE
#' component (e.g. FEM-based rSPDE models), it uses the standard
#' `inlabru::bru_set_missing()` + `bru_rerun()` + `inlabru::generate()` path
#' shared with [rSPDE::cross_validation()].
#'
#' @param models A fitted model from `inlabru::bru()` or a list of such models.
#'   All models must have the same number of likelihoods and be fitted to
#'   identical datasets.
#' @param model_names Character vector of model names. Defaults to `names(models)`
#'   or `Model 1`, `Model 2`, ... if absent.
#' @param scores Subset of `c("mae", "mse", "crps", "scrps", "dss", "wcrps", "swcrps")`.
#' @param cv_type One of `"k-fold"`, `"loo"`, `"lpo"`.
#' @param weight_thr Threshold for `wcrps` / `swcrps`. Required if either is requested.
#' @param k Number of folds for `k-fold`.
#' @param percentage Train percentage (1-99) for `lpo`.
#' @param number_folds Number of folds for `lpo`.
#' @param n_samples Number of posterior samples for sample-based scoring.
#' @param return_scores_folds If `TRUE`, return per-fold score matrices.
#' @param orientation_results One of `"negative"` (lower is better) or
#'   `"positive"` (higher is better).
#' @param include_best Add a `Best` row indicating the best model per score.
#' @param train_test_indexes Optional pre-built fold list. If supplied,
#'   `cv_type`, `k`, `percentage`, `number_folds` are ignored.
#' @param return_train_test Return the train/test indices used.
#' @param return_post_samples Return posterior response samples (forces
#'   `return_scores_folds = TRUE`).
#' @param return_true_test_values Return the true response values at test points.
#' @param parallelize_RP Parallelize CRPS/SCRPS pairwise integrals.
#' @param n_cores_RP Number of cores for `parallelize_RP`.
#' @param true_CV If `TRUE`, refit the model on each training fold; if `FALSE`,
#'   sample directly from the supplied fit without refitting.
#' @param save_settings If `TRUE`, return the CV settings used.
#' @param model_options_bru List of options passed to inlabru.
#' @param print Print partial progress.
#' @param fit_verbose Pass `verbose` through to INLA when refitting.
#'
#' @return Either a `data.frame` (default), or a list with `scores_df` plus the
#'   optional components requested by `return_*` / `save_settings` arguments.
#' @export
cross_validation <- function(models, model_names = NULL,
                             scores = c("mae", "mse", "crps", "scrps", "dss"),
                             cv_type = c("k-fold", "loo", "lpo"),
                             weight_thr = NULL,
                             k = 5, percentage = 20, number_folds = 10,
                             n_samples = 1000,
                             return_scores_folds = FALSE,
                             orientation_results = c("negative", "positive"),
                             include_best = TRUE,
                             train_test_indexes = NULL,
                             return_train_test = FALSE,
                             return_post_samples = FALSE,
                             return_true_test_values = FALSE,
                             parallelize_RP = FALSE,
                             n_cores_RP = parallel::detectCores() - 1,
                             true_CV = TRUE, save_settings = FALSE,
                             model_options_bru = list(),
                             print = TRUE,
                             fit_verbose = FALSE) {
  if (!requireNamespace("inlabru", quietly = TRUE)) {
    stop("Package 'inlabru' is required for cross_validation().")
  }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required for cross_validation().")
  }

  orientation_results <- orientation_results[[1]]
  if (!(orientation_results %in% c("positive", "negative"))) {
    stop("orientation_results must be either 'positive' or 'negative'!")
  }

  scores <- tolower(scores)
  if (any(scores %in% c("wcrps", "swcrps")) && is.null(weight_thr)) {
    stop("weight_thr must be supplied if 'wcrps' or 'swcrps' are requested!")
  }
  scores <- intersect(scores, c("mae", "mse", "crps", "scrps", "dss",
    "wcrps", "swcrps"))

  cv_type <- cv_type[[1]]
  if (!(cv_type %in% c("k-fold", "loo", "lpo"))) {
    stop("cv_type must be one of 'k-fold', 'loo', or 'lpo'!")
  }
  if (!is.numeric(percentage)) stop("percentage must be numeric!")
  if (percentage %% 1 != 0) {
    warning("Non-integer percentage given, rounding.")
    percentage <- round(percentage)
  }
  if (percentage <= 0 || percentage >= 100) {
    stop("percentage must be a number between 1 and 99!")
  }
  if (!is.numeric(number_folds)) stop("number_folds must be numeric!")
  if (number_folds %% 1 != 0) {
    warning("Non-integer number_folds given, rounding.")
    number_folds <- round(number_folds)
  }
  if (number_folds <= 0) stop("number_folds must be positive!")

  if (inherits(models, "bru")) {
    models <- list(models)
  } else {
    for (i in seq_along(models)) {
      if (!inherits(models[[i]], "bru")) {
        stop("models must be a bru fit or a list of bru fits.")
      }
    }
  }

  lhoods <- inlabru::as_bru_obs_list(models[[1]])
  n_likelihoods <- length(lhoods)
  for (mn in seq_along(models)) {
    if (length(inlabru::as_bru_obs_list(models[[mn]])) != n_likelihoods) {
      stop(paste("Model", mn,
        "does not have the same number of likelihoods as the first model."))
    }
  }

  if (is.null(model_names) && is.list(models)) {
    model_names <- names(models)
  }
  if (!is.null(model_names)) {
    if (!is.character(model_names)) stop("model_names must be a character vector!")
    if (length(models) != length(model_names)) {
      stop("model_names must contain one name per model!")
    }
  } else {
    model_names <- paste("Model", seq_along(models))
  }

  if (!is.numeric(n_samples)) stop("n_samples must be numeric!")
  if (n_samples %% 1 != 0) {
    warning("Non-integer n_samples given, rounding.")
    n_samples <- round(n_samples)
  }
  if (n_samples <= 0) stop("n_samples must be positive!")

  cluster_tmp <- NULL
  if (parallelize_RP) {
    if (!requireNamespace("parallel", quietly = TRUE) ||
      !requireNamespace("doParallel", quietly = TRUE) ||
      !requireNamespace("foreach", quietly = TRUE)) {
      warning("'parallel', 'doParallel' and 'foreach' are required for parallelize_RP. Disabling.")
      parallelize_RP <- FALSE
    } else {
      cluster_tmp <- parallel::makeCluster(n_cores_RP)
      doParallel::registerDoParallel(cluster_tmp)
    }
  }

  if (is.null(train_test_indexes)) {
    lhoods <- inlabru::as_bru_obs_list(models[[1]])
    data_list <- lapply(seq_len(n_likelihoods), function(i) lhoods[[i]]$data)
    train_test_indexes <- rSPDE::create_train_test_indices(data_list,
      cv_type = cv_type, k = k,
      percentage = percentage, number_folds = number_folds
    )
  } else {
    if (!is.list(train_test_indexes)) stop("train_test_indexes must be a list.")
    for (i in seq_along(train_test_indexes)) {
      fe <- train_test_indexes[[i]]
      if (!is.list(fe) || is.null(fe$train) || is.null(fe$test)) {
        stop(sprintf("train_test_indexes[[%d]] must be a list with 'train' and 'test'.", i))
      }
      if (!is.list(fe$train) || !is.list(fe$test)) {
        stop(sprintf("train_test_indexes[[%d]]$train and $test must be lists.", i))
      }
      if (length(fe$train) != n_likelihoods || length(fe$test) != n_likelihoods) {
        stop(sprintf("train/test entries in fold %d must have length %d.", i, n_likelihoods))
      }
    }
  }

  n_folds <- length(train_test_indexes)
  n_models <- length(models)

  post_samples <- list()
  true_test_values <- list()
  for (mn in seq_len(n_models)) {
    post_samples[[model_names[[mn]]]] <- vector("list", length = n_folds)
    if (return_true_test_values) {
      true_test_values[[model_names[[mn]]]] <- vector("list", length = n_folds)
    }
    for (j in seq_len(n_folds)) {
      post_samples[[model_names[[mn]]]][[j]] <- vector("list", length = n_likelihoods)
      if (return_true_test_values) {
        true_test_values[[model_names[[mn]]]][[j]] <- vector("list", length = n_likelihoods)
      }
    }
  }

  init_score_mat <- function() {
    lapply(seq_len(n_likelihoods), function(i) {
      matrix(numeric(n_folds * n_models), ncol = n_models)
    })
  }
  dss <- init_score_mat(); mse <- init_score_mat(); mae <- init_score_mat()
  crps <- init_score_mat(); scrps <- init_score_mat()
  wcrps <- init_score_mat(); swcrps <- init_score_mat()

  needs_paired_samples <- any(c("crps", "scrps", "dss", "wcrps", "swcrps") %in% scores)
  new_n_samples <- if (needs_paired_samples) 2 * n_samples else n_samples

  # Per-model dispatch flags.
  is_mg <- vapply(models, function(m) {
    !is.null(.cv_find_mg_component(m))
  }, logical(1))

  for (fold in seq_len(n_folds)) {
    train_list <- train_test_indexes[[fold]][["train"]]
    test_list <- train_test_indexes[[fold]][["test"]]
    for (mn in seq_len(n_models)) {
      if (print) {
        cat(sprintf("Fold: %d / %d\n", fold, n_folds))
        cat(sprintf("Model: %s\n", model_names[[mn]]))
      }
      cur_fit <- models[[mn]]
      mg_info <- if (is_mg[mn]) .cv_find_mg_component(cur_fit) else NULL

      if (true_CV) {
        if (is_mg[mn]) {
          new_fit <- .cv_bru_rerun_metric_graph(cur_fit,
            idx_data_per_lik = train_list,
            true_CV = TRUE, fit_verbose = fit_verbose,
            model_options_bru = model_options_bru,
            mg_info = mg_info
          )
        } else {
          new_fit <- .cv_bru_rerun_generic(cur_fit,
            idx_data = train_list,
            true_CV = TRUE, fit_verbose = fit_verbose,
            model_options_bru = model_options_bru
          )
        }
      } else {
        new_fit <- cur_fit
      }
      fit_for_hyper <- if (true_CV) new_fit else cur_fit

      for (i_lik in seq_len(n_likelihoods)) {
        if (is_mg[mn] && true_CV) {
          post_lp <- .cv_sample_post_lp_metric_graph(new_fit, i_lik, test_list,
            new_n_samples, print, full_bru_fit = cur_fit
          )
        } else {
          post_lp <- .cv_sample_post_lp_generic(new_fit, i_lik, test_list,
            new_n_samples, print
          )
        }

        Y_samples <- .cv_get_response_samples(post_lp, fit_for_hyper, i_lik,
          new_n_samples, print
        )
        Y_samples <- do.call(rbind, Y_samples)
        post_samples[[model_names[[mn]]]][[fold]][[i_lik]] <- Y_samples

        lhoods_full <- inlabru::as_bru_obs_list(cur_fit)
        test_data <- lhoods_full[[i_lik]]$response_data$BRU_response[test_list[[i_lik]]]
        if (return_true_test_values) {
          true_test_values[[model_names[[mn]]]][[fold]][[i_lik]] <- test_data
        }

        n_test <- length(test_data)
        if (n_test == 0) next

        post_mean <- rowMeans(Y_samples)
        if ("mse" %in% scores) {
          val <- mean((test_data - post_mean)^2, na.rm = TRUE)
          if (orientation_results == "positive") val <- -val
          mse[[i_lik]][fold, mn] <- val
          if (print) cat(sprintf("MSE - Likelihood %d: %g\n", i_lik, val))
        }
        if ("mae" %in% scores) {
          val <- mean(abs(test_data - post_mean), na.rm = TRUE)
          if (orientation_results == "positive") val <- -val
          mae[[i_lik]][fold, mn] <- val
          if (print) cat(sprintf("MAE - Likelihood %d: %g\n", i_lik, val))
        }
        if ("dss" %in% scores) {
          Y1 <- Y_samples[, seq_len(n_samples), drop = FALSE]
          Y2 <- Y_samples[, (n_samples + 1):(2 * n_samples), drop = FALSE]
          post_var <- rowMeans(Y1^2) - (rowMeans(Y1))^2
          val <- mean((test_data - rowMeans(Y2))^2 / post_var + log(post_var))
          dss[[i_lik]][fold, mn] <- val
          if (print) cat(sprintf("DSS - Likelihood %d: %g\n", i_lik, val))
        }

        if (any(c("crps", "scrps", "wcrps", "swcrps") %in% scores)) {
          Y1 <- Y_samples[, seq_len(n_samples), drop = FALSE]
          Y2 <- Y_samples[, (n_samples + 1):(2 * n_samples), drop = FALSE]
        }

        if (any(c("crps", "scrps") %in% scores)) {
          if (parallelize_RP) {
            E1 <- foreach::`%dopar%`(foreach::foreach(i = seq_len(n_test)),
              { mean(abs(Y1[i, ] - test_data[i])) })
            E2 <- foreach::`%dopar%`(foreach::foreach(i = seq_len(n_test)),
              { mean(abs(Y1[i, ] - Y2[i, ])) })
          } else {
            E1 <- lapply(seq_len(n_test), function(i) mean(abs(Y1[i, ] - test_data[i])))
            E2 <- lapply(seq_len(n_test), function(i) mean(abs(Y1[i, ] - Y2[i, ])))
          }
        }

        if (any(c("wcrps", "swcrps") %in% scores)) {
          if (parallelize_RP) {
            E1t <- foreach::`%dopar%`(foreach::foreach(i = seq_len(n_test)), {
              mean(abs((Y1[i, ] > weight_thr) * (Y1[i, ] - weight_thr) -
                (test_data[i] > weight_thr) * (test_data[i] - weight_thr)))
            })
            E2t <- foreach::`%dopar%`(foreach::foreach(i = seq_len(n_test)), {
              mean(abs((Y1[i, ] > weight_thr) * (Y1[i, ] - weight_thr) -
                (Y2[i, ] > weight_thr) * (Y2[i, ] - weight_thr)))
            })
          } else {
            E1t <- lapply(seq_len(n_test), function(i) {
              mean(abs((Y1[i, ] > weight_thr) * (Y1[i, ] - weight_thr) -
                (test_data[i] > weight_thr) * (test_data[i] - weight_thr)))
            })
            E2t <- lapply(seq_len(n_test), function(i) {
              mean(abs((Y1[i, ] > weight_thr) * (Y1[i, ] - weight_thr) -
                (Y2[i, ] > weight_thr) * (Y2[i, ] - weight_thr)))
            })
          }
        }

        if ("crps" %in% scores) {
          tmp <- mean(unlist(lapply(seq_len(n_test),
            function(i) -E1[[i]] + 0.5 * E2[[i]])))
          if (orientation_results == "negative") tmp <- -tmp
          crps[[i_lik]][fold, mn] <- tmp
          if (print) cat(sprintf("CRPS - Likelihood %d: %g\n", i_lik, tmp))
        }
        if ("scrps" %in% scores) {
          tmp <- mean(unlist(lapply(seq_len(n_test),
            function(i) -E1[[i]] / E2[[i]] - 0.5 * log(E2[[i]]))))
          if (orientation_results == "negative") tmp <- -tmp
          scrps[[i_lik]][fold, mn] <- tmp
          if (print) cat(sprintf("SCRPS - Likelihood %d: %g\n", i_lik, tmp))
        }
        if ("wcrps" %in% scores) {
          tmp <- mean(unlist(lapply(seq_len(n_test),
            function(i) 0.5 * E2t[[i]] - E1t[[i]])))
          if (orientation_results == "negative") tmp <- -tmp
          wcrps[[i_lik]][fold, mn] <- tmp
          if (print) cat(sprintf("wCRPS - Likelihood %d: %g\n", i_lik, tmp))
        }
        if ("swcrps" %in% scores) {
          if (any(unlist(E2t) == 0)) {
            warning("swCRPS cannot be computed; lower weight_thr or raise n_samples.")
            swcrps[[i_lik]][fold, mn] <- NA
          } else {
            tmp <- mean(unlist(lapply(seq_len(n_test),
              function(i) -E1t[[i]] / E2t[[i]] - 0.5 * log(E2t[[i]]))))
            if (orientation_results == "negative") tmp <- -tmp
            swcrps[[i_lik]][fold, mn] <- tmp
            if (print) cat(sprintf("swCRPS - Likelihood %d: %g\n", i_lik, tmp))
          }
        }
      }
    }
  }

  result_df <- data.frame(Model = model_names, stringsAsFactors = FALSE)
  add_score_cols <- function(score_list, score_name) {
    if (n_likelihoods > 1) {
      for (i in seq_len(n_likelihoods)) {
        result_df[[paste0(score_name, "_lik", i)]] <<- colMeans(score_list[[i]])
      }
      weighted <- matrix(0,
        nrow = nrow(score_list[[1]]),
        ncol = ncol(score_list[[1]])
      )
      total_w <- matrix(0,
        nrow = nrow(score_list[[1]]),
        ncol = ncol(score_list[[1]])
      )
      for (i in seq_len(n_likelihoods)) {
        w <- vapply(seq_len(nrow(score_list[[i]])), function(f) {
          length(train_test_indexes[[f]][["test"]][[i]])
        }, numeric(1))
        weighted <- weighted + w * score_list[[i]]
        total_w <- total_w + w
      }
      result_df[[paste0(score_name, "_total")]] <<- colMeans(weighted / total_w)
    } else {
      result_df[[score_name]] <<- colMeans(score_list[[1]])
    }
  }

  if ("mse" %in% scores) add_score_cols(mse, "mse")
  if ("mae" %in% scores) add_score_cols(mae, "mae")
  if ("dss" %in% scores) add_score_cols(dss, "dss")
  if ("crps" %in% scores) add_score_cols(crps, "crps")
  if ("scrps" %in% scores) add_score_cols(scrps, "scrps")
  if ("wcrps" %in% scores) add_score_cols(wcrps, "wcrps")
  if ("swcrps" %in% scores) add_score_cols(swcrps, "swcrps")

  settings_list <- NULL
  if (save_settings) {
    settings_list <- list(
      n_samples = n_samples, cv_type = cv_type, true_CV = true_CV,
      orientation_results = orientation_results
    )
    if (cv_type == "k-fold") settings_list$k <- k
    if (cv_type == "lpo") {
      settings_list$percentage <- percentage
      settings_list$number_folds <- number_folds
    }
  }

  if (include_best) {
    score_prefixes <- c("mse", "mae", "dss", "crps", "scrps", "wcrps", "swcrps")
    final_row <- c("Best")
    for (j in 2:ncol(result_df)) {
      colname <- names(result_df)[j]
      is_metric <- any(vapply(score_prefixes, function(x) {
        startsWith(colname, x)
      }, logical(1)))
      if (!is_metric) {
        final_row <- c(final_row, "")
        next
      }
      vals <- result_df[, j]
      if (all(is.na(vals))) {
        final_row <- c(final_row, NA)
        next
      }
      idx <- if (orientation_results == "negative") {
        which.min(vals)
      } else {
        which.max(vals)
      }
      final_row <- c(final_row, model_names[idx])
    }
    result_df <- rbind(result_df, final_row)
    rownames(result_df)[nrow(result_df)] <- ""
  }

  if (!is.null(cluster_tmp)) parallel::stopCluster(cluster_tmp)

  if (return_post_samples) return_scores_folds <- TRUE

  if (!return_scores_folds) {
    if (save_settings) {
      out <- list(scores_df = result_df, settings = settings_list)
      if (return_train_test) out$train_test <- train_test_indexes
      if (return_true_test_values) out$true_test_values <- true_test_values
    } else if (return_train_test) {
      out <- list(scores_df = result_df, train_test = train_test_indexes)
      if (return_true_test_values) out$true_test_values <- true_test_values
    } else if (return_true_test_values) {
      out <- list(scores_df = result_df, true_test_values = true_test_values)
    } else {
      out <- result_df
    }
  } else {
    add_model_names <- function(score_list) {
      lapply(score_list, function(mat) {
        colnames(mat) <- model_names
        mat
      })
    }
    scores_folds <- list()
    if ("dss" %in% scores) scores_folds$dss <- add_model_names(dss)
    if ("mse" %in% scores) scores_folds$mse <- add_model_names(mse)
    if ("mae" %in% scores) scores_folds$mae <- add_model_names(mae)
    if ("crps" %in% scores) scores_folds$crps <- add_model_names(crps)
    if ("scrps" %in% scores) scores_folds$scrps <- add_model_names(scrps)
    if ("wcrps" %in% scores) scores_folds$wcrps <- add_model_names(wcrps)
    if ("swcrps" %in% scores) scores_folds$swcrps <- add_model_names(swcrps)
    out <- list(scores_df = result_df, scores_folds = scores_folds)
    if (save_settings) out$settings <- settings_list
    if (return_train_test) out$train_test <- train_test_indexes
    if (return_true_test_values) out$true_test_values <- true_test_values
    if (return_post_samples) out$post_samples <- post_samples
  }
  out
}
