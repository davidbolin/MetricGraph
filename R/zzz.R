# from: https://github.com/tidyverse/hms/blob/master/R/zzz.R
# Thu Apr 19 10:53:24 CEST 2018
register_s3_method <- function(pkg, generic, class, fun = NULL) {
  stopifnot(is.character(pkg), length(pkg) == 1)
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  if (is.null(fun)) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
  } else {
    stopifnot(is.function(fun))
  }

  if (pkg %in% loadedNamespaces()) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }

  # Always register hook in case package is later unloaded & reloaded
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  )
}

#' @importFrom utils installed.packages
#' @importFrom utils packageVersion
register_all_s3_methods = function() {
  register_s3_method("base", "summary", "metric_graph_spde_result")
  
  inlabru_installed <- "inlabru" %in% rownames(installed.packages())
  if(inlabru_installed){
    ## Delayed registration of these are handled
    ## by S3Method NAMESPACE directives instead
    register_s3_method("inlabru", "bru_mapper", "inla_metric_graph_spde")
    register_s3_method("inlabru", "ibm_n", "bru_mapper_inla_metric_graph_spde")
    register_s3_method("inlabru", "ibm_values", "bru_mapper_inla_metric_graph_spde")
    register_s3_method("inlabru", "ibm_jacobian", "bru_mapper_inla_metric_graph_spde")

    ## After inlabru > 2.5.3 is released on cran, change this to an
    ## S3Method directive in inlabru_rspde.R instead.
    inlabru_version <- utils::packageVersion("inlabru")
    if(inlabru_version >= "2.5.3.9002"){
      register_s3_method("inlabru", "bru_get_mapper", "inla_metric_graph_spde")
    }
  }
}

.onLoad = function(libname, pkgname) {
  register_all_s3_methods()
}