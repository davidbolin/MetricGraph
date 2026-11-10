# Helpers that make the INLA and inlabru tests safe, following rSPDE.

# Skip the test unless INLA is available and can be used safely. INLA is set
# to one thread while the calling test runs, and its options are restored when
# the test ends. This uses rSPDE's checks when the installed rSPDE has them.
local_metricgraph_safe_inla <- function(envir = parent.frame()) {
  if ("local_rspde_safe_inla" %in% getNamespaceExports("rSPDE")) {
    rSPDE::local_rspde_safe_inla(envir = envir)
  } else {
    testthat::skip_if_not_installed("INLA")
  }
  invisible(TRUE)
}

# Path to the INLA executable, or NULL if it cannot be found. INLA versions
# before 26.11 store it in the 'inla.call' option. Newer versions leave that
# option unset and locate the executable with inla.binary.path().
metricgraph_inla_call <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    return(NULL)
  }
  inla_call <- tryCatch(
    INLA::inla.getOption("inla.call"),
    error = function(e) NULL
  )
  is_path <- is.character(inla_call) && length(inla_call) >= 1L &&
    all(nzchar(inla_call)) &&
    !any(tolower(inla_call) %in% c("remote", "inla.remote"))
  if (is_path) {
    return(inla_call)
  }
  inla_ns <- asNamespace("INLA")
  if (!exists("inla.binary.path", envir = inla_ns, inherits = FALSE)) {
    return(NULL)
  }
  binary_path <- get("inla.binary.path", envir = inla_ns, inherits = FALSE)
  inla_call <- tryCatch(binary_path(check = TRUE), error = function(e) NULL)
  if (is.character(inla_call) && length(inla_call) == 1L &&
      nzchar(inla_call)) {
    return(inla_call)
  }
  NULL
}

# inlabru::bru() stops when inlabru cannot use the installed INLA. In
# particular, inlabru <= 2.15.0 expects INLA's 'inla.call' option to hold the
# path to the INLA executable, but INLA >= 26.11 leaves that option unset. In
# that case, point 'inla.call' to the INLA executable until the calling test
# ends. Skip the test if inlabru still cannot use INLA.
skip_if_inlabru_cannot_use_inla <- function(envir = parent.frame()) {
  testthat::skip_if_not_installed("inlabru")
  inlabru_can_use_inla <- function() {
    tryCatch(
      isTRUE(inlabru::bru_safe_inla(multicore = TRUE, quietly = TRUE)),
      error = function(e) FALSE
    )
  }
  if (!inlabru_can_use_inla() &&
      requireNamespace("INLA", quietly = TRUE) &&
      is.null(INLA::inla.getOption("inla.call"))) {
    inla_call <- metricgraph_inla_call()
    if (!is.null(inla_call)) {
      INLA::inla.setOption(inla.call = inla_call)
      withr::defer(INLA::inla.setOption(inla.call = NULL), envir = envir)
    }
  }
  testthat::skip_if_not(
    inlabru_can_use_inla(),
    "inlabru cannot use the installed INLA."
  )
}

# Skip unless both INLA and inlabru can be used safely.
skip_if_no_safe_inla <- function(envir = parent.frame()) {
  local_metricgraph_safe_inla(envir = envir)
  skip_if_inlabru_cannot_use_inla(envir = envir)
}
