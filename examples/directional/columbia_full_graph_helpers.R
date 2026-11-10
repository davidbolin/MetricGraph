# Runtime helpers shared by the Columbia example scripts (crossval_methods_columbia.R,
# sim_speed_columbia.R). Data extraction lives in data-raw/columbia_main_component.R.


## Directional weight rule ---------------------------------------------------

# K2 directional-weight rule: sqrt(w / sum(w)), a shrunk-toward-uniform
# variant of the graph's default K1 rule (w / sum(w)). See the
# directional_rule comment above the methods list in
# crossval_methods_columbia.R for how K1/K2/reversed combine into the
# compared methods.
.columbia_k2_weights <- function(weight) {
  sqrt(weight / sum(weight))
}


## Data validation -----------------------------------------------------------

.columbia_component_fingerprint <- function(component) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required to validate the Columbia data.")
  }
  summary_without_fingerprint <- component$summary
  summary_without_fingerprint$fingerprint <- NULL
  payload <- list(
    weights_name = component$weights_name,
    edge_attributes = sf::st_drop_geometry(component$edges),
    edge_geometry = sf::st_as_binary(sf::st_geometry(component$edges)),
    observations = component$observations,
    summary = summary_without_fingerprint
  )

  paste0(
    "sha256:",
    digest::digest(
      payload,
      algo = "sha256",
      serialize = TRUE,
      # Pin the serialization format explicitly (3 has been R's own default
      # since R 3.6.0; format 2 was default before that -- see ?serialize)
      # so the hash stays reproducible even if a future R release changes
      # that default.
      serializeVersion = 3L
    )
  )
}

.columbia_require_component <- function(component) {
  if (!inherits(component, "columbia_main_component")) {
    stop("component must be the packaged columbia_main_component data set.")
  }
  required <- c(
    "edges", "observations", "summary", "weights_name", "fingerprint"
  )
  missing <- setdiff(required, names(component))
  if (length(missing)) {
    stop("Invalid Columbia component; missing: ", paste(missing, collapse = ", "))
  }
  required_observation_variables <- c(
    "columbia_obs_id", "STREAM_AUG", "ELEV", "SLOPE", "PRECIP",
    "edge_number", "distance_on_edge"
  )
  missing_observation_variables <- setdiff(
    required_observation_variables,
    names(component$observations)
  )
  if (length(missing_observation_variables)) {
    stop(
      "Invalid Columbia observations; missing: ",
      paste(missing_observation_variables, collapse = ", ")
    )
  }
  expected_fingerprint <- .columbia_component_fingerprint(component)
  if (!identical(component$fingerprint, expected_fingerprint) ||
      !identical(component$summary$fingerprint, expected_fingerprint)) {
    stop("The Columbia component fingerprint does not match its contents.")
  }
  invisible(component)
}


## Graph construction --------------------------------------------------------

.columbia_graph_summary <- function(graph, reversed) {
  indegree <- graph$get_degrees("indegree")
  outdegree <- graph$get_degrees("outdegree")

  list(
    reversed = reversed,
    vertices = graph$nV,
    edges = graph$nE,
    observations = nrow(graph$get_PtE()),
    sources = sum(indegree == 0L),
    outlets = sum(outdegree == 0L),
    maximum_indegree = max(indegree),
    maximum_outdegree = max(outdegree),
    replicates = length(graph$get_groups())
  )
}

columbia_make_graph <- function(
    component,
    reversed = FALSE,
    # Callers pass their own weights_name (see crossval_methods_columbia.R's
    # Settings block); it must still equal the packaged value checked below,
    # so a caller variable that has drifted from the packaged data is caught
    # rather than silently used.
    weights_name = component$weights_name) {
  .columbia_require_component(component)
  if (!is.logical(reversed) || length(reversed) != 1L || is.na(reversed)) {
    stop("reversed must be TRUE or FALSE.")
  }
  if (!identical(weights_name, component$weights_name)) {
    stop("weights_name must match the extracted Columbia component.")
  }

  edges <- component$edges
  observations <- component$observations

  if (isTRUE(reversed)) {
    edges <- sf::st_reverse(edges)
    observations$distance_on_edge <- 1 - observations$distance_on_edge
  }

  graph <- MetricGraph::metric_graph$new(
    edges = edges,
    include_edge_weights = TRUE,
    include_obs = FALSE,
    check_connected = TRUE
  )
  graph$set_edge_weights(directional_weights = weights_name)
  graph$add_observations(data = observations, normalized = TRUE)

  # Weights and connectivity must survive reconstruction unchanged.
  graph_weights <- graph$get_edge_weights(format = "tibble")[[weights_name]]
  if (length(graph_weights) != graph$nE ||
      anyNA(graph_weights) ||
      any(!is.finite(graph_weights)) ||
      any(graph_weights <= 0)) {
    stop(weights_name, " changed or became invalid during graph reconstruction.")
  }
  undirected <- igraph::make_graph(
    edges = c(t(graph$E)),
    directed = FALSE,
    n = graph$nV
  )
  if (igraph::components(undirected, mode = "weak")$no != 1L) {
    stop("The reconstructed Columbia component is disconnected.")
  }

  # Edge count must also be unchanged: a duplicated edge between two
  # already-connected vertices would leave the graph connected (so the
  # check above wouldn't catch it) while still breaking the tree property.
  if (graph$nE != nrow(component$edges) || graph$nE != graph$nV - 1L) {
    stop("The extracted Columbia component must be a connected tree.")
  }
  if (nrow(graph$get_PtE()) != nrow(component$observations)) {
    stop("Columbia observation count changed during graph reconstruction.")
  }

  # Observation identity must survive reconstruction.
  data <- graph$get_data(format = "tibble", drop_na = FALSE)
  if (!("columbia_obs_id" %in% names(data)) ||
      anyNA(data$columbia_obs_id) ||
      anyDuplicated(data$columbia_obs_id) ||
      !setequal(data$columbia_obs_id, component$observations$columbia_obs_id)) {
    stop("Stable Columbia observation identifiers were not preserved.")
  }

  # The graph must have the expected up/downstream tree shape.
  summary <- .columbia_graph_summary(graph, reversed)
  if (summary$replicates != 1L) {
    stop("The Columbia component must contain exactly one replicate.")
  }
  if (!isTRUE(reversed)) {
    if (summary$outlets != 1L || summary$maximum_outdegree > 1L) {
      stop("The original Columbia component must be a one-outlet river tree.")
    }
  } else {
    if (summary$sources != 1L || summary$maximum_indegree > 1L) {
      stop("The reversed Columbia component must be a one-source tree.")
    }
  }

  graph
}


## Checkpoint saving ---------------------------------------------------------

# Write to a temp file in the same directory, then rename into place: a kill
# mid-write can't corrupt or truncate a previously-good checkpoint at path.
# file.copy() is the fallback for a rename across filesystems, where
# file.rename() can fail.
columbia_save_checkpoint <- function(result, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "-"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(result, temporary_path)

  if (!file.rename(temporary_path, path)) {
    if (!file.copy(temporary_path, path, overwrite = TRUE)) {
      stop("Could not move the completed checkpoint into place: ", path)
    }
  }
  invisible(path)
}
