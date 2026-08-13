# Rebuild data/columbia_main_component.rda from the canonical Figshare data.
# Run from the MetricGraph package root. This script is not run during package
# installation or checks.

pkgload::load_all(".", reset = TRUE, export_all = TRUE)

.columbia_component_fingerprint <- function(component) {
  summary_without_fingerprint <- component$summary
  summary_without_fingerprint$fingerprint <- NULL
  payload <- list(
    weights_name = component$weights_name,
    edge_attributes = sf::st_drop_geometry(component$edges),
    edge_geometry = sf::st_as_binary(sf::st_geometry(component$edges)),
    observations = component$observations,
    summary = summary_without_fingerprint
  )

  if (requireNamespace("digest", quietly = TRUE)) {
    return(paste0(
      "sha256:",
      digest::digest(
        payload,
        algo = "sha256",
        serialize = TRUE,
        serializeVersion = 3L
      )
    ))
  }

  temporary_path <- tempfile(pattern = "columbia-fingerprint-", fileext = ".rds")
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(payload, temporary_path, version = 3L)
  paste0("md5:", unname(tools::md5sum(temporary_path)))
}
columbia_extract_main_component <- function(
    ssn,
    weights_name = "h2oAreaKm2") {
  cat("Preparing the largest Mid-Columbia component...\n")

  if (!inherits(ssn, "SSN") || is.null(ssn$obs) ||
      !inherits(ssn$obs, "data.frame")) {
    stop("ssn must be an SSN object containing an observation data frame.")
  }
  if ("columbia_obs_id" %in% names(ssn$obs)) {
    stop("The source data already contains the reserved columbia_obs_id column.")
  }
  source_observation_count <- nrow(ssn$obs)
  ssn_with_ids <- ssn
  ssn_with_ids$obs$columbia_obs_id <- seq_len(source_observation_count)

  graph_all <- MetricGraph::metric_graph$new(
    ssn_with_ids,
    check_connected = FALSE
  )
  rm(ssn_with_ids)
  graph_all$set_edge_weights(directional_weights = weights_name)

  weights <- graph_all$get_edge_weights(format = "tibble")[[weights_name]]
  if (any(graph_all$E[, 1L] == graph_all$E[, 2L])) {
    stop("The Columbia graph contains loop edges.")
  }

  undirected <- igraph::make_graph(
    edges = c(t(graph_all$E)),
    directed = FALSE,
    n = graph_all$nV
  )
  components <- igraph::components(undirected, mode = "weak")
  edge_component <- components$membership[graph_all$E[, 1L]]
  length_by_component <- tapply(
    graph_all$edge_lengths,
    edge_component,
    sum
  )
  component_id <- as.integer(names(which.max(length_by_component)))
  keep_edges <- edge_component == component_id
  kept_edge_indices <- which(keep_edges)
  component_weights <- weights[keep_edges]
  if (length(weights) != graph_all$nE ||
      anyNA(component_weights) ||
      any(!is.finite(component_weights)) ||
      any(component_weights <= 0)) {
    stop(
      weights_name,
      " must be finite and strictly positive on every retained edge."
    )
  }

  full_to_component <- integer(graph_all$nE)
  full_to_component[kept_edge_indices] <- seq_along(kept_edge_indices)

  component_edges <- graph_all$get_edges(format = "sf")[keep_edges, ]
  if (!all(sf::st_geometry_type(component_edges) == "LINESTRING")) {
    stop("The Columbia component must contain only LINESTRING edges.")
  }

  raw_data <- graph_all$.__enclos_env__$private$data
  full_observation_edges <- raw_data[[".edge_number"]]
  normalized_distance <- raw_data[[".distance_on_edge"]]
  source_observation_ids <- raw_data[["columbia_obs_id"]]
  if (length(source_observation_ids) != length(full_observation_edges) ||
      anyNA(source_observation_ids) ||
      anyDuplicated(source_observation_ids) ||
      any(source_observation_ids < 1L) ||
      any(source_observation_ids > source_observation_count)) {
    stop("Source Columbia observation identifiers were not preserved on import.")
  }
  if (anyNA(full_observation_edges) ||
      any(full_observation_edges < 1L) ||
      any(full_observation_edges > graph_all$nE)) {
    stop("The source Columbia graph has invalid observation edge numbers.")
  }
  if (anyNA(normalized_distance) ||
      any(!is.finite(normalized_distance)) ||
      any(normalized_distance < -1e-12) ||
      any(normalized_distance > 1 + 1e-12)) {
    stop("The source Columbia graph has invalid normalized observation positions.")
  }
  normalized_distance <- pmin(1, pmax(0, normalized_distance))

  keep_observations <- full_to_component[full_observation_edges] > 0L
  data_columns <- grep("^\\.", names(raw_data), invert = TRUE, value = TRUE)
  observations <- as.data.frame(
    lapply(raw_data[data_columns], function(x) x[keep_observations]),
    stringsAsFactors = FALSE
  )
  observations$edge_number <-
    full_to_component[full_observation_edges[keep_observations]]
  observations$distance_on_edge <-
    normalized_distance[keep_observations]

  summary <- list(
    full_vertices = graph_all$nV,
    full_edges = graph_all$nE,
    full_observations = length(full_observation_edges),
    source_observations = source_observation_count,
    post_merge_observations = length(full_observation_edges),
    merged_source_observations =
      source_observation_count - length(full_observation_edges),
    component_id = component_id,
    component_edges = nrow(component_edges),
    component_observations = nrow(observations),
    component_total_length = sum(graph_all$edge_lengths[keep_edges]),
    dropped_edges = sum(!keep_edges),
    dropped_observations = sum(!keep_observations)
  )

  cat(sprintf(
    "  kept %d / %d edges and %d / %d observations\n",
    summary$component_edges,
    summary$full_edges,
    summary$component_observations,
    summary$full_observations
  ))
  if (summary$merged_source_observations > 0L) {
    cat(sprintf(
      "  source rows: %d; unique ungrouped graph locations: %d\n",
      summary$source_observations,
      summary$post_merge_observations
    ))
  }

  rm(graph_all, undirected, components, raw_data)
  invisible(gc(verbose = FALSE))

  component <- structure(
    list(
      edges = component_edges,
      observations = observations,
      summary = summary,
      weights_name = weights_name,
      fingerprint = NULL
    ),
    class = "columbia_main_component"
  )
  component$fingerprint <- .columbia_component_fingerprint(component)
  component$summary$fingerprint <- component$fingerprint
  component
}


archive_url <- "https://ndownloader.figshare.com/files/42339810"
archive_md5 <- "022dd02ba2feea7c384f7deda63b2cbc"
archive_path <- tempfile(fileext = ".zip")
extract_path <- tempfile(pattern = "metricgraph-columbia-")
dir.create(extract_path)
on.exit(unlink(c(archive_path, extract_path), recursive = TRUE), add = TRUE)

utils::download.file(archive_url, archive_path, mode = "wb", quiet = FALSE)
actual_md5 <- unname(tools::md5sum(archive_path))
if (!identical(actual_md5, archive_md5)) {
  stop("Mid-Columbia archive checksum mismatch: ", actual_md5)
}
utils::unzip(archive_path, exdir = extract_path)

ssn_candidates <- list.dirs(extract_path, recursive = TRUE, full.names = TRUE)
ssn_candidates <- ssn_candidates[grepl("\\.ssn$", ssn_candidates)]
if (length(ssn_candidates) != 1L) {
  stop("Expected exactly one .ssn directory; found ", length(ssn_candidates), ".")
}

ssn <- SSN2::ssn_import(ssn_candidates, overwrite = FALSE)
columbia_main_component <- columbia_extract_main_component(
  ssn,
  weights_name = "h2oAreaKm2"
)

expected_summary <- c(
  component_edges = 18668,
  component_observations = 2080,
  source_observations = 9521,
  post_merge_observations = 2758
)
actual_summary <- unlist(columbia_main_component$summary[names(expected_summary)])
if (!identical(unname(actual_summary), unname(expected_summary))) {
  stop(
    "Regenerated component counts differ from the audited data: ",
    paste(names(actual_summary), actual_summary, collapse = ", ")
  )
}

expected_fingerprint <-
  "sha256:eaf68cffcdd7e33dbbbc44744df3054a0fcfef1dac4990021b02fa3f24a6d2f7"
if (!identical(columbia_main_component$fingerprint, expected_fingerprint)) {
  stop(
    "Regenerated component fingerprint differs from the audited data: ",
    columbia_main_component$fingerprint
  )
}

save(
  columbia_main_component,
  file = "data/columbia_main_component.rda",
  compress = "xz",
  version = 3L
)
