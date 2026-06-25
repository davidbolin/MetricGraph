# test_that("Inlabru predict method works for alpha = 1", {
#     skip_if_not_installed("INLA")
#     skip_if_not_installed("inlabru")

#     set.seed(1)
#     library(INLA)
#     library(inlabru)

#     # Create a simple graph with two edges
#     edge <- rbind(c(0, 0), c(1, 0))
#     edges <- list(edge)
#     graph <- metric_graph$new(edges = edges)

#     # Add observations at vertices only (2 points)
#     # Use fixed values instead of sampling
#     y <- c(2.0, 4.0) # Observations at the two vertices
#     obs_loc <- rbind(c(1, 0), c(1, 1)) # Edge 1, positions 0 and 1

#     # Add observations to graph
#     df_graph <- data.frame(
#         y = y, edge_number = obs_loc[, 1],
#         distance_on_edge = obs_loc[, 2]
#     )
#     graph$add_observations(data = df_graph, normalized = TRUE)

#     # Create SPDE model
#     spde_model <- graph_spde(graph, alpha = 1)

#     # Create component
#     cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_model)

#     # Create data for fitting
#     data_spde <- graph_data_spde(spde_model, loc_name = "loc")

#     # Fit the model
#     spde_fit <- bru(cmp,
#         data = data_spde[["data"]],
#         options = list(
#             num.threads = "1:1",
#             verbose = FALSE,
#             control.inla = list(int.strategy = "eb")
#         ),
#         allow_combine = FALSE
#     )

#     # Predict at midpoint 
#     pred_loc <- list(loc = data.frame(.edge_number = 1, .distance_on_edge = 0.5), Intercept = 1) 

#     field_pred <- predict(spde_model, cmp, spde_fit,
#         newdata = pred_loc,
#         formula = ~ Intercept + field
#     )

#     # Get predicted mean
#     pred_mean <- field_pred$pred$mean

#     # Get average of observations
#     obs_mean <- mean(y)

#     # Check that prediction is close to average (with reasonable tolerance)
#     expect_equal(pred_mean, obs_mean, tolerance = 0.5)
# })


# Helper: build a small graph with two edges and `n` grouped observations,
# fit an SPDE model with inlabru, and predict the field at a single location.
# `group` is the vector of grouping variable names passed to add_observations()
# (NULL for no grouping). Returns the predicted field mean.
.fit_and_predict_grouped <- function(group, alpha = 1, n = 16) {
    set.seed(3)
    edge1 <- rbind(c(0, 0), c(1, 0))
    edge2 <- rbind(c(1, 0), c(1, 1))
    graph <- metric_graph$new(edges = list(edge1, edge2), verbose = 0)
    df <- data.frame(
        y = rnorm(n),
        edge_number = rep(c(1, 2), each = n / 2),
        distance_on_edge = rep(c(0.1, 0.3, 0.6, 0.9), times = n / 4),
        g1 = rep(c("a", "b"), each = n / 2),
        g2 = rep(c("x", "y"), times = n / 2)
    )
    graph$add_observations(
        data = df, normalized = TRUE, group = group,
        verbose = 0, suppress_warnings = TRUE
    )
    spde_model <- graph_spde(graph, alpha = alpha)
    data_spde <- graph_data_spde(spde_model, loc_name = "loc")
    cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_model)
    # Specify is_rowwise explicitly so inlabru does not emit a "guessing
    # is_rowwise" warning for the list-like graph data.
    fit <- bru(cmp,
        data = data_spde[["data"]],
        is_rowwise = TRUE,
        options = list(
            num.threads = "1:1", verbose = FALSE,
            control.inla = list(int.strategy = "eb")
        )
    )
    pred_loc <- list(
        loc = data.frame(.edge_number = 1, .distance_on_edge = 0.5),
        Intercept = 1
    )
    field_pred <- predict(spde_model, cmp, fit,
        newdata = pred_loc, formula = ~ Intercept + field
    )
    field_pred$pred$mean
}

test_that("add_observations stores one attribute name per group variable", {
    edge1 <- rbind(c(0, 0), c(1, 0))
    edge2 <- rbind(c(1, 0), c(1, 1))
    graph <- metric_graph$new(edges = list(edge1, edge2), verbose = 0)
    df <- data.frame(
        y = 1:8,
        edge_number = rep(c(1, 2), each = 4),
        distance_on_edge = rep(c(0.2, 0.8), times = 4),
        g1 = rep(c("a", "b"), each = 4),
        g2 = rep(c("x", "y"), times = 4)
    )

    # No grouping -> sentinel ".none"
    graph$add_observations(data = df, normalized = TRUE, verbose = 0,
                           suppress_warnings = TRUE)
    gv <- attr(graph$.__enclos_env__$private$data, "group_variables")
    expect_identical(gv, ".none")

    # Two grouping variables -> length-2 attribute (the case that used to break
    # predict's `if (group_variables == ".none")` check).
    graph$clear_observations()
    graph$add_observations(data = df, normalized = TRUE, group = c("g1", "g2"),
                           verbose = 0, suppress_warnings = TRUE)
    gv2 <- attr(graph$.__enclos_env__$private$data, "group_variables")
    expect_length(gv2, 2)
    expect_setequal(gv2, c("g1", "g2"))
})

test_that("Inlabru predict works without grouping", {
    skip_if_not_installed("INLA")
    skip_if_not_installed("inlabru")
    library(INLA)
    library(inlabru)
    # No grouping -> no "first replicate" warning should be emitted.
    expect_no_warning(m <- .fit_and_predict_grouped(group = NULL))
    expect_true(is.finite(m))
})

test_that("Inlabru predict works with a single group variable", {
    skip_if_not_installed("INLA")
    skip_if_not_installed("inlabru")
    library(INLA)
    library(inlabru)
    # newdata has no replicate column, so predict falls back to the first
    # replicate and warns about it (intended behaviour).
    expect_warning(
        m <- .fit_and_predict_grouped(group = "g1"),
        "first replicate"
    )
    expect_true(is.finite(m))
})

test_that("Inlabru predict works with multiple group variables (regression: condition has length > 1)", {
    skip_if_not_installed("INLA")
    skip_if_not_installed("inlabru")
    library(INLA)
    library(inlabru)
    # Before the fix this raised:
    #   Error in if (group_variables == ".none") : the condition has length > 1
    # Now it predicts the first replicate and warns about it instead.
    expect_warning(
        m <- .fit_and_predict_grouped(group = c("g1", "g2")),
        "first replicate"
    )
    expect_true(is.finite(m))
})


# test_that("Inlabru predict method works for alpha = 2", {
#     skip_if_not_installed("INLA")
#     skip_if_not_installed("inlabru")

#     set.seed(2)
#     library(INLA)
#     library(inlabru)

#     # Create a simple graph with two edges
#     edge1 <- rbind(c(0, 0), c(1, 0))
#     edges <- list(edge1)
#     graph <- metric_graph$new(edges = edges)

#     # Add observations at vertices only (2 points)
#     # Use fixed values instead of sampling
#     y <- c(10.5, 20.5) # Observations at the two vertices
#     obs_loc <- rbind(c(1, 0), c(1, 1)) # Edge 1, positions 0 and 1

#     # Add observations to graph
#     df_graph <- data.frame(
#         y = y, edge_number = obs_loc[, 1],
#         distance_on_edge = obs_loc[, 2]
#     )
#     graph$add_observations(data = df_graph, normalized = TRUE)

#     # Create SPDE model with alpha = 2
#     spde_model <- graph_spde(graph, alpha = 2)

#     # Create component
#     cmp <- y ~ -1 + field(loc, model = spde_model)

#     # Create data for fitting
#     data_spde <- graph_data_spde(spde_model, loc_name = "loc")

#     # Fit the model
#     spde_fit <- bru(cmp,
#         data = data_spde[["data"]],
#         options = list(
#             num.threads = "1:1",
#             verbose = FALSE
#         ),
#         allow_combine = FALSE
#     )

#     # Predict at midpoint - use matrix format
#     pred_loc <- list(loc = data.frame(.edge_number = 1, .distance_on_edge = 0.5)) 
    
#     field_pred <- predict(spde_model, cmp, spde_fit,
#         newdata = pred_loc,
#         formula = ~ field
#     )

#     # Get predicted mean
#     pred_mean <- field_pred$pred$mean

#     # Get average of observations
#     obs_mean <- mean(y)

#     # Check that prediction is close to average (with reasonable tolerance)
#     expect_equal(pred_mean, obs_mean, tolerance = 5.0)
# })
