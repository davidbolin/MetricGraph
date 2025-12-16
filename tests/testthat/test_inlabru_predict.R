test_that("Inlabru predict method works for alpha = 1", {
    skip_if_not_installed("INLA")
    skip_if_not_installed("inlabru")

    set.seed(1)
    library(INLA)
    library(inlabru)

    # Create a simple graph with two edges
    edge <- rbind(c(0, 0), c(1, 0))
    edges <- list(edge)
    graph <- metric_graph$new(edges = edges)

    # Add observations at vertices only (2 points)
    # Use fixed values instead of sampling
    y <- c(2.0, 4.0) # Observations at the two vertices
    obs_loc <- rbind(c(1, 0), c(1, 1)) # Edge 1, positions 0 and 1

    # Add observations to graph
    df_graph <- data.frame(
        y = y, edge_number = obs_loc[, 1],
        distance_on_edge = obs_loc[, 2]
    )
    graph$add_observations(data = df_graph, normalized = TRUE)

    # Create SPDE model
    spde_model <- graph_spde(graph, alpha = 1)

    # Create component
    cmp <- y ~ -1 + Intercept(1) + field(loc, model = spde_model)

    # Create data for fitting
    data_spde <- graph_data_spde(spde_model, loc_name = "loc")

    # Fit the model
    spde_fit <- bru(cmp,
        data = data_spde[["data"]],
        options = list(
            num.threads = "1:1",
            verbose = FALSE,
            control.inla = list(int.strategy = "eb")
        ),
        allow_combine = FALSE
    )

    # Predict at midpoint 
    pred_loc <- list(loc = data.frame(.edge_number = 1, .distance_on_edge = 0.5), Intercept = 1) 

    field_pred <- predict(spde_model, cmp, spde_fit,
        newdata = pred_loc,
        formula = ~ Intercept + field
    )

    # Get predicted mean
    pred_mean <- field_pred$pred$mean

    # Get average of observations
    obs_mean <- mean(y)

    # Check that prediction is close to average (with reasonable tolerance)
    expect_equal(pred_mean, obs_mean, tolerance = 0.5)
})


test_that("Inlabru predict method works for alpha = 2", {
    skip_if_not_installed("INLA")
    skip_if_not_installed("inlabru")

    set.seed(2)
    library(INLA)
    library(inlabru)

    # Create a simple graph with two edges
    edge1 <- rbind(c(0, 0), c(1, 0))
    edges <- list(edge1)
    graph <- metric_graph$new(edges = edges)

    # Add observations at vertices only (2 points)
    # Use fixed values instead of sampling
    y <- c(10.5, 20.5) # Observations at the two vertices
    obs_loc <- rbind(c(1, 0), c(1, 1)) # Edge 1, positions 0 and 1

    # Add observations to graph
    df_graph <- data.frame(
        y = y, edge_number = obs_loc[, 1],
        distance_on_edge = obs_loc[, 2]
    )
    graph$add_observations(data = df_graph, normalized = TRUE)

    # Create SPDE model with alpha = 2
    spde_model <- graph_spde(graph, alpha = 2)

    # Create component
    cmp <- y ~ -1 + field(loc, model = spde_model)

    # Create data for fitting
    data_spde <- graph_data_spde(spde_model, loc_name = "loc")

    # Fit the model
    spde_fit <- bru(cmp,
        data = data_spde[["data"]],
        options = list(
            num.threads = "1:1",
            verbose = FALSE
        ),
        allow_combine = FALSE
    )

    # Predict at midpoint - use matrix format
    pred_loc <- list(loc = data.frame(.edge_number = 1, .distance_on_edge = 0.5)) 
    
    field_pred <- predict(spde_model, cmp, spde_fit,
        newdata = pred_loc,
        formula = ~ field
    )

    # Get predicted mean
    pred_mean <- field_pred$pred$mean

    # Get average of observations
    obs_mean <- mean(y)

    # Check that prediction is close to average (with reasonable tolerance)
    expect_equal(pred_mean, obs_mean, tolerance = 5.0)
})
