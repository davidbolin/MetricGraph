# Example demonstrating precompute_lgcp_graph functionality
# This example shows how to use precompute_lgcp_graph for efficient LGCP model fitting

library(MetricGraph)

# This is a demonstration of the API - actual usage would require real data
cat("Demonstrating precompute_lgcp_graph API...\n")

# Example usage pattern:
# 
# # Step 1: Precompute expensive quantities once
# precomputed_data <- precompute_lgcp_graph(
#   graph = my_graph,
#   covariates = c("covariate1", "covariate2", "covariate3"),
#   spde_model = my_spde_model,  # optional
#   use_current_mesh = TRUE
# )
# 
# # For maximum performance (modifies the original graph)
# precomputed_fast <- precompute_lgcp_graph(
#   graph = my_graph,
#   covariates = c("covariate1", "covariate2", "covariate3"),
#   clone_graph = FALSE  # Works directly on graph (faster)
# )
# 
# # Step 2: Fit multiple models efficiently using precomputed data
# fit1 <- lgcp_graph(
#   y ~ covariate1,
#   graph = my_graph,
#   precomputed_data = precomputed_data
# )
# 
# fit2 <- lgcp_graph(
#   y ~ covariate1 + covariate2,
#   graph = my_graph,
#   precomputed_data = precomputed_data
# )
# 
# fit3 <- lgcp_graph(
#   y ~ covariate3,
#   graph = my_graph,
#   precomputed_data = precomputed_data
# )
# 
# # For maximum performance (when preserving graph state is not needed)
# fit_fast <- lgcp_graph(
#   y ~ covariate1,
#   graph = my_graph,
#   clone_graph = FALSE  # Works directly on graph (faster)
# )

cat("Benefits of using precompute_lgcp_graph:\n")
cat("- Expensive integration point creation done once\n")
cat("- Graph processing and SPDE setup done once\n") 
cat("- Multiple model fits with different formulas are much faster\n")
cat("- Only covariates included in precomputation can be used\n")
cat("- Backwards compatible - existing code works unchanged\n")
cat("- Use clone_graph = FALSE for maximum performance when graph preservation not needed\n") 