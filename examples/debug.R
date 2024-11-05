graph$build_mesh(h = 0.1)
get_edge_cov <- graph$get_edge_weights()
df_pred <- data.frame(edge_number = graph$mesh$PtE[,1],
                      distance_on_edge = graph$mesh$PtE[,2],
                      SLOPE = c(get_edge_cov[graph$mesh$PtE[,1],"SLOPE"]$SLOPE),
                      netID = factor(get_edge_cov[graph$mesh$PtE[,1],"netID"]$netID))
pred_alpha1 <- augment(model.wm1 , newdata = df_pred,
                       normalized = TRUE)
pred_alpha2 <- augment(res.wm1.dir , newdata = df_pred,
                       normalized = TRUE)
pred_alpha1$.fitted = pred_alpha1$.fitted - pred_alpha2$.fitted
p <- graph$plot_function(
  newdata = pred_alpha1,
  data = ".fitted",
  improve_plot = TRUE,
  vertex_size = 0,
  edge_width = 0.5,
  mapview_caption = "Pred alpha1"
)
