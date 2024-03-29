rm(list=ls())
library(MetricGraph)
library(ggplot2)
library(hexSticker)
library(viridis)
graph <- metric_graph$new(lines = logo_lines(), merge_intersections = "all_intersections",
                          tolerance_intersections = 1e-3,
                          tolerance_overlapping = 1e-1, tolerance = 1e-3)
graph$plot(degree = TRUE)
graph$build_mesh(h=0.01)
if(1){
  package_brown = rgb(104/255,92/255,54/255)
  package_dark = rgb(77/255,76/255,139/255)
  package_light = rgb(1,1,1)
  package_bg = rgb(0.5*77/255,0.5*76/255,0.5*139/255)
  graph$build_mesh(h=0.01)
  set.seed(10)
  u <- sample_spde(kappa = 0.2, sigma = 1, alpha=2, type = "mesh", graph = graph)
  p <- graph$plot_function(u, vertex_size = 0.2, vertex_color = package_brown,
                           line_width = 0.75)
  p

  p <- p + theme_void() + theme_transparent() +
    scale_colour_gradient(low = package_dark,
                          high = package_light, guide = "none")
  #print(p)
  sticker(p, package=NULL, p_size=20, s_x=1, s_y=1, s_width=1.3, s_height=3,
               h_fill = package_bg,
               h_color = package_brown,
               filename="metricgraph.png")
  print(p)

}
