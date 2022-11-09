
library(osmdata)
set_overpass_url("https://maps.mail.ru/osm/tools/overpass/api/interpreter")
call <- opq(bbox = c(12.4,55.5,12.8,55.9))
call <- add_osm_feature(call, key = "highway",value=c("motorway",
                                                        "primary","secondary"))
data_sp <- osmdata_sp(call)

city_lines <- SpatialLines(mydata$osm_lines@lines)
#tmp2 <- line_to_vertex(tmp)
graph <- metric_graph$new(city_lines)
#graph$plot()
graph$build_mesh(h = 0.01)
u <- sample_spde_mesh(kappa = 10, sigma = 2, alpha = 1, graph = graph)
graph$plot_function(u)
