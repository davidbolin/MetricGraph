library(osmdata)
call <- opq(bbox = c("Stockholm")) 
call <- add_osm_feature(call, key = "highway",value=c("motorway",
"primary","secondary"))
mydata <- osmdata_sf(call)
mydata