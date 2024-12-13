 
 library(sf)
library(leaflet)

# Step 1: Retrieve edges in the original CRS and transform to WGS84
edges_sf <- graph$get_edges(format = "sf")
edges_sf_wgs84 <- st_transform(edges_sf, crs = 4326)

# Step 2: Calculate the bounding box using st_coordinates to ensure accuracy
coords <- st_coordinates(edges_sf_wgs84)
lng1 <- -115.43
lng2 <- -114.98
lat1 <- 44.28
lat2 <- 44.5

# Step 3: Create mapview object and set boundaries
m <- graph$plot(data = "Temperature", vertex_size = 0, type = "mapview", data_size = 2)
leaflet_map <- m@map

# Clear controls, add legend, and apply fitBounds using calculated boundaries
leaflet_map <- leaflet_map %>%
  clearControls() %>%
  addLegend(
    position = "bottomright",
    pal = colorNumeric("viridis", graph$get_data()$Temperature),
    values = graph$get_data()[["Temperature"]],
    title = "Temperature"
  ) %>%
  fitBounds(
    lng1 = lng1,
    lat1 = lat1,
    lng2 = lng2,
    lat2 = lat2
  )

# Display the customized map or save it with mapshot2
mapshot2(leaflet_map, file = "test.png", vwidth = 700, vheight = 500)
