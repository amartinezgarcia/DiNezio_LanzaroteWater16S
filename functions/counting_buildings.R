
latitudes <- stations$

count_buildings_osm <- function(latitudes, longitudes, radius = 1000) {
  
  library(osmdata)
  library(sf)
  library(dplyr)
  
  
  # Check inputs
  if (length(latitudes) != length(longitudes)) {
    stop("latitudes and longitudes must have the same length")
  }
  
  # Initialize results
  results <- data.frame(
    latitude = latitudes,
    longitude = longitudes,
    n_buildings = NA_integer_
  )
  
  # Loop over each coordinate pair
  for (i in seq_along(latitudes)) {
    cat("Processing point", i, "of", length(latitudes), "...\n")
    
    # Create point and buffer
    point <- st_sfc(st_point(c(longitudes[i], latitudes[i])), crs = 4326)
    point_proj <- st_transform(point, 3857)  # project to meters
    buffer <- st_buffer(point_proj, radius)
    buffer_wgs84 <- st_transform(buffer, 4326)
    
    # Get bounding box
    bbox <- st_bbox(buffer_wgs84)
    
    # Query OSM for buildings
    query <- opq(bbox = bbox) %>%
      add_osm_feature(key = "building")
    
    # Try fetching OSM data
    buildings <- tryCatch({
      osmdata_sf(query)$osm_polygons
    }, error = function(e) {
      warning(paste("Failed to fetch OSM data for point", i))
      return(NULL)
    })
    
    # Count buildings
    if (!is.null(buildings) && nrow(buildings) > 0) {
      buildings <- st_transform(buildings, st_crs(buffer_wgs84))
      in_buffer <- buildings[st_intersects(buildings, buffer_wgs84, sparse = FALSE), ]
      results$n_buildings[i] <- nrow(in_buffer)
    } else {
      results$n_buildings[i] <- 0
    }
  }
  
  return(results)
}