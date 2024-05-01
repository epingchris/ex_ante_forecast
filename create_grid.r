rm(list = ls())

library(sf)
library(readr)

# Function to create a grid of polygons within a given bounding box
create_grid = function(bbox, width) {
  xmin = bbox["xmin"]
  ymin = bbox["ymin"]
  xmax = bbox["xmax"]
  ymax = bbox["ymax"]

  rows = floor((ymax - ymin) / width)
  cols = floor((xmax - xmin) / width)
  polygons = vector("list", (rows + 1) * (cols + 1))

  for (i in 0:rows) {
    for (j in 0:cols) {
      poly = st_polygon(list(rbind(c(xmin + j * width, ymin + i * width),
                                   c(xmin + (j + 1) * width, ymin + i * width),
                                   c(xmin + (j + 1) * width, ymin + (i + 1) * width),
                                   c(xmin + j * width, ymin + (i + 1) * width),
                                   c(xmin + j * width, ymin + i * width))))
      polygons[[j + 1 + i * (cols + 1)]] = poly
    }
  }

  grid = st_sfc(polygons)
  st_crs(grid) = st_crs(3857)
  return(grid)
}

# Set project
proj = "1201"
dir_path = paste0("/maps/epr26/tmf-data-grid/", proj)
dir.create(dir_path)

# Define grid width (e.g., 5 km, 25 km2)
grid_width = 5000  # in meters

# Load polygon
shp = st_read(paste0("/maps/epr26/tmf-data/projects/", proj, ".geojson"))

# Reproject to a projected CRS (e.g., EPSG:3857 for Web Mercator)
shp = st_transform(shp, crs = st_crs(3857))

# Create a grid within the extent of the polygon
grid = create_grid(st_bbox(shp), grid_width)

# Intersect the grid with the polygon
intersected_grid = st_intersection(grid, shp)
write_rds(x = intersected_grid, file = paste0(dir_path, "/", proj, "_intersected_grid.rds"))

# Split multipart polygons
splitted_grid = vector("list", sum(sapply(intersected_grid, function(x) length(unique(st_coordinates(x)[, "L2"])))))
ind = 1
for (i in seq_along(intersected_grid)) {
  coords = st_coordinates(intersected_grid[i])
  n_parts = length(unique(coords[, "L2"]))
  if (n_parts > 1) {
  #alternatively: if (class(intersected_grid[[i]])[2] == "MULTIPOLYGON") {
    for(j in seq_len(n_parts)) {
      part_coords = coords[coords[, "L2"] == j, c("X", "Y")]
      part_geometry = st_polygon(list(part_coords)) #create an sfg object from the coordinates
      splitted_grid[[ind]] = part_geometry #store it in the list
      ind = ind + 1
    }
  } else {
    splitted_grid[[ind]] = intersected_grid[[i]] #get the unchanged sfg and store as an item in the list
    ind = ind + 1
  }
}
splitted_grid = st_sfc(splitted_grid) #convert the list to an sfc object
st_crs(splitted_grid) = st_crs(intersected_grid)
write_rds(x = splitted_grid, file = paste0(dir_path, "/", proj, "_splitted_grid.rds"))

# Save each splitted grid polygon as a separate GeoJSON file before merging
for (i in seq_along(splitted_grid)) {
  cell = splitted_grid[i]
  st_write(cell, paste0(dir_path, "/splitted_", proj, "_", i, ".geojson"), driver = "GeoJSON")
}

# Find grids smaller than a certain area
#(4 km2 = 400 ha: this means a minimum of 100 points in samples K to ensure representativeness)
remerged_grid = splitted_grid
small_i = which(as.numeric(st_area(remerged_grid)) < 4e6)
union_dst = rep(NA, length(remerged_grid))
for (i in seq_along(remerged_grid)) {
  if (i %in% small_i) {
      #find all polygons that are not iself ("F") and with touching sides ("1")
      touching_i = st_relate(remerged_grid, pattern = "F***1****")[i][[1]]
      union_i = sample(touching_i %>% as.character(), size = 1) %>% as.numeric() #randomly pick one
      union_dst[i] = union_i
      remerged_grid[union_i] = st_union(remerged_grid[i], remerged_grid[union_i])
  }
}

# Re-number
new_id = rep(NA, length(remerged_grid))
no_union = which(is.na(union_dst))
new_id[no_union] = seq_along(no_union)

grid_info_df = data.frame(old_id = seq_along(remerged_grid),
                          small = seq_along(remerged_grid) %in% small_i,
                          union_to = union_dst,
                          new_id = new_id)
remerged_grid = remerged_grid[-small_i]

# Add area attribute
area_ha = st_area(remerged_grid) / 10000
grid_info_df$area = NA
grid_info_df$area[no_union] = area_ha
write_rds(x = grid_info_df, file = paste0(dir_path, "/", proj, "_grid_info_df.rds"))

# Reproject to geodetic CRS (WGS 84)
remerged_grid = st_transform(remerged_grid, crs = st_crs(4326))
write_rds(x = remerged_grid, file = paste0(dir_path, "/", proj, "_remerged_grid.rds"))

# Save each grid cell as a separate GeoJSON file
for (i in 1:length(remerged_grid)) {
  cell = remerged_grid[i]
  st_write(cell, paste0(dir_path, "/", proj, "_", i, ".geojson"), driver = "GeoJSON")
}


# Replicate previous merges from intersected_grid ----
grid_info_df = read_rds(paste0(dir_path, "/", proj, "_grid_info_df.rds"))

for (i in seq_along(intersected_grid)) {
  if (grid_info_df$small[i]) {
      union_i = grid_info_df$union_to[i]
      intersected_grid[union_i] = st_union(intersected_grid[i], intersected_grid[union_i])
  }
}
intersected_grid = intersected_grid[-which(grid_info_df$small)]