

a = Sys.time()
ecoregion_control = lapply(seq_along(control$geometry), function(i) {
  pt = control$geometry[i]
  coord = sf::st_coordinates(pt)
  WorE = ifelse(coord[1] > 0, "E", "W")
  NorS = ifelse(coord[2] > 0, "N", "S")
  lng_orig = abs(floor(coord[1] / 10) * 10)
  lat_orig = abs(ceiling(coord[2] / 10) * 10)

  ecoregion_filename = paste0("ecoregion_", NorS, lat_orig, "_", WorE, lng_orig)
  ecoregion_rast = terra::rast(paste0("/maps/4C/ecoregions/", ecoregion_filename, ".tif")) #each file indicates its top-left corner

  sf::st_crs(pt) = sf::st_crs(ecoregion_rast)
  val = terra::extract(ecoregion_rast, coord) %>% as.vector()

  return(val)
}) %>%
  unlist()
b = Sys.time()
b - a #30 seconds per pair, too slow