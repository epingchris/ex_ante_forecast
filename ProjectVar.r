ProjectVar = function(proj_id) { #extract project-level variables
  #get t0 and country
  if(analysis_type == "old_source") {
    cat("Use new results from epr26 instead.\n")
    return(NULL)
  } else if(analysis_type == "epr26") {
    proj_id_bare = proj_id %>% str_replace("a", "")
    myproj = proj_meta %>% filter(ID == proj_id_bare)
    t0 = myproj$t0
    country = myproj$COUNTRY
  } else if(analysis_type == "control") {
    t0 = "2011"
    country = "Brazil"
  }

  #get area_ha
#   if(analysis_type == "old_source") {
#     #previous: from supplier path (NEED TO CONSISTENTLY SET AOI NAME)
#     site = paste("VCS", proj_id_bare, sep = "_")
#     if(!str_detect(supplier_path, ".shp")) {
#       aoi_path = file.path(supplier_path, site, "GIS", "aoi.shp")
#       if(!file.exists(aoi_path)) aoi_path = file.path(supplier_path, site, paste0(site, ".shp"))
#     }
#     aoi_project = read_sf(aoi_path)
  if(analysis_type == "epr26") {
    aoi_project = st_read(paste0("/maps/epr26/tmf-data/projects/", proj_id, ".geojson"))
  } else if(analysis_type == "control") {
    aoi_project = st_read(paste0("/maps/epr26/tmf-data-grid/0000/", proj_id, ".geojson"))
  }
  area_ha = aoi_project %>%
    st_make_valid() %>%
    st_union() %>%
    st_transform(4326) %>%
    st_area_ha() #find area in hectares

  #get ACD of undisturbed forest
  if(analysis_type == "epr26") {
    acd = read.csv(paste0("/maps/epr26/tmf_pipe_out/", proj_id, "/", proj_id, "carbon-density.csv"))
  } else if(analysis_type == "control") {
    acd = read.csv(paste0("/maps/epr26/tmf_pipe_out/0000_grid/", proj_id, "/", "0000_", proj_id, "carbon-density.csv"))
  }
  acd_u = acd %>% filter(land.use.class == 1) %>% pull(carbon.density)
  if(length(acd_u) == 0) acd_u = NA

  return(project_var = data.frame(t0 = t0, country = country, acd_u = acd_u, area_ha = area_ha))
}