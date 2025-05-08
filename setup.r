# This script generates a dataframe containing the information needed for subsequent analysis and saves it in out_path
# The dataframe contains the following columns:
# ID: project ID
# COUNTRY: project country
# t0: project start year
# area_ha: project area in hectares
# cdens_1 to cdens_6: carbon density (MgC ha-1) per land class
# code (for ongoing projects): anonymised project code

# It requires the following inputs:
#1. analysis_type: "ongoing" for ongoing projects; "placebo" for placebo projects
#2. project_dir: absolute path of the directory containing implementation output
#3. polygon_dir: absolute path of the directory containing project shapefiles
#4. proj_info: data frame containing project ID, country and start year
#5. include_strings: vector containing strings to include when searching for directories containing implementation output
#6. exclude_strings: vector containing strings to exclude when searching for directories containing implementation output
#7. out_path: absolute path where the output data frame is stored

rm(list = ls())

#Load packages
library(tidyverse) #ggplot2, dplyr, stringr, plotPlacebo/plotBaseline.r: tibble to store labels with bquote()
library(magrittr) #pipe operators
library(units) #units::set_units
library(sf) #sf::st_area; runs on GDAL 3.10

#Load wrapper function
FindFiles = function(dir, include = NULL, exclude = NULL, full = F) {
  files = list.files(dir, full = full)

  if (!is.null(include)) {
    include_pattern = paste(include, collapse = "|")
    files = files %>% str_subset(include_pattern)
  }

  if (!is.null(exclude)) {
    exclude_pattern = paste(exclude, collapse = "|")
    files = files %>% str_subset(exclude_pattern, negate = T)
  }

  if(length(files) == 0) {
    return(NA)
  } else {
    return(files)
  }
}

#define input variables
analysis_type = "ongoing" #analysis type
polygon_dir = "/maps/epr26/tmf-data/projects/" #where polygons are stored
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored

#define input variables that depend on analysis type
if(analysis_type == "ongoing") {
  #define where to look for directories containing implementation outputs
  project_dir = "/maps/epr26/tmf_pipe_out/"

  #load basic information (csv file copied from Tom's directory)
  proj_info = read.csv("proj_meta.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #define include and exclude strings
  include_strings = NULL
  exclude_strings = c("archive", "slopes", "elevation", "srtm", "\\.", "\\_")

} else if(analysis_type == "placebo") {
  #define where to look for directories containing implementation outputs
  project_dir = "/maps/epr26/tmf_pipe_out_luc_t/"

  #load basic information
  proj_info = read.csv("proj_meta_placebo.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #define include and exclude strings
  include_strings = c("asn", "af", "sa")
  exclude_strings = "\\."
}

#Find directories containing implementation outputs and save project names in vector "projects"
projects = FindFiles(project_dir, include = include_strings, exclude = exclude_strings)


#Retrieve data frames containing carbon density (MgC/ha) per LUC
cdens_paths = rep(NA, length(projects))
cdens_list = vector("list", length(projects))
for(i in seq_along(projects)) {
  cdens_path = FindFiles(paste0(project_dir, projects[i]), "carbon-density", full = T)
  if(!is.na(cdens_path)) {
    cdens = read.csv(cdens_path)[, 1:2] #use only the second column as carbon density value
    colnames(cdens) = c("land.use.class", "carbon.density")
    for(class in 1:6) {
      if(class %in% cdens$land.use.class == F) cdens = rbind(cdens, c(class, NA))
    }
    cdens_list[[i]] = cdens[order(cdens$land.use.class), ] #order land class from 1 to 6
    cdens_paths[i] = cdens_path
  }
}
names(cdens_list) = projects

#Check if carbon density values for LUC 1, 2, 3, and 4 are available
is_carbon_complete = sapply(cdens_list, function(x) !is.na(sum(x[1:4, ]$carbon.density)))

#Check if pairs parquet files are present (indicating complete output)
is_done = sapply(projects, function(x) FindFiles(paste0(project_dir, x, "/pairs"), ".parquet") %>% length() == 200)

#Select projects with complete output
projects_status = data.frame(project = projects, carbon_complete = is_carbon_complete, done = is_done)
projects = subset(projects_status, carbon_complete & done)$project

#Retrieve project variables
project_var = proj_info[match(projects, proj_info$ID), ]
t0_vec = project_var$t0

#Retrieve project areas (ha)
polygon_paths = paste0(polygon_dir, projects, ".geojson")
area_ha_vec = sapply(seq_along(projects), function(i) {
  area_ha_i = st_read(polygon_paths[i]) %>%
    st_make_valid() %>%
    st_union() %>%
    st_transform(4326) %>%
    st_area() %>% #area in m^2
    set_units(ha) #convert into hectares
  return(area_ha_i)
})
project_var$area_ha = area_ha_vec

#Retrieve carbon density
cdens_df = cdens_list %>%
  imap(function(.x, .y) { #map list element name to new column
    .x %>%
      mutate(project = .y) %>%
      filter(land.use.class != 0) %>% #land use class 0 doesn't mean anything
      pivot_wider(names_from = land.use.class, values_from = carbon.density, names_prefix = "cdens_")
  }) %>%
  list_rbind()
project_var = merge(project_var, cdens_df, by.x = "ID", by.y = "project")

#Sort rows by project ID
if(analysis_type == "ongoing") {
  project_var = project_var %>%
    arrange(as.numeric(ID)) %>%
    mutate(code = LETTERS[1:nrow(project_var)])
} else {
  project_var = project_var %>%
    arrange(ID)
}

#Output: project status check results
write.csv(projects_status, paste0(out_path, "_project_status.csv"), row.names = F)

#Output: project-level variables
write.csv(project_var, paste0(out_path, "_project_var.csv"), row.names = F)
write.csv(subset(project_var, select = -ID), paste0(out_path, "_project_var_anon.csv"), row.names = F)
