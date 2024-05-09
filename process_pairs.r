# set out a table of parameters you want to change
# read the config
# write the config with the parameters
# Run the script

rm(list = ls())

# install.packages(c("arrow","configr", "tidyverse", "magrittr", "sf", "magrittr", "MatchIt",
#                    "rnaturalearthdata", "configr", "terra", "pbapply", "cleangeo", "doParallel",
#                    "foreach", "readr", "lwgeom", "rnaturalearth", "stars"), depends = TRUE)

library(tidyverse)
library(configr)
library(magrittr)
library(readr)
library(sf)
library(stars)
library(arrow) #arrow::read_parquet
library(MatchIt) #MatchIt::matchit
library(parallel) #parallel::mclapply
# library(rnaturalearthdata)
# library(terra)
# library(pbapply)
# library(cleangeo)
# library(foreach)
# library(lwgeom)
# library(countrycode)

source("functions.r") #cpc_rename, tmfemi_reformat
source("ProcessPairs.r") #RetrievePoints, ProcessPairs

orig_dir = getwd()
setwd("/home/tws36/4c_evaluations")

# The list of projects to be run in this evaluation:
proj_meta = read.csv("./data/project_metadata/proj_meta.csv")

config = read.config("./config/fixed_config_sherwood.ini")
config$USERPARAMS$data_path = "/maps/pf341/tom"
#write.config(config, "./config/fixed_config_tmp.ini") #error: permission denied

# Load user-defined functions that Tom wrote
sapply(list.files("./R", full = T, pattern = ".R$"), source)

# Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = FALSE)

source("./R/scripts/setup_and_reusable/load_config.R")
# source("./R/scripts/0.2_load_project_details.R")

setwd(orig_dir)

class_prefix = "JRC"
match_years = c(0, -5, -10)
match_classes = c(1, 3)

# Load the IDs of projects to process ----
#input: project_dir, exclude_id, acd_id
#output: projects, pair_dirs, acd_dir
analysis_type = "control" #"old_source", "epr26", "grid", "ac", "control"


if(analysis_type == "old_source") {
  #previous: from /maps/pf341
  project_dir = "/maps/pf341/results/2024-january-pipeline"
  projects = list.files(project_dir, full = T) %>%
    str_subset("pairs") %>%
    basename() %>%
    str_replace("_pairs", "")
  #select those with carbon density data, unselect problematic projects (1566, 1067, 958, 1133) and projects not in TMF extent (562)
  exclude_id = c("1566", "1067", "958", "1133", "562")
  #select projects with a non-empty carbon_density.csv
  acd_dirs = "/maps/pf341/results/live-pipeline/"
  acd_id = list.files(acd_dir, full = T) %>%
    str_subset("carbon-density") %>%
    basename() %>%
    str_replace("-carbon-density.csv", "")
  projects = projects[which(projects %in% acd_id) & !(projects %in% exclude_id)] %>%
    as.numeric() %>%
    sort() %>%
    as.character()
  pair_dirs = paste0("/maps/pf341/results/2024-january-pipeline/", projects, "_pairs/")
} else if(analysis_type == "epr26") {
  #now: from /maps/epr26
  project_dir = "/maps/epr26/tmf_pipe_out/"
  projects = list.files(project_dir, full = T) %>%
    str_subset("\\.", negate = T) %>%
    str_subset("\\_grid", negate = T) %>%
    str_subset("fit\\_distribution", negate = T) %>%
    str_replace("/maps/epr26/tmf_pipe_out//", "")
  #remove: 0000, 9999 are controls and test; 562 and 1202 have incomplete ACD; 1399 and 1408 anomalous
  exclude_id = c("0000", "9999", "562", "1202", "1399", "1408")
  projects = projects[!(projects %in% exclude_id)] %>%
    str_replace("a", "") %>%
    as.numeric() %>%
    sort() %>%
    as.character()
  projects[which(projects %in% c("612", "1340"))] %<>% str_c("a") #612, 1340 replaced by 612a, 1340a
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")
} else if(analysis_type == "grid") {
  project_dir = "/maps/epr26/tmf_pipe_out/1201_grid/"
  projects = 1:10
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")
} else if(analysis_type == "control") {
  project_dir = "/maps/epr26/tmf_pipe_out/0000_grid/"
  projects = c(2, 3, 4, 5, 7, 8)
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")
} else if(analysis_type == "ac") {
  project_dir = "/maps/epr26/tmf_pipe_out/"
  projects = list.files(project_dir, full = T) %>%
    str_subset("ac") %>%
    str_replace("/maps/epr26/tmf_pipe_out//", "")
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")
}

# Loop through all projects ----
proj_list = lapply(seq_along(projects), function(i) { #mclapply() does not work on Windows
  a = Sys.time()

  proj_id = projects[i]
  if(analysis_type == "old_source") {
    cat("Use new results from epr26 instead.\n")
    return(NULL)
  } else if(analysis_type == "epr26") {
    proj_id_bare = proj_id %>% str_replace("a", "")
    myproj = proj_meta %>% filter(ID == proj_id_bare)
    t0 = myproj$t0
    country = myproj$COUNTRY

    file_prefix = paste0(project_dir, proj_id, "/", proj_id)

    aoi_project = st_read(paste0("/maps/epr26/tmf-data/projects/", proj_id, ".geojson"))
    acd = read.csv(paste0("/maps/epr26/tmf_pipe_out/", proj_id, "/", proj_id, "carbon-density.csv"))
  } else if(analysis_type == "grid") {
    myproj = proj_meta %>% filter(ID == "1201")
    t0 = myproj$t0
    country = myproj$COUNTRY

    file_prefix = paste0(project_dir, proj_id, "/1201_", proj_id)

    aoi_project = st_read(paste0("/maps/epr26/tmf-data-grid/1201/1201_", proj_id, ".geojson"))
    acd = read.csv(paste0("/maps/epr26/tmf_pipe_out/1201_grid/", proj_id, "/", "1201_", proj_id, "carbon-density.csv"))
  } else if(analysis_type == "control") {
    t0 = 2011
    country = "Brazil"

    file_prefix = paste0(project_dir, proj_id, "/0000_", proj_id)

    aoi_project = st_read(paste0("/maps/epr26/tmf-data-grid/0000/0000_", proj_id, ".geojson"))
    acd = read.csv(paste0("/maps/epr26/tmf_pipe_out/0000_grid/", proj_id, "/", "0000_", proj_id, "carbon-density.csv"))
  }

  #get area_ha
  area_ha = aoi_project %>%
    st_make_valid() %>%
    st_union() %>%
    st_transform(4326) %>%
    st_area_ha() #find area in hectares

  #get ACD of undisturbed forest, gather project-level variables
  acd_u = acd %>% filter(land.use.class == 1) %>% pull(carbon.density)
  if(length(acd_u) == 0) acd_u = NA

  #find paths to match and unmatached points in each sampled pairs
  pair_paths = list.files(pair_dirs[i], full = T) %>% str_subset(".parquet")
  matchless_ind = pair_paths %>% str_detect("matchless")
  matchless_paths = pair_paths[matchless_ind]
  matched_paths = pair_paths[!matchless_ind]

  #find biome of project pixels
  k = read_parquet(paste0(file_prefix, "k.parquet")) %>%
    dplyr::select(c("lat", "lng", "ecoregion"))

  #find biome of matched pixels
  matches = read_parquet(paste0(file_prefix, "matches.parquet")) %>%
    dplyr::select(c("lat", "lng", "ecoregion"))

  #loop through all sampled pairs
  pairs_out = lapply(seq_along(matched_paths), function(j) {
    ProcessPairs(matched_path = matched_paths[j], matchless_path = matchless_paths[j],
                 k = k, matches = matches, t0 = t0, area_ha = area_ha, acd = acd, pair_id = j)
  })

  #combine pair-level variables and project-level variables
  pair_var_summary = lapply(pairs_out, function(x) x$pair_var) %>%
    do.call(rbind, .) %>%
    group_by(var) %>%
    summarise(min = min(val), median = median(val), max = max(val)) %>%
    pivot_longer(cols = min:max, names_to = "stat", values_to = "val") %>%
    mutate(var = paste0(var, "_", stat)) %>%
    dplyr::select(c(var, val)) %>%
    pivot_wider(names_from = "var", values_from = "val")

  project_var = data.frame(t0 = t0, country = country, acd_u = acd_u, area_ha = area_ha) %>%
    cbind(., pair_var_summary) %>%
    mutate(project = paste0("0000_", proj_id))

  project_estimates = lapply(pairs_out, function(x) x$out_df) %>%
    do.call(rbind, .) %>%
    mutate(started = ifelse(year > t0, T, F))

  b = Sys.time()
  cat(proj_id, ":", b - a, "\n")
  return(list(project_var = project_var, project_estimates = project_estimates))
})


# Save outputs ----
out_path = "/maps/epr26/tmf_pipe_out/"

project_var = lapply(proj_list, function(x) x$project_var) %>% do.call(rbind, .)

project_estimates = lapply(proj_list, function(x) x$project_estimates)
names(project_estimates) = projects

if(analysis_type == "ac") {
  saveRDS(project_var, file.path(paste0(out_path, "ac_project_var.rds")))
  saveRDS(project_estimates, file.path(paste0(out_path, "ac_project_estimates.rds")))
} else if(analysis_type == "grid") {
  saveRDS(project_var, file.path(paste0(out_path, "grid_1201_project_var.rds")))
  saveRDS(project_estimates, file.path(paste0(out_path, "grid_1201_project_estimates.rds")))
} else if(analysis_type == "control") {
  saveRDS(project_var, file.path(paste0(out_path, "control_project_var.rds")))
  saveRDS(project_estimates, file.path(paste0(out_path, "control_project_estimates.rds")))
} else {
  saveRDS(project_var, file.path(paste0(out_path, "project_var.rds")))
  saveRDS(project_estimates, file.path(paste0(out_path, "project_estimates.rds")))
}