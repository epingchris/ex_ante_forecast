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
analysis_type = "ac" #"old_source", "epr26", "ac"
process_grid = F
grid_id = "1201"

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
  acd_dir = "/maps/pf341/results/live-pipeline/"
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
  acd_dir = paste0(project_dir, projects, "/")
} else if(analysis_type == "ac") {
  project_dir = "/maps/epr26/tmf_pipe_out/"
  projects = list.files(project_dir, full = T) %>%
    str_subset("ac") %>%
    str_replace("/maps/epr26/tmf_pipe_out//", "")
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dir = paste0(project_dir, projects, "/")
}

# Loop through all projects ----
proj_list = mclapply(seq_along(projects), mc.cores = 15, function(i) { #mclapply() does not work on Windows
  a = Sys.time()
  proj_id = projects[i]
  proj_id_bare = proj_id %>% str_replace("a", "")
  pair_dir = pair_dirs[i]

  #extract project-level variables
  myproj = proj_meta %>% filter(ID == proj_id_bare)
  t0 = myproj$t0
  country = myproj$COUNTRY

  #extract the area of the region
  if(previous) {
    #previous: from supplier path (NEED TO CONSISTENTLY SET AOI NAME)
    site = paste("VCS", proj_id_bare, sep = "_")
    if(!str_detect(supplier_path, ".shp")) {
      aoi_path = file.path(supplier_path, site, "GIS", "aoi.shp")
      if(!file.exists(aoi_path)) aoi_path = file.path(supplier_path, site, paste0(site, ".shp"))
    }
    aoi_project = read_sf(aoi_path)
  } else {
    #now: from geojson
    aoi_project = st_read(paste0("/maps/epr26/tmf-data/projects/", proj_id, ".geojson"))
  }
  area_ha = aoi_project %>%
    st_make_valid() %>%
    st_union() %>%
    st_transform(4326) %>%
    st_area_ha() #find area in hectares

  #extract ACD per LUC
  acd_path = ifelse(length(acd_dir) > 1, acd_dir[i], acd_dir)
  acd = read.csv(paste0(acd_path, proj_id, "carbon-density.csv"))
  acd_u = acd %>% filter(land.use.class == 1) %>% pull(carbon.density)
  if(length(acd_u) == 0) acd_u = NA

  #gather project-level independent variables: area, ACD of undisturbed forest, country
  project_var = data.frame(acd_u = acd_u,
                           area_ha = area_ha,
                           country = country)

  #find paths to match and unmatached points in each sampled pairs
  pair_paths = list.files(pair_dir, full = T) %>% str_subset(".parquet")
  matchless_ind = pair_paths %>% str_detect("matchless")
  matchless_paths = pair_paths[matchless_ind]
  matched_paths = pair_paths[!matchless_ind]

  #find biome of project pixels
  k = read_parquet(paste0(project_dir, proj_id, "/", proj_id, "k.parquet")) %>%
    dplyr::select(c("lat", "lng", "ecoregion"))

  #find biome of matched pixels
  matches = read_parquet(paste0(project_dir, proj_id, "/", proj_id, "matches.parquet")) %>%
    dplyr::select(c("lat", "lng", "ecoregion"))

  #loop through all sampled pairs
  project_estimates = lapply(seq_along(matched_paths), function(j) {
    pairs = read_parquet(matched_paths[j]) %>%
      dplyr::left_join(., k, by = join_by(k_lat == lat, k_lng == lng)) %>%
      rename(k_ecoregion = ecoregion) %>%
      dplyr::left_join(., matches, by = join_by(s_lat == lat, s_lng == lng)) %>%
      rename(s_ecoregion = ecoregion) %>%
      mutate(s_id = 1:n(), k_id = 1:n())

    unmatched_pairs = read_parquet(matchless_paths[j])

    control = pairs %>%
      dplyr::select(starts_with("s_")) %>%
      rename_with(~str_replace(.x, "s_", "")) %>%
      mutate(treatment = "control") %>%
      tmfemi_reformat(t0 = t0)

    treat = pairs %>%
      dplyr::select(starts_with("k_")) %>%
      rename_with(~str_replace(.x, "k_", "")) %>%
      mutate(treatment = "treatment") %>%
      tmfemi_reformat(t0 = t0)

    exp_n_pairs = nrow(treat) + nrow(unmatched_pairs)

    pts_matched = rbind(treat, control)

    # Pair-level independent variables: median of all pixels in each pair (control + treat), then min/median/max across 100 pairs
    # elevation, slope, accessibility, cpc0/5/10_u, cpc0/5/10_d, defor_5_0 = cpc5_u - cpc0_u, defor_10_5 = cpc10_u - cpc5_u
    pair_var = pts_matched %>%
      dplyr::select(elevation:cpc10_d) %>%
      mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5) %>%
      reframe(elevation = median(elevation),
              slope = median(slope),
              accessibility = median(accessibility),
              cpc0_u = median(cpc0_u),
              cpc0_d = median(cpc0_d),
              cpc5_u = median(cpc5_u),
              cpc5_d = median(cpc5_d),
              cpc10_u = median(cpc10_u),
              cpc10_d = median(cpc10_d),
              defor_5_0 = median(defor_5_0),
              defor_10_5 = median(defor_10_5)) %>%
      pivot_longer(cols = elevation:defor_10_5, names_to = "var", values_to = "val") %>%
      mutate(pair = j)

    #calculate annual LUC change, carbon flux and additionality
    luc_series = simulate_area_series(pts_matched,
                                      class_prefix, t0 = t0, match_years, match_classes,
                                      exp_n_pairs, area_ha,
                                      verbose = F)

    carbon_series = luc_series$series %>%
      merge(., acd, by.x = "class", by.y = "land.use.class", all.x = T) %>%
      mutate(carbon_content = class_area * carbon.density) %>%
      group_by(treatment, year) %>%
      summarise(carbon_content = sum(carbon_content, na.rm = T)) %>%
      ungroup()

    carbon_wide = pivot_wider(carbon_series, names_from = "treatment", values_from = "carbon_content")
    out_df = data.frame(pair = j,
                        year = carbon_wide$year[-1],
                        c_loss = -diff(carbon_wide$control),
                        t_loss = -diff(carbon_wide$treatment)) %>%
      mutate(additionality = c_loss - t_loss)

    if(process_grid) {
      pts_matched = pts_matched %>% mutate(pair = j, pair_id = paste0(pair, "_", id))
      return(list(pair_var = pair_var, out_df = out_df, pts_matched = pts_matched))

    } else {
      return(list(pair_var = pair_var, out_df = out_df))
    }
  })

  pair_var_df = lapply(project_estimates, function(x) x$pair_var) %>% do.call(rbind, .)
  pair_var_summary = pair_var_df %>%
    group_by(var) %>%
    summarise(min = min(val), median = median(val), max = max(val)) %>%
    pivot_longer(cols = min:max, names_to = "stat", values_to = "val") %>%
    mutate(var = paste0(var, "_", stat)) %>%
    dplyr::select(c(var, val)) %>%
    pivot_wider(names_from = "var", values_from = "val")

  project_var_all = cbind(project_var, pair_var_summary) %>%
    mutate(project = proj_id)

  project_estimates = lapply(project_estimates, function(x) x$out_df) %>%
    do.call(rbind, .) %>%
    mutate(started = ifelse(year > t0, T, F))


  # Process gridded subplots ----
  if(process_grid & proj_id == grid_id) {
    pts_matched = lapply(project_estimates, function(x) x$pts_matched) %>% do.call(rbind, .)

    grid_path = paste0("/maps/epr26/tmf-data-grid/", proj_id)
    remerged_grid = read_rds(paste0(grid_path, "/", proj_id, "_remerged_grid.rds")) %>% st_as_sf()
    grid_n = nrow(remerged_grid)
    remerged_grid$name = 1:grid_n
    remerged_grid$area_ha = st_area_ha(remerged_grid)
    st_crs(pts_matched) = st_crs(remerged_grid)
    pts_matched_grid = st_join(pts_matched, remerged_grid, left = T)

    #balance check
    # table(pts_matched_grid$name, pts_matched_grid$treatment) #check balance: 31 is 0, 27 is 48 for Gola

    grid_list = lapply(seq_len(grid_n), function(i) {
      #calculate observed additionality in each grid
      pts_treatment_grid = pts_matched_grid %>% filter(treatment == "treatment", name == i)
      pts_control_grid = pts_matched_grid %>% filter(treatment == "control", pair_id %in% pts_treatment_grid$pair_id)
      pts_grid = rbind(pts_treatment_grid, pts_control_grid) %>% dplyr::select(-c("name", "area_ha"))

      if(nrow(pts_grid) == 0) return(list(out_df = NULL, carbon_loss = NULL, vicinity_area = NA, grid_var = NULL))

      luc_series = simulate_area_series(pts_grid,
                                        class_prefix, t0 = t0, match_years, match_classes,
                                        nrow(pts_grid) / 2, remerged_grid$area_ha[i],
                                        verbose = F)
      carbon_series = luc_series$series %>%
        merge(., acd, by.x = "class", by.y = "land.use.class", all.x = T) %>%
        mutate(carbon_content = class_area * carbon.density) %>%
        group_by(treatment, year) %>%
        summarise(carbon_content = sum(carbon_content, na.rm = T)) %>%
        ungroup()

      carbon_wide = pivot_wider(carbon_series, names_from = "treatment", values_from = "carbon_content")
      out_df = data.frame(grid = i,
                          year = carbon_wide$year[-1],
                          c_loss = -diff(carbon_wide$control),
                          t_loss = -diff(carbon_wide$treatment)) %>%
        mutate(additionality = c_loss - t_loss)

      grid_dir = paste0(project_dir, proj_id, "_grid/", i, "/")
      #ACD change associated with LUC change from undisturbed to deforested for each grid
      acd = read.csv(paste0(grid_dir, proj_id, "_", i, "carbon-density.csv")) #for vicinity baseline carbon loss rate (MgC / ha)
      acd_change = ifelse(sum(is.na(acd$carbon.density[c(1, 3)])) == 0, acd$carbon.density[1] - acd$carbon.density[3], NA)

      #calculate vicinity baseline carbon loss rate for each grid
      carbon_loss = read_parquet(paste0(grid_dir, proj_id, "_", i, "matches.parquet")) %>%
        dplyr::select(access, cpc0_u:cpc10_d) %>%
        mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
        mutate(carbon_loss_5_0 = defor_5_0 * acd_change, carbon_loss_10_5 = defor_10_5 * acd_change, carbon_loss_10_0 = defor_10_0 * acd_change) %>%
        mutate(grid = i)

      vicinity_area = nrow(carbon_loss) * 900 / 10000 #convert from number of 30x30m2 pixels to hectare

      # Grid-level independent variables: min/median/max  of all pixels in each grid (control + treat)
      # elevation, slope, accessibility, cpc0/5/10_u, cpc0/5/10_d, defor_5_0 = cpc5_u - cpc0_u, defor_10_5 = cpc10_u - cpc5_u
      grid_var = pts_grid %>%
        dplyr::select(elevation:biome) %>%
        mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
        reframe(elevation = summary(elevation),
                slope = summary(slope),
                accessibility = summary(accessibility),
                cpc0_u = summary(cpc0_u),
                cpc0_d = summary(cpc0_d),
                cpc5_u = summary(cpc5_u),
                cpc5_d = summary(cpc5_d),
                cpc10_u = summary(cpc10_u),
                cpc10_d = summary(cpc10_d),
                defor_5_0 = summary(defor_5_0),
                defor_10_5 = summary(defor_10_5),
                defor_10_0 = summary(defor_10_0)) %>%
        mutate(type = c("min", "q1", "median", "mean", "q3", "max")) %>%
        pivot_longer(cols = elevation:defor_10_0, names_to = "var", values_to = "val") %>%
        mutate(var = paste0(var, "_", type)) %>%
        select(-type) %>%
        pivot_wider(names_from = "var", values_from = "val") %>%
        mutate(biome = names(table(pts_grid$biome))[which.max(table(pts_grid$biome))], grid = i)

      return(list(out_df = out_df, carbon_loss = carbon_loss, vicinity_area = vicinity_area, grid_var = grid_var))
    })

    obs_add_grid = lapply(grid_list, function(x) x$out_df) %>% do.call(rbind, .)
    carbon_loss_grid = lapply(grid_list, function(x) x$carbon_loss) %>% do.call(rbind, .)
    vicinity_area_grid = sapply(grid_list, function(x) x$vicinity_area)
    grid_var = lapply(grid_list, function(x) x$grid_var) %>% do.call(rbind, .)

    grid_out_list = list(obs_add = obs_add_grid, vicinity_area = vicinity_area_grid, c_loss = c_loss_grid)
  }

  b = Sys.time()
  cat(proj_id, ":", b - a, "\n")
  if(process_grid & proj_id == grid_id) {
    return(list(project_estimates = project_estimates, project_var = project_var_all, grid_out = grid_out_list))
  } else {
    return(list(project_estimates = project_estimates, project_var = project_var_all))
  }
})

# Save outputs ----
out_path = paste0("/maps/epr26/tmf_pipe_out/")

project_var = lapply(proj_list, function(x) x$project_var) %>% do.call(rbind, .)
project_estimates = lapply(proj_list, function(x) x$project_estimates)
names(project_estimates) = projects

if(analysis_type == "ac") {
  saveRDS(project_var, file.path(paste0(out_path, "ac_project_var.rds")))
  saveRDS(project_estimates, file.path(paste0(out_path, "ac_project_estimates.rds")))
} else {
  saveRDS(project_var, file.path(paste0(out_path, "project_var.rds")))
  saveRDS(project_estimates, file.path(paste0(out_path, "project_estimates.rds")))
}




#@@@to be sorted out: plot to compare addditionaly whole or by grid@@@
add_by_grid = obs_add_grid %>%
  do.call(rbind, .) %>%
  group_by(year) %>%
  summarise(additionality = sum(additionality)) %>%
  ungroup() %>%
  filter(!is.na(year))

add_whole = project_estimates %>%
  group_by(year) %>%
  summarise(additionality = mean(additionality),
            add.max = max(additionality),
            add.min = min(additionality)) %>%
  ungroup() %>%
  filter(!is.na(year))

p_add_compare = ggplot(data = add_whole, aes(x = year, y = additionality)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = add.min, ymax = add.max), color = "grey", alpha = 0.5) +
  geom_line(data = add_by_grid, color = "red") +
  theme_bw()