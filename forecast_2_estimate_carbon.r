# 0. Setup ----
rm(list = ls())

#Load packages
library(tidyverse) #ggplot2, dplyr, and stringr used in plotPlacebo/plotBaseline.r: tibble to store labels with bquote()
library(magrittr) #pipe operators
library(sf) #st_drop_geometry() used in GetCarbonLoss.r (runs on GDAL 3.10)
library(arrow) #read_parquet()
library(MatchIt) #matchit(), used in the customised function AssessBalance()

options(dplyr.summarise.inform = F) #remove dplyr summarise grouping message because it prints a lot

#Load pre-defined functions
source("FindFiles.r") #wrapper function to search files or folders based on inclusion/exclusion keywords
source("ReformatPixels.r") #wrapper function to reformat column names of data frame of matched pixels
source("AssessBalance.r") #function to assess matching balance
source("GetCarbonLoss.r") #function to generate time series of annual change in project-level average carbon density
source("BootOut.r")

#Define input variables needed to read TMF implementation output and other data

# It requires the following input variables to read TMF implementation output and other data.
# All variables are vectors containing one value for each project to be analysed:

#1. projects: an index of all projects to be analysed
# This should correspond to the filenames of the shapefiles and to the -p argument in the implementation code
# It is usually be the ongoing projects' VCS ID or customised (e.g. prefixed series of integers)

#2. project_out_dirs: absolute paths of the directories containing all implementation outputs
# The directory should containing pairs of parquet files with the same file name, with and without the "_matchless" suffix.
# (typically "/pairs/xxx.parquet" and  "/pairs/xxx_matchless.parquet")
# This is used to calculate estimated observed additionality.

#3. k_paths: absolute paths of the set K (typically "k.parquet")
#4. m_paths: absolute paths of the set M (typically "matches.parquet")
# Both should be in parquet format, containing the following columns:
# "lat", "lng" (degrees), "slope" (degrees), "elevation" (metres), "access" (remoteness, in minutes), "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d" (from 0 to 1), "luc_[t-10]" to "luc_2021" (categorical, 1-6, based on the JRC-TMF dataset, Vancutsem et al. 2021), "ecoregion" (categorical, based on the RESOLVE dataset, Dinerstein et al. 2017)

#5. acd_paths: absolute paths of the carbon density per LUC (typically "carbon-density.csv")
# This should be an csv file (although the script could be modiified in the future to support txt format) containing columns "land.use.class" (categorical, 1-6) and "carbon.density" (MgC/ha) for all six LUCs, although the script checks and fill missing LUC with NAs

#6. polygon_paths: absolute paths of the shapefile of the project extent
# This should be a geojson file containing valid geometries in WGS84 (EPSG: 4326), although the script checks for both conditions.
# This is currently only used to calculate project area (ha), but could be useful for other purposes in the future.

#7. country: country of the project
#8. t0: year of start of the project (real or hypothetical)
#9. OPTIONAL: proj_ID: this is only used to remove the trailing "a" that I added to the filename of some shapefiles that needed fixing,
# So that I can retrieve country and t0 information from the the proj_info data frame
#10. out_path: absolute paths of the directory where outputs are to be saved; include file prefix if desired


# A. Read input ----

analysis_type = "ongoing"
in_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #path of setup script output
project_dir = "/maps/epr26/tmf_pipe_out/" #path to directories containing implementation outputs
lagged_dir = "/maps/epr26/tmf_pipe_out_lagged/" #path to directories containing implementation outputs (lagged matching)
project_var = read.csv(paste0(in_path, "_project_var.csv"), header = T)
projects = project_var$ID
t0_vec = project_var$t0
area_ha_vec = project_var$area_ha
cdens_list = project_var %>%
  dplyr::select(ID, cdens_1:cdens_6) %>%
  pivot_longer(cdens_1:cdens_6, names_to = "land.use.class", names_prefix = "cdens_", values_to = "carbon.density") %>%
  split(f = .$ID)

#Retrieve input paths needed for the analysis
project_out_dirs = paste0(project_dir, projects)
k_paths = rep(NA, length(projects))
m_paths = rep(NA, length(projects))
for(i in seq_along(projects)) {
  k_paths[i] = FindFiles(project_out_dirs[i], "k.parquet", full = T)
  m_paths[i] = FindFiles(project_out_dirs[i], "matches.parquet", full = T)
}

fig_path = paste0("/maps/epr26/ex_ante_forecast_out/out_") #where figures are stored
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored


# B. Get ex post additionality, project and regional carbon loss rates ----
for(i in seq_along(projects)) {
  t0 = t0_vec[i]
  luc_t_20 = paste0("luc_", t0_vec[i] - 20)
  luc_t_10 = paste0("luc_", t0_vec[i] - 10)
  luc_t0 = paste0("luc_", t0_vec[i])
  area_ha = area_ha_vec[i]
  cdens = cdens_list[[i]]
  pair_dir = paste0(project_out_dirs[i], "/pairs/")

  a = Sys.time()
  #find paths to match and unmatached points in each sampled pairs
  pair_paths = FindFiles(pair_dir, include = ".parquet", full = T)
  matched_paths = pair_paths %>% str_subset("matchless", negate = T)
  matchless_paths = pair_paths %>% str_subset("matchless")

  if(length(matched_paths) == 0) {
    cat("No matches\n")
    closs_df = NULL
    expost_add = NULL
    closs_project = NULL
  } else {
    pairs_out = lapply(seq_along(matched_paths), function(j) {
      pair_start = Sys.time()

      matched_path = matched_paths[j]
      matchless_path = matchless_paths[j]

      pairs = read_parquet(matched_path) %>%
        dplyr::select(-dplyr::ends_with(c("_x", "_y", "_trt", "_cluster")))
      unmatched_pairs = read_parquet(matchless_path) %>%
        dplyr::select(-dplyr::ends_with(c("_x", "_y", "_trt", "_cluster")))

      #calculate proportion of unmatched pixels, used to adjust the area represented by sample: is this still necessary?
      area_adj_ratio = nrow(pairs) / (nrow(pairs) + nrow(unmatched_pairs))

      pixels_cf = ReformatPixels(in_df = pairs, prefix = "s_", t0 = t0, treatment = "counterfactual", pair = j)
      pixels_p = ReformatPixels(in_df = pairs, prefix = "k_", t0 = t0, treatment = "project", pair = j)
      pixels_matched = rbind(pixels_cf, pixels_p)

      #assess matching balance
      balance_assessment = AssessBalance(pixels_matched, t0 = t0)

      #retrieve carbon time series for project and matched counterfactual
      carbon_cf = GetCarbonLoss(pixels_cf, t0, area_ha, area_adj_ratio, cdens, pair = j)
      carbon_p = GetCarbonLoss(pixels_p, t0, area_ha, area_adj_ratio, cdens, pair = j)

      pair_end = Sys.time()
      cat(j, ":", pair_end - pair_start, "\n")

      return(list(carbon_cf = carbon_cf, carbon_p = carbon_p,
                  pixels_matched = pixels_matched, area_adj_ratio = area_adj_ratio,
                  balance_assessment = balance_assessment))
    })
  }

  b = Sys.time()
  cat("Project", i, "/", length(projects), "-", projects[i], "- project carbon loss :", b - a, "\n")

  a = Sys.time()
  setM = read_parquet(m_paths[i])
  if(nrow(setM) > 250000) setM = setM[sample(nrow(setM), 250000), ]
  pixels_region = ReformatPixels(in_df = setM, prefix = "", t0 = t0, treatment = "region", pair = 1)

  #retrieve carbon time series for surrounding region
  closs_region = GetCarbonLoss(pixels_region, t0, area_ha, area_adj_ratio = 1, cdens, pair = 1)
  b = Sys.time()
  cat("Project", i, "/", length(projects), "-", projects[i], "- regional carbon loss :", b - a, "\n")

  carbon_cf = lapply(pairs_out, function(x) x$carbon_cf) %>%
    list_rbind() %>%
    mutate(treatment = "counterfactual")

  carbon_p = lapply(pairs_out, function(x) x$carbon_p) %>%
    list_rbind() %>%
    mutate(treatment = "project")

  carbon_matched = rbind(carbon_cf, carbon_p) %>%
    pivot_wider(values_from = "carbon_density", names_from = "treatment")
  tmax = max(carbon_matched$year)

  #calculate ex post annual per-area additionality
  additionality = carbon_matched %>%
    filter(year >= t0) %>%
    group_by(pair) %>%
    mutate(diff_cf = first(counterfactual) - counterfactual,
           diff_p = first(project) - project,
           additionality_whole = diff_cf - diff_p,
           additionality_annual = exp(log(additionality_whole) / (year - t0)))

  #calculate ex ante within-project annual per-area carbon loss rates
  ante_project = carbon_matched %>%
    filter(year <= t0) %>%
    group_by(pair) %>%
    mutate(closs = project - last(project),
           closs_annual = exp(log(closs) / (t0 - year)))

  #calculate ex ante regional annual per-area carbon loss rates
  ante_region = closs_region %>%
    filter(year <= t0) %>%
    group_by(pair) %>%
    mutate(closs = carbon_density - last(carbon_density),
           closs_annual = exp(log(closs) / (t0 - year))) %>%
    mutate(closs_annual = ifelse(closs_annual >= 0, closs_annual, NA))

  write.csv(additionality, paste0(out_path, "_additionality_", projects[i], ".csv"), row.names = F)
  write.csv(ante_project, paste0(out_path, "_project_closs_rate_", projects[i], ".csv"), row.names = F)
  write.csv(ante_region, paste0(out_path, "_regional_closs_rate_", projects[i], ".csv"), row.names = F)
}


# C. Bootstrap outcomes ----
additionality_boot_list = vector("list", length(projects))
ante_project_boot_list = vector("list", length(projects))
ante_region_boot_list = vector("list", length(projects))

for(i in seq_along(projects)) {
  t0 = t0_vec[i]
  project_i = projects[i]
  additionality = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T)
  ante_project = read.csv(paste0(out_path, "_project_closs_rate_", projects[i], ".csv"), header = T)
  ante_region = read.csv(paste0(out_path, "_regional_closs_rate_", projects[i], ".csv"), header = T)
  tmax = max(additionality$year)

  #bootstrap ex post additionality
  a = Sys.time()
  additionality_boot_list[[i]] = BootOut(in_df = additionality, column = "additionality_annual", from = t0 + 1, to = tmax) %>%
    mutate(project = project_i)
  b = Sys.time()
  cat("Project", i, "/", length(projects), "-", projects[i], "- additionality bootstrapped:", b - a, "\n")

  #bootstrap ex ante project carbon loss rate
  a = Sys.time()
  ante_project_boot_list[[i]] = BootOut(in_df = ante_project, column = "closs_annual", from = t0 - 10, to = t0 - 1) %>%
    mutate(project = project_i)
  b = Sys.time()
  cat("Project", i, "/", length(projects), "-", projects[i], "- project rate bootstrapped:", b - a, "\n")

  #bootstrap ex ante regional carbon loss rate
  a = Sys.time()
  ante_region_boot_list[[i]] = BootOut(in_df = ante_region, column = "closs_annual", from = t0 - 10, to = t0 - 1) %>%
    mutate(project = project_i)
  b = Sys.time()
  cat("Project", i, "/", length(projects), "-", projects[i], "- project rate bootstrapped:", b - a, "\n")

  # effectiveness = observed_add_boot_out$t / baseline_best_boot_out$t
  # effectiveness_list[[i]] = data.frame(project = projects[i],
  #                                       eff = mean(effectiveness, na.rm = T),
  #                                       eff_lower = quantile(effectiveness, probs = 0.025, na.rm = T),
  #                                       eff_upper = quantile(effectiveness, probs = 0.975, na.rm = T))
}

additionality_boot_df = list_rbind(additionality_boot_list)
ante_project_boot_df = list_rbind(ante_project_boot_list)
ante_region_boot_df = list_rbind(ante_region_boot_list)

write.csv(additionality_boot_df, paste0(out_path, "_boot_additionality.csv"), row.names = F)
write.csv(ante_project_boot_df, paste0(out_path, "_boot_project_closs_rate.csv"), row.names = F)
write.csv(ante_region_boot_df, paste0(out_path, "_boot_regional_closs_rate.csv"), row.names = F)