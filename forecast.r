# 0. Setup ----
rm(list = ls())

#Load packages
# install.packages(c("arrow","configr", "tidyverse", "magrittr", "sf", "magrittr", "MatchIt",
#                    "rnaturalearthdata", "configr", "terra", "pbapply", "cleangeo", "doParallel",
#                    "foreach", "readr", "lwgeom", "rnaturalearth", "stars", "Metrics", "patchwork"), depends = TRUE)

library(tidyverse) #ggplot2, dplyr, stringr, plotPlacebo/plotBaseline.r: tibble to store labels with bquote()
library(magrittr) #pipe operators
library(units) #units::set_units
library(sf) #sf::st_area; runs on GDAL 3.10
library(arrow) #arrow::read_parquet
library(MatchIt) #MatchIt::matchit
library(boot) #boot::boot
library(scales) #scales::trans_break
library(Metrics) #CalcError.r: rmse, mae
library(patchwork)
#library(pryr) #pryr::object_size
#library(parallel) #parallel::mclapply

#Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = F)

#Load pre-defined functions
source("functions.r") #cpc_rename, tmfemi_reformat, simulate_area_series, make_area_series, assess_balance, make_match_formula
source("AdditionalityPair.r")
source("CalcError.r")
source("plotPlacebo.r")
source("plotBaseline.r")

FindFiles = function(dir, pattern, full = F, negate = F) {
  file_matched = list.files(dir, full = full) %>% str_subset(pattern, negate = negate)
  if(length(file_matched) == 0) {
    return(NA)
  } else {
    return(file_matched)
  }
}

BootOut = function(type, in_df, boot_n = 1000) {
  boot_out = boot::boot(data = in_df,
                        statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                        R = boot_n)
  boot_ci = boot::boot.ci(boot.out = boot_out, type = "perc")
  boot_summ = data.frame(type = type,
                          mean = mean(boot_out$t),
                          ci_lower = boot_ci$percent[4],
                          ci_upper = boot_ci$percent[5])
  return(list(t = boot_out$t, summ = boot_summ))
}


#Define input variables needed to read TMF implementation output and other data

# It requires the following input variables to read TMF implementation output and other data.
# All variables are vectors containing one value for each project to be analysed:

#1. projects: an index of all projects to be analysed
# This should correspond to the filenames of the shapefiles and to the -p argument in the implementation code
# It is usually be the ongoing projects' VCS ID or customised (e.g. prefixed series of integers)

#2. pair_dirs: absolute paths of the directories containing all matched pair pixel sets (typically "/pairs/xxx.parquet" and  "/pairs/xxx_matchless.parquet")
# The directory should containing pairs of parquet files with the same file name, with and without the "_matchless" suffix.
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


# A. Read input (E-Ping's workflow) ----

#Define analysis type
analysis_type = "control"
#"control": placebo "projects"
#"ongoing": ongoing REDD+ projects (best-matched and loosely-matched baselines)
#"ac": Amazonian Collective polygons
ofir = F

polygon_dir = "/maps/epr26/tmf-data/projects/" #where polygons are stored
lagged_dir = "/maps/epr26/tmf_pipe_out_lagged/" #where the results for the lagged baselines are stored
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored
fig_path = paste0("/maps/epr26/ex_ante_forecast_out/out_") #where figures are stored
projects_to_exclude = c("674", "934", "2502", "1408", "sa11", "1686", "1122", "1650") #which projects to exclude manually; 1122 not done yet

#Based on analysis type, find project names and store in vector "projects", and project basic info and store in dataframe "proj_info"
if(analysis_type == "ongoing") {
  project_dir = "/maps/epr26/tmf_pipe_out/" #define where implementation code results are stored

  #Load basic information (csv file copied from Tom's directory)
  proj_info = read.csv("proj_meta.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #Find all project IDs
  exclude_strings = c("slopes", "elevation", "srtm", "ac", "as", "\\.", "\\_", "0000", "9999")
  projects = map(exclude_strings, function(x) FindFiles(project_dir, x, negate = T)) %>%
    reduce(intersect)

} else if(analysis_type == "control") {
  project_dir = "/maps/epr26/tmf_pipe_out_luc_t/" #define where implementation code results are stored

  #Load basic information
  proj_info = read.csv("proj_meta_control.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #Find all project IDs
  projects = map(c("asn", "af", "sa"), function(x) FindFiles(project_dir, x)) %>%
    reduce(union) %>%
    str_subset("\\.", negate = T)

} else if(analysis_type == "ac") {
  project_dir = "/maps/epr26/tmf_pipe_out/" #define where implementation code results are stored

  #Load basic information
  proj_info = data.frame(ID = paste0("ac", sprintf("%02d", c(1:6))), COUNTRY = "Brazil", t0 = 2021)

  #Amazonian Collective polygons
  projects = paste0("ac", sprintf("%02d", c(1:6)))

}

done_vec = sapply(projects, function(x) FindFiles(paste0(project_dir, x, "/pairs"), ".parquet") %>% length() == 200)

#only keep projects with complete ACD values for LUC 1, 2, 3, and 4
full_acd_vec = sapply(projects, function(x) {
  acd_path = FindFiles(paste0(project_dir, x), "carbon-density", full = T)
  if(!is.na(acd_path)) acd = read.csv(acd_path)
  acd = acd[, 1:2]
  colnames(acd) = c("land.use.class", "carbon.density")
  return(Reduce("&", 1:4 %in% acd$land.use.class))
})

projects_df = data.frame(project = projects, done = done_vec, full_acd = full_acd_vec)
projects_df$to_exclude = projects_df$project %in% projects_to_exclude
write.csv(projects_df, paste0(out_path, "_project_status.csv"), row.names = F)
projects = subset(projects_df, done & full_acd & !to_exclude)$project

#Produce input variables needed for the analysis
pair_dirs = paste0(project_dir, projects, "/pairs/")
k_paths = rep(NA, length(projects))
m_paths = rep(NA, length(projects))
acd_paths = rep(NA, length(projects))
for(i in seq_along(projects)) {
  project_out_dir = paste0(project_dir, projects[i])
  k_paths[i] = FindFiles(project_out_dir, "k.parquet", full = T)
  m_paths[i] = FindFiles(project_out_dir, "matches.parquet", full = T)
  acd_paths[i] = FindFiles(project_out_dir, "carbon-density.csv", full = T)
}

pair_dirs_lagged = paste0(lagged_dir, projects, "/pairs/")
k_paths_lagged = rep(NA, length(projects))
m_paths_lagged = rep(NA, length(projects))
for(i in seq_along(projects)) {
  project_out_dir_lagged = paste0(lagged_dir, projects[i])
  k_paths_lagged[i] = FindFiles(project_out_dir_lagged, "k.parquet", full = T)
  m_paths_lagged[i] = FindFiles(project_out_dir_lagged, "matches.parquet", full = T)
}

polygon_paths = paste0(polygon_dir, projects, ".geojson")
proj_ID = gsub("(?<=[0-9])a$", "", projects, perl = T)
project_var = proj_info[match(proj_ID, proj_info$ID), ]
t0_vec = project_var$t0

#vector containing area (ha) of every project
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

#list containing data frame of ACD (MgC/ha) per LUC of every project
acd_list = lapply(seq_along(projects), function(i) {
  acd_i = read.csv(acd_paths[i])
  acd_i = acd_i[, 1:2]
  colnames(acd_i) = c("land.use.class", "carbon.density")
  for(class in 1:6) {
    if(class %in% acd_i$land.use.class == F) acd_i = rbind(acd_i, c(class, NA))
  }
  return(acd_i)
})
names(acd_list) = projects
project_var$acd_undisturbed = sapply(acd_list, function(x) filter(x, land.use.class == 1)$carbon.density)

if(analysis_type == "ongoing") {
  project_var = project_var %>%
    arrange(ID) %>%
    mutate(code = LETTERS[1:nrow(project_var)])
  project_var = project_var[match(projects, project_var$ID), ]
}

#Output: project-level variables
write.csv(project_var, paste0(out_path, "_project_var.csv"), row.names = F)
#project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)

#Output: carbon density per land class for S1 in manuscript
lapply(acd_list, function(x) {
  x = x %>%
    sort(land.use.class)
})

acd_df = acd_list %>%
  imap(function(.x, .y) {
    .x %>%
      mutate(project = .y) %>%
      filter(land.use.class != 0) %>% #land use class 0 doesn't mean anything
      pivot_wider(names_from = land.use.class, values_from = carbon.density, names_prefix = "class_")
  }) %>%
  list_rbind() %>%
  relocate(project, sort(tidyselect::peek_vars())) #sort columns by land class

if(analysis_type == "ongoing") {
  acd_df = acd_df %>%
    arrange(as.numeric(project)) #sort rows by project ID
} else {
  acd_df = acd_df %>%
    arrange(project) #sort rows by project ID
}
write.csv(acd_df, paste0(out_path, "_carbon_density_per_project.csv"), row.names = F)
#project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)


# A2. Read input (user-defined) ----

#projects = NULL
#pair_dirs = NULL
#k_paths = NULL
#m_paths = NULL
#pair_dirs_lagged = NULL
#k_paths_lagged = NULL
#m_paths_lagged = NULL
#polygon_paths = NULL
#t0_vec = NULL
#area_ha_vec = NULL
#acd_list = NULL
#out_path = NULL
#fig_path = NULL

pairs_lagged_list = vector("list", length(projects))

# B. Get additionality and baseline ----
for(i in seq_along(projects)) {

    t0 = t0_vec[i]
    luc_t_20 = paste0("luc_", t0_vec[i] - 20)
    luc_t_10 = paste0("luc_", t0_vec[i] - 10)
    luc_t0 = paste0("luc_", t0_vec[i])

    area_ha = area_ha_vec[i]
    acd = acd_list[[i]]
    pair_dir = pair_dirs[i]

    k = read_parquet(k_paths[i]) %>%
      rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
      as.data.frame() %>%
      dplyr::select(lat, lng, k_ecoregion)

    setM = read_parquet(m_paths[i]) %>%
      rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
      as.data.frame()

    matches = setM %>%
      dplyr::select(lat, lng, s_ecoregion)

    # a = Sys.time()
    # pairs_best = AdditionalityPair(pair_dir = pair_dir, t0 = t0, area_ha = area_ha, acd = acd,
    #                                k = k, matches = matches, lagged = F)
    # b = Sys.time()
    # cat("Project", i, "/", length(projects), "-", projects[i], "- pairs_best :", b - a, "\n")

    # additionality = lapply(pairs_best, function(x) x$out_df) %>%
    #     list_rbind() %>%
    #     mutate(started = ifelse(year > t0, T, F)) %>%
    #     mutate(project = projects[i])

    # baseline_best = lapply(pairs_best, function(x) {
    #     filter(x$out_df, year <= t0 & year > t0 - 10) %>%
    #     dplyr::select(c_loss) %>%
    #     mutate(c_loss = c_loss / area_ha)
    # }) %>%
    #     list_rbind() %>%
    #     mutate(project = projects[i])

    # write.csv(additionality, paste0(out_path, "_additionality_", projects[i], ".csv"), row.names = F)
    # write.csv(baseline_best, paste0(out_path, "_baseline_best_", projects[i], ".csv"), row.names = F)

    # setM = setM %>%
    #   dplyr::select(-starts_with("luc_")) %>%
    #   dplyr::select(-starts_with("cpc"))
    # if(nrow(setM) > 250000) setM = setM[sample(nrow(setM), 250000), ]

    # baseline_loose = setM %>%
    #     mutate(acd10 = acd$carbon.density[match(luc10, acd$land.use.class)],
    #            acd0 = acd$carbon.density[match(luc0, acd$land.use.class)],
    #            c_loss = (acd10 - acd0) / 10) %>%
    #     dplyr::select(c_loss) %>%
    #     mutate(project = projects[i])

    # write.csv(baseline_loose, paste0(out_path, "_baseline_loose_", projects[i], ".csv"), row.names = F)

    baseline_lagged = data.frame(c_loss = numeric())
    if(!(is.na(k_paths_lagged[i]) | is.na(m_paths_lagged[i]))) {
        pair_dir_lagged = pair_dirs_lagged[i]
        k_lagged = read_parquet(k_paths_lagged[i]) %>%
          rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
          as.data.frame() %>%
          dplyr::select(lat, lng, k_ecoregion)
        matches_lagged = read_parquet(m_paths_lagged[i]) %>%
          rename(luc10 = all_of(luc_t_20), luc0 = all_of(luc_t_10), s_ecoregion = ecoregion) %>%
          as.data.frame() %>%
          dplyr::select(lat, lng, s_ecoregion)

        a = Sys.time()
        pairs_lagged_list[[i]] = AdditionalityPair(pair_dir = pair_dir_lagged, t0 = t0, area_ha = area_ha, acd = acd,
                                         k = k_lagged, matches = matches_lagged, lagged = T)
        b = Sys.time()
        cat("Project", i, "/", length(projects), "-", projects[i], "- pairs_lagged :", b - a, "\n")
}

for(i in seq_along(projects)) {
        pairs_lagged = pairs_lagged_list[[i]]
        baseline_lagged = lapply(pairs_lagged, function(x) {
            filter(x$out_df, year < 10 & year > 0) %>%
            dplyr::select(c_loss) %>%
            mutate(c_loss = c_loss / area_ha)
        }) %>%
            list_rbind() %>%
            mutate(project = projects[i])
    write.csv(baseline_lagged, paste0(out_path, "_baseline_lagged_new_", projects[i], ".csv"), row.names = F)
}


# C. Bootstrap baselines ----
pre_cf_c_loss_boot_list = vector("list", length(projects))
pre_p_c_loss_boot_list = vector("list", length(projects))
post_cf_c_loss_boot_list = vector("list", length(projects))
post_p_c_loss_boot_list = vector("list", length(projects))
observed_add_boot_list = vector("list", length(projects))
baseline_best_boot_list = vector("list", length(projects))
baseline_loose_boot_list = vector("list", length(projects))
baseline_lagged_boot_list = vector("list", length(projects))

for(i in seq_along(projects)) {
  # area_i = area_ha_vec[i]
  # observed = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T) %>%
  #   mutate(t_loss = t_loss / area_i,
  #          c_loss = c_loss / area_i,
  #          additionality = additionality / area_i)
  # obs_pre = observed %>% filter(started == F)
  # obs_post = observed %>% filter(started == T)

  # baseline_best = read.csv(paste0(out_path, "_baseline_best_", projects[i], ".csv"), header = T)
  # baseline_loose = read.csv(paste0(out_path, "_baseline_loose_", projects[i], ".csv"), header = T)
  baseline_lagged = read.csv(paste0(out_path, "_baseline_lagged_new_", projects[i], ".csv"), header = T)

  # time_a = Sys.time()
  # pre_cf_c_loss_boot_list[[i]] = BootOut(type = "cf_c_loss", in_df = dplyr::select(obs_pre, c_loss))$summ %>%
  #   mutate(project = projects[i])
  # time_b = Sys.time()
  # cat("Project", i, "/", length(projects), "-", projects[i], "- pre_cf_c_loss_boot :", time_b - time_a, "\n")

  # pre_p_c_loss_boot_list[[i]] = BootOut(type = "p_c_loss", in_df = dplyr::select(obs_pre, t_loss))$summ %>%
  #   mutate(project = projects[i])
  # time_c = Sys.time()
  # cat("Project", i, "/", length(projects), "-", projects[i], "- pre_p_c_loss_boot :", time_c - time_b, "\n")

  # post_cf_c_loss_boot_list[[i]] = BootOut(type = "cf_c_loss", in_df = dplyr::select(obs_post, c_loss))$summ %>%
  #   mutate(project = projects[i])
  # time_d = Sys.time()
  # cat("Project", i, "/", length(projects), "-", projects[i], "- post_cf_c_loss_boot :", time_d - time_c, "\n")

  # post_p_c_loss_boot_list[[i]] = BootOut(type = "p_c_loss", in_df = dplyr::select(obs_post, t_loss))$summ %>%
  #   mutate(project = projects[i])
  # time_e = Sys.time()
  # cat("Project", i, "/", length(projects), "-", projects[i], "- post_p_c_loss_boot :", time_e - time_d, "\n")

  # observed_add_boot_out = BootOut(type = "additionality", in_df = dplyr::select(obs_post, additionality))
  # observed_add_boot_list[[i]] = observed_add_boot_out$summ %>%
  #   mutate(project = projects[i])
  # time_f = Sys.time()
  # cat("Project", i, "/", length(projects), "-", projects[i], "- add_boot :", time_f - time_e, "\n")

  # baseline_best_boot_out = BootOut(type = "best", in_df = dplyr::select(baseline_best, c_loss))
  # baseline_best_boot_list[[i]] = baseline_best_boot_out$summ %>%
  #     mutate(project = projects[i])
  # time_g = Sys.time()
  # cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_best_boot :", time_g - time_f, "\n")

  # baseline_loose_boot_list[[i]] = BootOut(type = "loose", in_df = dplyr::select(baseline_loose, c_loss))$summ %>%
  #     mutate(project = projects[i])
   time_h = Sys.time()
  # cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_loose_boot :", time_h - time_g, "\n")

  if(nrow(baseline_lagged) > 0) {
    baseline_lagged_boot_list[[i]] = BootOut(type = "lagged", in_df = dplyr::select(baseline_lagged, c_loss))$summ %>%
      mutate(project = projects[i])
  } else {
    baseline_lagged_boot_list[[i]] = data.frame(type = character(), mean = numeric(), ci_lower = numeric(), ci_upper = numeric(), project = character())
  }
  time_i = Sys.time()
  cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_lagged_boot :", time_i - time_h, "\n")

  # effectiveness = observed_add_boot_out$t / baseline_best_boot_out$t
  # effectiveness_list[[i]] = data.frame(project = projects[i],
  #                                       eff = mean(effectiveness, na.rm = T),
  #                                       eff_lower = quantile(effectiveness, probs = 0.025, na.rm = T),
  #                                       eff_upper = quantile(effectiveness, probs = 0.975, na.rm = T))
}

pre_cf_c_loss_boot_df = list_rbind(pre_cf_c_loss_boot_list)
pre_p_c_loss_boot_df = list_rbind(pre_p_c_loss_boot_list)
post_cf_c_loss_boot_df = list_rbind(post_cf_c_loss_boot_list)
post_p_c_loss_boot_df = list_rbind(post_p_c_loss_boot_list)
observed_add_boot_df = list_rbind(observed_add_boot_list)
baseline_best_boot_df = list_rbind(baseline_best_boot_list)
baseline_loose_boot_df = list_rbind(baseline_loose_boot_list)
baseline_lagged_boot_df = list_rbind(baseline_lagged_boot_list)

append_result = T
write.table(pre_cf_c_loss_boot_df, paste0(out_path, "_pre_cf_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_pre_cf_c_loss.csv")), row.names = F, append = append_result)
write.table(pre_p_c_loss_boot_df, paste0(out_path, "_pre_p_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_pre_p_c_loss.csv")), row.names = F, append = append_result)
write.table(post_cf_c_loss_boot_df, paste0(out_path, "_post_cf_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_post_cf_c_loss.csv")), row.names = F, append = append_result)
write.table(post_p_c_loss_boot_df, paste0(out_path, "_post_p_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_post_p_c_loss.csv")), row.names = F, append = append_result)
write.table(observed_add_boot_df, paste0(out_path, "_observed_add.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_observed_add.csv")), row.names = F, append = append_result)
write.table(baseline_best_boot_df, paste0(out_path, "_baseline_best_boot.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_baseline_best_boot.csv")), row.names = F, append = append_result)
write.table(baseline_loose_boot_df, paste0(out_path, "_baseline_loose_boot.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_baseline_loose_boot.csv")), row.names = F, append = append_result)
write.table(baseline_lagged_boot_df, paste0(out_path, "_baseline_lagged_new_boot.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_baseline_lagged_new_boot.csv")), row.names = F, append = append_result)


# D. Generate results ----

if(analysis_type == "control") {
  closs_placebo = read.csv(paste0(out_path, "_post_p_c_loss.csv"), header = T)
  closs_ante_regional = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
  closs_ante_project = read.csv(paste0(out_path, "_baseline_best_boot.csv"), header = T)
  closs_ante_historical = read.csv(paste0(out_path, "_baseline_lagged_new_boot.csv"), header = T)
  closs_post_counterfactual = read.csv(paste0(out_path, "_post_cf_c_loss.csv"), header = T)

  # Create plots
  out_regional = plotPlacebo(dat = rbind(closs_placebo, closs_ante_regional),
                             label_to_x = "p_c_loss", col = "#006CD1")
  out_project = plotPlacebo(dat = rbind(closs_placebo, closs_ante_project),
                            label_to_x = "p_c_loss", col = "#40B0A6")
  out_historical = plotPlacebo(dat = rbind(closs_placebo, closs_ante_historical),
                               label_to_x = "p_c_loss", col = "#CDAC60")
  out_post = plotPlacebo(dat = rbind(closs_placebo, closs_post_counterfactual),
                         label_to_x = "p_c_loss")

  p_regional = out_regional$p
  p_project = out_project$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_historical = out_historical$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_post = out_post$p

  # Create column labels
  col_ante = ggplot() +
    ggtitle(bquote("A." ~ italic("Ex ante") ~ "forecasts")) +
    theme_void() + theme(plot.title = element_text(size = 64, hjust = 0.5, margin = margin(t = 10)))
  col_post = ggplot() +
    ggtitle(bquote("B." ~ italic("Ex post") ~ "match")) +
    theme_void() + theme(plot.title = element_text(size = 64, hjust = 0.5, margin = margin(t = 10)))

  col_regional = ggplot() +
    ggtitle("Recent regional") +
    theme_void() + theme(plot.title = element_text(size = 60, hjust = 0.5, margin = margin(t = 10)))
  col_project = ggplot() +
    ggtitle("Recent project") +
    theme_void() + theme(plot.title = element_text(size = 60, hjust = 0.5, margin = margin(t = 10)))
  col_historical = ggplot() +
    ggtitle("Time-shifted historical match") +
    theme_void() + theme(plot.title = element_text(size = 60, hjust = 0.5, margin = margin(t = 10)))

# Figure 3a. ex ante forecasts with placebo areas
  plots = (p_regional + p_project + p_historical) +
    plot_layout(nrow = 1, guides = "collect", axis_titles = "collect")
  cols = (col_regional + col_project + col_historical) +
    plot_layout(nrow = 1)
  plot_complete =
    col_ante / cols / plots +
    plot_layout(nrow = 3, heights = c(0.02, 0.01, 1))
  ggsave(paste0(fig_path, "figure3a_placebo_ex_ante_new.png"), width = 48, height = 20, units = "in")

# Figure 3b. ex post estimation with placebo areas
  plot_complete =
    col_post / p_post +
    plot_layout(nrow = 2, heights = c(0.02, 1))
  ggsave(paste0(fig_path, "figure3b_placebo_ex_post.png"), width = 18, height = 20, units = "in")

# Figure 4. t-test result

  # Create a named vector for labels of different types
  label_types = c(
    "Recent regional" = "Recent regional",
    "Recent project" = "Recent project",
    "Time-shifted historical match" = "Time-shifted\nhistorical match",
    "Ex post match" = expression(italic("Ex post") ~ "\nmatch")
  )
  colors = c("#006CD1", "#40B0A6", "#CDAC60", "black")

  t_out_list = list(out_regional$t, out_project$t, out_historical$t, out_post$t)
  t_out_df = lapply(t_out_list, function(x) {
    data.frame(estimate = as.numeric(x$estimate),
               ci_lower = as.numeric(x$conf.int[1]),
               ci_upper = as.numeric(x$conf.int[2]),
               pval = as.numeric(x$p.value))
    }) %>%
    list_rbind() %>%
    mutate(type = factor(names(label_types), levels = names(label_types)))

  p_t_out = ggplot(data = t_out_df, aes(x = type)) +
    geom_col(aes(y = estimate, fill = type)) +
    geom_hline(yintercept = 0, linewidth = 1, linetype = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.5) +
    geom_text(aes(y = pmax(0, ci_upper) + 0.05, label = ifelse(round(pval, 2) == 0, "< 0.01", round(pval, 2))), size = 5) +
    scale_fill_manual(values = colors, guide = NULL) +
    scale_x_discrete(labels = label_types) +
    labs(x = "Estimate type", y = "Estimate-to-observed difference") +
    theme_bw() +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.text.x = element_text(margin = margin(t = 30), hjust = 0.5, vjust = 0.5))
  ggsave(paste0(fig_path, "figure4_t_test_out.png"), width = 10, height = 10, units = "in")
  }


scale_color_manual(values = c(best = "#40B0A6", loose = "#006CD1", lagged = "#CDAC60")) +

# Figure in SI: ongoing projects
if(analysis_type == "ongoing") {
  closs_ante_regional = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
  closs_ante_project = read.csv(paste0(out_path, "_baseline_best_boot.csv"), header = T)
  closs_ante_historical = read.csv(paste0(out_path, "_baseline_lagged_boot.csv"), header = T)
  closs_post_counterfactual = read.csv(paste0(out_path, "_post_cf_c_loss.csv"), header = T)

  # Create plots
  out_regional = plotPlacebo(dat = rbind(closs_post_counterfactual, closs_ante_regional),
                             label_to_x = "cf_c_loss", col = "#006CD1")
  out_project = plotPlacebo(dat = rbind(closs_post_counterfactual, closs_ante_project),
                            label_to_x = "cf_c_loss", col = "#40B0A6")
  out_historical = plotPlacebo(dat = rbind(closs_post_counterfactual, closs_ante_historical),
                               label_to_x = "cf_c_loss", col = "#CDAC60")

  p_regional = out_regional$p
  p_project = out_project$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_historical = out_historical$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # Create column labels
  col_ante = ggplot() +
    ggtitle(bquote(italic("Ex ante") ~ "forecasts")) +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))

  col_regional = ggplot() +
    ggtitle("Recent regional") +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))
  col_project = ggplot() +
    ggtitle("Recent project") +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))
  col_historical = ggplot() +
    ggtitle("Time-shifted historical match") +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))

  plots = (p_regional + p_project + p_historical) +
    plot_layout(nrow = 1, guides = "collect", axis_titles = "collect")
  cols = (col_regional + col_project + col_historical) +
    plot_layout(nrow = 1)
  plot_complete =
    col_ante / cols / plots +
    plot_layout(nrow = 3, heights = c(0.02, 0.01, 1))
  ggsave(paste0(fig_path, "figure_s4_ongoing_projects_new.png"), width = 48, height = 16, units = "in")

# ---TO BE REMOVED AFTER THIS LINE--- #

  summ = list_rbind(list(post_cf_c_loss_boot_df,
                         baseline_best_boot_df,
                         baseline_loose_boot_df,
                         baseline_lagged_boot_df)) %>%
    left_join(project_var["ID"], by = join_by(project == ID))

  #fit LM while forcing through the origin (0, 0)
  summ_wide = summ %>%
    dplyr::select(type, mean, project) %>%
    pivot_wider(names_from = "type", values_from = "mean", id_expand = T)
  lm_best = lm(cf_c_loss ~ best - 1, data = summ_wide)
  lm_loose = lm(cf_c_loss ~ loose - 1, data = summ_wide)
  lm_lagged = lm(cf_c_loss ~ lagged - 1, data = summ_wide)

  #calculate R2
  r2 = c(summary(lm_best)$r.squared, summary(lm_loose)$r.squared, summary(lm_lagged)$r.squared) %>% round(., 3)

  #calculate MAE
  mae_val = c(mean(abs(summ_wide$cf_c_loss - summ_wide$best)),
              mean(abs(summ_wide$cf_c_loss - summ_wide$loose)),
              mean(abs(summ_wide$cf_c_loss - summ_wide$lagged), na.rm = T)) %>% round(., 2)

  #non-stationarity correction factor
  corr_fact = c(summary(lm_best)$coefficients[1],
               summary(lm_loose)$coefficients[1],
               summary(lm_lagged)$coefficients[1]) %>% round(., 2)
  corr_confint = data.frame(best = as.numeric(confint(lm_best)),
                            loose = as.numeric(confint(lm_loose)),
                            lagged = as.numeric(confint(lm_lagged))) %>% round(., 2)
  write.csv(data.frame(type = c("best", "loose", "lagged"), correction_factor = corr_fact), paste0(out_path, "_corr_fact.csv"))

  summ_corr = list_rbind(list(post_cf_c_loss_boot_df,
                              baseline_best_boot_df %>% mutate(across(mean:ci_upper, function(x) x * corr_fact[1])),
                              baseline_loose_boot_df %>% mutate(across(mean:ci_upper, function(x) x * corr_fact[2])),
                              baseline_lagged_boot_df %>% mutate(across(mean:ci_upper, function(x) x * corr_fact[3]))))
  summ_corr = summ_corr %>%
    left_join(project_var[c("ID")], by = join_by(project == ID))

  #fit LM while forcing through the origin (0, 0)
  summ_corr_wide = summ_corr %>%
    dplyr::select(type, mean, project) %>%
    pivot_wider(names_from = "type", values_from = "mean", id_expand = T)
  lm_best_corr = lm(cf_c_loss ~ best - 1, data = summ_corr_wide)
  lm_loose_corr = lm(cf_c_loss ~ loose - 1, data = summ_corr_wide)
  lm_lagged_corr = lm(cf_c_loss ~ lagged - 1, data = summ_corr_wide)

  #calculate MAE
  mae_corr = c(mean(abs(summ_corr_wide$cf_c_loss - summ_corr_wide$best)),
               mean(abs(summ_corr_wide$cf_c_loss - summ_corr_wide$loose)),
               mean(abs(summ_corr_wide$cf_c_loss - summ_corr_wide$lagged), na.rm = T)) %>% round(., 2)

  # Figure 5. show how baseline compares to counterfactual carbon loss in ongoing projects (before vs after correction)
  p1 = plotBaseline(dat = summ, baseline_used = "best", metrics = mae_val[1])
  p2 = plotBaseline(dat = summ, baseline_used = "loose", metrics = mae_val[2])
  p3 = plotBaseline(dat = summ, baseline_used = "lagged", metrics = mae_val[3])
  p4 = plotBaseline(dat = summ_corr, baseline_used = "best", metrics = c(mae_corr[1], corr_fact[1], corr_confint$best), corr = T)
  p5 = plotBaseline(dat = summ_corr, baseline_used = "loose", metrics = c(mae_corr[2], corr_fact[2], corr_confint$loose), corr = T)
  p6 = plotBaseline(dat = summ_corr, baseline_used = "lagged", metrics = c(mae_corr[3], corr_fact[3], corr_confint$lagged), corr = T)
  # Create column labels
  col1 = ggplot() + ggtitle("A. Close matching") + theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))
  col2 = ggplot() + ggtitle("B. Loose matching") + theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))
  col3 = ggplot() + ggtitle("C. Time-lagged matching") + theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))

  # Create row labels
  row1 = ggplot() + annotate("text", x = 1, y = 0.5, label = "Before correction", angle = 270, size = 15) + theme_void()
  row2 = ggplot() + annotate("text", x = 1, y = 0.5, label = "After correction", angle = 270, size = 15) + theme_void()

  # Create r2 labels
  r2_1 = ggplot() + annotate("text", x = 1, y = 0.5, label = bquote(R^2 * ": " * .(r2[1])), size = 10) + theme_void()
  r2_2 = ggplot() + annotate("text", x = 1, y = 0.5, label = bquote(R^2 * ": " * .(r2[2])), size = 10) + theme_void()
  r2_3 = ggplot() + annotate("text", x = 1, y = 0.5, label = bquote(R^2 * ": " * .(r2[3])), size = 10) + theme_void()


  # Combine plots with labels
  plots = (p1 + p2 + p3 + row1 + p4 + p5 + p6 + row2) +
    plot_layout(nrow = 2, axes = "collect", axis_titles = "collect", widths = c(1, 1, 1, 0.2)) #row first
  cols = (col1 + col2 + col3 + plot_spacer()) +
    plot_layout(nrow = 1, widths = c(1, 1, 1, 0.2))
  r2_lab = (r2_1 + r2_2 + r2_3 + plot_spacer()) +
    plot_layout(nrow = 1, widths = c(1, 1, 1, 0.2))
  plot_complete = #then column
    cols / r2_lab / plots +
    plot_layout(nrow = 3, heights = c(0.01, 0.05, 1))
  ggsave(paste0(fig_path, "figure3c_control.png"), width = 7500, height = 5500, units = "px")
  ggsave(paste0(fig_path, "figure4_ongoing.png"), width = 7500, height = 5500, units = "px")


  # Calculate project performance ratio
  eff_list = vector("list", length(projects))

  corr_fact_mean = summary(lm_best)$coefficients[1]
  corr_fact_se = summary(lm_best)$coefficients[2]

  for(i in seq_along(projects)) {
    area_i = area_ha_vec[i]
    obs_post = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T) %>%
      mutate(t_loss = t_loss / area_i,
             c_loss = c_loss / area_i,
             additionality = additionality / area_i) %>%
      filter(started == T)
    baseline_best = read.csv(paste0(out_path, "_baseline_best_", projects[i], ".csv"), header = T)

    observed_add_boot = BootOut(type = "additionality", in_df = dplyr::select(obs_post, additionality))$t
    baseline_best_boot = BootOut(type = "best", in_df = dplyr::select(baseline_best, c_loss))$t
    corr_fact_boot = rnorm(1000, corr_fact_mean, corr_fact_se)
    baseline_best_boot_corr = baseline_best_boot * corr_fact_boot

    eff_boot = observed_add_boot / baseline_best_boot
    eff_corr_boot = observed_add_boot / baseline_best_boot_corr
    eff_list[[i]] = data.frame(project = projects[i],
                               eff = mean(eff_boot, na.rm = T),
                               eff_lower = quantile(eff_boot, 0.025, na.rm = T),
                               eff_upper = quantile(eff_boot, 0.975, na.rm = T),
                               eff_corr = mean(eff_corr_boot, na.rm = T),
                               eff_corr_lower = quantile(eff_corr_boot, 0.025, na.rm = T),
                               eff_corr_upper = quantile(eff_corr_boot, 0.975, na.rm = T))
  }

  eff_df = list_rbind(eff_list)
  write.table(eff_df, paste0(out_path, "_effectiveness.csv"), sep = ",", row.names = F)

  eff_df = read.csv(paste0(out_path, "_effectiveness.csv"), header = T) %>%
    left_join(project_var[c("ID")], by = join_by(project == ID))

  eff_plot = eff_df %>%
    filter(eff > 0) %>%
    arrange(eff) %>%
    mutate(code = factor(code, levels = code))

  eff_med = median(eff_plot$eff_corr)
  median(eff_plot$eff_corr_lower)
  median(eff_plot$eff_corr_upper)

  eff_5perc = quantile(eff_plot$eff_corr, 0.05)
  quantile(eff_plot$eff_corr_lower, 0.05)
  quantile(eff_plot$eff_corr_upper, 0.05)

  p8 = ggplot(data = eff_plot) +
    geom_segment(aes(x = code, y = eff_corr_lower, yend = eff_corr_upper), color = "blue", linewidth = 2) +
    geom_rect(aes(xmin = as.numeric(code) - 0.4,
                  xmax = as.numeric(code) + 0.4,
                  ymin = 0.1, ymax = eff_corr), fill = "lightblue") +
    geom_text(aes(x = code, y = eff_corr * 0.85, label = code), color = "darkblue", size = 10) +
    geom_hline(yintercept = c(eff_5perc, eff_med, 1), linetype = 3, linewidth = 1.2) +
    scale_x_discrete(name = "", labels = NULL) +
    scale_y_continuous(limits = c(0.1, 20),
                       breaks = c(0.1, 0.2, 0.5, 1, 1.5, 2, 3, 5, 10, 15, 20),
                       expand = c(0, 0),
                       transform = scales::transform_log10()) +
    labs(title = "",
         x = "Project code",
         y = "Project performance ratio") +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_blank(),
          plot.title = element_text(size = 48, hjust = 0.5, margin = margin(b = 10)),
          axis.title = element_text(size = 40),
          axis.text = element_text(size = 36),
          axis.title.y = element_text(margin = margin(r = 10)),
          axis.text.y = element_text(margin = margin(r = 10)),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth = 2),
          axis.ticks.length.y = unit(.5, "cm"),
          legend.position = "none")

  (p7 + p8) +
    plot_layout(axes = "collect", axis_titles = "collect")
  p8
  ggsave(paste0(fig_path, "figure5_effectiveness_after.png"), width = 4000, height = 4000, units = "px")


}


#OPTIONAL output: only basic variables, additionality distribution data to send to Ofir
if(ofir) {
    write.table(project_var %>% dplyr::select(project, t0, country, area_ha),
                paste0(out_path, "_project_var_basic.csv"), sep = ",", row.names = F)

    additionality_distribution = lapply(seq_along(projects), function(i) {
        additionality = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T)

        additionality %>%
            filter(started) %>%
            dplyr::select(year, additionality, pair) %>%
            mutate(project = projects[i])
    }) %>%
        list_rbind()
    write.csv(additionality_distribution, paste0("/maps/epr26/tmf_pipe_out/additionality_distribution.csv"), row.names = F)
}


#Compare with Jody's output
analysis_type = "ongoing"
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored
observed_add_boot_df = read.csv(paste0(out_path, "_observed_add.csv"), header = T)

observed_add_jody = read.csv("/maps/epr26/ex_ante_forecast_out/merged_additionality_by_jody.csv", header = T)
year_max = observed_add_jody[nrow(observed_add_jody), ]$year

observed_add_jody_mean = observed_add_jody %>%
  filter(year == year_max) %>%
  rename_with(function(x) gsub("X", "", gsub("_additionality", "", x))) %>%
  pivot_longer(cols = !year, names_to = "ID", values_to = "cumul_add") %>%
  mutate(ID = as.numeric(ID)) %>%
  left_join(x = ., y = project_var %>% dplyr::select("ID", "t0", "area_ha"), by = "ID") %>%
  filter(!is.na(t0)) %>%
  mutate(additionality_jody = cumul_add / ((year - t0) * area_ha))

observed_add_boot_compare = observed_add_boot_df %>%
  left_join(x = ., y = project_var %>% dplyr::select("ID", "t0", "area_ha"), by = join_by(project == ID)) %>%
  left_join(x = ., y = observed_add_jody_mean %>% dplyr::select("ID", "additionality_jody"), by = join_by(project == ID))
write.csv(observed_add_boot_compare, paste0("/maps/epr26/ex_ante_forecast_out/additionality_compare.csv"), row.names = F)
