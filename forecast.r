# 0. Setup ----
rm(list = ls())

#Load packages
# install.packages(c("arrow","configr", "tidyverse", "magrittr", "sf", "magrittr", "MatchIt",
#                    "rnaturalearthdata", "configr", "terra", "pbapply", "cleangeo", "doParallel",
#                    "foreach", "readr", "lwgeom", "rnaturalearth", "stars"), depends = TRUE)

library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators
library(MatchIt) #MatchIt::matchit
library(sf) #sf::st_area
library(ggpubr) #ggpubr::ggarrange
library(units) #units::set_units
library(arrow) #arrow::read_parquet
library(pryr) #pryr::object_size
library(cowplot)
library(boot) #boot::boot
#library(parallel) #parallel::mclapply

#Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = FALSE)

#Load pre-defined functions
source("functions.r") #cpc_rename, tmfemi_reformat, simulate_area_series, make_area_series, assess_balance, make_match_formula
source("AdditionalityPair.r")
source("PredictDefor.r")
BootDF = function(type, in_df, boot_n = 1000) {
  data.frame(type = type,
             val = as.vector(boot::boot(data = in_df,
                                        statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                                        R = boot_n)$t)) #R: number of bootstrapping repetitions
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
#9. OPTIONAL: proj_ID: ID used to identify country t0 from the proj_info_all data frame
# This is only used to remove the trailing "a" that I added to the filename of some shapefiles that needed fixing
# If unspecified, the projects variable will be used
#10. out_path: absolute paths of the directory where outputs are to be saved; include file prefix if desired


# 0a. E-Ping's workflow to obtain input variables ----

#Define analysis type
analysis_type = "ongoing"
#"ongoing": ongoing REDD+ projects (best-matched and loosely-matched baselines)
#"offset": ongoing REDD+ projects (time-offsetted best-matched baseline)
#"control": non-project polygons
#"ac": Amazonian Collective polygons

polygon_dir = "/maps/epr26/tmf-data/projects/" #where polygons are stored
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored
projects_to_exclude = c("674", "934", "2502", "1408") #which projects to exclude manually

#Based on analysis type, find project names and store in vector "projects", and project basic info and store in dataframe "proj_info"
if(analysis_type == "ongoing") {
  project_dir = "/maps/epr26/tmf_pipe_out/" #define where implementation code results are stored

  #Load basic information (csv file copied from Tom's directory)
  proj_info = read.csv("proj_meta.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #Find all project IDs
  exclude_strings = c("slopes", "elevation", "srtm", "ac", "as", "\\.", "\\_", "0000", "9999")
  projects = map(exclude_strings, function(x) str_subset(string = list.files(project_dir), pattern = x, negate = T)) %>%
    reduce(intersect)

} else if(analysis_type == "offset") {
  project_dir = "/maps/epr26/tmf_pipe_out_offset/" #define where implementation code results are stored

  #Load basic information (csv file copied from Tom's directory)
  proj_info = read.csv("proj_meta.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #Find all project IDs
  exclude_strings = c("slopes", "elevation", "srtm", "ac", "as", "\\.", "\\_", "0000", "9999")
  projects = map(exclude_strings, function(x) str_subset(string = list.files(project_dir), pattern = x, negate = T)) %>%
    reduce(intersect)

} else if(analysis_type == "control") {
  project_dir = "/maps/epr26/tmf_pipe_out_luc_t/" #define where implementation code results are stored

  #Load basic information
  asn_info = read.csv("asian_tropics_controls.csv")
  sa_info = read.csv("neotropics_controls.csv")
  af_info = read.csv("afrotropics_controls.csv")
  proj_info = do.call(bind_rows, list(asn_info, sa_info, af_info)) %>%
    mutate(ID = proj_name, COUNTRY = name, t0 = 2011) %>%
    filter(ID != "") %>%
    dplyr::select(ID, COUNTRY, t0)

  #Find all project IDs
  projects = list.files(project_dir) %>%
    map(c("asn", "af", "sa"), str_subset, string = .) %>%
    reduce(union) %>%
    str_subset("\\.", negate = T)

} else if(analysis_type == "ac") {
  project_dir = "/maps/epr26/tmf_pipe_out/" #define where implementation code results are stored

  #Load basic information
  proj_info = data.frame(ID = paste0("ac", sprintf("%02d", c(1:6))), COUNTRY = "Brazil", t0 = 2021)

  #Amazonian Collective polygons
  projects = paste0("ac", sprintf("%02d", c(1:6)))

}

#only keep projects who have finished running ("additionality.csv" exists)
if(analysis_type == "offset") {
  done_vec = sapply(projects, function(x) list.files(paste0(project_dir, x, "/pairs")) %>% str_subset(".parquet") %>% length() == 200)
} else {
  done_vec = sapply(projects, function(x) list.files(paste0(project_dir, x)) %>% str_subset("additionality.csv") %>% length() > 0)
}

#only keep projects with complete ACD values for LUC 1, 2, 3, and 4
full_acd_vec = sapply(projects, function(x) {
  acd = read.csv(list.files(paste0(project_dir, x), full = T) %>% str_subset("carbon-density"))
  acd = acd[, 1:2]
  colnames(acd) = c("land.use.class", "carbon.density")
  Reduce("&", 1:4 %in% acd$land.use.class)
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
  k_paths[i] = list.files(paste0(project_dir, projects[i]), full = T) %>% str_subset("k.parquet")
  m_paths[i] = list.files(paste0(project_dir, projects[i]), full = T) %>% str_subset("matches.parquet")
  acd_paths[i] = list.files(paste0(project_dir, projects[i]), full = T) %>% str_subset("carbon-density.csv")
}

polygon_paths = paste0(polygon_dir, projects, ".geojson")
proj_ID = gsub("(?<=[0-9])a$", "", projects, perl = T)
proj_info_selected = proj_info[match(proj_ID, proj_info$ID), ]
country = proj_info_selected$COUNTRY
t0_vec = proj_info_selected$t0
visualise = T #define whether one wants to generate plots (for the manuscript) or not


# 0b. User-defined input variables ----

#projects = NULL
#pair_dirs = NULL
#k_paths = NULL
#m_paths = NULL
#acd_paths = NULL
#polygon_paths = NULL
#country = NULL
#t0 = NULL
#proj_ID = NULL
#out_path = NULL
#visualise = FALSE

# A. Read data ----

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
proj_info_selected$area_ha = area_ha_vec
write.csv(proj_info_selected, paste0(out_path, "_project_selected_info.csv"), row.names = F)

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
acd_undisturbed = sapply(acd_list, function(x) filter(x, land.use.class == 1)$carbon.density)

#list containing set K of every project
setK = lapply(seq_along(projects), function(i) {
  luc_t_10 = paste0("luc_", t0_vec[i] - 10)
  luc_t0 = paste0("luc_", t0_vec[i])
  read_parquet(k_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
    as.data.frame()
})

#list containing set M of every project
setM = lapply(seq_along(projects), function(i) {
  if(analysis_type == "offset") {
    luc_t_10 = paste0("luc_", t0_vec[i] - 20)
    luc_t0 = paste0("luc_", t0_vec[i] - 10)
  } else {
    luc_t_10 = paste0("luc_", t0_vec[i] - 10)
    luc_t0 = paste0("luc_", t0_vec[i])
  }
  read_parquet(m_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
    as.data.frame()
})

#data frame of project-level variables
project_var = data.frame(project = proj_ID, t0 = t0_vec, country = country, area_ha = area_ha_vec, acd_undisturbed = acd_undisturbed)


# B. Get observed additionality ----
#additionality_df = read.csv(paste0(out_path, "_additionality_estimates.csv"), header = T)
#existing_projects = unique(additionality_df$project)
#projects = projects[projects %in% existing_projects == F]

#additionality_out = mclapply(seq_along(projects), mc.cores = 15, function(i) {
additionality_out = lapply(seq_along(projects), function(i) { #mclapply() does not work on Windows
  a = Sys.time()

  t0 = t0_vec[i]
  area_ha = area_ha_vec[i]
  offset = (analysis_type == "offset")

  #find paths to match and unmatached points in each sampled pairs
  pair_paths = list.files(pair_dirs[i], full = T) %>% str_subset(".parquet")
  matched_paths = pair_paths %>% str_subset("matchless", negate = T)
  matchless_paths = pair_paths %>% str_subset("matchless")

  #exit if no matches
  if(length(matched_paths) == 0) {
    return(list(pair_var = NULL, additionality_estimates = NULL))
  }

  #loop through all sampled pairs, get matched points and additionality series in each pair - this is the bottleneck
  pairs_out = lapply(seq_along(matched_paths), function(j) {
    AdditionalityPair(matched_path = matched_paths[j], matchless_path = matchless_paths[j],
                      k = dplyr::select(setK[[i]], c("lat", "lng", "k_ecoregion")),
                      matches = dplyr::select(setM[[i]], c("lat", "lng", "s_ecoregion")),
                      t0 = t0, area_ha = area_ha, acd = acd_list[[i]], pair_id = j, offset = offset)
  })

  if(offset) {
    baseline_best = lapply(pairs_out, function(x) {
      filter(x$out_df, year <= 0 & year > -10) %>%
        dplyr::select(c_loss) %>%
        mutate(c_loss = c_loss / area_ha)
    }) %>%
      do.call(rbind, .)
  } else {
    baseline_best = lapply(pairs_out, function(x) {
      filter(x$out_df, year <= t0 & year > t0 - 10) %>%
        dplyr::select(c_loss) %>%
        mutate(c_loss = c_loss / area_ha)
    }) %>%
      do.call(rbind, .)
  }

  #gather pair-level summary variables
  pair_ecoregion = lapply(pairs_out, function(x) x$pts_matched$ecoregion) %>%
    unlist() %>%
    table() %>%
    which.max() %>%
    names()

  pair_var = lapply(pairs_out, function(x) x$pair_var) %>%
    do.call(dplyr::bind_rows, .) %>%
    group_by(var) %>%
    summarise(min = min(val), median = median(val), max = max(val)) %>%
    pivot_longer(cols = min:max, names_to = "stat", values_to = "val") %>%
    mutate(var = paste0(var, "_", stat)) %>%
    dplyr::select(c(var, val)) %>%
    pivot_wider(names_from = "var", values_from = "val") %>%
    mutate(ecoregion = pair_ecoregion)

  additionality_estimates = lapply(pairs_out, function(x) x$out_df) %>%
    do.call(dplyr::bind_rows, .) %>%
    mutate(started = ifelse(year > t0, T, F)) %>%
    mutate(project = projects[i])

  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")
  return(list(pair_var = pair_var, additionality_estimates = additionality_estimates, baseline_best = baseline_best))
})

#Output: project-level variables
project_var = cbind(project_var,
                    lapply(additionality_out, function(x) x$pair_var) %>% do.call(dplyr::bind_rows, .))
#use bind_rows because of potentially different numbers of columns
write.csv(project_var, paste0(out_path, "_project_var.csv"), row.names = F)
#project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)

#Output: additionality time series
if(analysis_type != "offset") {
  additionality_estimates = lapply(additionality_out, function(x) x$additionality_estimates)
  additionality_df = do.call(rbind, additionality_estimates)
  write.csv(additionality_df, paste0(out_path, "_additionality_estimates.csv"), row.names = F)
  #additionality_df = read.csv(paste0(out_path, "_additionality_estimates.csv"), header = T)
}

#OPTIONAL output: only basic variables, additionality distribution data to send to Ofir
write.table(project_var %>% dplyr::select(project, t0, country, area_ha),
             paste0(out_path, "_project_var_basic.csv"), sep = ",", row.names = F)

additionality_distribution = lapply(seq_along(projects), function(i) {
  additionality_estimates[[i]] %>%
    filter(started) %>%
    dplyr::select(year, additionality, pair) %>%
    mutate(project = projects[i])
}) %>%
  do.call(rbind, .)
write.csv(additionality_distribution, paste0("/maps/epr26/tmf_pipe_out/additionality_distribution.csv"), row.names = F)


# C. Get best-matched baseline ----
#get baseline
baseline_best = lapply(additionality_out, function(x) x$baseline_best)
saveRDS(baseline_best, paste0(out_path, "_baseline_best.rds"))
#baseline_best = read_rds(paste0(out_path, "_baseline_best.rds"))

#bootstrap baselines
baseline_best_boot_df = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  baseline_best_i = baseline_best[[i]] %>% dplyr::select(c_loss)
  baseline_best_boot_i = boot::boot(data = baseline_best_i,
                                    statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                                    R = 1000)
  baseline_best_boot_ci = boot::boot.ci(boot.out = baseline_best_boot_i, type = "norm")
  baseline_best_boot_df_i = data.frame(project = projects[i],
                                       type = "best",
                                       mean = mean(as.vector(baseline_best_boot_i$t)),
                                       ci_lower = baseline_best_boot_ci$normal[2],
                                       ci_upper = baseline_best_boot_ci$normal[3])
  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")

  return(baseline_best_boot_df_i)
}) %>%
  do.call(rbind, .)
write.csv(baseline_best_boot_df, paste0(out_path, "_baseline_best_boot.csv"), row.names = F)
#baseline_best_boot_df = read.csv(paste0(out_path, "_baseline_best_boot.csv"), header = T)

# D. Get loosely-matched baseline (if analysis_type == "ongoing") ----
if(analysis_type != "offset") {
   baseline_loose = lapply(seq_along(projects), function(i) {
    a = Sys.time()
    acd_i = acd_list[[i]]
    M = setM[[i]] %>%
      dplyr::select(-starts_with("luc_")) %>%
      dplyr::select(-starts_with("cpc")) %>%
      as.data.frame()
    if(nrow(M) > 250000) M = M[sample(nrow(M), 250000), ]

    baseline_i = M %>%
      mutate(acd10 = acd_i$carbon.density[match(luc10, acd_i$land.use.class)],
            acd0 = acd_i$carbon.density[match(luc0, acd_i$land.use.class)],
            c_loss = (acd10 - acd0) / 10)
    b = Sys.time()
    cat(projects[i], ":", b - a, "\n")
    return(baseline_i)
  })
  names(baseline_loose) = projects
  saveRDS(baseline_loose, paste0(out_path, "_baseline_loose.rds"))
  #baseline_loose = read_rds(paste0(out_path, "_baseline_loose.rds"))

  #bootstrap baselines
  baseline_loose_boot_df = lapply(seq_along(projects), function(i) {
    a = Sys.time()
    baseline_loose_i = baseline_loose[[i]] %>% dplyr::select(c_loss)
    baseline_loose_boot_i = boot::boot(data = baseline_loose_i,
                                      statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                                      R = 1000)
    baseline_loose_boot_ci = boot::boot.ci(boot.out = baseline_loose_boot_i, type = "norm")
    baseline_loose_boot_df_i = data.frame(project = projects[i],
                                          type = "loose",
                                          mean = mean(as.vector(baseline_loose_boot_i$t)),
                                          ci_lower = baseline_loose_boot_ci$normal[2],
                                          ci_upper = baseline_loose_boot_ci$normal[3])
    b = Sys.time()
    cat(projects[i], ":", b - a, "\n")

    return(baseline_loose_boot_df_i)
  }) %>%
    do.call(rbind, .)
  write.csv(baseline_loose_boot_df, paste0(out_path, "_baseline_loose_boot.csv"), row.names = F)
  #baseline_loose_boot_df = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
}
