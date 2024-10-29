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
FindFiles = function(dir, pattern, full = F, negate = F) {
  file_matched = list.files(dir, full = full) %>% str_subset(pattern, negate = negate)
  if(length(file_matched) == 0) {
    return(NA)
  } else {
    return(file_matched)
  }
}

BootSumm = function(type, in_df, boot_n = 1000) {
  boot_out = boot::boot(data = in_df,
                        statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                        R = boot_n)
  boot_ci = boot::boot.ci(boot.out = boot_out, type = "norm")
  boot_summ = data.frame(type = type,
                         mean = mean(as.vector(boot_out$t)),
                         ci_lower = boot_ci$normal[2],
                         ci_upper = boot_ci$normal[3])
  return(boot_summ)
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
analysis_type = "offset"
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
  projects = map(exclude_strings, function(x) FindFiles(project_dir, x, negate = T)) %>%
    reduce(intersect)

} else if(analysis_type == "offset") {
  project_dir = "/maps/epr26/tmf_pipe_out_offset/" #define where implementation code results are stored

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
  asn_info = read.csv("asian_tropics_controls.csv")
  sa_info = read.csv("neotropics_controls.csv")
  af_info = read.csv("afrotropics_controls.csv")
  proj_info = do.call(bind_rows, list(asn_info, sa_info, af_info)) %>%
    mutate(ID = proj_name, COUNTRY = name, t0 = 2011) %>%
    filter(ID != "") %>%
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

#only keep projects who have finished running ("additionality.csv" exists)
if(analysis_type == "offset") {
  done_vec = sapply(projects, function(x) FindFiles(paste0(project_dir, x, "/pairs"), ".parquet") %>% length() == 200)
} else {
  done_vec = sapply(projects, function(x) FindFiles(paste0(project_dir, x), "additionality.csv") %>% length() > 0)
}

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
#block_paths = rep(NA, length(projects))
acd_paths = rep(NA, length(projects))
for(i in seq_along(projects)) {
  project_out_dir = paste0(project_dir, projects[i])
  k_paths[i] = FindFiles(project_out_dir, "k.parquet", full = T)
  m_paths[i] = FindFiles(project_out_dir, "matches.parquet", full = T)
#  block_paths[i] = FindFiles(project_out_dir, "block_baseline.parquet", full = T)
  acd_paths[i] = FindFiles(project_out_dir, "carbon-density.csv", full = T)
}

polygon_paths = paste0(polygon_dir, projects, ".geojson")
proj_ID = gsub("(?<=[0-9])a$", "", projects, perl = T)
project_var = proj_info[match(proj_ID, proj_info$ID), ]
country = project_var$COUNTRY
t0_vec = project_var$t0


# 0b. User-defined input variables ----

#projects = NULL
#pair_dirs = NULL
#k_paths = NULL
#m_paths = NULL
#block_paths = NULL
#acd_paths = NULL
#polygon_paths = NULL
#country = NULL
#t0 = NULL
#proj_ID = NULL
#out_path = NULL


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
project_var$acd_undisturbed = sapply(acd_list, function(x) filter(x, land.use.class == 1)$carbon.density)

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

# #list containing block baseline of every project
# set_block = lapply(seq_along(projects), function(i) {
#   if(analysis_type == "offset") {
#     luc_t_10 = paste0("luc_", t0_vec[i] - 20)
#     luc_t0 = paste0("luc_", t0_vec[i] - 10)
#   } else {
#     luc_t_10 = paste0("luc_", t0_vec[i] - 10)
#     luc_t0 = paste0("luc_", t0_vec[i])
#   }
#   if(!is.na(block_paths[i])) {
#       read_parquet(block_paths[i]) %>%
#         rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
#         as.data.frame()
#   } else {
#     return(NULL)
#   }
# })


#Output: project-level variables
write.csv(project_var, paste0(out_path, "_project_var.csv"), row.names = F)
#project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)


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
  pair_paths = FindFiles(pair_dirs[i], ".parquet", full = T)
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

  additionality_estimates = lapply(pairs_out, function(x) x$out_df) %>%
    do.call(dplyr::bind_rows, .) %>%
    mutate(started = ifelse(year > t0, T, F)) %>%
    mutate(project = projects[i])

  b = Sys.time()
  cat("Project", i, "/", length(projects), "-", projects[i], ":", b - a, "\n")
  return(list(additionality_estimates = additionality_estimates, baseline_best = baseline_best))
})

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
  baseline_best_boot_df_i = BootSumm(type = "best", in_df = baseline_best_i) %>%
    mutate(project = projects[i])
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
    baseline_loose_boot_df_i = BootSumm(type = "loose", in_df = baseline_loose_i) %>%
      mutate(project = projects[i])
    b = Sys.time()
    cat(projects[i], ":", b - a, "\n")

    return(baseline_loose_boot_df_i)
  }) %>%
    do.call(rbind, .)
  write.csv(baseline_loose_boot_df, paste0(out_path, "_baseline_loose_boot.csv"), row.names = F)
  #baseline_loose_boot_df = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
}

# E. Get block-based baseline (if analysis_type == "ongoing") ----
if(analysis_type != "offset") {
  #  baseline_block = lapply(seq_along(projects), function(i) {
  #   a = Sys.time()
  #   if(is.null(set_block[[i]])) return(NULL)
    
  #   acd_i = acd_list[[i]]
  #   block = set_block[[i]] %>%
  #     filter(elevation >= 0 & access >= 0 & slope >= 0) %>%
  #     dplyr::select(-starts_with("luc_")) %>%
  #     dplyr::select(-starts_with("cpc")) %>%
  #     as.data.frame()
  #   if(nrow(block) > 250000) block = block[sample(nrow(block), 250000), ]

  #   baseline_i = block %>%
  #     mutate(acd10 = acd_i$carbon.density[match(luc10, acd_i$land.use.class)],
  #           acd0 = acd_i$carbon.density[match(luc0, acd_i$land.use.class)],
  #           c_loss = (acd10 - acd0) / 10)
  #   b = Sys.time()
  #   cat(projects[i], ":", b - a, "\n")
  #   return(baseline_i)
  # })
  # names(baseline_block) = projects
  # saveRDS(baseline_block, paste0(out_path, "_baseline_block.rds"))
  # #baseline_loose = read_rds(paste0(out_path, "_baseline_loose.rds"))

  # #bootstrap baselines
  # baseline_block_boot_df = lapply(seq_along(projects), function(i) {
  #   a = Sys.time()
  #   baseline_block_i = baseline_block[[i]] %>% dplyr::select(c_loss)
  #   baseline_block_boot_df_i = BootSumm(type = "block", in_df = baseline_block_i) %>%
  #     mutate(project = projects[i])
  #   b = Sys.time()
  #   cat(projects[i], ":", b - a, "\n")

  #   return(baseline_block_boot_df_i)
  # }) %>%
  #   do.call(rbind, .)
  # write.csv(baseline_block_boot_df, paste0(out_path, "_baseline_block_boot.csv"), row.names = F)
  # #baseline_block_boot_df = read.csv(paste0(out_path, "_baseline_block_boot.csv"), header = T)
}
