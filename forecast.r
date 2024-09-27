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

#Load country and t0 info (csv file copied from Tom's directory)
proj_info = read.csv("proj_meta.csv") %>%
  dplyr::select(ID, COUNTRY, t0)
ac_info = data.frame(ID = paste0("ac", sprintf("%02d", c(1:6))), COUNTRY = "Brazil", t0 = 2021)
control_as_info = read.csv("as_info.csv") %>%
  mutate(ID = paste0("as", sprintf("%02d", index)), t0 = 2011) %>%
  dplyr::select(ID, COUNTRY, t0)
proj_info_all = do.call(rbind, list(proj_info, ac_info, control_as_info))

#Based on analysis type, find project names (file prefixes) and store in vector "projects"
project_dir = "/maps/epr26/tmf_pipe_out/" #define where implementation code results are stored
polygon_dir = "/maps/epr26/tmf-data/projects/" #define where polygons are stored
analysis_type = "ongoing" #"ongoing": ongoing REDD+ projects; "control": control areas; "ac": Amazonian Collective
if(analysis_type == "ongoing") {

  exclude_strings = c("slopes", "elevation", "srtm", "ac", "as", "\\.", "\\_", "0000", "9999")
  projects = map(exclude_strings, function(x) str_subset(string = list.files(project_dir), pattern = x, negate = T)) %>%
    reduce(intersect)

  #only keep projects who have finished running ("additionality.csv" exists)
  done_vec = sapply(projects, function(x) list.files(paste0(project_dir, x)) %>% str_subset("additionality.csv") %>% length() > 0)

  #only keep projects with complete ACD values for LUC 1, 2, 3, and 4
  full_acd_vec = sapply(projects, function(x) {
    acd = read.csv(list.files(paste0(project_dir, x), full = T) %>% str_subset("carbon-density"))
    Reduce("&", 1:4 %in% acd$land.use.class)
  })

  projects_df = data.frame(project = projects, done = done_vec, full_acd = full_acd_vec)
  write.csv(projects_df, paste0(project_dir, "project_status.csv"), row.names = F)
  projects = subset(projects_df, done & full_acd)$project

  projects_to_exclude = c("674", "934", "2502", "1408")

  projects = projects[which(projects %in% projects_to_exclude == F)]

} else if(analysis_type == "control") {

  #control polygons in Asia
  projects = list.files(project_dir) %>% str_subset("as") %>% str_subset("\\.", negate = T)

  #only keep projects who have finished running ("additionality.csv" exists)
  done_vec = sapply(projects, function(x) list.files(paste0(project_dir, x)) %>% str_subset("additionality.csv") %>% length() > 0)

  #only keep projects with complete ACD values for LUC 1, 2, 3, and 4
  full_acd_vec = sapply(projects, function(x) {
    acd = read.csv(list.files(paste0(project_dir, x), full = T) %>% str_subset("carbon-density"))
    Reduce("&", 1:4 %in% acd$land.use.class)
  })

  projects_df = data.frame(project = projects, done = done_vec, full_acd = full_acd_vec)
  projects = subset(projects_df, done & full_acd)$project

} else if(analysis_type == "ac") {

  #Amazonian Collective polygons
  projects = paste0("ac", sprintf("%02d", c(1:6)))

}

#Produce input variables needed for the analysis
pair_dirs = paste0(project_dir, projects, "/pairs/")
k_paths = list.files(paste0(project_dir, projects), full = T) %>% str_subset("k.parquet")
m_paths = list.files(paste0(project_dir, projects), full = T) %>% str_subset("matches.parquet")
acd_paths = list.files(paste0(project_dir, projects), full = T) %>% str_subset("carbon-density.csv")
polygon_paths = paste0(polygon_dir, projects, ".geojson")
proj_ID = gsub("(?<=[0-9])a$", "", projects, perl = T)
country = proj_info_all[match(proj_ID, proj_info_all$ID), ]$COUNTRY
t0 = proj_info_all[match(proj_ID, proj_info_all$ID), ]$t0
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type)
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
area_ha = sapply(seq_along(projects), function(i) {
  area_ha = st_read(polygon_paths[i]) %>%
    st_make_valid() %>%
    st_union() %>%
    st_transform(4326) %>%
    st_area() %>% #area in m^2
    set_units(ha) #convert into hectares
  return(area_ha)
})

#list containing data frame of ACD (MgC/ha) per LUC of every project
acd = lapply(seq_along(projects), function(i) {
  acd_i = read.csv(acd_paths[i])
  for(class in 1:6) {
    if(class %in% acd_i$land.use.class == F) acd_i = rbind(acd_i, c(class, NA))
  }
  return(acd_i)
})
acd_undisturbed = sapply(acd, function(x) filter(x, land.use.class == 1)$carbon.density)

#list containing set K of every project
setK = lapply(seq_along(projects), function(i) {
  luc_t_10 = paste0("luc_", t0[i] - 10)
  luc_t0 = paste0("luc_", t0[i])
  read_parquet(k_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
    as.data.frame()
})

#list containing set M of every project
setM = lapply(seq_along(projects), function(i) {
  luc_t_10 = paste0("luc_", t0[i] - 10)
  luc_t0 = paste0("luc_", t0[i])
  read_parquet(m_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
    as.data.frame()
})

#data frame of project-level variables
project_var = data.frame(project = proj_ID, t0 = t0, country = country, area_ha = area_ha, acd_undisturbed = acd_undisturbed)


# B. Get observed additionality ----
#additionality_out = mclapply(seq_along(projects), mc.cores = 15, function(i) {
additionality_out = lapply(seq_along(projects), function(i) { #mclapply() does not work on Windows
  a = Sys.time()

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
                      t0 = t0[i], area_ha = area_ha[i], acd = acd[[i]], pair_id = j)
  })

  baseline_best = lapply(pairs_out, function(x) {
    filter(x$out_df, year <= t0[i] & year > t0[i] - 10) %>%
      dplyr::select(c_loss) %>%
      mutate(c_loss = c_loss / area_ha[1])
  }) %>%
    do.call(rbind, .)

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
    mutate(started = ifelse(year > t0[i], T, F)) %>%
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


# OPTIONAL output: only basic variables
write.table(project_var %>% dplyr::select(project, t0, country, area_ha),
             paste0(out_path, "_project_var_basic.csv"), sep = ",", row.names = F)

#Output: additionality time series
additionality_estimates = lapply(additionality_out, function(x) x$additionality_estimates)
saveRDS(additionality_estimates, paste0(out_path, "_additionality_estimates.rds"))
#additionality_estimates = read_rds(paste0(out_path, "_additionality_estimates.rds"))

additionality_df = do.call(rbind, additionality_estimates)
write.csv(additionality_df, paste0(out_path, "_additionality_estimates.csv"), row.names = F)
#additionality_df = read.csv(paste0(out_path, "_additionality_estimates.csv"), header = T)


#OPTIONAL output: additionality distribution data to send to Ofir
additionality_distribution = lapply(seq_along(projects), function(i) {
  additionality_estimates[[i]] %>%
    filter(started) %>%
    dplyr::select(year, additionality, pair) %>%
    mutate(project = projects[i])
}) %>%
  do.call(rbind, .)
write.csv(additionality_distribution, paste0("/maps/epr26/tmf_pipe_out/additionality_distribution.csv"), row.names = F)


# C1. Get best-matched baseline ----
#get baseline
baseline_best = lapply(additionality_out, function(x) x$baseline_best)
saveRDS(baseline_best, paste0(out_path, "_baseline_best.rds"))
#baseline_best = read_rds(paste0(out_path, "_baseline_best.rds"))

#bootstrap C loss values
baseline_best_boot = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  baseline_i = baseline_best[[i]] %>% dplyr::select(c_loss)
  bootstrap_i = BootDF(type = projects[i], in_df = baseline_i)

  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")
  return(bootstrap_i)
})
names(baseline_best_boot) = projects
saveRDS(baseline_best_boot, paste0(out_path, "_baseline_best_boot.rds"))
#baseline_best_boot = read_rds(paste0(out_path, "_baseline_best_boot.rds"))

#get mean bootstrapped C loss values
baseline_best_boot_df = baseline_best_boot %>%
  do.call(rbind, .) %>%
  group_by(type) %>%
  summarise(mean = mean(val), .groups = "drop")


# C2. Get loosely-matched baseline ----
#get baseline
baseline_loose = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  acd_i = acd[[i]]
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

#bootstrap C loss values
baseline_loose_boot = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  baseline_i = baseline_loose[[i]] %>% dplyr::select(c_loss)
  bootstrap_i = BootDF(type = projects[i], in_df = baseline_i)

  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")
  return(bootstrap_i)
})
names(baseline_loose_boot) = projects
saveRDS(baseline_loose_boot, paste0(out_path, "_baseline_loose_boot.rds"))
#baseline_loose_boot = read_rds(paste0(out_path, "_baseline_loose_boot.rds"))

#get mean bootstrapped C loss values
baseline_loose_boot_df = baseline_loose_boot %>%
  do.call(rbind, .) %>%
  group_by(type) %>%
  summarise(mean = mean(val), .groups = "drop")


# D. Non-project control areas ----

#bootstrapping mean of observed C loss values for control areas (target and counterfactual pixels)
c_loss_obs_control_boot = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  area_i = project_var$area_ha[i]
  dat_i = additionality_estimates[[i]] %>%
    mutate(t_loss = t_loss / area_i,
           c_loss = t_loss / area_i,
           additionality = additionality / area_i)

  p_loss_pre_i = subset(dat_i, started == F) %>% dplyr::select(t_loss)
  p_loss_pre_boot_i = BootDF(type = "t_loss_pre", in_df = p_loss_pre_i)

  cf_loss_pre_i = subset(dat_i, started == F) %>% dplyr::select(c_loss)
  cf_loss_pre_boot_i = BootDF(type = "cf_loss_pre", in_df = cf_loss_pre_i)

  p_loss_post_i = subset(dat_i, started == T) %>% dplyr::select(t_loss)
  p_loss_post_boot_i = BootDF(type = "p_loss_post", in_df = p_loss_post_i)

  cf_loss_post_i = subset(dat_i, started == T) %>% dplyr::select(c_loss)
  cf_loss_post_boot_i = BootDF(type = "cf_loss_post", in_df = cf_loss_post_i)

  add_i = subset(dat_i, started == T) %>% dplyr::select(additionality)
  add_boot_i = BootDF(type = "additionality", in_df = add_i)

  b = Sys.time() # less than one second per run

  bootstrap_i = rbind(p_loss_pre_boot_i, cf_loss_pre_boot_i, p_loss_post_boot_i, cf_loss_post_boot_i, add_boot_i) %>%
    mutate(project = projects[i])
  cat(projects[i], ":", b - a, "\n")
  return(bootstrap_i)
}) %>%
  do.call(rbind, .)

#Output: bootstrapped observed carbon loss values for control areas (project and counterfactual pixels)
saveRDS(c_loss_obs_control_boot, paste0(out_path, "_c_loss_obs_control_boot.rds"))
#c_loss_obs_control_boot = read_rds(paste0(out_path, "_c_loss_obs_control_boot.rds"))

#Output: bootstrapped mean of C loss values of control area (target/counterfactual pixels) and baseline
c_loss_obs_control_boot_df = c_loss_obs_control_boot %>%
  group_by(type, project) %>%
  summarise(mean = mean(val), .groups = "drop") %>%
  pivot_wider(names_from = "type", values_from = "mean", id_expand = T) %>%
  mutate(baseline_best = baseline_best_boot_df$mean, #baseline carbon loss
         baseline_loose = baseline_loose_boot_df$mean)

#Figure 4. show that matching works
ggplot(data = c_loss_control_df, aes(x = t_loss_pre, y = cf_loss_pre)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_x_continuous(limits = c(0, 2.1), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(0, 2.1), expand = c(0, 0)) +
  labs(x = "Mean carbon loss in project pixels (MgC/ha/yr)",
       y = "Mean carbon loss in counterfactual pixels (MgC/ha/yr)",
       title = "Before t0 in control areas") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure4_c_loss_t_vs_cf_pre.png"), width = 3000, height = 3000, units = "px")

#expression(paste0("Mean carbon loss in target pixels before ", t[0], "(MgC/ha/yr)"))

#Figure 5. show that there is no bias in additionality estimation
ggplot(data = c_loss_control_df, aes(x = t_loss_post, y = cf_loss_post)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_x_continuous(limits = c(0, 2.1), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(0, 2.1), expand = c(0, 0)) +
  labs(x = "Mean carbon loss in project pixels (MgC/ha/yr)",
       y = "Mean carbon loss in counterfactual pixels (MgC/ha/yr)",
       title = "After t0 in control areas") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure5_c_loss_t_vs_cf_post.png"), width = 2500, height = 2500, units = "px")

#Figure 6. show how baseline compares to additionality
ggplot(data = c_loss_control_df, aes(x = baseline_best, y = additionality)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  labs(x = "Mean best-matched baseline carbon loss (MgC/ha/yr)",
       y = "Mean additionality (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure6a_baseline_best_vs_add.png"), width = 2500, height = 2500, units = "px")

ggplot(data = c_loss_control_df, aes(x = baseline_loose, y = additionality)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  labs(x = "Mean loosely-matched baseline carbon loss (MgC/ha/yr)",
       y = "Mean additionality (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure6a_baseline_loose_vs_add.png"), width = 2500, height = 2500, units = "px")


# E. Ongoing project areas ----

#bootstrapped mean counterfactual C loss and mean additionality for ongoing projects
c_loss_ongoing_boot = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  area_i = project_var$area_ha[i]
  cf_loss_obs_i = additionality_estimates[[i]] %>%
    dplyr::select(c_loss) %>%
    mutate(c_loss = c_loss / area_i)
  cf_loss_obs_boot_i = BootDF(type = "cf_c_loss_obs", in_df = cf_loss_obs_i)


  add_i = additionality_estimates[[i]] %>%
    dplyr::select(additionality) %>%
    mutate(additionality = additionality / area_i)
  add_boot_i = BootDF(type = "additionality", in_df = add_i)
  b = Sys.time()

  bootstrap_i = rbind(cf_loss_obs_boot_i, add_boot_i) %>%
    mutate(project = projects[i])
  cat(projects[i], ":", b - a, "\n")
  return(bootstrap_i)
}) %>%
  do.call(rbind, .)

#Output: bootstrapped mean counterfactual C loss for ongoing projects
write.csv(c_loss_ongoing_boot, paste0(out_path, "_c_loss_ongoing_boot.csv"), row.names = F)
#c_loss_ongoing_boot = read.csv(paste0(out_path, "_c_loss_ongoing_boot.csv"), header = T)

#Output: bootstrapped mean observed counterfactual C loss and additionality and baseline
c_loss_ongoing_df = c_loss_ongoing_boot %>%
  group_by(type, project) %>%
  summarise(mean = mean(val), .groups = "drop") %>%
  pivot_wider(names_from = "type", values_from = "mean", id_expand = T) %>%
  mutate(baseline_best = baseline_best_boot_df$mean,
         baseline_loose = baseline_loose_boot_df$mean)
write.csv(c_loss_ongoing_df, paste0(out_path, "_c_loss_ongoing_df.csv"), row.names = F)
#c_loss_ongoing_df = read.csv(paste0(out_path, "_c_loss_ongoing_df.csv"), header = T)

#Figure 7. show how best-matched baseline compares to counterfactual carbon loss
ggplot(data = c_loss_ongoing_df, aes(x = baseline_best, y = cf_c_loss_obs)) +
  geom_point() +
  geom_text(aes(label = project), hjust = -0.1, size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-0.05, 1.5), expand = c(0.01, 0.01)) +
  labs(x = "Mean best-matched baseline carbon loss (MgC/ha/yr)",
       y = "Mean observed counterfactual carbon loss (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure7a_baseline_best_vs_c_loss_cf.png"), width = 3500, height = 3500, units = "px")

ggplot(data = c_loss_ongoing_df, aes(x = baseline_loose, y = cf_c_loss_obs)) +
  geom_point() +
  geom_text(aes(label = project), hjust = -0.1, size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-0.05, 1.5), expand = c(0.01, 0.01)) +
  labs(x = "Mean loosely-matched baseline carbon loss (MgC/ha/yr)",
       y = "Mean observed counterfactual carbon loss (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure7b_baseline_loose_vs_c_loss_cf.png"), width = 3500, height = 3500, units = "px")


#Figure 8. show how baseline compares to additionality
ggplot(data = c_loss_ongoing_df, aes(x = baseline_best, y = additionality)) +
  geom_point() +
  geom_text(aes(label = project), hjust = -0.1, size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-0.2, 0.4), expand = c(0.01, 0.01)) +
  labs(x = "Mean best-matched baseline carbon loss (MgC/ha/yr)",
       y = "Mean additionality (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure8a_baseline_best_vs_add.png"), width = 3500, height = 3500, units = "px")

ggplot(data = c_loss_ongoing_df, aes(x = baseline_loose, y = additionality)) +
  geom_point() +
  geom_text(aes(label = project), hjust = -0.1, size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-0.2, 0.4), expand = c(0.01, 0.01)) +
  labs(x = "Mean loosely-matched baseline carbon loss (MgC/ha/yr)",
       y = "Mean additionality (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_figure8b_baseline_loose_vs_add.png"), width = 3500, height = 3500, units = "px")
