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
#library(microbenchmark) #microbenchmark::microbenchmark

#Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = FALSE)

#Load pre-defined functions
#fillNA = function(x) if (length(x) == 0) return(NA) else return(x) #turn empty object as NA
#findLUCC = function(x) sum(x %in% c("1_2", "1_3", "1_4")) / length(x) #find LUC transitions from undisturbed to disturbed

source("functions.r") #cpc_rename, tmfemi_reformat, simulate_area_series, make_area_series, assess_balance, make_match_formula
source("AdditionalityPair.r")
source("PredictDefor.r")

#Load user-defined functions that Tom wrote
#sapply(list.files(paste0("/home/tws36/4c_evaluations/R"), full = T, pattern = ".R$"), source)

#library(configr) #configr::read.config
# config = configr::read.config(paste0("/home/tws36/4c_evaluations/config/fixed_config_sherwood.ini"))
# config$USERPARAMS$data_path = "/maps/pf341/tom"
# write.config(config, "./config/fixed_config_tmp.ini") #error: permission denied
# source(paste0("/home/tws36/4c_evaluations/R/scripts/setup_and_reusable/load_config.R"))
# source("./R/scripts/0.2_load_project_details.R")


#Define input variables needed to read TMF implementation output and other data

# It requires the following input variables to read TMF implementation output and other data.
# All variables are vectors containing one value for each project to be analysed:

#1. projects: an index of all projects to be analysed; it could be the projects' VCS ID or customised (e.g. simply a series of integers)

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
#9. OPTIONAL: proj_name: full name of the project for readability (if unspecified, the projects variable will be used)
#10. out_path: absolute paths of the directory where outputs are to be saved; include file prefix if desired


# 0a. E-Ping's workflow to obtain input variables ----

#Pre-defined settings for input variables based on analysis type
analysis_type = "full" #"full", "grid", "ac", "control"
forecast = (analysis_type == "ac")
visualise = T #generate plots or not

#Load basic info (csv file copied from Tom's directory) for country and t0 input
proj_meta = read.csv(paste0("proj_meta.csv"))

if(analysis_type == "full") {

  project_dir = "/maps/epr26/tmf_pipe_out/" #new results from E-Ping's pipeline run
  projects = list.files(project_dir) %>% #full = T and basename() negates one another
    str_subset("\\.", negate = T) %>%
    str_subset("\\_", negate = T) %>%
    str_subset("ac", negate = T) %>%
    setdiff(c("0000", "9999")) #reserved for control and grid

  #only keep projects who have finished running ("additionality.csv" exists)
  done_id = sapply(projects, function(x) list.files(paste0(project_dir, x)) %>% str_subset("additionality.csv") %>% length() > 0)
  projects = projects[done_id]

  #only keep projects with complete ACD values for LUC 1, 2, 3, and 4
  full_acd_id = sapply(projects, function(x) {
    acd = read.csv(paste0(project_dir, x, "/", x, "carbon-density.csv"))
    Reduce("&", 1:4 %in% acd$land.use.class)
  })
  projects = projects[full_acd_id]
  #1399 and 1408 might potentially be excluded due to weird results: to be investigated

  in_paths = paste0(project_dir, projects, "/", projects)

  pair_dirs = paste0(project_dir, projects, "/pairs/")
  k_paths = paste0(in_paths, "k.parquet")
  m_paths = paste0(in_paths, "matches.parquet")
  acd_paths = paste0(in_paths, "carbon-density.csv")
  polygon_paths = paste0("/maps/epr26/tmf-data/projects/", projects, ".geojson")
  country = proj_meta[match(str_replace(projects, "a", ""), proj_meta$ID), ]$COUNTRY
  t0 = proj_meta[match(str_replace(projects, "a", ""), proj_meta$ID), ]$t0
  proj_name = str_replace(projects, "a", "")

} else if(analysis_type == "grid") {

  project_dir = "/maps/epr26/tmf_pipe_out/1201_grid/"
  projects = 1:49 #27, 31 with no matches
  in_paths = paste0(project_dir, projects, "/1201_", projects)

  pair_dirs = paste0(project_dir, projects, "/pairs/")
  k_paths = paste0(in_paths, "k.parquet")
  m_paths = paste0(in_paths, "matches.parquet")
  acd_paths = paste0(in_paths, "carbon-density.csv")
  polygon_paths = paste0("/maps/epr26/tmf-data-grid/1201/1201_", projects, ".geojson")
  country = filter(proj_meta, ID == "1201")$COUNTRY
  t0 = filter(proj_meta, ID == "1201")$t0
  proj_name = paste0("1201_", projects)

} else if(analysis_type == "control") {

  project_dir = "/maps/epr26/tmf_pipe_out/0000_grid/"
  projects = c(2:5, 7, 8, 10)
  in_paths = paste0(project_dir, projects, "/0000_", projects)

  pair_dirs = paste0(project_dir, projects, "/pairs/")
  k_paths = paste0(in_paths, "k.parquet")
  m_paths = paste0(in_paths, "matches.parquet")
  acd_paths = paste0(in_paths, "carbon-density.csv")
  polygon_paths = paste0("/maps/epr26/tmf-data-grid/0000/0000_", projects, ".geojson")
  country = "Brazil"
  t0 = 2011
  proj_name = paste0("0000_", projects)

} else if(analysis_type == "ac") {

  project_dir = "/maps/epr26/tmf_pipe_out/"
  projects = list.files(project_dir) %>%
    str_subset("ac\\d\\d")
  in_paths = paste0(project_dir, projects, "/", projects)

  pair_dirs = paste0(project_dir, projects, "/pairs/")
  k_paths = paste0(in_paths, "k.parquet")
  m_paths = paste0(in_paths, "matches.parquet")
  acd_paths = paste0(in_paths, "carbon-density.csv")
  polygon_paths = paste0("/maps/epr26/tmf-data/projects/", projects, ".geojson")
  country = "Brazil"
  t0 = 2021
  proj_name = projects

}

out_path = paste0("/maps/epr26/tmf_pipe_out/out_",
                  ifelse(analysis_type == "grid", "grid_1201", analysis_type))


# 0b. User-defined input variables ----

#projects = NULL
#pair_dirs = NULL
#k_paths = NULL
#m_paths = NULL
#acd_paths = NULL
#polygon_paths = NULL
#country = NULL
#t0 = NULL
#proj_name = NULL
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
project_var = data.frame(project = proj_name, t0 = t0, country = country, area_ha = area_ha, acd_undisturbed = acd_undisturbed)


# B. Obtain observed additionality ----
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

  #loop through all sampled pairs, get matched points and additionality series in each pair
  pairs_out = lapply(seq_along(matched_paths), function(j) {
    AdditionalityPair(matched_path = matched_paths[j], matchless_path = matchless_paths[j],
                      k = dplyr::select(setK[[i]], c("lat", "lng", "k_ecoregion")),
                      matches = dplyr::select(setM[[i]], c("lat", "lng", "s_ecoregion")),
                      t0 = t0[i], area_ha = area_ha[i], acd = acd[[i]], pair_id = j)
  })

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
    mutate(started = ifelse(year > t0[i], T, F))

  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")
  return(list(pair_var = pair_var, additionality_estimates = additionality_estimates))
})

#Output: project-level variables
project_var = cbind(project_var,
                    lapply(additionality_out, function(x) x$pair_var) %>% do.call(dplyr::bind_rows, .))
#use bind_rows because of potentially different numbers of columns
write.table(project_var, paste0(out_path, "_project_var.csv"), sep = ",", row.names = F)
#project_var = read.table(paste0(out_path, "_project_var.csv"), header = T, sep = ",")

# OPTIONAL output: only basic variables
write.table(project_var %>% dplyr::select(project, t0, country, area_ha),
             paste0(out_path, "_project_var_basic.csv"), sep = ",", row.names = F)

#Output: additionality time series
additionality_estimates = lapply(additionality_out, function(x) x$additionality_estimates)
names(additionality_estimates) = projects
saveRDS(additionality_estimates, paste0(out_path, "_additionality_estimates.rds"))
#additionality_estimates = read_rds(paste0(out_path, "_additionality_estimates.rds"))

#OPTIONAL output: additionality distribution data to send to Ofir
additionality_distribution = lapply(seq_along(projects), function(i) {
  additionality_estimates[[i]] %>%
    filter(started) %>%
    dplyr::select(year, additionality, pair) %>%
    mutate(project = projects[i])
}) %>%
  do.call(rbind, .)
write.csv(additionality_distribution, paste0("/maps/epr26/tmf_pipe_out/additionality_distribution.csv"), row.names = F)


# C. Get baseline ----
baseline = lapply(seq_along(projects), function(i) {
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
names(baseline) = projects
saveRDS(baseline, paste0(out_path, "_baseline.rds"))
#baseline = read_rds(paste0(out_path, "_baseline.rds"))


# D. Bootstrapping ----
BootMean = function(dat, ind) mean(dat[ind, ], na.rm = T)
boot_n = 1000

#bootstrapped mean of baseline C loss
c_loss_boot = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  baseline_i = baseline[[i]] %>%
    dplyr::select(c_loss)
  bootstrap_i = data.frame(project = projects[i],
                           val = as.vector(boot::boot(baseline_i, BootMean, R = boot_n)$t)) #around 23 seconds per run
  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")
  return(bootstrap_i)
})
names(c_loss_boot) = projects

#Output: bootstrapped baseline C loss values for all projects
saveRDS(c_loss_boot, paste0(out_path, "_bootstrapped_base_c_loss.rds"))
#c_loss_boot = read_rds(paste0(out_path, "_bootstrapped_base_c_loss.rds"))

#bootstrapped mean of observed additionality
obs_add_boot = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  area_i = project_var$area_ha[i]
  obs_i = additionality_estimates[[i]] %>%
    dplyr::select(additionality) %>%
    mutate(additionality = additionality / area_i)
  bootstrap_i = data.frame(project = projects[i],
                           val = as.vector(boot::boot(obs_i, BootMean, R = boot_n)$t)) #around 23 seconds per run
  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")
  return(bootstrap_i)
})
names(obs_add_boot) = projects

#Output: bootstrapped observed addditionality values for all projects
saveRDS(obs_add_boot, paste0(out_path, "_bootstrapped_obs_add.rds"))
#c_loss_boot = read_rds(paste0(out_path, "_bootstrapped_obs_add.rds"))

# D. Generate additionality forecast and estimate project effectiveness ----
forecast_summ = lapply(seq_along(projects), function(i) {
  c_loss_i = c_loss_boot[[i]]$val
  c_loss_mean = mean(c_loss_i, na.rm = T)
  c_loss_ci = (qt(p = 0.975, df = boot_n - 1) * sd(c_loss_i, na.rm = T) / sqrt(boot_n))

  obs_add_i = obs_add_boot[[i]]$val
  obs_add_mean = mean(obs_add_i, na.rm = T)
  obs_add_ci = (qt(p = 0.975, df = boot_n - 1) * sd(obs_add_i, na.rm = T) / sqrt(boot_n))

  out_df = data.frame(c_loss_mean = c_loss_mean, c_loss_ci = c_loss_ci,
                      obs_add_mean = obs_add_mean, obs_add_ci = obs_add_ci) %>%
    mutate(c_loss_upper = c_loss_mean + c_loss_ci,
           c_loss_lower = c_loss_mean - c_loss_ci,
           obs_add_upper = obs_add_mean + obs_add_ci,
           obs_add_lower = obs_add_mean - obs_add_ci)
  return(out_df)
}) %>%
  do.call(rbind, .) %>%
  mutate_all(function(x) signif(x, 3)) %>%
  mutate(project = projects) %>%
  mutate(effectiveness_mean = obs_add_mean / c_loss_mean, #calculate project effectiveness
         effectiveness_upper = obs_add_upper / c_loss_lower,
         effectiveness_lower = obs_add_lower / c_loss_upper)


#Output: additionality forecast summaries
write.table(forecast_summ, paste0(out_path, "_forecast_summ.csv"), sep = ",", col.names = NA, row.names = T)

if(visualise) {
  #Visualisation: Figure S1: side-by-side comparison of forecast vs. observed additionality distributions for each project
  df_forecast_obs = lapply(seq_along(projects), function(i) {
    rbind(data.frame(type = "Forecast",
                     value = c_loss_boot[[i]]$val),
          data.frame(type = "Observed",
                     value = obs_add_boot[[i]]$val)) %>%
      mutate(project = projects[i])
  }) %>% do.call(rbind, .)

  p_forecast_obs = ggplot(data = df_forecast_obs, aes(x = type, y = value)) +
    geom_boxplot(aes(color = type)) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_wrap(vars(project), ncol = 5) +
    scale_x_discrete(labels = c("Forecast", "Observed")) +
    scale_color_manual(values = c("red", "black"),
                      labels = c("Forecast", "Observed")) +
    labs(x = "", y = "Annual additionality (MgC/ha/yr)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          strip.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.text.x = element_text())
  ggsave(paste0(out_path, "_forecast_obs.png"), width = 4000, height = 4000, units = "px")

  # p_forecast_obs_list = lapply(seq_along(projects), function(i) {
  #   ggplot(data = subset(df_forecast_obs, project == projects[i]), aes(x = type, y = value)) +
  #   geom_boxplot(aes(color = type)) +
  #   geom_hline(yintercept = 0, linetype = 3) +
  #   scale_x_discrete(labels = c("Forecast", "Observed")) +
  #   scale_color_manual(values = c("red", "black"),
  #                      labels = c("Forecast", "Observed")) +
  #   labs(x = "", y = "", title = projects[i]) +
  #   theme_bw() +
  #   theme(panel.grid = element_blank(),
  #         legend.position = "bottom",
  #         strip.text = element_text(size = 16),
  #         axis.title = element_text(size = 16),
  #         axis.text = element_text(size = 14),
  #         axis.text.x = element_blank(),
  #         legend.title = element_blank(),
  #         legend.text = element_text(size = 14))
  # }) %>%
  #   ggpubr::ggarrange(plotlist = ., ncol = 5, nrow = 4, common.legend = TRUE,
  #                     legend = "bottom", align = "hv") + bgcolor("white")
  # annotate_figure(p_forecast_obs_list,
  #                 left = textGrob("Annual additionality (MgC/ha/yr)", rot = 90, vjust = 2, gp = gpar(cex = 2)))
  # ggsave(paste0(out_path, "_forecast_obs_ggarrange.png"), width = 5000, height = 4000, units = "px")

  #Visualisation: Figure 1: scatterplot of mean forecast vs mean observed
  p_forecast_obs_mean = ggplot(data = forecast_summ, aes(x = c_loss_mean, y = obs_add_mean)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_text(aes(x = c_loss_mean + 0.005, y = obs_add_mean + 0.005, label = project)) +
      scale_x_continuous(limits = c(0, 0.55), expand = c(0, 0)) + #ensures no padding
      scale_y_continuous(limits = c(-0.4, 0.4), expand = c(0, 0)) +
      labs(x = "Mean forecasted annual additionality (MgC/ha/yr)", y = "Mean observed annual additionality (MgC/ha/yr)") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
    ggsave(paste0(out_path, "_forecast_obs_mean.png"), width = 4000, height = 4000, units = "px")

  #Visualisation: Figure 2: project effectiveness
  p_effectiveness = ggplot(data = forecast_summ, aes(x = c_loss_mean, y = effectiveness_mean)) +
      geom_point() +
      geom_hline(yintercept = 1, linetype = 2) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_text(aes(x = c_loss_mean + 0.01, y = effectiveness_mean, label = project)) +
      scale_x_continuous(limits = c(0, 0.55), expand = c(0, 0)) + #ensures no padding
      scale_y_continuous(limits = c(-10, 30), expand = c(0, 0)) +
      labs(x = "Mean forecasted annual additionality (MgC/ha/yr)", y = "Project effectiveness") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  ggsave(paste0(out_path, "_effectiveness.png"), width = 4000, height = 4000, units = "px")

  #Visualisation: Figure 3: over-crediting risk if the forecast is used to issue credits ex ante
  df_overcredit = lapply(seq_along(projects), function(i) {
    forecast_i = subset(forecast_summ, project == projects[i])$c_loss_mean
    obs_i = additionality_estimates[[i]] %>%
      mutate(additionality = additionality / area_i) %>%
      pull(additionality)
    prob_overcredit = length(obs_i[obs_i < forecast_i]) / length(obs_i)
    mean_overcredit = sum(forecast_i - obs_i[obs_i < forecast_i]) / length(obs_i)

    df_overcredit_i = data.frame(project = projects[i],
                                 forecast = forecast_i,
                                 prob_overcredit = prob_overcredit,
                                 mean_overcredit = mean_overcredit)
    return(df_overcredit_i)
    }) %>%
      do.call(rbind, .) %>%
      mutate(project = factor(project, levels = projects))

    p_prob_overcredit = ggplot(data = df_overcredit, aes(x = forecast, y = prob_overcredit)) +
      geom_point() +
      geom_text(aes(x = forecast, y = prob_overcredit + 0.01, label = project)) +
      labs(x = "Mean forecasted annual additionality (MgC/ha/yr)", y = "Over-crediting risk") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
    ggsave(paste0(out_path, "_forecast_overcrediting_prob.png"), width = 3000, height = 3000, units = "px")

    # p_mean_overcredit = ggplot(data = df_overcredit, aes(x = forecast, y = mean_overcredit)) +
    #   geom_point() +
    #   geom_abline(intercept = 0, slope = 1, linetype = 2) +
    #   geom_text(aes(x = forecast, y = mean_overcredit + 0.005, label = project)) +
    #   scale_x_continuous(limits = c(0, 0.75)) + #ensures no padding
    #   scale_y_continuous(limits = c(0, 0.75)) +
    #   labs(x = "Mean forecasted annual additionality (MgC/ha/yr)", y = "Expected amount of over-crediting (MgC/ha/yr)") +
    #   theme_bw() +
    #   theme(panel.grid = element_blank(),
    #         axis.title = element_text(size = 16),
    #         axis.text = element_text(size = 14))
    # ggsave(paste0(out_path, "_forecast_overcrediting_mean.png"), width = 3000, height = 3000, units = "px")
}