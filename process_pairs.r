



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
library(pryr) #pryr::object_size
library(ggpubr) #ggpubr::ggarrange
# library(rnaturalearthdata)
# library(terra)
# library(pbapply)
# library(cleangeo)
# library(foreach)
# library(lwgeom)
# library(countrycode)

source("functions.r") #cpc_rename, tmfemi_reformat
source("ProcessPairs.r") #RetrievePoints, ProcessPairs
source("PlotExAnte.r")

# Load Tom's script ----
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

# Load parameters ----
#input: project_dir, exclude_id, acd_id
#output: projects, pair_dirs, acd_dir
analysis_type = "full" #"old_source", "full", "grid", "ac", "control"
forecast = (analysis_type == "ac")
stratified = T #stratify project pixels by accessibility
pr_vec = seq(0.01, 0.99, by = 0.01) #different quantiles of baseline carbon loss to test

if(analysis_type == "old_source") {

  project_dir = "/maps/pf341/results/2024-january-pipeline" #old results from Patrick's pipeline run
  projects = list.files(project_dir, full = T) %>%
    str_subset("pairs") %>%
    basename() %>%
    str_replace("_pairs", "")

  #find projects with a non-empty carbon_density.csv
  acd_dirs = "/maps/pf341/results/live-pipeline/"
  acd_id = list.files(acd_dir, full = T) %>%
    str_subset("carbon-density") %>%
    basename() %>%
    str_replace("-carbon-density.csv", "")

  #remove: projects 1566, 1067, 958, 1133 are anomalous, 562 is not in TMF extent
  #remove anomalous projects and projects with incomplete carbon density data
  exclude_id = c("1566", "1067", "958", "1133", "562")
  projects = projects[which(projects %in% acd_id) & !(projects %in% exclude_id)] %>%
    as.numeric() %>%
    sort() %>%
    as.character()
  pair_dirs = paste0("/maps/pf341/results/2024-january-pipeline/", projects, "_pairs/")

} else if(analysis_type == "full") {

  project_dir = "/maps/epr26/tmf_pipe_out/" #new results from E-Ping's pipeline run
  projects = list.files(project_dir, full = T) %>%
    str_subset("\\.", negate = T) %>%
    str_subset("\\_grid", negate = T) %>%
    str_subset("ac\\_", negate = T) %>%
    str_subset("fit\\_distribution", negate = T) %>%
    basename()

  #remove: 0000, 9999 are controls and test; 562, 1202, 1340 have incomplete ACD; 1399 and 1408 anomalous
  #remove anomalous projects
  exclude_id = c("0000", "9999", "562", "1202", "1340a", "1399", "1408")
  projects = projects[!(projects %in% exclude_id)] %>%
    str_replace("a", "") %>%
    as.numeric() %>%
    sort() %>%
    as.character()
  projects[which(projects == "612")] %<>% str_c("a") #612, 1340 replaced by 612a, 1340a
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")

} else if(analysis_type == "grid") {

  project_dir = "/maps/epr26/tmf_pipe_out/1201_grid/"
  projects = 1:49 #27, 31 with no matches
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")

} else if(analysis_type == "control") {

  project_dir = "/maps/epr26/tmf_pipe_out/0000_grid/"
  projects = c(2:5, 7, 8, 10)
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")

} else if(analysis_type == "ac") {

  project_dir = "/maps/epr26/tmf_pipe_out/"
  projects = list.files(project_dir, full = T) %>%
    str_subset("ac\\d\\d") %>%
    basename()
  pair_dirs = paste0(project_dir, projects, "/pairs/")
  acd_dirs = paste0(project_dir, projects, "/")

}

out_prefix = ifelse(analysis_type == "grid", "grid_1201", analysis_type)
out_path = paste0("/maps/epr26/tmf_pipe_out/out_", out_prefix)

# Loop through all projects ----
#additionality_out = mclapply(seq_along(projects), mc.cores = 7, function(i) { #mclapply() does not work on Windows
k_elevation_plot = vector("list", length(projects))
k_slope_plot = vector("list", length(projects))
k_access_plot = vector("list", length(projects))

additionality_out = lapply(seq_along(projects), function(i) { #mclapply() does not work on Windows
  a = Sys.time()

  if(analysis_type == "old_source") {
    cat("Use new results from epr26 instead.\n")
    return(NULL)
  }

  proj_id = projects[i]

  proj_name = switch(analysis_type,
              "full" = proj_id,
              "grid" = paste0("1201_", proj_id),
              "control" = paste0("0000_", proj_id),
              "ac" = proj_id)

  t0 = switch(analysis_type,
              "full" = filter(proj_meta, ID == str_replace(proj_id, "a", ""))$t0,
              "grid" = filter(proj_meta, ID == "1201")$t0,
              "control" = 2011,
              "ac" = 2021)

  #column names to select from vicinity
  luc_t_10 = paste0("luc_", t0 - 10)
  luc_t0 = paste0("luc_", t0)

  country = switch(analysis_type,
                   "full" = filter(proj_meta, ID == str_replace(proj_id, "a", ""))$COUNTRY,
                   "grid" = filter(proj_meta, ID == "1201")$COUNTRY,
                   "control" = "Brazil",
                   "ac" = "Brazil")

  #get area_ha
  aoi_path = switch(analysis_type,
                    "full" = "/maps/epr26/tmf-data/projects/",
                    "grid" = "/maps/epr26/tmf-data-grid/1201/1201_",
                    "control" = "/maps/epr26/tmf-data-grid/0000/0000_",
                    "ac" = "/maps/epr26/tmf-data/projects/")
  area_ha = st_read(paste0(aoi_path, proj_id, ".geojson")) %>%
    st_make_valid() %>%
    st_union() %>%
    st_transform(4326) %>%
    st_area_ha() #area in hectares

  #get biome variables and project vicinity
  file_prefix = paste0(project_dir, proj_id, "/",
                       switch(analysis_type,
                              "full" = "",
                              "grid" = "1201_",
                              "control" = "0000_",
                              "ac" = ""),
                       proj_id)

  k = read_parquet(paste0(file_prefix, "k.parquet")) %>%
    dplyr::select(c("lat", "lng", "ecoregion", "elevation", "slope", "access", luc_t_10, luc_t0)) %>%
    mutate(lucc = paste(.data[[luc_t_10]], .data[[luc_t0]], sep = "_"),
           defor = ifelse(lucc %in% c("1_2", "1_3", "1_4"), 1, 0)) %>%
    rename(k_ecoregion = ecoregion) #for biome of project pixels
  logit_k = glm(defor ~ elevation + slope + access, data = k, family = "binomial")
  summary(logit_k)
  confint(logit_k)
  logit_k_coef = summary(logit_k)$coefficients

  if(logit_k_coef["elevation", "Pr(>|z|)"] < 0.05) {
    k_pred = data.frame(elevation = seq(min(k$elevation), max(k$elevation), len = 10000),
                        slope = mean(k$slope),
                        access = mean(k$access)) %>%
      mutate(defor = predict(logit_k, newdata = ., type = "response"))
    thres = k_pred$elevation[which(k_pred$defor < 0.01)[1]]

  k_elevation_plot[[i]] = ggplot(data = k_pred, aes(x = elevation, y = defor)) +
    geom_line() +
    geom_point(data = k, aes(x = elevation, y = defor * max(k_pred$defor))) +
    geom_hline(yintercept = 0.01, linetype = 3) +
    geom_vline(xintercept = thres, color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank())
  }
  if(logit_k_coef["slope", "Pr(>|z|)"] < 0.05) {
    k_pred = data.frame(elevation = mean(k$elevation),
                        slope = seq(min(k$slope), max(k$slope), len = 10000),
                        access = mean(k$access)) %>%
      mutate(defor = predict(logit_k, newdata = ., type = "response"))
    thres = k_pred$slope[which(k_pred$defor < 0.01)[1]]

  k_slope_plot[[i]] = ggplot(data = k_pred, aes(x = slope, y = defor)) +
    geom_line() +
    geom_point(data = k, aes(x = slope, y = defor * max(k_pred$defor))) +
    geom_hline(yintercept = 0.01, linetype = 3) +
    geom_vline(xintercept = thres, color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank())
  }
  if(logit_k_coef["access", "Pr(>|z|)"] < 0.05) {
    k_pred = data.frame(elevation = mean(k$elevation),
                        slope = mean(k$slope),
                        access = seq(min(k$access), max(k$access), len = 10000)) %>%
      mutate(defor = predict(logit_k, newdata = ., type = "response"))
    thres = k_pred$access[which(k_pred$defor < 0.01)[1]]

  k_access_plot[[i]] = ggplot(data = k_pred, aes(x = access, y = defor)) +
    geom_line() +
    geom_point(data = k, aes(x = access, y = defor * max(k_pred$defor))) +
    geom_hline(yintercept = 0.01, linetype = 3) +
    geom_vline(xintercept = thres, color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank())
  }

#1 / (1 + exp(-(logit_k_coef[1, 1] + logit_k_coef[2, 1] * 242.6669 + logit_k_coef[3, 1] * 6.0885 + logit_k_coef[4, 1] * 22)))

  matches = read_parquet(paste0(file_prefix, "matches.parquet"))
  vicinity_area = nrow(matches) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares
  matches = matches %>%
    dplyr::select(c("lat", "lng", "ecoregion", "elevation", "slope", "access", luc_t_10, luc_t0)) %>%
    mutate(lucc_10_0 = paste(.data[[luc_t_10]], .data[[luc_t0]], sep = "_"),
           defor_10_0 = ifelse(lucc_10_0 %in% c("1_2", "1_3", "1_4"), 1, 0)) %>%
    rename(s_ecoregion = ecoregion) #for biome of matched pixels
  vicinity = matches[sample(nrow(matches), 2500000), ] %>% #sub-sample project vicinity down to around the smallest vicinity size of the 15 projects (2559309)
    dplyr::select(c(luc_t_10, luc_t0, "access"))

  #get ACD of undisturbed forest
  acd = read.csv(paste0(file_prefix, "carbon-density.csv"))
  acd_1 = filter(acd, land.use.class == 1)$carbon.density %>% fillNA()

  #gather common project-level variables
  project_var = data.frame(t0 = t0, country = country, acd_1 = acd_1, area_ha = area_ha, vicinity_area = vicinity_area, project = proj_name)

  #find paths to match and unmatached points in each sampled pairs
  pair_paths = list.files(pair_dirs[i], full = T) %>% str_subset(".parquet")
  matchless_ind = pair_paths %>% str_detect("matchless")
  matchless_paths = pair_paths[matchless_ind]
  matched_paths = pair_paths[!matchless_ind]

  #exit if no matches
  if(length(matched_paths) == 0) return(list(project_var = project_var, additionality_estimates = NULL))

  #loop through all sampled pairs, get matched points and additionality series in each pair
  pairs_out = lapply(seq_along(matched_paths), function(j) {
    ProcessPairs(matched_path = matched_paths[j], matchless_path = matchless_paths[j],
                  k = k, matches = matches %>% dplyr::select(c("lat", "lng", "s_ecoregion")),
                  t0 = t0, area_ha = area_ha, acd = acd, pair_id = j)
  })

  #stratify by accessibility
  if(stratified) {
    pts_matched = lapply(pairs_out, function(x) x$pts_matched) %>%
      do.call(dplyr::bind_rows, .) %>%
      mutate(transition_10_0 = paste(.data[[paste0("JRC", t0 - 10)]], .data[[paste0("JRC", t0)]], sep = "_"))
    access_breaks = c(0, quantile(c(vicinity$access, pts_matched$accessibility), seq(0.1, 1, by = 0.1)))
    vicinity$access_class = cut(vicinity$access, breaks = access_breaks, labels = 1:10, include.highest = T, right = T)
    pts_matched$access_class = cut(pts_matched$accessibility, breaks = access_breaks, labels = 1:10, include.highest = T, right = T)
  }

  #project pixel and matched pixel's pre-project CPC change
  cpcc = lapply(pairs_out, function(x) {
    x$pts_matched %>%
      dplyr::select(c("treatment", "defor_10_0")) %>%
      st_drop_geometry()
  }) %>%
    do.call(dplyr::bind_rows, .) #use bind_rows() for speed

  #project pixel and matched pixel's pre-project and during-project LUC change
  lucc = lapply(pairs_out, function(x) x$lucc) %>% do.call(dplyr::bind_rows, .)

  #combine pair-level variables and project-level variables
  pair_biome = lapply(pairs_out, function(x) x$pts_matched$biome) %>%
    unlist() %>%
    table() %>%
    which.max() %>%
    names()

  pair_var_summary = lapply(pairs_out, function(x) x$pair_var) %>%
    do.call(dplyr::bind_rows, .) %>%
    group_by(var) %>%
    summarise(min = min(val), median = median(val), max = max(val)) %>%
    pivot_longer(cols = min:max, names_to = "stat", values_to = "val") %>%
    mutate(var = paste0(var, "_", stat)) %>%
    dplyr::select(c(var, val)) %>%
    pivot_wider(names_from = "var", values_from = "val") %>%
    mutate(biome = pair_biome)

  project_var = cbind(project_var, pair_var_summary)

  additionality_estimates = lapply(pairs_out, function(x) x$out_df) %>%
    do.call(dplyr::bind_rows, .) %>%
    mutate(started = ifelse(year > t0, T, F))

  b = Sys.time()
  cat(proj_id, ":", b - a, "\n")
  if(stratified) {
    return(list(project_var = project_var, additionality_estimates = additionality_estimates,
                cpcc = cpcc, lucc = lucc, vicinity = vicinity, pts_matched = pts_matched))
  } else {
    return(list(project_var = project_var, additionality_estimates = additionality_estimates,
            cpcc = cpcc, lucc = lucc, vicinity = vicinity))
  }
})
names(additionality_out) = projects


# Stratified analysis ----
plot_stratified = lapply(seq_along(projects), function(i) {
  pts_matched = additionality_out[[i]]$pts_matched %>%
    dplyr::select(c("defor_10_0", "transition_10_0", "accessibility", "access_class")) %>%
    st_drop_geometry()
  vicinity = additionality_out[[i]]$vicinity %>%
    mutate(defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
    dplyr::select(c("defor_10_0", "transition_10_0", "access", "access_class"))

  counterfactual_luc = data.frame(
    access_class = 1:10,
    lucc = tapply(pts_matched$transition_10_0, pts_matched$access_class, findLUCC),
    source = "Counterfactual"
  )


# lucc = rep(NA, 10)
# for(i in 1:10) {
#   a = Sys.time()
#   lucc = boot(data = pts_matched, statistic = function(x) findLUCC(subset(x, access_class == i)$transition_10_0), R = 100)
#   b = Sys.time()
#   b - a
# }


  vicinity_luc = data.frame(
    access_class = 1:10,
    lucc = tapply(vicinity$transition_10_0, vicinity$access_class, findLUCC),
    source = "Vicinity"
  )

  data_stratified = rbind(counterfactual_luc, vicinity_luc)

  x_labels = paste0("(", seq(1, 91, by = 10), ", ", seq(10, 100, by = 10), "]")

  p = ggplot(data = data_stratified, aes(x = access_class, y = lucc * 100, group = source)) +
    geom_line(aes(color = source)) +
    scale_x_continuous(limits = c(1, 10), breaks = 1:10, labels = x_labels) +
    scale_y_continuous() +
#    scale_y_sqrt() +
    scale_color_manual(values = c("red", "blue")) +
    labs(title = projects[i], x = "Accessibility class (percentile)", y = "Annual forest loss (%)\n(square root transformed)") +
    theme_bw() +
    theme(title = element_text(size = 24),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.5),
)

  return(p)
})

ggpubr::ggarrange(plotlist = plot_stratified, ncol = 3, nrow = 5)

# Save additionality outputs ----
#project-level variables
project_var = lapply(additionality_out, function(x) x$project_var) %>% do.call(dplyr::bind_rows, .)
#use bind_rows because of potentially different numbers of columns
saveRDS(project_var, paste0(out_path, "_project_var.rds"))
#project_var = read_rds(paste0(out_path, "_project_var.rds"))

#basic information
basic_df = project_var %>% dplyr::select(t0, country, area_ha, vicinity_area, project)
write.table(basic_df, paste0(out_path, "_basic_info.csv"), sep = ",", row.names = F)

#additionality time series
additionality_estimates = lapply(additionality_out, function(x) x$additionality_estimates)
names(additionality_estimates) = projects
saveRDS(additionality_estimates, paste0(out_path, "_additionality_estimates.rds"))
#additionality_estimates = read_rds(paste0(out_path, "_additionality_estimates.rds"))

#deforestation rates
cpcc_list = lapply(additionality_out, function(x) x$cpcc)
names(cpcc_list) = projects
saveRDS(cpcc_list, paste0(out_path, "_cpcc_list.rds"))
#cpcc_list = read_rds(paste0(out_path, "_cpcc_list.rds"))

lucc_list = lapply(additionality_out, function(x) x$lucc)
names(lucc_list) = projects
saveRDS(lucc_list, paste0(out_path, "_lucc_list.rds"))
#lucc_list = read_rds(paste0(out_path, "_lucc_list.rds"))

vicinity_list = lapply(additionality_out, function(x) x$vicinity)
names(vicinity_list) = projects
saveRDS(vicinity_list, paste0(out_path, "_vicinity_list.rds"))
#vicinity_list = read_rds(paste0(out_path, "_vicinity_list.rds"))


# Compare CPC- and LUC- based deforestation rates ----
cpcc_lucc_list = lapply(seq_along(projects), function(i) {
  cpcc = cpcc_list[[i]]
  lucc = lucc_list[[i]]
  t0 = project_var$t0[i]
  luc_t_10 = paste0("luc_", t0 - 10)
  luc_t0 = paste0("luc_", t0)
  vicinity = vicinity_list[[i]] %>% mutate(transition_10_0 = paste(.data[[luc_t_10]], .data[[luc_t0]], sep = "_"))

  #pre-project CPC change of project and matched pixels (in each year in each pair)
  project_cpc = subset(cpcc, treatment == "treatment") %>% pull(defor_10_0) %>% fillNA()
  counterfactual_cpc = subset(cpcc, treatment == "control") %>% pull(defor_10_0) %>% fillNA()

  #pre-project CPC change of project vicinity (in each pixel)
  vicinity_cpc = (vicinity$cpc10_u - vicinity$cpc0_u) / 10

  #pre-project LUC change of project and matched pixels (in each year in each pair)
  project_luc = subset(lucc, treatment == "treatment" & !started) %>% pull(prop_df) %>% fillNA()
  counterfactual_luc = subset(lucc, treatment == "control" & !started) %>% pull(prop_df) %>% fillNA()

  #pre-project LUC change of project vicinity
  vicinity_luc = findLUCC(vicinity$transition_10_0)
  subsamp_size = nrow(cpcc) / (2 * 100)
  vicinity_luc = rep(NA, 100)
  for(i in 1:100) {
    vicinity_luc[i] = findLUCC(slice_sample(vicinity, n = subsamp_size)$transition_10_0)
  }

  #during-project LUC change of project and matched pixels (in each year in each pair)
  post_t0_project_luc = subset(lucc, treatment == "treatment" & started) %>% pull(prop_df) %>% fillNA()
  post_t0_counterfactual_luc = subset(lucc, treatment == "control" & started) %>% pull(prop_df) %>% fillNA()

  cpcc_lucc_df = do.call(rbind,
    list(data.frame(Type = "CPC", Period = "Pre-t0", var = "project_cpc", val = project_cpc),
         data.frame(Type = "CPC", Period = "Pre-t0", var = "counterfactual_cpc", val = counterfactual_cpc),
         data.frame(Type = "CPC", Period = "Pre-t0", var = "vicinity_cpc", val = vicinity_cpc),
         data.frame(Type = "LUC", Period = "Pre-t0", var = "project_luc", val = project_luc),
         data.frame(Type = "LUC", Period = "Pre-t0", var = "counterfactual_luc", val = counterfactual_luc),
         data.frame(Type = "LUC", Period = "Pre-t0", var = "vicinity_luc", val = vicinity_luc),
         data.frame(Type = "LUC", Period = "Post-t0", var = "post_t0_project_luc", val = post_t0_project_luc),
         data.frame(Type = "LUC", Period = "Post-t0", var = "post_t0_counterfactual_luc", val = post_t0_counterfactual_luc)))
  cpcc_lucc_df$Type = factor(cpcc_lucc_df$Type, levels = c("CPC", "LUC"))
  cpcc_lucc_df$Period = factor(cpcc_lucc_df$Period, levels = c("Pre-t0", "Post-t0"))
  cpcc_lucc_df$var = factor(cpcc_lucc_df$var, levels = c("project_cpc", "counterfactual_cpc", "vicinity_cpc",
                                                         "project_luc", "counterfactual_luc", "vicinity_luc",
                                                         "post_t0_project_luc", "post_t0_counterfactual_luc"))

  cpcc_lucc_summary = cpcc_lucc_df %>%
    group_by(var) %>%
    reframe(min = min(val, na.rm = T),
            q1 = quantile(val, 0.25, na.rm = T),
            median = median(val, na.rm = T),
            mean = mean(val, na.rm = T),
            q3 = quantile(val, 0.75, na.rm = T),
            max = max(val, na.rm = T))

  return(list(cpcc_lucc_df = cpcc_lucc_df, cpcc_lucc_summary = cpcc_lucc_summary))
})

cpcc_lucc_summary_list = lapply(cpcc_lucc_list, function(x) x$cpcc_lucc_summary)
names(cpcc_lucc_summary_list) = projects
saveRDS(cpcc_lucc_summary_list, paste0(out_path, "_cpcc_lucc_summary.rds"))

for(i in seq_along(cpcc_lucc_summary_list)) {
  write.csv(cpcc_lucc_summary_list[[i]] %>% mutate(project = projects[i]),
            paste0("/maps/epr26/tmf_pipe_out/out_", out_prefix, "_cpcc_lucc_summ_", i, ".csv"))
}

#CPC-based and LUC-based deforestation rates
cpcc_lucc_max = sapply(cpcc_lucc_list, function(x) range(x$cpcc_lucc_df$val * 100)) %>% max()
plot_cpcc_lucc = lapply(seq_along(cpcc_lucc_list), function(i) {
  ggplot(data = cpcc_lucc_list[[i]]$cpcc_lucc_df, aes(x = var, y = val * 100)) +
    geom_boxplot(aes(color = Type, linetype = Period)) +
    scale_x_discrete(labels = c("Project-CPC", "Counterfactual-CPC", "Vicinity-CPC",
                                "Project-LUC", "Counterfactual-LUC", "Vicinity-LUC",
                                "Post-t0 project-LUC", "Post-t0 counterfactual-LUC")) +
    scale_y_sqrt(limits = c(0, ceiling(cpcc_lucc_max * 2) / 2)) +
    scale_color_manual(values = c("red", "blue")) +
    scale_linetype_manual(values = 1:2) +
    labs(title = projects[i], x = NULL, y = "Annual forest loss (%)\n(square root transformed)") +
    theme_bw() +
    theme(title = element_text(size = 24),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(size = 18, angle = 90),
          axis.text.y = element_text(size = 18),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 20))
}) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 4, nrow = 2, common.legend = T, legend = "bottom")
SaveMultiPagePlot(plot_cpcc_lucc, "cpcc_lucc", width = 6000, height = 6000)

ggsave(paste0(out_path, "cpcc_lucc"), plot_cpcc_lucc, width = 6000, height = 6000, units = "px", bg = "white")
#  ggsave("cpcc_lucc_summ.png", width = 1000, height = 1000, units = "px")


# Run function for ex ante analysis ----
ex_ante_out = lapply(seq_along(projects), function(i) {
  proj_id = projects[i]
  area_ha = project_var$area_ha[i]
  if(is.null(additionality_estimates[[i]])) {
    return(list(plot_df = NULL, forecast_df = NULL,
                p0 = NULL, p1 = NULL, p2 = NULL, p_legend_grob = NULL, p_perc = NULL, p_overcredit = NULL))
  }

  obs_val = additionality_estimates[[i]] %>%
    filter(started) %>%
    mutate(c_loss = c_loss / area_ha, additionality = additionality / area_ha)

  path = switch(analysis_type,
          "full" = paste0(project_dir, proj_id, "/", proj_id),
          "ac" = paste0(project_dir, proj_id, "/", proj_id),
          "grid" = paste0(project_dir, proj_id, "/1201_", proj_id),
          "control" = paste0(project_dir, proj_id, "/0000_", proj_id))

  acd = read.csv(paste0(path, "carbon-density.csv"))

  PlotExAnte(proj_id = proj_id, area_ha = area_ha, obs_val = obs_val, path = path, acd = acd, forecast = forecast)
})
names(ex_ante_out) = projects


# Save ex ante outputs ----

#ex ante analysis output
ex_ante_plot_df = lapply(ex_ante_out, function(x) x$plot_df)
ex_ante_forecast_df = lapply(ex_ante_out, function(x) x$forecast_df)
names(ex_ante_plot_df) = projects
names(ex_ante_forecast_df) = projects
write_rds(ex_ante_plot_df, paste0(out_path, "_ex_ante_plot_df.rds"))
write_rds(ex_ante_forecast_df, paste0(out_path, "_ex_ante_forecast_df.rds"))

#summary statistics of observed values
obs_c_loss_summ = lapply(ex_ante_out, function(x) {
  if(is.null(x$plot_df)) return(NA)
  quantile(filter(x$plot_df, Type == "obs_c_loss")$Value, c(0.05, 0.1, 0.25, 0.5, 0.75), na.rm = T)
}) %>% do.call(rbind, .)
rownames(obs_c_loss_summ) = projects
write.table(obs_c_loss_summ, paste0(out_path, "_obs_c_loss.csv"), sep = ",", row.names = F)

obs_add_summ = lapply(ex_ante_out, function(x) {
  if(is.null(x$plot_df)) return(NA)
  quantile(filter(x$plot_df, Type == "obs_add")$Value, c(0.05, 0.1, 0.25, 0.5, 0.75), na.rm = T)
}) %>% do.call(rbind, .)
rownames(obs_add_summ) = projects
write.table(obs_add_summ, paste0(out_path, "_obs_add.csv"), sep = ",", row.names = F)


# Make ex ante plots ----
#plot of different baseline periods
x_max = sapply(ex_ante_out, function(x) {
  if(is.null(x$plot_df)) return(NA)
  x$plot_df %>%
    filter(Period != "after") %>%
    pull(Value) %>%
    max()
})
plot_period = lapply(ex_ante_out, function(x) {
  if(is.null(x$p0)) return(NULL)
  x$p0 + scale_x_continuous(limits = c(0, max(x_max, na.rm = T)))
}) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = T, legend = "bottom")
SaveMultiPagePlot(plot_period, "c_loss_periods")

#plot of distributions
x_range = sapply(ex_ante_out, function(x) {
  if(is.null(x$plot_df)) return(NA)
  x$plot_df %>%
    filter(str_detect(Period, "10_0") | str_detect(Type, "obs")) %>%
    pull(Value) %>%
    range()
}) %>% unlist()
x_legend = ex_ante_out[[1]]$p_legend_grob
plot_distr = lapply(seq_along(ex_ante_out), function(i) {
  x = ex_ante_out[[i]]
  if(is.null(x$p1)) return(NULL)
  x1 = x$p1 + ggtitle("") + scale_x_continuous(limits = c(min(x_range, na.rm = T), max(x_range, na.rm = T)))
  x2 = x$p2 + ggtitle("") + scale_x_continuous(limits = c(min(x_range, na.rm = T), max(x_range, na.rm = T)))
  p_1_2 = ggpubr::ggarrange(x1, x2, ncol = 2, nrow = 1, legend = "none")
  ggpubr::annotate_figure(p_1_2, top = ggpubr::text_grob(projects[i], face = "bold", size = 14))
})
plot_distr_all = plot_distr %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = T, legend = "bottom", legend.grob = x_legend)
SaveMultiPagePlot(plot_distr_all, "distribution", width = 4000, height = 4000)

#plot of forecast and overcrediting risk
plot_forecast = lapply(seq_along(ex_ante_out), function(i) {
  x = ex_ante_out[[i]]
  if(is.null(x$p_perc)) return(NULL)
  x1 = x$p_perc + ggtitle("")
  x2 = x$p_overcredit + ggtitle("")
  p_1_2 = ggpubr::ggarrange(x1, x2, ncol = 2, nrow = 1)
  ggpubr::annotate_figure(p_1_2, top = ggpubr::text_grob(projects[i], face = "bold", size = 14))
}) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = T, legend = "bottom")
SaveMultiPagePlot(plot_forecast, "forecast")