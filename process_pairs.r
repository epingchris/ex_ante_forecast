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
library(boot) #boot::boot
library(microbenchmark) #microbenchmark::microbenchmark
# library(rnaturalearthdata)
# library(terra)
# library(pbapply)
# library(cleangeo)
# library(foreach)
# library(lwgeom)
# library(countrycode)

source("functions.r") #cpc_rename, tmfemi_reformat
source("FilterBaseline.r")
source("ProcessPairs.r") #RetrievePoints, ProcessPairs
source("CalcExAnte.r")


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
  exclude_id = c("0000", "9999", "562", "612a", "1202", "1340a", "1399", "1408")
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


# Get additionality of all projects ----
#additionality_out = mclapply(seq_along(projects), mc.cores = 15, function(i) { #mclapply() does not work on Windows
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

  #get ACD and ACD of undisturbed forest
  acd = read.csv(paste0(file_prefix, "carbon-density.csv"))
  for(class in 1:6) {
    if(class %in% acd$land.use.class == F) acd = rbind(acd, c(class, NA))
  }
  acd_1 = filter(acd, land.use.class == 1)$carbon.density

  k = read_parquet(paste0(file_prefix, "k.parquet")) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion)
  matches = read_parquet(paste0(file_prefix, "matches.parquet")) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion)

  #gather common project-level variables
  project_var = data.frame(project = proj_name, t0 = t0, country = country, acd_1 = acd_1, area_ha = area_ha)

  #find paths to match and unmatached points in each sampled pairs
  pair_paths = list.files(pair_dirs[i], full = T) %>% str_subset(".parquet")
  matchless_ind = pair_paths %>% str_detect("matchless")
  matchless_paths = pair_paths[matchless_ind]
  matched_paths = pair_paths[!matchless_ind]

  #exit if no matches
  if(length(matched_paths) == 0) {
    return(list(project_var = project_var, additionality_estimates = NULL))
  }

  #loop through all sampled pairs, get matched points and additionality series in each pair
  pairs_out = lapply(seq_along(matched_paths), function(j) {
    ProcessPairs(matched_path = matched_paths[j], matchless_path = matchless_paths[j],
                  k = dplyr::select(k, c("lat", "lng", "k_ecoregion")),
                  matches = dplyr::select(matches, c("lat", "lng", "s_ecoregion")),
                  t0 = t0, area_ha = area_ha, acd = acd, pair_id = j)
  })

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
  return(list(project_var = project_var, additionality_estimates = additionality_estimates))
})
names(additionality_out) = projects


# Save outputs ----
#project-level variables
project_var = lapply(additionality_out, function(x) x$project_var) %>% do.call(dplyr::bind_rows, .)
#use bind_rows because of potentially different numbers of columns
saveRDS(project_var, paste0(out_path, "_project_var.rds"))
#project_var = read_rds(paste0(out_path, "_project_var.rds"))

#basic information
basic_df = project_var %>% dplyr::select(t0, country, area_ha, project)
write.table(basic_df, paste0(out_path, "_basic_info.csv"), sep = ",", row.names = F)

#additionality time series
additionality_estimates = lapply(additionality_out, function(x) x$additionality_estimates)
names(additionality_estimates) = projects
saveRDS(additionality_estimates, paste0(out_path, "_additionality_estimates.rds"))
#additionality_estimates = read_rds(paste0(out_path, "_additionality_estimates.rds"))


# Find filtering thresholds ----
var_vec = c("slope", "elevation", "access")
var_label = c("Slope (dgree)", "Elevation (meter)", "Remoteness (minutes)")
filter_out = lapply(projects, function(x) FilterBaseline(analysis_type = analysis_type, proj_id = x))

#visualise and compare results when filtering is done with logistic fitted to K or to M: the results are similar, so only results using M is used
p_k_m_list = lapply(filter_out, function(x) x$plotlist %>% cowplot::plot_grid(plotlist = ., byrow = F, nrow = 3, ncol = 2))

df_effect = lapply(filter_out, function(x) x$df_summary$by_m[1:3]) %>% do.call(rbind, .)
colnames(df_effect) = var_vec

var_max = c(25, 3000, 2000)
text_adjust = c(1, 200, 200)
for(i in seq_along(var_vec)) {
  var_i = var_vec[i]
  var_label_i = var_label[i]
  df_filtering = lapply(seq_along(filter_out), function(j) {
    out_j = filter_out[[j]]$baseline
    out_df = data.frame(all = mean(out_j %>% pull(all_of(var_i)), na.rm = T),
                        low = mean(filter(out_j, risk_m == "low") %>% pull(all_of(var_i)), na.rm = T),
                        high = mean(filter(out_j, risk_m == "high") %>% pull(all_of(var_i)), na.rm = T)) %>%
      mutate(max = apply(., 1, max),
             project = projects[j])
    return(out_df)
  }) %>%
    do.call(rbind, .) %>%
    mutate(project = factor(project, levels = projects),
           rank = rank(all))

ggplot(data = df_filtering, aes(x = rank)) +
    geom_point(aes(y = low, color = "blue")) +
    geom_point(aes(y = high, color = "red")) +
    geom_point(aes(y = all, color = "black")) +
    geom_segment(aes(y = all, xend = rank, yend = low), color = "blue") +
    geom_segment(aes(y = all, xend = rank, yend = high), color = "red") +
    geom_text(aes(x = rank, y = max + text_adjust[i], label = project), size = 4) +
    scale_color_manual(name = "Pixel type", values = c("black", "blue", "red"),
                       labels = c("All", "Low deforestation probability", "High deforestation probability")) +
    labs(x = "Project", y = var_label_i) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 14))
  ggsave(paste0(out_path, "_filtering_comparison_new_", var_i, ".png"), width = 2000, height = 2000, unit = "px")
}

df_slope_filtering = lapply(seq_along(filter_out), function(i) {
  data.frame(project = projects[i],
             all = mean(filter_out[[i]]$baseline$slope, na.rm = T),
             low = mean(filter(filter_out[[i]]$baseline, risk_m == "low")$slope, na.rm = T),
             high = mean(filter(filter_out[[i]]$baseline, risk_m == "high")$slope, na.rm = T))
  }) %>%
  do.call(rbind, .)
df_elevation_filtering = lapply(seq_along(filter_out), function(i) {
  data.frame(project = projects[i],
             all = mean(filter_out[[i]]$baseline$elevation, na.rm = T),
             low = mean(filter(filter_out[[i]]$baseline, risk_m == "low")$elevation, na.rm = T),
             high = mean(filter(filter_out[[i]]$baseline, risk_m == "high")$elevation, na.rm = T))
  }) %>%
  do.call(rbind, .)
df_access_filtering = lapply(seq_along(filter_out), function(i) {
  data.frame(project = projects[i],
             all = mean(filter_out[[i]]$baseline$access, na.rm = T),
             low = mean(filter(filter_out[[i]]$baseline, risk_m == "low")$access, na.rm = T),
             high = mean(filter(filter_out[[i]]$baseline, risk_m == "high")$access, na.rm = T))
  }) %>%
  do.call(rbind, .)


ggplot(data = df_slope_filtering, aes(x = all)) +
  geom_point(aes(y = low, color = "blue")) +
  geom_point(aes(y = high, color = "red")) +
  geom_segment(aes(y = all, xend = all, yend = low), arrow = arrow(length = unit(0.01, "npc")), color = "blue") +
  geom_segment(aes(y = all, xend = all, yend = high), arrow = arrow(length = unit(0.01, "npc")), color = "red") +
  geom_text(aes(x = all, y = low * 1.1, label = project), size = 6) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = c("blue", "red"), labels = c("Low-risk", "High-risk")) +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "Slope (degree) (before filtering)", y = "Slope (degree) (after filtering)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16))
ggsave("filtering_comparison_slope.png", width = 800, height = 800)

effect_pval_df = lapply(filter_out, function(x) x$pval) %>% do.call(rbind, .) %>% signif(., 2)
write.table(effect_pval_df, paste0(out_path, "_effect_pval_by_", by_source, ".csv"), sep = ",", row.names = F)
effect_labels_df = lapply(filter_out, function(x) x$effect_labels) %>% do.call(rbind, .)
write.table(effect_labels_df, paste0(out_path, "_effect_labels_by_", by_source, ".csv"), sep = ",", row.names = F)
exclude_ratio_df = lapply(filter_out, function(x) x$exclude_ratio) %>% do.call(rbind, .)
write.table(exclude_ratio_df %>% round(., 3), paste0(out_path, "_exclude_ratio_by_", by_source, ".csv"), sep = ",", row.names = F)
thres_df = lapply(filter_out, function(x) x$thres) %>% do.call(rbind, .) %>% as.data.frame()
var_vec = c("slope", "elevation", "access")
colnames(thres_df) = var_vec
thres_df$project = projects
write.table(thres_df %>% mutate_at(var_vec, function(x) round(x, 1)),
            paste0(out_path, "_thres_by_", by_source, ".csv"), sep = ",", row.names = F)

n = switch(analysis_type,
           "full" = 4,
           "control" = 4,
           "ac" = 3)

vicinity_slope_plot = lapply(filter_out, function(x) x$plotlist$slope)
SaveMultiPagePlot(vicinity_slope_plot, paste0("vicinity_slope_by_", by_source), n = n, width = 4000, height = 4000)
vicinity_elev_plot = lapply(filter_out, function(x) x$plotlist$elevation)
SaveMultiPagePlot(vicinity_elev_plot, paste0("vicinity_elevation_by_", by_source), n = n, width = 4000, height = 4000)
vicinity_access_plot = lapply(filter_out, function(x) x$plotlist$access)
SaveMultiPagePlot(vicinity_access_plot, paste0("vicinity_accessibility_by_", by_source), n = n, width = 4000, height = 4000)

#project vicinity
vicinity_list = lapply(filter_out, function(x) x$vicinity)
names(vicinity_list) = projects
saveRDS(vicinity_list, paste0(out_path, "_vicinity_by_", by_source, ".rds"))
#vicinity_list = read_rds(paste0(out_path, "vicinity_by_", by_source, ".rds"))

#basic information about the vicinity
vicinity_df = lapply(seq_along(projects), function(i) {
  data.frame(project = projects[i],
             vicinity_area = filter_out[[i]]$vicinity_area,
             vicinity_area_filtered = filter_out[[i]]$vicinity_area_filtered) %>%
    mutate(exclusion_rate = 1 - (vicinity_area_filtered / vicinity_area))
}) %>%
  do.call(rbind, .)

#calculate vicinity C loss before and after filtering
meanCLoss = function(dat, ind) mean(dat[ind, ]$c_loss, na.rm = T)
boot_n = 1000

filter_diff_out = lapply(seq_along(projects), function(i) {
  vicinity = vicinity_list[[i]]
  baseline_c_loss = boot::boot(vicinity, meanCLoss, R = boot_n) #around 700 seconds

  vicinity_filtered = subset(vicinity, !exclude)
  if(nrow(vicinity_filtered) == 0) {
    return(list(baseline_c_loss = baseline_c_loss, baseline_c_loss_filtered = NULL,
                ses = NULL, plot_df = NULL, risk_overcrediting = NULL, risk_reversal = NULL))
  }
  baseline_c_loss_filtered = boot::boot(vicinity_filtered, meanCLoss, R = boot_n) #around 700 seconds

  #collate data
  baseline_val = data.frame(base_c_loss = as.vector(baseline_c_loss$t),
                            base_c_loss_filt = as.vector(baseline_c_loss_filtered$t))

  #calculate standardized effect size
  ses = (mean(baseline_val$base_c_loss_filt) - mean(baseline_val$base_c_loss)) / sd(baseline_val$base_c_loss)

  #long data for plotting
  baseline_long = baseline_val %>%
    pivot_longer(everything(), names_to = "Type", values_to = "Value") %>%
    mutate(project = projects[i])

  return(list(ses = ses, baseline_long = baseline_long))
})
saveRDS(filter_diff_out, paste0(out_path, "vicinity_filtering_by_", by_source, ".rds"))

p_filter_df = lapply(filter_diff_out, function(x) x$baseline_long) %>%
  do.call(rbind, .) %>%
  mutate(project = factor(project, levels = projects))
y_range = range(p_filter_df$Value)
ses_plot_df = data.frame(project = projects, ses = sapply(filter_diff_out, function(x) x$ses)) %>%
  mutate(ses_text = paste0("SES: ", round(ses, 2)))

p_filter = ggplot(data = p_filter_df, aes(x = Type, y = Value)) +
  geom_boxplot(aes(color = Type)) +
  facet_wrap(vars(project), ncol = 5) +
  ggpubr::stat_compare_means(method = "t.test", aes(label = ..p.signif..),
                             label.x = 1.5, label.y = y_range[2] * 1.1, size = 5) +
  geom_text(data = ses_plot_df,
            mapping = aes(x = 1.5, y = y_range[2] * 1.2, label = ses_text), size = 5) +
  scale_x_discrete(labels = c("Unfiltered", "Filtered")) +
  scale_y_continuous(limits = c(0, y_range[2] * 1.3)) +
  scale_color_manual(values = c("red", "blue"),
                      labels = c("Unfiltered", "Filtered")) +
  labs(x = "", y = "Annual carbon loss (Mg/ha)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_vicinity_filtering_by_", by_source, ".png"), width = 4000, height = 4000, units = "px")


# Ex ante analysis ----

#generate additionality forecast
vicinity_summ = lapply(filter_diff_out, function(x) {
  val = subset(x$baseline_long, Type == "base_c_loss_filt")$Value
  scenarios = c(100, 75, 50, 25)
  val_mean = mean(val, na.rm = T) * scenarios / 100
  val_ci = (qt(p = 0.975, df = boot_n - 1) * sd(val, na.rm = T) / sqrt(boot_n)) * scenarios / 100
  column_order = c(rbind(paste0("mean_", scenarios), paste0("ci_", scenarios)))
  out = data.frame(mean = val_mean,
                   ci = val_ci,
                   scenario = scenarios) %>%
    pivot_wider(names_from = scenario, values_from = c(mean, ci), names_sep = "_") %>%
    dplyr::select(all_of(column_order))
}) %>%
  do.call(rbind, .) %>%
  mutate_all(function(x) signif(x, 3)) %>%
  mutate(project_area = round(project_var$area_ha),
         k_filtered_ratio = sapply(filter_out, function(x) x$k_area_ratio) %>% round(., 3),
         project_area_filtered = round(project_area * k_filtered_ratio))
rownames(vicinity_summ) = projects
write.table(vicinity_summ, paste0(out_path, "_vicinity_summ_matches.csv"), sep = ",", col.names = NA, row.names = T)

ex_ante_out = lapply(seq_along(projects), function(i) {
  proj_id = projects[i]
  area_ha = project_var$area_ha[i]
  vicinity = vicinity_list[[i]]
  if(forecast) {
    obs_val = data.frame(c_loss = NA, t_loss = NA, additionality = NA)
  } else {
    obs_val = additionality_estimates[[i]] %>%
      filter(started) %>%
      dplyr::select(c("c_loss", "t_loss", "additionality")) %>%
      mutate_all(function(x) x = x / area_ha)
  }

  CalcExAnte(proj_id = proj_id, vicinity = vicinity, obs_val = obs_val)
})
names(ex_ante_out) = projects


# Save ex ante outputs ----
#before filtering
baseline_c_loss = lapply(ex_ante_out, function(x) x$baseline_c_loss)
names(baseline_c_loss) = projects
write_rds(baseline_c_loss, paste0(out_path, "_baseline_c_loss_", by_source, ".rds"))
#baseline_c_loss = read_rds(paste0(out_path, "_baseline_c_loss.rds"))

#after filtering
baseline_c_loss_filtered = lapply(ex_ante_out, function(x) x$baseline_c_loss_filtered)
names(baseline_c_loss_filtered) = projects
write_rds(baseline_c_loss_filtered, paste0(out_path, "_baseline_c_loss_filtered_", by_source, ".rds"))
#baseline_c_loss_filtered = read_rds(paste0(out_path, "_baseline_c_loss_filtered.rds"))

#standard effect size of filtering
ses_df = data.frame(project = projects, ses = sapply(ex_ante_out, function(x) x$ses)) %>%
  mutate(project = factor(project, levels = projects))
write_rds(ses_df, paste0(out_path, "_ses_df_", by_source, ".rds"))
#ses_df = read_rds(paste0(out_path, "_ses_df.rds"))

#data for plotting
plot_df = lapply(ex_ante_out, function(x) x$plot_df)
names(plot_df) = projects
write_rds(plot_df, paste0(out_path, "_plot_df_", by_source, ".rds"))
#plot_df = read_rds(paste0(out_path, "_plot_df.rds"))

#overcrediting and reversal risk
risk_summary = lapply(ex_ante_out, function(x) {
  x$risk_overcrediting %>%
    pivot_wider(names_from = "scenario", values_from = c("forecast", "risk"))
}) %>%
  do.call(rbind, .) %>%
  mutate(risk_reversal = sapply(ex_ante_out, function(x) x$risk_reversal))
write.table(risk_summary, paste0(out_path, "_risk_summary_", by_source, ".csv"), sep = ",", row.names = F)

#summary statistics of observed values
obs_c_loss_summ = lapply(ex_ante_out, function(x) {
  quantile(filter(x$plot_df, Type == "obs_c_loss")$Value, c(0.05, 0.1, 0.25, 0.5, 0.75), na.rm = T)
}) %>% do.call(rbind, .)
rownames(obs_c_loss_summ) = projects
write.table(obs_c_loss_summ, paste0(out_path, "_obs_c_loss_", by_source, ".csv"), sep = ",", row.names = F)


# Plot ex ante outputs ----
#plot C loss distributions
p_c_loss_df = lapply(seq_along(plot_df), function(i) {
  plot_df[[i]] %>%
    filter(Type != "base_c_loss") %>%
    mutate(project = projects[i])
}) %>%
  do.call(rbind, .) %>%
  mutate(project = factor(project, levels = projects),
         Type = factor(Type, levels = c("base_c_loss_filt", "obs_c_loss", "obs_t_loss", "additionality")))

p_c_loss = ggplot(data = p_c_loss_df, aes(x = Type, y = Value)) +
  geom_boxplot(aes(color = Type)) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(vars(project), ncol = 5) +
  scale_x_discrete(labels = c("Baseline\nC loss", "Observed\nCounterfactual\nC loss", "Observed\nProject\nC loss", "Additionality")) +
  scale_color_manual(values = c("red", "black", "black", "blue"),
                     labels = c("Baseline\nC loss", "Observed\nCounterfactual\nC loss", "Observed\nProject\nC loss", "Additionality")) +
  labs(x = "", y = "Annual carbon flux (Mg/ha)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
ggsave(paste0(out_path, "_c_loss_disribution_", by_source, ".png"), width = 4000, height = 6000, units = "px")

#plot over-claiming risks
p_risk_df = lapply(seq_along(projects), function(i) {
  ex_ante_out[[i]]$risk_overcrediting %>%
    dplyr::select(c("scenario", "risk")) %>%
    mutate(project = projects[i],
           scenario = scenario * 100)
}) %>%
  do.call(rbind, .) %>%
  mutate(scenario = factor(scenario, levels = c(25, 50, 75, 100)))

project_label = data.frame(project = projects,
                           y = filter(p_risk_df, scenario == 100)$risk)

p_risk = ggplot(data = p_risk_df, aes(x = scenario, y = risk, group = project)) +
  geom_line(aes(color = project), linewidth = 2) +
  geom_text(data = project_label, aes(x = 4.1, y = y, label = project)) +
  labs(x = "Assumed scenario (% effectiveness)", y = "Over-claiming risk") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_risk_overclaiming_", by_source, ".png"), width = 4000, height = 4000, units = "px")