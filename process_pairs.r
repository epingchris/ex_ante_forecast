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
source("FilterPoints.r")
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


# Find filtering thresholds ----
filter_out = lapply(projects, function(x) FilterPoints(analysis_type = analysis_type, proj_id = x))

effect_pval_df = lapply(filter_out, function(x) x$pval) %>% do.call(rbind, .) %>% signif(., 2)
write.table(effect_pval_df, paste0(out_path, "_effect_pval.csv"), sep = ",", row.names = F)
effect_labels_df = lapply(filter_out, function(x) x$effect_labels) %>% do.call(rbind, .)
write.table(effect_labels_df, paste0(out_path, "_effect_labels.csv"), sep = ",", row.names = F)
exclude_ratio_df = lapply(filter_out, function(x) x$exclude_ratio) %>% do.call(rbind, .)
write.table(exclude_ratio_df %>% round(., 3), paste0(out_path, "_exclude_ratio.csv"), sep = ",", row.names = F)
thres_df = lapply(filter_out, function(x) x$thres) %>% do.call(rbind, .)
write.table(thres_df %>% round(., 1), paste0(out_path, "_thres.csv"), sep = ",", row.names = F)
n = switch(analysis_type,
           "full" = 4,
           "ac" = 3)
vicinity_slope_plot = lapply(filter_out, function(x) x$plotlist$slope)
SaveMultiPagePlot(vicinity_slope_plot, "vicinity_slope", n = n, width = 4000, height = 4000)
vicinity_elev_plot = lapply(filter_out, function(x) x$plotlist$elevation)
SaveMultiPagePlot(vicinity_elev_plot, "vicinity_elevation", n = n, width = 4000, height = 4000)
vicinity_access_plot = lapply(filter_out, function(x) x$plotlist$access)
SaveMultiPagePlot(vicinity_access_plot, "vicinity_accessibility", n = n, width = 4000, height = 4000)


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

  #load filtering threshold and find vicinity with filtering
  effect_labels = effect_labels_df[i, ]
  thres = thres_df[i, ]

  slope_exclude = switch(effect_labels$slope,
                         "Neg." = (matches$slope >= thres$slope),
                         "Pos." = (matches$slope <= thres$slope),
                         "N.S." = rep(F, nrow(matches)))
  elevation_exclude = switch(effect_labels$elevation,
                             "Neg." = (matches$elevation >= thres$elevation),
                             "Pos." = (matches$elevation <= thres$elevation),
                             "N.S." = rep(F, nrow(matches)))
  access_exclude = switch(effect_labels$access,
                             "Neg." = (matches$access >= thres$access),
                             "Pos." = (matches$access <= thres$access),
                             "N.S." = rep(F, nrow(matches)))

  vicinity = matches %>%
    dplyr::select(c("elevation", "slope", "access", "luc10", "luc0")) %>%
    mutate(exclude = (slope_exclude | elevation_exclude | access_exclude),
           lucc = paste(luc10, luc0, sep = "_"),
           defor = ifelse(lucc %in% c("1_2", "1_3", "1_4"), 1, 0),
           acd_t_10 = acd$carbon.density[match(luc10, acd$land.use.class)],
           acd_t0 = acd$carbon.density[match(luc0, acd$land.use.class)],
           c_loss = (acd_t_10 - acd_t0) / 10)

  vicinity_area = nrow(vicinity) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares
  vicinity_area_filtered = nrow(subset(vicinity, !exclude)) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares

  #sub-sample project vicinity down to around 1/10 of the smallest vicinity size of the 15 projects (2559309)
  if(nrow(vicinity) > 250000) vicinity = vicinity[sample(nrow(vicinity), 250000), ]

  #gather common project-level variables
  project_var = data.frame(project = proj_name, t0 = t0, country = country, acd_1 = acd_1, area_ha = area_ha,
                           vicinity_area = vicinity_area, vicinity_area_filtered = vicinity_area_filtered)

  #find paths to match and unmatached points in each sampled pairs
  pair_paths = list.files(pair_dirs[i], full = T) %>% str_subset(".parquet")
  matchless_ind = pair_paths %>% str_detect("matchless")
  matchless_paths = pair_paths[matchless_ind]
  matched_paths = pair_paths[!matchless_ind]

  #exit if no matches
  if(length(matched_paths) == 0) {
    return(list(vicinity = vicinity, project_var = project_var, additionality_estimates = NULL))
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
  return(list(vicinity = vicinity, project_var = project_var, additionality_estimates = additionality_estimates))
})
names(additionality_out) = projects


# Save outputs ----
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

#project vicinity
vicinity_list = lapply(additionality_out, function(x) x$vicinity)
names(vicinity_list) = projects
saveRDS(vicinity_list, paste0(out_path, "_vicinity.rds"))
#vicinity_list = read_rds(paste0(out_path, "_vicinity.rds"))


# Ex ante analysis ----
ex_ante_out = lapply(seq_along(projects), function(i) {
  proj_id = projects[i]
  area_ha = project_var$area_ha[i]
  vicinity = vicinity_list[[i]]
  if(forecast) {
    obs_val = NULL
  } else {
    obs_val = additionality_estimates[[i]] %>%
      filter(started) %>%
      mutate(c_loss = c_loss / area_ha)
  }

  CalcExAnte(proj_id = proj_id, vicinity = vicinity, obs_val = obs_val)
})
names(ex_ante_out) = projects


# Save ex ante outputs ----
#before filtering
baseline_c_loss = lapply(ex_ante_out, function(x) x$baseline_c_loss)
names(baseline_c_loss) = projects
write_rds(baseline_c_loss, paste0(out_path, "_baseline_c_loss.rds"))

#after filtering
baseline_c_loss_filtered = lapply(ex_ante_out, function(x) x$baseline_c_loss_filtered)
names(baseline_c_loss_filtered) = projects
write_rds(baseline_c_loss_filtered, paste0(out_path, "_baseline_c_loss_filtered.rds"))

plot_df = lapply(ex_ante_out, function(x) x$plot_df)
names(plot_df) = projects
write_rds(plot_df, paste0(out_path, "_plot_df.rds"))
#plot_df = read_rds(paste0(out_path, "_plot_df.rds"))

p_filter = lapply(ex_ante_out, function(x) x$p_filter)
names(p_filter) = projects
write_rds(p_filter, paste0(out_path, "_p_filter.rds"))

p_c_loss = lapply(ex_ante_out, function(x) x$p_c_loss)
names(p_c_loss) = projects
write_rds(p_c_loss, paste0(out_path, "_p_c_loss.rds"))
#p_c_loss = read_rds(paste0(out_path, "_p_c_loss.rds"))

risk_summary = data.frame(project = projects,
                          risk_overcrediting = sapply(ex_ante_out, function(x) x$risk_overcrediting),
                          risk_reversal = sapply(ex_ante_out, function(x) x$risk_reversal))
write.table(risk_summary, paste0(out_path, "_risk_summary.csv"), sep = ",", row.names = F)

#summary statistics of observed values
obs_c_loss_summ = lapply(ex_ante_out, function(x) {
  quantile(filter(x$plot_df, Type == "obs_c_loss")$Value, c(0.05, 0.1, 0.25, 0.5, 0.75), na.rm = T)
}) %>% do.call(rbind, .)
rownames(obs_c_loss_summ) = projects
write.table(obs_c_loss_summ, paste0(out_path, "_obs_c_loss.csv"), sep = ",", row.names = F)


# Plot ex ante outputs ----

#calculate SES of filtering
baseline_c_loss_ses = sapply(ex_ante_out, function(x) {
  (mean(x$baseline_c_loss_filtered$t) - mean(x$baseline_c_loss$t)) / sd(x$baseline_c_loss$t)
})

sapply(ex_ante_out, function(x) {
  c(mean(x$baseline_c_loss_filtered$t), mean(x$baseline_c_loss$t))
})


sapply(ex_ante_out, function(x) {
  c(sd(x$baseline_c_loss_filtered$t), sd(x$baseline_c_loss$t))
})

#plot difference before/after filtering
y_range = sapply(plot_df, function(x) {
  if(is.null(x)) return(NA)
  x %>%
    filter(Type %in% c("base_c_loss", "base_c_loss_filt")) %>%
    pull(Value) %>%
    range()
}) %>%
  unlist() %>%
  range(., na.rm = T)
p_filter_adj = lapply(seq_along(p_filter), function(i) {
  if(is.null(p_filter[[i]])) return(NULL)
  p_filter[[i]] +
    scale_y_continuous(limits = c(y_range[1], y_range[2] * 1.5)) +
    ggpubr::stat_compare_means(method = "t.test", aes(label = ..p.signif..),
                                 label.x = 1.5, label.y = y_range[2] * 1.1, size = 5) +
    annotate("text", x = 1.5, y = y_range[2] * 1.3,
             label = paste0("SES = ", round(baseline_c_loss_ses[i], 2)))
})
SaveMultiPagePlot(p_filter_adj, "vicinity_filtering", width = 4000, height = 4000)

#plot of distributions
y_range = sapply(plot_df, function(x) {
  if(is.null(x)) return(NA)
  x %>%
    filter(Type != "base_c_loss") %>%
    pull(Value) %>%
    range()
}) %>%
  unlist() %>%
  range(., na.rm = T)
p_c_loss_adj = lapply(p_c_loss, function(x) {
  if(is.null(x)) return(NULL)
  x + scale_y_continuous(limits = c(y_range[1], y_range[2])) +
      scale_x_discrete(labels = c("Baseline", "Observed")) +
      scale_color_manual(values = c("red", "black"),
                         labels = c("Baseline", "Observed")) +
      scale_linetype_manual(values = c(2, 1),
                            labels = c("Baseline", "Observed")) +
      labs(x = "", y = "Annual carbon loss (Mg/ha)")
})
SaveMultiPagePlot(p_c_loss_adj, n = 5, "c_loss_disribution_test", width = 4000, height = 4000)

p_c_loss_df = lapply(seq_along(plot_df), function(i) {
  plot_df[[i]] %>%
    filter(Type != "base_c_loss") %>%
    mutate(project = projects[i])
}) %>%
  do.call(rbind, .) %>%
  mutate(project = factor(project, levels = projects))

    p_c_loss = ggplot(data = p_c_loss_df, aes(x = Type, y = Value)) +
      geom_boxplot(aes(color = Type)) +
      facet_wrap(vars(project), ncol = 5) +
      scale_x_discrete(labels = c("Baseline", "Observed")) +
      scale_color_manual(values = c("red", "black"),
                         labels = c("Baseline", "Observed")) +
      scale_linetype_manual(values = c(2, 1),
                            labels = c("Baseline", "Observed")) +
      labs(x = "", y = "Annual carbon loss (Mg/ha)") +
#      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
ggsave(paste0(out_path, "_c_loss_disribution_test.png"), width = 5000, height = 5000, units = "px")

