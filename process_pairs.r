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
analysis_type = "control" #"old_source", "full", "grid", "ac", "control"
stratified = F #stratify project pixels by accessibility
forecast = (analysis_type == "ac") #forecast; do not print observed values in output of ex ante analysis
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

# Loop through all projects ----
additionality_out = lapply(seq_along(projects), function(i) { #mclapply() does not work on Windows
  a = Sys.time()

  if(analysis_type == "old_source") {
    cat("Use new results from epr26 instead.\n")
    return(NULL)
  }

  proj_id = projects[i]
  t0 = switch(analysis_type,
              "full" = filter(proj_meta, ID == str_replace(proj_id, "a", ""))$t0,
              "grid" = filter(proj_meta, ID == "1201")$t0,
              "control" = 2011,
              "ac" = 2021)

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

  #get ACD and biome variables
  file_prefix = paste0(project_dir, proj_id, "/",
                       switch(analysis_type,
                              "full" = "",
                              "grid" = "1201_",
                              "control" = "0000_",
                              "ac" = ""),
                       proj_id)
  acd = read.csv(paste0(file_prefix, "carbon-density.csv"))
  k = read_parquet(paste0(file_prefix, "k.parquet")) #for biome of project pixels
  matches = read_parquet(paste0(file_prefix, "matches.parquet")) #for biome of matched pixels
  vicinity_area = nrow(matches) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares

  #get ACD of undisturbed forest, gather project-level variables
  acd_u = filter(acd, land.use.class == 1)$carbon.density
  if(length(acd_u) == 0) acd_u = NA

  #find paths to match and unmatached points in each sampled pairs
  pair_paths = list.files(pair_dirs[i], full = T) %>% str_subset(".parquet")
  matchless_ind = pair_paths %>% str_detect("matchless")
  matchless_paths = pair_paths[matchless_ind]
  matched_paths = pair_paths[!matchless_ind]

  #exit if no matches
  if(length(matched_paths) == 0) {
    project_var = data.frame(t0 = t0, country = country, acd_u = acd_u, area_ha = area_ha, vicinity_area = vicinity_area) %>%
      mutate(project = paste0("0000_", proj_id))
    return(list(project_var = project_var, additionality_estimates = NULL))
  }

  if(stratified) {
    #loop through all sampled pairs and only extract points for stratification
    out = lapply(seq_along(matched_paths), function(j) {
      RetrievePoints(matched_path = matched_paths[j], matchless_path = matchless_paths[j],
                           k = k, matches = matches, t0 = t0)
    })

    exp_n_pairs_out = lapply(out, function(x) x$exp_n_pairs)
    points_out = lapply(out, function(x) x$pts_matched)

    points_out_df = points_out %>% do.call(rbind, .)
  } else {
    #loop through all sampled pairs and get additionality series in each pair
    pairs_out = lapply(seq_along(matched_paths), function(j) {
      ProcessPairs(matched_path = matched_paths[j], matchless_path = matchless_paths[j],
                   k = k, matches = matches, t0 = t0, area_ha = area_ha, acd = acd, pair_id = j)
    })
  }

  #compare deforestation rate based on proportion of pixel-level LUC change and based on CPC change
  transition_10_0 = as.data.frame(table(paste(matches$luc_2001, matches$luc_2011, sep = "_")))
  #pre-project CPC change of project and matched pixels (in each year in each pair)
  #pre-project CPC change of project vicinity (in each pixel)
  #pre-project LUC change of project and matched pixels (in each year in each pair)
  #pre-project CPC change of project vicinity (in each pixel)
  #during-project LUC change of project and matched pixels (in each year in each pair)
  cpcc_lucc_long = do.call(rbind,
    list(data.frame(Type = "CPC", Period = "Pre-project", var = "obs_cpcc_treatment", val = lapply(pairs_out, function(x) x$cpcc$treatment) %>% unlist()),
         data.frame(Type = "CPC", Period = "Pre-project", var = "obs_cpcc_control", val = lapply(pairs_out, function(x) x$cpcc$control) %>% unlist()),
         data.frame(Type = "CPC", Period = "Pre-project", var = "vicinity_cpcc", val = (matches$cpc10_u - matches$cpc0_u) / 10),
         data.frame(Type = "LUC", Period = "Pre-project", var = "obs_lucc_treatment_pre", val = lapply(pairs_out, function(x) x$lucc$treatment_pre) %>% unlist()),
         data.frame(Type = "LUC", Period = "Pre-project", var = "obs_lucc_control_pre", val = lapply(pairs_out, function(x) x$lucc$control_pre) %>% unlist()),
         data.frame(Type = "LUC", Period = "Pre-project", var = "vicinity_lucc", val = sum(filter(transition_10_0, Var1 %in% c("1_2", "1_3", "1_4"))$Freq) / nrow(matches)),
         data.frame(Type = "LUC", Period = "During-project", var = "obs_lucc_treatment_during", val = lapply(pairs_out, function(x) x$lucc$treatment_during) %>% unlist()),
         data.frame(Type = "LUC", Period = "During-project", var = "obs_lucc_control_during", val = lapply(pairs_out, function(x) x$lucc$control_during) %>% unlist())))
  cpcc_lucc_long$Type = factor(cpcc_lucc_long$Type, levels = c("CPC", "LUC"))
  cpcc_lucc_long$Period = factor(cpcc_lucc_long$Period, levels = c("Pre-project", "During-project"))
  cpcc_lucc_long$var = factor(cpcc_lucc_long$var, levels = c("obs_cpcc_treatment", "obs_cpcc_control", "vicinity_cpcc",
                                                             "obs_lucc_treatment_pre", "obs_lucc_control_pre", "vicinity_lucc",
                                                             "obs_lucc_treatment_during", "obs_lucc_control_during"))

  cpcc_lucc_summ = cpcc_lucc_long %>%
    group_by(var) %>%
    reframe(min = min(val),
            q1 = quantile(val, 0.25),
            median = median(val),
            mean = mean(val),
            q3 = quantile(val, 0.75),
            max = max(val)) %>%
    mutate(project = proj_id)

  p = ggplot(data = cpcc_lucc_long, aes(x = var, y = sqrt(val * 100))) +
    geom_boxplot(aes(color = Type, linetype = Period)) +
    scale_x_discrete(labels = c("Project", "Counterfactual", "Vicinity",
                                "Project", "Counterfactual", "Vicinity",
                                "Project", "Counterfactual")) +
    scale_color_manual(values = c("red", "blue")) +
    scale_linetype_manual(values = 1:2) +
    labs(x = "Variable", y = "Sqrt(Annual forest loss(%))") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12))

  #combine pair-level variables and project-level variables
  pair_biome = lapply(pairs_out, function(x) x$biome) %>%
    unlist() %>%
    table() %>%
    which.max() %>%
    names()

  pair_var_summary = lapply(pairs_out, function(x) x$pair_var) %>%
    do.call(rbind, .) %>%
    group_by(var) %>%
    summarise(min = min(val), median = median(val), max = max(val)) %>%
    pivot_longer(cols = min:max, names_to = "stat", values_to = "val") %>%
    mutate(var = paste0(var, "_", stat)) %>%
    dplyr::select(c(var, val)) %>%
    pivot_wider(names_from = "var", values_from = "val") %>%
    mutate(biome = pair_biome)

  project_var = data.frame(t0 = t0, country = country, acd_u = acd_u, area_ha = area_ha, vicinity_area = vicinity_area) %>%
    cbind(., pair_var_summary) %>%
    mutate(project = proj_id)

  additionality_estimates = lapply(pairs_out, function(x) x$out_df) %>%
    do.call(rbind, .) %>%
    mutate(started = ifelse(year > t0, T, F))

  b = Sys.time()
  cat(proj_id, ":", b - a, "\n")
  return(list(project_var = project_var, additionality_estimates = additionality_estimates,
              plot_cpcc_lucc = p, cpcc_lucc_summ = cpcc_lucc_summ))
})


# Save additionality outputs ----
out_prefix = ifelse(analysis_type == "grid", "grid_1201", analysis_type)
out_path = paste0("/maps/epr26/tmf_pipe_out/out_", out_prefix)

#project-level variables
project_var = lapply(additionality_out, function(x) x$project_var) %>% do.call(dplyr::bind_rows, .)
saveRDS(project_var, paste0(out_path, "_project_var.rds"))

#additionality time series
additionality_estimates = lapply(additionality_out, function(x) x$additionality_estimates)
names(additionality_estimates) = projects
saveRDS(additionality_estimates, paste0(out_path, "_additionality_estimates.rds"))

#CPC-based and LUC-based deforestation rates
plot_cpcc_lucc = lapply(additionality_out, function(x) x$plot_cpcc_lucc) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = T, legend = "bottom")
ggsave(paste0(out_path, "_cpcc_lucc.png"), plot_cpcc_lucc, width = 6000, height = 3000, units = "px", bg = "white")


cpcc_lucc_summary = lapply(additionality_out, function(x) x$cpcc_lucc_summ)
names(cpcc_lucc_summary) = projects
saveRDS(cpcc_lucc_summary, paste0(out_path, "_cpcc_lucc_summary.rds"))
for(i in seq_along(cpcc_lucc_summary)) {
  write.csv(cpcc_lucc_summary[[i]], paste0("/maps/epr26/tmf_pipe_out/out_", out_prefix, "_cpcc_lucc_summ_", i, ".csv"))
}


# Run function for ex ante analysis ----

ex_ante_out = lapply(seq_along(projects), function(i) {
  proj_id = projects[i]
  area_ha = project_var$area_ha[i]
  if(is.null(additionality_estimates[[i]])) {
    return(list(vicinity_area = NA, plot_df = NULL, forecast_df = NULL,
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

  CalcExAnte(proj_id = proj_id, area_ha = area_ha, obs_val = obs_val, path = path, forecast = forecast)
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

#basic information
basic_df = project_var %>%
  dplyr::select(t0, country, area_ha, project) %>%
  filter(project %in% paste0("0000_", projects)) %>%
  mutate(vicinity_area = sapply(ex_ante_out, function(x) x$vicinity_area))
write.table(basic_df, paste0(out_path, "_basic_info.csv"), sep = ",", row.names = F)

#summary statistics of observed values
obs_c_loss_summ = lapply(ex_ante_out, function(x) {
  if(is.null(x$plot_df)) return(NA)
   summary(filter(x$plot_df, Type == "obs_c_loss")$Value, na.rm = T)
}) %>% do.call(rbind, .)
rownames(obs_c_loss_summ) = projects
write.table(obs_c_loss_summ, paste0(out_path, "_obs_c_loss.csv"), sep = ",", row.names = F)

obs_add_summ = lapply(ex_ante_out, function(x) {
  if(is.null(x$plot_df)) return(NA)
   summary(filter(x$plot_df, Type == "obs_add")$Value, na.rm = T)
}) %>% do.call(rbind, .)
rownames(obs_add_summ) = projects
write.table(obs_add_summ, paste0(out_path, "_obs_add.csv"), sep = ",", row.names = F)

#functions for multi-page plots
SaveMultiPagePlot = function(plots, suffix) {
  if(class(plots)[1] == "list") {
    lapply(seq_along(plots), function(i) {
      ggsave(paste0(out_path, "_", suffix, "_", i, ".png"), plots[[i]],
             width = 5000, height = 3000, units = "px", bg = "white")
    })
  } else {
    ggsave(paste0(out_path, "_", suffix, ".png"), plots,
           width = 5000, height = 3000, units = "px", bg = "white")
  }
}

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
plot_distr = lapply(seq_along(ex_ante_out), function(i) {
  x = ex_ante_out[[i]]
  if(is.null(x$p1)) return(NULL)
  x1 = x$p1 + ggtitle("") + scale_x_continuous(limits = c(min(x_range, na.rm = T), max(x_range, na.rm = T)))
  x2 = x$p2 + ggtitle("") + scale_x_continuous(limits = c(min(x_range, na.rm = T), max(x_range, na.rm = T)))
  x_legend = x$p_legend_grob
  p_1_2 = ggpubr::ggarrange(x1, x2, ncol = 2, nrow = 1, common.legend = T, legend = "bottom", legend.grob = x_legend)
  ggpubr::annotate_figure(p_1_2, top = ggpubr::text_grob(projects[i], face = "bold", size = 14))
}) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = 3, common.legend = T, legend = "bottom")
SaveMultiPagePlot(plot_distr, "distribution")

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