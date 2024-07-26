# 0. Setup ----
rm(list = ls())

#Load packages
# install.packages(c("arrow","configr", "tidyverse", "magrittr", "sf", "magrittr", "MatchIt",
#                    "rnaturalearthdata", "configr", "terra", "pbapply", "cleangeo", "doParallel",
#                    "foreach", "readr", "lwgeom", "rnaturalearth", "stars"), depends = TRUE)

library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators
library(sf) #sf::st_area
library(ggpubr) #ggpubr::ggarrange
library(units) #units::set_units
library(arrow) #arrow::read_parquet
library(pryr) #pryr::object_size
library(cowplot)
#library(MatchIt) #MatchIt::matchit
#library(boot) #boot::boot
#library(parallel) #parallel::mclapply
#library(microbenchmark) #microbenchmark::microbenchmark

#Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = FALSE)

#source("functions.r") #cpc_rename, tmfemi_reformat, simulate_area_series, make_area_series, assess_balance, make_match_formula
#source("AdditionalityPair.r")
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

  #only keep projects  who have finished running ("additionality.csv" exists)
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
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion)
})

#list containing set M of every project
setM = lapply(seq_along(projects), function(i) {
  luc_t_10 = paste0("luc_", t0[i] - 10)
  luc_t0 = paste0("luc_", t0[i])
  read_parquet(m_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion)
})

#data frame of project-level variables
project_var = data.frame(project = proj_name, t0 = t0, country = country, area_ha = area_ha, acd_undisturbed = acd_undisturbed)


# B. Predict deforestation probability of baseline pixels using logistic regression ----
#models fitted to K or to M are similar, so only results fitted to M are used (model_by = "M")
var_vec = c("slope", "elevation", "access")
var_label = c("Slope (dgree)", "Elevation (meter)", "Remoteness (minutes)")
predict_defor_out = lapply(seq_along(projects), function(i) PredictDefor(proj_id = projects[i], t0 = t0[i], acd = acd[[i]], K = setK[[i]], M = setM[[i]]))
names(predict_defor_out) = projects
saveRDS(predict_defor_out, paste0(out_path, "_predict_defor_out.rds"))
#predict_defor_out = read_rds(paste0(out_path, "_predict_defor_out.rds"))

#Output: predicted baseline deforestation probability
baseline = lapply(predict_defor_out, function(x) x$baseline)
names(baseline) = projects
saveRDS(baseline, paste0(out_path, "_baseline.rds"))
#baseline = read_rds(paste0(out_path, "_baseline.rds"))

#Output: predicted project deforestation probability
project_defor_prob = lapply(predict_defor_out, function(x) x$project_defor_prob)
names(project_defor_prob) = projects
saveRDS(project_defor_prob, paste0(out_path, "_project_defor_prob.rds"))
#project_defor_prob = read_rds(paste0(out_path, "_project_defor_prob.rds"))

#Output: total range of predicted baseline deforestation probability
range_defor_prob = sapply(baseline, function(x) range(x$defor_prob)) #0 - 0.58

#Output: sensitivity and specificity of each logistic regression model
threshold = seq(0, 1, by = 0.01)
roc_out = lapply(seq_along(predict_defor_out), function(i) {
  # specificity = rep(NA, length(threshold))
  # sensitivity = rep(NA, length(threshold))
  baseline = predict_defor_out[[i]]$baseline
  roc_baseline = roc(data = baseline, response = "defor", predictor = "defor_prob")
  best_threshold = coords(roc_baseline, x = "best")
  # p_roc = plot(roc_baseline)
  # for(j in seq_along(threshold)) {
  #   baseline$defor_pred = ifelse(baseline$defor_prob > threshold[j], 1, 0)
  #   specificity[j] = nrow(subset(baseline, defor_pred == 0 & defor == 0)) / nrow(subset(baseline, defor == 0))
  #   sensitivity[j] = nrow(subset(baseline, defor_pred == 1 & defor == 1)) / nrow(subset(baseline, defor == 1))
  # }
  df_roc = data.frame(threshold = roc_baseline$thresholds,
                      specificity = roc_baseline$specificities,
                      sensitivity = roc_baseline$sensitivities) %>%
    mutate(inv_specificity = 1 - specificity)
  return(list(roc_baseline = roc_baseline, best_threshold = best_threshold, df_roc = df_roc))
})
names(roc_out) = projects
saveRDS(roc_out, paste0(out_path, "_roc_out.rds"))
#roc_out = read_rds(paste0(out_path, "_roc_out.rds"))

#Visualisation: ROC curves for each logistic model
p_roc = lapply(seq_along(roc_out), function(i) {
  best_threshold = roc_out[[i]]$best_threshold
  ggplot(data = roc_out[[i]]$df_roc, aes(x = inv_specificity, y = sensitivity)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = 1 - best_threshold$specificity) +
    geom_hline(yintercept = best_threshold$sensitivity) +
    annotate("text", x = 0.75, y = 0.5,
             label = paste0("AUC: ", round(auc(roc_out[[i]]$roc_baseline), 2))) +
    annotate("text", x = 0.75, y = 0.4,
             label = paste0("Opt. thres.:\n", signif(best_threshold$threshold, 3))) +
    labs(x = "1 - Specificity", y = "Sensitivity", title = projects[i]) +
    theme_classic()
}) %>%
  cowplot::plot_grid(plotlist = ., ncol = 4) #ROC plots
p_roc
ggsave(paste0(out_path, "_roc_curves.png"), p_roc,
        width = 4000, height = 4000, units = "px", bg = "white")

#Visualisation: plots of predicted deforestation probability across each environmental variable for all projects
p_model = lapply(seq_along(predict_defor_out), function(i) {
  glm_out = predict_defor_out[[i]]$glm_out
  baseline = predict_defor_out[[i]]$baseline

  p_model_list = vector("list", 3)
  for(j in seq_along(var_vec)) {
    p_model_list[[j]] = ggplot(data = baseline, aes_string(x = var_vec[j], y = defor_prob)) +
    geom_point() +
    labs(x = var_label[j], y = "Predicted deforestation probability")
    theme_classic() +
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 14))
  }
  return(cowplot::plot_grid(plotlist = p_model_list, ncol = 3))
}) %>%
  cowplot::plot_grid(plotlist = ., nrow = 5, ncol = 1)
#@@@

#Visualisation: plots of difference in environmental variables between low vs high-risk pixels
if(visualise) {
  p_slope_by_risk = lapply(predict_defor_out, function(x) x$plotlist$slope) %>%
    cowplot::plot_grid(plotlist = ., ncol = 4)
  ggsave(paste0(out_path, "_var_by_risk_slope.png"), p_slope_by_risk,
        width = 4000, height = 4000, units = "px", bg = "white")
  p_elevation_by_risk = lapply(predict_defor_out, function(x) x$plotlist$elevation) %>%
    cowplot::plot_grid(plotlist = ., ncol = 4)
  ggsave(paste0(out_path, "_var_by_risk_elevation.png"), p_elevation_by_risk,
        width = 4000, height = 4000, units = "px", bg = "white")
  p_access_by_risk = lapply(predict_defor_out, function(x) x$plotlist$access) %>%
    cowplot::plot_grid(plotlist = ., ncol = 4)
  ggsave(paste0(out_path, "_var_by_risk_access.png"), p_access_by_risk,
        width = 4000, height = 4000, units = "px", bg = "white")
}

#Visualisation: Figure 1: plot of difference in environmental variables between low vs high-risk pixels for every project
if(visualise) {
  var_max = c(25, 3000, 2000)
  text_adjust = c(1, 100, 100)
  for(i in seq_along(var_vec)) {
    var_i = var_vec[i]
    var_label_i = var_label[i]
    df_defor_prob = lapply(seq_along(predict_defor_out), function(j) {
      out_j = predict_defor_out[[j]]$baseline
      out_df = data.frame(all = mean(out_j %>% pull(all_of(var_i)), na.rm = T),
                          low = mean(filter(out_j, risk == "low") %>% pull(all_of(var_i)), na.rm = T),
                          high = mean(filter(out_j, risk == "high") %>% pull(all_of(var_i)), na.rm = T)) %>%
        mutate(max = apply(., 1, max),
              project = projects[j])
      return(out_df)
    }) %>%
      do.call(rbind, .) %>%
      mutate(project = factor(project, levels = projects),
            rank = rank(all))

    ggplot(data = df_defor_prob, aes(x = rank)) +
      geom_point(aes(y = low, color = "blue"), size = 2) +
      geom_point(aes(y = high, color = "red"), size = 2) +
      geom_point(aes(y = all, color = "black"), size = 1) +
      geom_segment(aes(y = all, xend = rank, yend = low), color = "blue") +
      geom_segment(aes(y = all, xend = rank, yend = high), color = "red") +
      geom_text(aes(x = rank, y = max + text_adjust[i], label = project), size = 4) +
      scale_color_manual(name = "Pixel type", values = c("black", "blue", "red"),
                        labels = c("All", "Low-risk (< 1%)", "High-risk (>= 1%)")) +
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
    ggsave(paste0(out_path, "_all_var_by_risk_", var_i, ".png"), width = 2000, height = 2000, unit = "px")
  }
}

#Output: basic information about the baseline
baseline_summary = lapply(seq_along(projects), function(i) {
  data.frame(project = projects[i],
             baseline_area = nrow(predict_defor_out[[i]]$baseline),
             low_risk_ratio = predict_defor_out[[i]]$low_risk_ratio[1]) %>%
    cbind(., t(data.frame(predict_defor_out[[i]]$effects)))
}) %>%
  do.call(rbind, .)
rownames(baseline_summary) = NULL
write.table(baseline_summary, paste0(out_path, "_baseline_summary.csv"), sep = ",", row.names = F)
#baseline_summary = read.table(paste0(out_path, "_baseline_summary.csv"), header = T, sep = ",")



#@@@After bootstrapping@@@#


#Output: standardised effect size by using only high-risk pixels instead of all pixels
df_ses = data.frame(project = projects, ses = sapply(c_loss_out, function(x) x$ses)) %>%
  mutate(ses_text = paste0("SES: ", round(ses, 2)))
rownames(df_ses) = NULL
write.table(df_ses, paste0(out_path, "baseline_ses_high_vs_all.csv"), sep = ",", row.names = F)
#df_ses = read.table(paste0(out_path, "baseline_ses_high_vs_all.csv"), header = T, sep = ",")


#Visualisation: Figure 2: ratio between C loss of high-risk pixels in baseline vs. in the entire baseline
if(visualise) {
  df_c_loss_ratio = c_loss_boot %>%
    group_by(project, type) %>%
    summarise(mean = mean(val)) %>%
    ungroup() %>%
    pivot_wider(names_from = "type", values_from = "mean") %>%
    mutate(ratio = high_risk / all)

  p_c_loss_ratio = ggplot(data = df_c_loss_ratio, aes(x = all, y = ratio)) +
    geom_point() +
    geom_text(aes(x = all, y = ratio + 0.2, label = project)) +
    labs(x = "Baseline annual C loss rate (Mg/ha/yr)", y = "Ratio between C loss rate in high-risk pixels vs. all pixels") +
    theme_classic()
  ggsave(paste0(out_path, "_c_loss_ratio.png"), width = 2000, height = 2000, unit = "px")
}

#Visualisation: Figure S3: baseline C loss in all vs high-risk pixels
if(visualise) {
  y_max = max(c_loss_boot$val)
  p_c_loss_by_risk = ggplot(data = c_loss_boot, aes(x = type, y = val)) +
    geom_boxplot(aes(color = type)) +
    facet_wrap(vars(project), ncol = 5) +
    ggpubr::stat_compare_means(method = "t.test", aes(label = ..p.signif..),
                              label.x = 1.5, label.y = y_max * 1.1, size = 5) +
    geom_text(data = df_ses,
              mapping = aes(x = 1.5, y = y_max * 1.2, label = ses_text), size = 5) +
    scale_x_discrete(labels = c("All", "High-risk")) +
    scale_y_continuous(limits = c(0, y_max * 1.3)) +
    scale_color_manual(values = c("red", "blue"),
                        labels = c("All", "High-risk")) +
    labs(x = "", y = "Annual carbon loss (Mg/ha/yr)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          strip.text = element_text(size = 20),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))
  ggsave(paste0(out_path, "_baseline_c_loss_by_risk.png"), width = 4000, height = 4000, units = "px")
}