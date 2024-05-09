rm(list = ls())

library(tidyverse)
library(magrittr)
library(stars)
library(arrow)
library(parallel)
library(sf)
library(pryr) #pryr::object_size

#load functions
#source("./functions.r") #cpc_rename, tmfemi_reformat
source("./CalcExAnte.r")

#load parameters
pr_vec = seq(0.01, 0.99, by = 0.01) #different quantiles of baseline carbon loss to test
analysis_type = "grid" #epr26, grid, control

#load input data
if(analysis_type == "epr26") {
  project_dir = "/maps/epr26/tmf_pipe_out/"
  proj_estimate = readRDS(paste0(project_dir, "project_estimates.rds")) #for observed additionality
  projects = names(proj_estimate)
  projects_bare = projects %>% str_remove("a")

  proj_info = read.csv("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv") %>%
    filter(ID %in% projects_bare) %>%
    arrange(ID) %>%
    mutate(area_ha = readRDS(paste0(project_dir, "project_var.rds")) %>% pull(area_ha))
} else if(analysis_type == "grid") {
  project_dir = "/maps/epr26/tmf_pipe_out/"
  proj_estimate = readRDS(paste0(project_dir, "grid_1201_project_estimates.rds")) #for observed additionality
  projects = 1:10
  proj_info = readRDS("/maps/epr26/tmf_pipe_out/grid_1201_project_var.rds")
} else if(analysis_type == "control") { #controls
  projects = c(2, 3, 4, 5, 7)
  proj_area = sapply(seq_along(projects), function(i) {
    a = st_read(paste0("/maps/epr26/tmf-data-grid/0000/0000_", i, ".geojson")) %>%
      st_area_ha()
    return(a)
  })
  proj_info = data.frame(ID = projects,
                         area_ha = proj_area)
}
plot_nrow = ceiling(length(projects) / 3)

#run core function
if(analysis_type == "grid") {
  proj_list = lapply(seq_along(projects), function(i) {
    proj_id = projects[i]
    area_ha = proj_info$area_ha[i]
    obs_val = proj_estimate[[proj_id]] %>%
      filter(started) %>%
      mutate(c_loss = c_loss / area_ha, additionality = additionality / area_ha)
    path = paste0(project_dir, "1201_grid/", proj_id, "/1201_", proj_id)
    CalcExAnte(proj_id = proj_id, area_ha = area_ha, obs_val = obs_val, path = path)
  })
} else if(analysis_type == "control") {
  proj_list = lapply(seq_along(projects), function(i) {
    proj_id = paste0("0000_", projects[i])
    area_ha = proj_info$area_ha[i]
    obs_val = proj_estimate[[proj_id]] %>%
      filter(started) %>%
      mutate(c_loss = c_loss / area_ha, additionality = additionality / area_ha)
    path = paste0(project_dir, proj_id, "/", proj_id)
    CalcExAnte(proj_id = proj_id, area_ha = area_ha, obs_val = obs_val, path = path)
  })
} else if(analysis_type == "epr26") {
  proj_list = lapply(seq_along(projects), function(i) {
    proj_id = projects[i]
    area_ha = proj_info$area_ha[i]
    obs_val = proj_estimate[[proj_id]] %>%
      filter(started) %>%
      mutate(c_loss = c_loss / area_ha, additionality = additionality / area_ha)
    path = paste0(project_dir, proj_id, "/", proj_id)
    CalcExAnte(proj_id = proj_id, area_ha = area_ha, obs_val = obs_val, path = path)
  })
  names(proj_list) = projects
  write_rds(proj_list, "/maps/epr26/tmf_pipe_out/out_ex_ante.rds")
}

#summary statistics of observed additionality
obs_add_summ = lapply(proj_list, function(x) {
  if(is.null(x$plot_df)) return(NA)
  x$plot_df %>% filter(Type == "obs_add") %>% pull(Value) %>% summary()
}) %>% do.call(rbind, .)
rownames(obs_add_summ) = projects
write.table(obs_add_summ, paste0("out_obs_add.csv"), sep = ",", row.names = F)

basic_df = proj_info %>%
  mutate(vicinity_area = sapply(proj_list, function(x) x$vicinity_area))
write.table(basic_df, paste0("out_basic_info.csv"), sep = ",", row.names = F)

# carbon_loss_q_df = proj_info %>%
#   dplyr::select(c("NAME", "ID")) %>%
#   mutate(Variable = "Percentile") %>%
#   cbind(., lapply(proj_list, function(x) {
#     if(is.null(x$over_df)) return(NA)
#     return(x$over_df$carbon_loss_q)
#   }) %>%
#     do.call(rbind, .))

# p_over_df = proj_info %>%
#   dplyr::select(c("NAME", "ID")) %>%
#   mutate(Variable = "Overcrediting risk") %>%
#   cbind(., lapply(proj_list, function(x) {
#     if(is.null(x$over_df)) return(NA)
#     return(x$over_df$p_over)
#   }) %>%
#     do.call(rbind, .))
# write.table(forecast_df, paste0("out_forecast.csv"), sep = ",", row.names = F)

plot_period = lapply(proj_list, function(x) x$p0) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("plot_c_loss_periods.png"), plot_period, width = 5000, height = 5000, units = "px", bg = "white")

plot_distr = lapply(proj_list, function(x) ggpubr::ggarrange(x$p1, x$p2, ncol = 2, nrow = 1)) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("plot_distribution.png"), plot_distr, width = 5000, height = 5000, units = "px", bg = "white")

plot_forecast = lapply(proj_list, function(x) ggpubr::ggarrange(x$p_perc, x$p_overcredit, ncol = 2, nrow = 1)) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("plot_forecast.png"), plot_forecast, width = 5000, height = 5000, units = "px", bg = "white")
