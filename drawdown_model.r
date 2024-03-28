#multiple regression

rm(list = ls())

library(tidyverse)

path = paste0("/maps/epr26/tmf_pipe_out/")
fit_type = "normal"

fit_before_param = readRDS(file.path(paste0(path, fit_type, "/drawdown_fit_before_param.rds")))
fit_after_param = readRDS(file.path(paste0(path, fit_type, "/drawdown_fit_after_param.rds")))
project_var = readRDS(file.path(paste0(path, "project_var.rds")))

project_var = merge(fit_before_param, fit_after_param, by = "project") %>%
  merge(., project_var, by = "project") %>%
  dplyr::select(!c(mean.x)) %>%
  rename(sd.before = sd.x, mean.after = mean.y, sd.after = sd.y)

drawdown_lm = lm(data = project_var, mean.after + sd.after ~ sd.before + acd_u + area_ha + country + 
                           accessibility_median + cpc0_d_median + cpc0_u_median + cpc10_d_median +
                           cpc10_u_median + cpc5_d_median + cpc5_u_median + elevation_median + slope_median)
