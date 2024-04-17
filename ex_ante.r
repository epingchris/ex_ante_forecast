rm(list = ls())

library(tidyverse)
library(magrittr)
library(stars)
library(arrow)

#load functions
source("./functions.r") #cpc_rename, tmfemi_reformat

#load the list of projects to be run
projects = c("1112", "1113", "1201", "1396", "1399")

#load metadata
proj_meta = read.csv("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv")
proj_var = readRDS("/maps/epr26/tmf_pipe_out/project_var.rds")

#obtain baseline deforestation rate
proj_list = lapply(seq_along(projects), function(i) {
  c = Sys.time()
  proj = projects[i]
  matches = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "matches.parquet"))

  #get carbon loss associated with undisturbed -> deforested LUC change
  acd = read.csv(paste0('/maps/pf341/results/live-pipeline/', proj, '-carbon-density.csv'))
  acd_change = NA
  if (sum(is.na(acd$carbon.density[c(1, 3)])) == 0) acd_change = acd$carbon.density[1] - acd$carbon.density[3]

  #get project vicinity deforestation rate
  matches_defor = matches %>%
    mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
    mutate(acd_defor_5_0 = defor_5_0 * acd_change, acd_defor_10_5 = defor_10_5 * acd_change, acd_defor_10_0 = defor_10_0 * acd_change)

  df = matches_defor %>%
    dplyr::select(acd_defor_5_0:acd_defor_10_0) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Period", values_to = "value")

  p = ggplot(data = df, aes(value, after_stat(density), color = Period)) +
    geom_freqpoly() +
    scale_color_manual(values = c("blue", "black", "red"), guide = "none") +
    labs(x = "Annual C loss (Mg/ha)", y = "Count") +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(paste0("hist_", proj, ".png"), width = 840, height = 840, units = "px")

  d = Sys.time()
  cat(proj, ":", d - c, "\n") #7.111846 when mclapply outside, 2.133847 when mclapply inside
  return(matches_defor)
})

basic_df = proj_meta %>%
  filter(ID %in% projects) %>%
  mutate(vicinity_area = sapply(proj_list, nrow) * 900 / 10000, #convert from number of 30x30m2 pixels to hectare
         project_area = proj_var %>% filter(project %in% projects) %>% pull(area_ha),
         mean = sapply(proj_list, function(x) mean(x$acd_defor_5_0)),
         sd = sapply(proj_list, function(x) sd(x$acd_defor_5_0)))

stat_df = proj_meta %>%
  filter(ID %in% projects) %>%
  dplyr::select(c("NAME", "ID")) %>%
  mutate(mean = sapply(proj_list, function(x) mean(x$acd_defor_5_0)),
         q5 = sapply(proj_list, function(x) quantile(x$acd_defor_5_0, 0.05)),
         q25 = sapply(proj_list, function(x) quantile(x$acd_defor_5_0, 0.25)),
         median = sapply(proj_list, function(x) median(x$acd_defor_5_0)),
         q75 = sapply(proj_list, function(x) quantile(x$acd_defor_5_0, 0.75)),
         q95 = sapply(proj_list, function(x) quantile(x$acd_defor_5_0, 0.95)),
         sd = sapply(proj_list, function(x) sd(x$acd_defor_5_0)))

write.table(basic_df, "basic_df.csv", sep = ",", row.names = F)
write.table(stat_df, "stat_df.csv", sep = ",", row.names = F)
