rm(list = ls())

library(tidyverse)
library(configr)
library(magrittr)
library(MatchIt)
library(stars)
library(arrow)
library(parallel)
library(units)

#load functions
source("./functions.r") #cpc_rename, tmfemi_reformat

orig_dir = getwd()
setwd("/home/tws36/4c_evaluations")

# The list of projects to be run in this evaluation:
proj_meta = read.csv("./data/project_metadata/proj_meta.csv")
#proj_to_eval = read.table('./data/project_metadata/proj_to_eval.txt') %>% unlist() %>% as.numeric()
#projects_agb = read_csv('./data/GEDI/project_agb.csv')

# For testing on local
# data_suffix= '230313'
# proj_id = 1201
# site = 'Gola'

config = read.config("./config/fixed_config_sherwood.ini")
config$USERPARAMS$data_path = '/maps/pf341/tom'
#write.config(config, './config/fixed_config_tmp.ini') #error: permission denied

# Load user-defined functions that Tom wrote
sapply(list.files("./R", full.names = TRUE, pattern = '.R$'), source)

# Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = FALSE)

source("./R/scripts/setup_and_reusable/load_config.R")
# source('./R/scripts/0.2_load_project_details.R')

setwd(orig_dir)

#load the list of projects to be run
proj = "1201"
#grid_n = 13 #use nrow(remerged_grid)

#load metadata
proj_meta = read.csv("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv") #for t0
t0 = proj_meta %>% filter(ID == proj) %>% pull(t0)

proj_var = readRDS("/maps/epr26/tmf_pipe_out/project_var.rds") #for area
area_ha = proj_var %>% filter(project == proj) %>% pull(area_ha)

acd = read.csv(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "carbon-density.csv")) #for vicinity baseline carbon loss rate (MgC / ha)
acd_change = ifelse(sum(is.na(acd$carbon.density[c(1, 3)])) == 0, acd$carbon.density[1] - acd$carbon.density[3], NA)

proj_estimate = readRDS("/maps/epr26/tmf_pipe_out/project_estimates.rds") #for observed additionality
#obs_add = proj_estimate[[proj]] %>% filter(started) %>% pull(additionality) #observed additionality for entire project area


# Observed additionality for each grid ----

#find paths to match and unmatached points
myproject_path = "/maps/pf341/results/2024-january-pipeline/1201_pairs"
myproject_path = paste0("/maps/epr26/tmf_pipe_out/", proj, "/pairs")
pair_paths = list.files(myproject_path, full = TRUE)
parquet_ind = pair_paths %>% str_detect('.parquet')
matchless_ind = pair_paths %>% str_detect('matchless')
matchless_paths = pair_paths[matchless_ind & parquet_ind]
matched_paths = pair_paths[!matchless_ind & parquet_ind]

#read all points from all pairs
pts_matched = lapply(seq_along(matched_paths), function(j) {
  n_pair = matched_paths[j] %>% str_replace(paste0(myproject_path, "/"), "") %>% str_replace(".parquet", "")
  pts_pairs = read_parquet(matched_paths[j]) %>% mutate(s_id = seq_len(nrow(.)), k_id = seq_len(nrow(.)))
  unmatched_pairs = read_parquet(matchless_paths[j])

  control = pts_pairs %>%
    dplyr::select(starts_with('s_')) %>%
    rename_with(~str_replace(.x, 's_', '')) %>%
    mutate(treatment = 'control') %>%
    tmfemi_reformat(t0 = t0)

  treat = pts_pairs %>%
    dplyr::select(starts_with('k_')) %>%
    rename_with(~str_replace(.x, 'k_', '')) %>%
    mutate(treatment = 'treatment') %>%
    tmfemi_reformat(t0 = t0)

  exp_n_pairs = nrow(treat) + nrow(unmatched_pairs)
  pts_matched = rbind(treat, control) %>%
    mutate(pair = n_pair, pair_id = paste0(pair, "_", id))

  return(pts_matched)
}) %>%
  do.call(rbind, .)

#find which grid each point is located in
grid_path = paste0("/maps/epr26/tmf-data-grid/", proj)
remerged_grid = read_rds(paste0(grid_path, "/", proj, "_remerged_grid.rds")) %>% st_as_sf()
grid_n = nrow(remerged_grid)
remerged_grid$name = 1:grid_n
remerged_grid$area_ha = drop_units(st_area(remerged_grid) / 10000)
st_crs(pts_matched) = st_crs(remerged_grid)
project_pts_grid = st_join(pts_matched, remerged_grid, left = T)

#balance check
table(project_pts_grid$name, project_pts_grid$treatment) #check balance: 31 is 0, 27 is 48 for Gola

st_write(pts_matched %>% filter(treatment == "treatment"), paste0("/maps/epr26/tmf_pipe_out/orig_", proj, "pts_matched.geojson"), driver = "GeoJSON")
st_write(project_pts_grid %>% filter(treatment == "treatment"), paste0("/maps/epr26/tmf_pipe_out/orig_", proj, "pts_grid_join.geojson"), driver = "GeoJSON")

st_write(pts_matched %>% filter(treatment == "treatment"), paste0("/maps/epr26/tmf_pipe_out/", proj, "pts_matched.geojson"), driver = "GeoJSON")
st_write(project_pts_grid %>% filter(treatment == "treatment"), paste0("/maps/epr26/tmf_pipe_out/", proj, "pts_grid_join.geojson"), driver = "GeoJSON")

#calculate observed additionality in each grid
obs_add_grid = lapply(seq_len(grid_n), function(i) {
  pts_treatment_grid = project_pts_grid %>% filter(treatment == "treatment", name == i) %>% dplyr::select(-c("name", "area_ha"))
  pts_control_grid = project_pts_grid %>% filter(treatment == "control", pair_id %in% pts_treatment_grid$pair_id) %>% dplyr::select(-c("name", "area_ha"))
  pts_grid = rbind(pts_treatment_grid, pts_control_grid)

  if(nrow(pts_grid) == 0) return(NA)

  control_series = simulate_area_series(pts_grid,
                                        class_prefix, t0 = t0, match_years, match_classes,
                                        nrow(pts_grid) / 2, remerged_grid$area_ha[i],
                                        verbose = FALSE)
  y = control_series$series %>%
    merge(., acd, by.x = "class", by.y = "land.use.class", all.x = T) %>%
    mutate(carbon_content = class_area * carbon.density) %>%
    group_by(treatment, year) %>%
    summarise(carbon_content = sum(carbon_content, na.rm = T)) %>%
    ungroup()

  year = y %>% filter(treatment == 'control') %>% pull(year)
  yc = y %>% filter(treatment == 'control') %>% pull(carbon_content)
  yt = y %>% filter(treatment == 'treatment')  %>% pull(carbon_content)

  out_df = data.frame(pair = i, year = year[-1], c_loss = -diff(yc), t_loss = -diff(yt)) %>%
    mutate(additionality = c_loss - t_loss)

  return(out_df)
})

add_by_grid = obs_add_grid %>%
  do.call(rbind, .) %>%
  group_by(year) %>%
  summarise(additionality = sum(additionality)) %>%
  ungroup() %>%
  filter(!is.na(year))

add_whole = proj_estimate[["1201"]] %>%
  group_by(year) %>%
  summarise(additionality = mean(additionality),
            add.max = max(additionality),
            add.min = min(additionality)) %>%
  ungroup() %>%
  filter(!is.na(year))

ggplot(data = add_whole, aes(x = year, y = additionality)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = add.min, ymax = add.max), color = "grey", alpha = 0.5) +
  geom_line(data = add_by_grid, color = "red") +
  theme_bw()

# Vicinity baseline carbon loss rate for each grid ----
grid_list = lapply(seq_len(grid_n), function(i) {
  c = Sys.time()
  grid_path = paste0("/maps/epr26/tmf_pipe_out/", proj, "_grid/", i)

  #vicinity baseline carbon loss rate: LUC change from undisturbed to deforested
  defor = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "_grid/", i, "/", proj, "_", i, "matches.parquet")) %>%
    dplyr::select(access, cpc0_u:cpc10_d) %>%
    mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
    mutate(acd_defor_5_0 = defor_5_0 * acd_change, acd_defor_10_5 = defor_10_5 * acd_change, acd_defor_10_0 = defor_10_0 * acd_change)

  vicinity_area = nrow(defor) * 900 / 10000 #convert from number of 30x30m2 pixels to hectare

  #obtain different quantiles of baseline carbo loss
  q_vec = c(0.05, 0.25, 0.5, 0.75, 0.95)
  acd_defor_q = quantile(defor$acd_defor_5_0, q_vec)

  #calculate probability of overcrediting
  p_over = sapply(acd_defor_q, function(x) length(which(obs_add < x))) / length(obs_add)
  names(p_over) = q_vec

  d = Sys.time()
  cat(proj, ":", d - c, "\n")
  return(list(vicinity_area = vicinity_area, defor = defor, acd_defor_q = acd_defor_q, p_over = p_over))
})

defor_list = lapply(grid_list, function(x) x$defor)

defor_df = lapply(seq_along(defor_list), function(i) {
  defor = data.frame(Grid = i,
                     Value = defor_list[[i]] %>% pull(acd_defor_5_0),
                     Type = "Baseline carbon loss") %>%
          rbind(., data.frame(Grid = i,
                              Value = obs_add / area_ha,
                              Type = "Observed additionality"))

}) %>%
  do.call(rbind, .)

#@@@to be sorted out@@@#
ggplot(data = defor_df, aes(Value, after_stat(density))) +
  geom_freqpoly(aes(color = Type, linetype = Type)) +
  facet_wrap(vars(Grid)) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Baseline C loss", "Observed additionality")) +
  scale_linetype_manual(values = c(2, 1),
                     labels = c("Baseline C loss", "Observed additionality")) +
  xlab("Annual carbon flux (Mg/ha)") +
  ggtitle(proj) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical")
ggsave(paste0("hist_grid_", proj, ".png"), width = 1800, height = 2000, units = "px")

basic_df = proj_meta %>%
  filter(ID %in% projects) %>%
  mutate(project_area = proj_var %>% filter(project %in% projects) %>% pull(area_ha),
         vicinity_area = sapply(proj_list, function(x) x$vicinity_area)) #convert from number of 30x30m2 pixels to hectare
acd_defor_q_df = cbind(proj_meta %>% filter(ID %in% projects) %>% dplyr::select(c("NAME", "ID")),
                       lapply(proj_list, function(x) x$acd_defor_q) %>% do.call(rbind, .))
p_over_df = cbind(proj_meta %>% filter(ID %in% projects) %>% dplyr::select(c("NAME", "ID")),
                  lapply(proj_list, function(x) x$p_over) %>% do.call(rbind, .))

basic_df = basic_df[order(basic_df$ID), ]
acd_defor_q_df = acd_defor_q_df[order(acd_defor_q_df$ID), ]
p_over_df = p_over_df[order(p_over_df$ID), ]

thres_text = gsub("\\.", "_", thres)

write.table(basic_df, paste0("basic_df_", thres_text, ".csv"), sep = ",", row.names = F)
write.table(acd_defor_q_df, paste0("acd_defor_q_df_", thres_text, ".csv"), sep = ",", row.names = F)
write.table(p_over_df, paste0("p_over_df_", thres_text, ".csv"), sep = ",", row.names = F)

plot_list = lapply(proj_list, function(x) x$plot)
plot_nrow = ceiling(length(plot_list) / 3)
plot_all = ggpubr::ggarrange(plotlist = plot_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "right")
#ggpubr::ggexport(plot_all, filename = paste0("hist_", thres, ".png"), width = 1000, height = 800, units = "px")
ggsave(paste0("hist_", thres_text, ".png"), plot_all, width = 3000, height = 2000, units = "px", bg = "white")
#@@@to be sorted out@@@#





    control_series = simulate_area_series(pts_matched,
                                           class_prefix, t0 = t0, match_years, match_classes,
                                           exp_n_pairs, project_area_ha,
                                           verbose = FALSE)

    y = control_series$series %>%
      merge(., acd, by.x = "class", by.y = "land.use.class", all.x = T) %>%
      mutate(carbon_content = class_area * carbon.density) %>%
      group_by(treatment, year) %>%
      summarise(carbon_content = sum(carbon_content, na.rm = T)) %>%
      ungroup()

    year = y %>% filter(treatment == 'control') %>% pull(year)
    yc = y %>% filter(treatment == 'control') %>% pull(carbon_content)
    yt = y %>% filter(treatment == 'treatment')  %>% pull(carbon_content)

    out_df = data.frame(pair = j, year = year[-1], c_loss = -diff(yc), t_loss = -diff(yt)) %>%
      mutate(additionality = c_loss - t_loss)


  pair_var_df = lapply(project_estimates, function(x) x$pair_var) %>% do.call(rbind, .)

  #still need to process this to get project-level biome values
  #biome_df_list = lapply(project_estimates, function(x) x$biome_df)

  pair_var_summary = pair_var_df %>%
    group_by(var) %>%
    summarise(min = min(val), median = median(val), max = max(val)) %>%
    pivot_longer(cols = min:max, names_to = "stat", values_to = "val") %>%
    mutate(var = paste0(var, "_", stat)) %>%
    dplyr::select(c(var, val)) %>%
    pivot_wider(names_from = "var", values_from = "val")

  project_var_all = cbind(project_var, pair_var_summary) %>%
    mutate(project = proj_id)

  project_estimates = lapply(project_estimates, function(x) x$out_df) %>%
    do.call(rbind, .) %>%
    mutate(started = ifelse(year > t0, T, F))

  b = Sys.time()
  cat(b - a, "\n")
  return(list(project_estimates = project_estimates, project_var = project_var_all))







  

  # @@@grid-level independent variables: median of all pixels in each pair (control + treat), then min/median/max across 100 pairs
  # elevation, slope, accessibility, cpc0/5/10_u, cpc0/5/10_d, defor_5_0 = cpc5_u - cpc0_u, defor_10_5 = cpc10_u - cpc5_u@@@
  pair_var = rbind(control, treat) %>%
      dplyr::select(elevation:cpc10_d) %>%
      mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 2) %>%
      reframe(elevation = median(elevation),
              slope = median(slope),
              accessibility = median(accessibility),
              cpc0_u = median(cpc0_u),
              cpc0_d = median(cpc0_d),
              cpc5_u = median(cpc5_u),
              cpc5_d = median(cpc5_d),
              cpc10_u = median(cpc10_u),
              cpc10_d = median(cpc10_d),
              defor_5_0 = median(defor_5_0),
              defor_10_5 = median(defor_10_5)) %>%
      pivot_longer(cols = elevation:defor_10_5, names_to = "var", values_to = "val") %>%
      mutate(pair = j)

    biome_df = NULL
    if("biome" %in% colnames(control) & "biome" %in% colnames(treat)) {
      biome_df = rbind(control, treat) %>%
        pull(biome) %>%
        table() %>%
        as.data.frame()
    }