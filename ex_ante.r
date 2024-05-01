rm(list = ls())

library(tidyverse)
library(magrittr)
library(stars)
library(arrow)
library(parallel)
library(sf)

#load functions
source("./functions.r") #cpc_rename, tmfemi_reformat


# Load input data ----

#load the list of projects to be run
myproject_path = "/maps/epr26/tmf_pipe_out/"
exclude_id = c(562) # 562: not in TMF extent

pair_paths = list.files(myproject_path, full = TRUE)
matchless_ind = (pair_paths %>% str_detect("\\.") | pair_paths %>% str_detect("\\_grid") | pair_paths %>% str_detect("fit\\_distribution"))
project_dir = pair_paths[!matchless_ind]
projects = gsub(".*/(\\d+)$", "\\1", project_dir)
projects = projects[which(projects %in% c("0000", "9999") == F)] %>%
  as.numeric() %>%
  sort() %>%
  as.character()
plot_nrow = ceiling(length(projects) / 3)

#load metadata
proj_meta = read.csv("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv") #for t0
proj_var = readRDS("/maps/epr26/tmf_pipe_out/project_var.rds") #for area
proj_estimate = readRDS("/maps/epr26/tmf_pipe_out/project_estimates.rds") #for during-project additionality

proj_area = proj_var %>%
  filter(project %in% projects) %>%
  mutate(project = as.numeric(project)) %>%
  arrange(project) %>%
  pull(area_ha)
control = F

#load controls
control = T
t0 = 2011
projects = c(2, 3, 4, 5, 7)
proj_area = sapply(seq_along(projects), function(i) {
  a = st_read(paste0("/maps/epr26/tmf-data-grid/0000/0000_", i, ".geojson")) %>%
    st_area() / 10000
  return(a)
})

#load parameters
thres = 1 #what proportion of most inaccessibility pixels to filter out; 1 means no filtering
pr_vec = seq(0.01, 0.99, by = 0.01) #different quantiles of baseline carbon loss to test
make_plot = T

proj_list = lapply(seq_along(projects), function(i) {
  c = Sys.time()
  proj = projects[i]
  area_ha = proj_area[i]
  if(control) {
    #@@@obs_val =
  } else {
    obs_val = proj_estimate[[proj]] %>% filter(started)
  }

  #carbon loss associated with LUC change undisturbed -> deforested (MgC / ha)
  if(control) {
    acd = read.csv(paste0("/maps/epr26/tmf_pipe_out/0000_grid/", proj, "/0000_", proj, "-carbon-density.csv"))
  } else {
    acd = read.csv(paste0("/maps/pf341/results/live-pipeline/", proj, "-carbon-density.csv"))
  }
  acd_change = ifelse(sum(is.na(acd$carbon.density[c(1, 3)])) == 0,
                      acd$carbon.density[1] - acd$carbon.density[3], NA)

  #obtain vicinity baseline carbon loss rate
  if(control) {
    matches = read_parquet(paste0("/maps/epr26/tmf_pipe_out/0000_grid/", proj, "/0000_", proj, "matches.parquet"))
  } else {
    matches = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "matches.parquet"))
  }

  max_access = max(matches$access, na.rm = T)

  matches_filtered = matches %>%
    filter(access < max_access * thres) %>%
    dplyr::select(access, cpc0_u:cpc10_d) %>%
    mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
    mutate(carbon_loss_5_0 = defor_5_0 * acd_change, carbon_loss_10_5 = defor_10_5 * acd_change, carbon_loss_10_0 = defor_10_0 * acd_change)
#    mutate(defor = (cpc10_u - cpc0_u) / 10, carbon_loss = defor * acd_change)
    #convert from proportion of 30x30m2 pixels with change to MgC/ha
    #since previous data shows that the three time intervals are almost the same, only use one (t-10 to t0)

  vicinity_area = nrow(matches_filtered) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares
  obs_c_loss_ha = obs_val$c_loss / area_ha
  obs_add_ha = obs_val$additionality / area_ha

  c_loss_df = matches_filtered %>%
    dplyr::select(carbon_loss_5_0:carbon_loss_10_0) %>%
    pivot_longer(cols = everything(), values_to = "Value", names_to = "Period")

  plot_df = matches_filtered %>%
    dplyr::select(carbon_loss_10_0) %>%
    rename(Value = carbon_loss_10_0) %>%
    mutate(Type = "baseline_c_loss") %>% #Baseline carbon loss
    #add counterfactual carbon loss after the project starts
    rbind(., data.frame(Value = obs_c_loss_ha,
                        Type = "obs_c_loss")) %>% #Observed counterfactual carbon loss (MgC/ha)
    #add additionality after the project starts
    rbind(., data.frame(Value = obs_add_ha,
                        Type = "obs_add")) #Observed additionality (MgC/ha)

  if(make_plot) {
    plot_title = ifelse(control, paste("Control", proj), proj)

    p0 = ggplot(data = c_loss_df, aes(Value, after_stat(density))) +
      geom_freqpoly(aes(color = Period)) +
      scale_color_manual(values = c("red", "black", "blue"),
                         labels = c(expression(t [-10]*" to "*t [0]),
                                    expression(t [-10]*" to "*t [-5]),
                                    expression(t [-5]*" to "*t [0]))) +
      xlab("Annual carbon flux (Mg/ha)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "bottom",
            legend.direction = "vertical")

    p1 = ggplot(data = plot_df %>% filter(Type != "obs_add"), aes(Value, after_stat(density))) +
      geom_freqpoly(aes(color = Type, linetype = Type)) +
      scale_color_manual(values = c("red", "blue"),
                         labels = c("Baseline carbon loss",
                                  "Observed counterfactual carbon loss")) +
      scale_linetype_manual(values = c(2, 1),
                            labels = c("Baseline carbon loss",
                                       "Observed counterfactual carbon loss")) +
      xlab("Annual carbon flux (Mg/ha)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "bottom",
            legend.direction = "vertical")

    p2 = ggplot(data = plot_df %>% filter(Type != "obs_c_loss"), aes(Value, after_stat(density))) +
      geom_freqpoly(aes(color = Type, linetype = Type)) +
      scale_color_manual(values = c("red", "black"),
                         labels = c("Baseline carbon loss",
                                  "Observed additionality")) +
      scale_linetype_manual(values = c(2, 1),
                            labels = c("Baseline carbon loss",
                                     "Observed additionality")) +
      xlab("Annual carbon flux (Mg/ha)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "bottom",
            legend.direction = "vertical")
  }

  #obtain different quantiles of baseline carbon loss and calculate probability of overcrediting
  len_pr = length(pr_vec)
  carbon_loss_q = quantile(matches_filtered$carbon_loss_10_0, pr_vec)
  p_over = sapply(carbon_loss_q, function(x) length(which(obs_add_ha < x)) / length(obs_add_ha))
  over_df = data.frame(pr = pr_vec,
                       carbon_loss_q = carbon_loss_q,
                       p_over = p_over)
  over_long = data.frame(Variable = c(rep("Percentile score", len_pr), rep("Overcrediting risk", len_pr)),
                         Pr = pr_vec,
                         Value = c(carbon_loss_q, p_over),
                         ID = proj)

  d = Sys.time()
  cat(proj, ":", d - c, "\n")
  if(make_plot) {
    out_list = list(vicinity_area = vicinity_area, plot_df = plot_df, over_df = over_df, over_long = over_long, plot0 = p0, plot1 = p1, plot2 = p2)
  } else {
    out_list = list(vicinity_area = vicinity_area, plot_df = plot_df, over_df = over_df, over_long = over_long)
  }
  return(out_list)
})

#summary statistics of observed additionality
obs_add_summ = lapply(proj_list, function(x) x$plot_df %>% filter(Type == "obs_add") %>% pull(Value) %>% summary()) %>% do.call(rbind, .)
rownames(obs_add_summ) = projects
write.table(obs_add_summ, paste0("out_obs_add.csv"), sep = ",", row.names = F)

basic_df = proj_meta %>%
  filter(ID %in% projects) %>%
  arrange(ID) %>%
  mutate(project_area = proj_area,
         vicinity_area = sapply(proj_list, function(x) x$vicinity_area))

carbon_loss_q_df = proj_meta %>%
  filter(ID %in% projects) %>%
  arrange(ID) %>%
  dplyr::select(c("Name", "ID")) %>%
  mutate(Variable = "Percentile") %>%
  cbind(., lapply(proj_list, function(x) x$over_df$carbon_loss_q) %>% do.call(rbind, .))

p_over_df = proj_meta %>%
  filter(ID %in% projects) %>%
  arrange(ID) %>%
  dplyr::select(c("Name", "ID")) %>%
  mutate(Variable = "Overcrediting risk") %>%
  cbind(., lapply(proj_list, function(x) x$over_df$p_over) %>% do.call(rbind, .))

forecast_list = lapply(proj_list, function(x) x$over_long)
  mutate(Pr = as.numeric(Pr), ID = as.factor(ID))

forecast_p_list = lapply(seq_along(forecast_list), function(i) {
  x = forecast_list[[i]] %>% mutate(Pr = as.numeric(Pr))

  percentile_p = ggplot(data = x %>% filter(Variable == "Percentile")) +
    geom_line(aes(x = Pr, y = Value), color = "blue") +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(limits = c(0, 5)) +
    labs(title = projects[i], x = "Percentage", y = "C loss percentile") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))

  overcredit_p = ggplot(data = x %>% filter(Variable == "Overcrediting risk")) +
    geom_line(aes(x = Pr, y = Value), color = "red") +
    geom_hline(yintercept = c(0.05, 0.1), linetype = 2) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "", x = "Percentage", y = "Overcrediting risk") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))

  out_p = ggpubr::ggarrange(percentile_p, overcredit_p, ncol = 2, nrow = 1)
  return(out_p)
})
plot_over_all = ggpubr::ggarrange(plotlist = forecast_p_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("pr_overcrediting.png"), plot_over_all, width = 5000, height = 5000, units = "px", bg = "white")

thres_text = gsub("\\.", "_", thres)
write.table(basic_df, paste0("out_basic_info.csv"), sep = ",", row.names = F)
write.table(forecast_df, paste0("out_forecast_", thres_text, ".csv"), sep = ",", row.names = F)

plot0_list = lapply(proj_list, function(x) x$plot0)
plot0_all = ggpubr::ggarrange(plotlist = plot0_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("hist_c_loss_periods_", thres_text, ".png"), plot0_all, width = 3000, height = 5000, units = "px", bg = "white")
plot1_list = lapply(proj_list, function(x) x$plot1)
plot1_all = ggpubr::ggarrange(plotlist = plot1_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("hist_c_loss_", thres_text, ".png"), plot1_all, width = 3000, height = 5000, units = "px", bg = "white")
plot2_list = lapply(proj_list, function(x) x$plot2)
plot2_all = ggpubr::ggarrange(plotlist = plot2_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("hist_add_", thres_text, ".png"), plot2_all, width = 3000, height = 5000, units = "px", bg = "white")