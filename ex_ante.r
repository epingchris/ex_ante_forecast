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
analysis_type = "epr26"

#load data
if(analysis_type == "epr26") {
  project_dir = "/maps/epr26/tmf_pipe_out/"
  proj_meta = read.csv("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv") #for t0
  proj_var = readRDS(paste0(project_dir, "project_var.rds")) #for area
  proj_estimate = readRDS(paste0(project_dir, "project_estimates.rds")) #for observed additionality
  projects = names(proj_estimate)
  proj_area = proj_var$area_ha
  plot_nrow = ceiling(length(projects) / 3)
}

if(analysis_type == "control") {
  #load controls
  t0 = 2011
  projects = c(2, 3, 4, 5, 7)
  proj_area = sapply(seq_along(projects), function(i) {
    a = st_read(paste0("/maps/epr26/tmf-data-grid/0000/0000_", i, ".geojson")) %>%
      st_area_ha()
    return(a)
  })
}

#load parameters
thres = 1 #what proportion of most inaccessibility pixels to filter out; 1 means no filtering
pr_vec = seq(0.01, 0.99, by = 0.01) #different quantiles of baseline carbon loss to test
make_plot = T

proj_list = lapply(seq_along(projects), function(i) {
  c = Sys.time()
  proj_id = projects[i]
  area_ha = proj_area[i]

  #obtain observed counterfactual C loss and additionality
  if(analysis_type == "control") {
    #@@@obs_val =
  } else if(analysis_type == "epr26") {
    obs_val = proj_estimate[[proj_id]] %>% filter(started)
  }
  obs_c_loss_ha = obs_val$c_loss / area_ha
  obs_add_ha = obs_val$additionality / area_ha

  #obtain vicinity baseline carbon loss rate
  if(analysis_type == "control") {
    acd = read.csv(paste0(project_dir, "0000_grid/", proj_id, "/0000_", proj_id, "carbon-density.csv")) #@@@to be fixed
    matches = read_parquet(paste0(project_dir, "0000_grid/", proj_id, "/0000_", proj_id, "matches.parquet")) #@@@to be fixed
  } else if(analysis_type == "epr26") {
    acd = read.csv(paste0(project_dir, proj_id, "/", proj_id, "carbon-density.csv"))
    matches = read_parquet(paste0(project_dir, proj_id, "/", proj_id, "matches.parquet"))
  }
  #carbon loss associated with LUC change undisturbed -> deforested (MgC / ha)
  acd_1 = acd %>% filter(land.use.class == 1) %>% pull(carbon.density)
  acd_3 = acd %>% filter(land.use.class == 3) %>% pull(carbon.density)
  acd_change = ifelse(length(acd_1) + length(acd_3) == 2, acd_1 - acd_3, NA)

  max_access = max(matches$access, na.rm = T)
  matches_filtered = matches %>%
    filter(access < max_access * thres) %>%
    dplyr::select(access, cpc0_u:cpc10_d) %>%
    mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
    mutate(carbon_loss_5_0 = defor_5_0 * acd_change, carbon_loss_10_5 = defor_10_5 * acd_change, carbon_loss_10_0 = defor_10_0 * acd_change)
    #convert from proportion of 30x30m2 pixels with change to MgC/ha
    #since previous data shows that the three time intervals are almost the same, only use one (t-10 to t0)

  #get vicinity area: convert from numbers of 30x30m2 pixels to hectares
  vicinity_area = nrow(matches_filtered) * 900 / 10000

  #create data for plotting
  plot_df = matches_filtered %>%
    dplyr::select(carbon_loss_5_0:carbon_loss_10_0) %>%
    pivot_longer(cols = everything(), values_to = "Value", names_to = "Period") %>%
    mutate(Type = "baseline_c_loss") %>% #Baseline carbon loss in three periods
    #add counterfactual carbon loss after the project starts
    rbind(., data.frame(Value = obs_c_loss_ha,
                        Period = "after",
                        Type = "obs_c_loss")) %>% #Observed counterfactual carbon loss (MgC/ha)
    #add additionality after the project starts
    rbind(., data.frame(Value = obs_add_ha,
                        Period = "after",
                        Type = "obs_add")) #Observed additionality (MgC/ha)

  if(make_plot) {
    plot_title = ifelse(analysis_type == "control", paste("Control", proj_id), proj_id)

    p0 = ggplot(data = plot_df %>% filter(Period != "after"), aes(Value, after_stat(density))) +
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

    p1 = ggplot(data = plot_df %>% filter(Period == "carbon_loss_10_0" | Type == "obs_c_loss"), aes(Value, after_stat(density))) +
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

    p2 = ggplot(data = plot_df %>% filter(Period == "carbon_loss_10_0" | Type == "obs_add"), aes(Value, after_stat(density))) +
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


  distr_p = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1)
  #obtain different quantiles of baseline carbon loss and calculate probability of overcrediting
  len_pr = length(pr_vec)
  carbon_loss_q = quantile(matches_filtered$carbon_loss_10_0, pr_vec)
  p_over = sapply(carbon_loss_q, function(x) length(which(obs_add_ha < x)) / length(obs_add_ha))
  over_df = data.frame(Pr = pr_vec,
                       carbon_loss_q = carbon_loss_q,
                       p_over = p_over,
                       ID = proj_id)

  d = Sys.time()
  cat(proj, ":", d - c, "\n")
  if(make_plot) {
    out_list = list(vicinity_area = vicinity_area, plot_df = plot_df, over_df = over_df, plot0 = p0, plot_distr = distr_p)
  } else {
    out_list = list(vicinity_area = vicinity_area, plot_df = plot_df, over_df = over_df)
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
  dplyr::select(c("NAME", "ID")) %>%
  mutate(Variable = "Percentile") %>%
  cbind(., lapply(proj_list, function(x) x$over_df$carbon_loss_q) %>% do.call(rbind, .))

p_over_df = proj_meta %>%
  filter(ID %in% projects) %>%
  arrange(ID) %>%
  dplyr::select(c("NAME", "ID")) %>%
  mutate(Variable = "Overcrediting risk") %>%
  cbind(., lapply(proj_list, function(x) x$over_df$p_over) %>% do.call(rbind, .))

forecast_p_list = lapply(seq_along(proj_list), function(i) {
  x = proj_list[[i]]$over_df %>%
    pivot_longer(cols = carbon_loss_q:p_over, names_to = "Variable", values_to = "Value") %>%
    mutate(Pr = as.numeric(Pr), ID = as.factor(ID))

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

plot0_list = lapply(proj_list, function(x) x$plot0)
plot0_all = ggpubr::ggarrange(plotlist = plot0_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("hist_c_loss_periods_", thres_text, ".png"), plot0_all, width = 3000, height = 5000, units = "px", bg = "white")

plot_distr_list = lapply(proj_list, function(x) x$plot_distr)
plot_distr_all = ggpubr::ggarrange(plotlist = plot_distr_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("hist_c_loss_", thres_text, ".png"), plot_distr_all, width = 5000, height = 5000, units = "px", bg = "white")

plot_over_all = ggpubr::ggarrange(plotlist = forecast_p_list, ncol = 3, nrow = plot_nrow, common.legend = T, legend = "bottom")
ggsave(paste0("plot_distributions.png"), plot_over_all, width = 5000, height = 5000, units = "px", bg = "white")

thres_text = gsub("\\.", "_", thres)
write.table(basic_df, paste0("out_basic_info.csv"), sep = ",", row.names = F)
write.table(forecast_df, paste0("out_forecast_", thres_text, ".csv"), sep = ",", row.names = F)