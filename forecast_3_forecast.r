# 0. Setup ----
rm(list = ls())

#Load packages
library(tidyverse) #ggplot2, dplyr, and stringr used in plotPlacebo/plotBaseline.r: tibble to store labels with bquote()
library(magrittr) #pipe operators
#library(sf) #st_drop_geometry() used in GetCarbonLoss.r (runs on GDAL 3.10)
#library(arrow) #read_parquet()
#library(MatchIt) #matchit(), used in the customised function AssessBalance()
library(patchwork)

options(dplyr.summarise.inform = F) #remove dplyr summarise grouping message because it prints a lot

#Define input variables needed to read TMF implementation output and other data

# It requires the following input variables to read TMF implementation output and other data.
# All variables are vectors containing one value for each project to be analysed:

#1. projects: an index of all projects to be analysed
# This should correspond to the filenames of the shapefiles and to the -p argument in the implementation code
# It is usually be the ongoing projects' VCS ID or customised (e.g. prefixed series of integers)

#2. project_out_dirs: absolute paths of the directories containing all implementation outputs
# The directory should containing pairs of parquet files with the same file name, with and without the "_matchless" suffix.
# (typically "/pairs/xxx.parquet" and  "/pairs/xxx_matchless.parquet")
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
#9. OPTIONAL: proj_ID: this is only used to remove the trailing "a" that I added to the filename of some shapefiles that needed fixing,
# So that I can retrieve country and t0 information from the the proj_info data frame
#10. out_path: absolute paths of the directory where outputs are to be saved; include file prefix if desired


# Load data ----
analysis_type = "ongoing"
in_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored
project_var = read.csv(paste0(in_path, "_project_var.csv"), header = T)
projects = project_var$ID
t0_vec = project_var$t0
area_ha_vec = project_var$area_ha
cdens_list = project_var %>%
  dplyr::select(ID, cdens_1:cdens_6) %>%
  pivot_longer(cdens_1:cdens_6, names_to = "land.use.class", names_prefix = "cdens_", values_to = "carbon.density") %>%
  split(f = .$ID)

fig_path = paste0("/maps/epr26/ex_ante_forecast_out/out_") #where figures are stored
out_path = in_path

additionality_boot_df = read.csv(paste0(in_path, "_boot_additionality.csv"), header = T)
ante_project_boot_df = read.csv(paste0(in_path, "_boot_project_closs_rate.csv"), header = T)
ante_region_boot_df = read.csv(paste0(in_path, "_boot_regional_closs_rate.csv"), header = T)

additionality_boot_list = group_split(additionality_boot_df, project)
ante_project_boot_list = group_split(ante_project_boot_df, project)
ante_region_boot_list = group_split(ante_region_boot_df, project)

#express year in relative terms to project start
additionality_boot_rel = lapply(seq_along(projects), function(i) {
  additionality_boot_list[[i]] %<>%
    mutate(year = year - t0_vec[i])
}) %>%
  list_rbind()

ante_project_boot_rel = lapply(seq_along(projects), function(i) {
  ante_project_boot_list[[i]] %<>%
    mutate(year = year - t0_vec[i])
}) %>%
  list_rbind()

ante_region_boot_rel = lapply(seq_along(projects), function(i) {
  ante_region_boot_list[[i]] %<>%
    mutate(year = year - t0_vec[i])
}) %>%
  list_rbind()


# Generate forecasts from different intervals ----
#simple forecasts
forecast_prj = data.frame(project_used = numeric(), project = numeric(), year = numeric(),
                          forecast = numeric(), observed = numeric(), type = character())
forecast_reg = data.frame(region_used = numeric(), project = numeric(), year = numeric(),
                          forecast = numeric(), observed = numeric(), type = character())
for(i in -10:-1) {
  for(k in seq_along(projects)) {
    project_k = projects[k]

    #observed rates
    observed = additionality_boot_rel %>%
      filter(year <= 10 & project == project_k) %>%
      pull(mean)
    if(length(observed) < 10) observed = c(observed, rep(NA, 10 - length(observed)))

    #simple forecasts using project rates
    rate_project = ante_project_boot_rel %>%
      filter(year == i & project == project_k) %>%
      pull(mean)
    prj_i = data.frame(project_used = i, project = project_k, year = 1:10,
                       forecast = rate_project, observed = observed[1:10], type = "project")
    forecast_prj = rbind(forecast_prj, prj_i)

    #simple forecasts using regional rates
    rate_region = ante_region_boot_rel %>%
      filter(year == i & project == project_k) %>%
      pull(mean)
    reg_i = data.frame(region_used = i, project = project_k, year = 1:10,
                       forecast = rate_region, observed = observed[1:10], type = "region")
    forecast_reg = rbind(forecast_reg, reg_i)
  }
}

#hybrid forecasts
forecast_hyb = data.frame(project_used = numeric(), region_used = numeric(), project = numeric(), year = numeric(),
                          forecast = numeric(), observed = numeric(), type = character())
for(i in -10:-1) {
  for(j in -10:-1) {
    for(k in seq_along(projects)) {
      project_k = projects[k]

      #observed rates
      observed = additionality_boot_rel %>%
        filter(year <= 10 & project == project_k) %>%
        pull(mean)
      if(length(observed) < 10) observed = c(observed, rep(NA, 10 - length(observed)))

      #hybrid forecasts using project and regional rates
      rate_project = ante_project_boot_rel %>%
        filter(year == i & project == project_k) %>%
        pull(mean)
      rate_region = ante_region_boot_rel %>%
        filter(year == j & project == project_k) %>%
        pull(mean)
      rate_hybrid = ((rate_project ^ (9:0)) * (rate_region ^ (0:9))) ^ (1 / 10)

      hyb_i = data.frame(project_used = i, region_used = j, project = project_k, year = 1:10,
                         forecast = rate_hybrid, observed = observed[1:10], type = "hybrid")
      forecast_hyb = rbind(forecast_hyb, hyb_i)
    }
  }
}


# Summarise r2 across projects ----
forecast_prj_summ = forecast_prj %>%
  group_by(project_used, year) %>%
  summarise(r2 = summary(lm(observed ~ forecast))$r.squared) %>%
  ungroup()
forecast_reg_summ = forecast_reg %>%
  group_by(region_used, year) %>%
  summarise(r2 = summary(lm(observed ~ forecast))$r.squared) %>%
  ungroup()
forecast_hyb_summ = forecast_hyb %>%
  group_by(project_used, region_used, year) %>%
  summarise(r2 = summary(lm(observed ~ forecast))$r.squared) %>%
  ungroup()

write.csv(forecast_prj, paste0(out_path, "_forecast_1_prj.csv"), row.names = F)
write.csv(forecast_reg, paste0(out_path, "_forecast_2_reg.csv"), row.names = F)
write.csv(forecast_hyb, paste0(out_path, "_forecast_3_hyb.csv"), row.names = F)
write.csv(forecast_prj_summ, paste0(out_path, "_forecast_summ_1_prj.csv"), row.names = F)
write.csv(forecast_reg_summ, paste0(out_path, "_forecast_summ_2_reg.csv"), row.names = F)
write.csv(forecast_hyb_summ, paste0(out_path, "_forecast_summ_3_hyb.csv"), row.names = F)


# Plot forecasting performances ----
forecast_aggr = rbind(forecast_prj_summ %>% dplyr::select(year, r2) %>% mutate(type = "Project"),
                      forecast_reg_summ %>% dplyr::select(year, r2) %>% mutate(type = "Region"),
                      forecast_hyb_summ %>% dplyr::select(year, r2) %>% mutate(type = "Mixed"))

forecast_aggr_summ = forecast_aggr %>%
  group_by(type, year) %>%
  summarise(mean = mean(r2, na.rm = T),
            min = min(r2, na.rm = T),
            max = max(r2, na.rm = T),
            lower = quantile(r2, 0.025, na.rm = T),
            upper = quantile(r2, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(type = factor(type, levels = c("Project", "Region", "Mixed")))

axis_label_prj = expression(paste("Start of historical period for project rates (years before ", italic(t[0]), ")", sep = ""))
axis_label_reg = expression(paste("Start of historical period for regional rates (years before ", italic(t[0]), ")", sep = ""))
axis_label_target = expression(paste("End of target period (years after ", italic(t[0]), ")", sep = ""))

#Figure 5. r2 mean and range for each type of forecasts across all target periods
ggplot(data = forecast_aggr_summ, aes(x = year, y = mean)) +
  geom_line(aes(color = type), linewidth = 2) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#40B0A6", "#CDAC60", "#9467BD")) +
  scale_fill_manual(values = c("#40B0A6", "#CDAC60", "#9467BD")) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  scale_y_continuous(breaks = seq(0, 0.75, 0.25), labels = seq(0, 0.75, 0.25)) +
  labs(x = axis_label_target,
       y = expression(R^2),
       color = "Forecast type",
       fill = "Forecast type") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 26, hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.key.size = unit(1.5, "cm"))
ggsave(paste0(fig_path, "forecast_0_overall_r2.png"),
       width = 35, height = 30, unit = "cm")

#Supplementary: determine historical period for project-based forecasts
ggplot(data = forecast_prj_summ, aes(x = project_used, y = year)) +
  geom_tile(aes(fill = r2)) +
  geom_text(aes(label = round(r2, 3)), color = "white", size = 6) +
  scale_fill_gradient(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.25), low = "black", high = "#40B0A6") +
  scale_x_continuous(breaks = -10:-1, labels = 10:1, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:10, expand = c(0, 0)) +
  labs(title = "Project-based forecasts",
       x = expression(paste("Start of historical period (years before ", italic(t[0]), ")", sep = "")),
       y = expression(paste("End of forecasted period (years after ", italic(t[0]), ")", sep = "")),
       fill = expression(R^2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"))
ggsave(paste0(fig_path, "forecast_1_prj.png"),
       width = 30, height = 30, unit = "cm")

#Supplementary: determine historical period for region-based forecasts
ggplot(data = forecast_reg_summ, aes(x = region_used, y = year)) +
  geom_tile(aes(fill = r2)) +
  scale_fill_gradient(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.25), low = "black", high = "#CDAC60") +
  scale_x_continuous(breaks = -10:-1, labels = 10:1, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:10, expand = c(0, 0)) +
  geom_text(aes(label = round(r2, 3)), color = "white", size = 6) +
  labs(title = "Region-based forecasts",
       x = expression(paste("Start of historical period (years before ", italic(t[0]), ")", sep = "")),
       y = expression(paste("End of forecasted period (years after ", italic(t[0]), ")", sep = "")),
       fill = expression(R^2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"))
ggsave(paste0(fig_path, "forecast_2_reg.png"),
       width = 30, height = 30, unit = "cm")

#determine historical periods for mixed forecasts (using project and regional rates)

#average r2 over forecasted periods longer than 5 years
forecast_hyb_summ_yr = forecast_hyb_summ %>%
  filter(year >= 5) %>%
  group_by(project_used, region_used) %>%
  summarise(r2 = mean(r2, na.rm = T)) %>%
  ungroup() %>%
  mutate(rank = rank(-r2))

#find maximum of average r2
r2_max = forecast_hyb_summ_yr %>%
  filter(r2 == max(r2, na.rm = T))
#project_used = -4 and region_used = -1
#but actually as long as project_used is not -1 and region_used is 1, performance is good
r2_top10 = forecast_hyb_summ_yr %>%
  filter(rank <= 10)

#Figure 6. Average R2 for all combinations of historical periods for the mixed forecasting approach
ggplot(data = forecast_hyb_summ_yr, aes(x = project_used, y = region_used, fill = r2)) +
  geom_tile() +
  geom_tile(data = r2_top10, aes(x = project_used, y = region_used, fill = r2), color = "white", linewidth = 2, alpha = 0, width = 1, height = 1) +
  geom_tile(data = r2_max, aes(x = project_used, y = region_used, fill = r2), color = "red", linewidth = 2, alpha = 0, width = 1, height = 1) +
  scale_fill_gradient(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.25), low = "black", high = "#9467BD") +
  scale_x_continuous(breaks = -10:-1, labels = 10:1, expand = c(0, 0)) +
  scale_y_continuous(breaks = -10:-1, labels = 10:1, expand = c(0, 0)) +
  geom_text(aes(label = round(r2, 3)), color = "white", size = 6.5) +
  labs(title = "Mixed forecasts", x = axis_label_prj, y = axis_label_reg,
       fill = expression(paste("Average ", R^2, sep = ""))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"))
ggsave(paste0(fig_path, "forecast_3_hyb.png"),
        width = 35, height = 30, unit = "cm")

#observed vs forecasted for the best mixed forecast method
forecast_best = forecast_hyb %>%
  filter(project_used == max_r2$project_used & region_used == max_r2$region_used) %>%
  mutate(ratio = observed / forecast)

forecast_best_summ = matrix(NA, ncol = 5, nrow = length(5:10))
forecast_best_plot = vector("list", length(5:10))
for(i in seq_along(5:10)) {
  yr_i = seq(5, 10)[i]
  dat_i = filter(forecast_best, year == yr_i)
  lm_out = lm(observed ~ forecast, data = dat_i)
  lm_r2 = summary(lm_out)$r.squared
  lm_coef = coef(lm_out)
  ratio_p5 = quantile(dat_i$ratio, 0.05, na.rm = T)
  ratio_min = min(dat_i$ratio, na.rm = T)
  dat_i = dat_i %>%
    dplyr::select(-type) %>%
    rename(forecast_Original = forecast) %>%
    mutate(forecast_Adjusted = forecast_Original * ratio_min)

  dat_i_comp = dat_i %>%
    filter(!is.na(forecast_Original) & !is.na(observed))
  mae_orig = Metrics::mae(dat_i_comp$observed, dat_i_comp$forecast_Original)
  mae_adj = Metrics::mae(dat_i_comp$observed, dat_i_comp$forecast_Adjusted)
  forecast_best_summ[i, ] = c(yr_i, lm_r2, mae_orig, ratio_min, mae_adj)

  dat_long = dat_i %>%
    pivot_longer(all_of(c("forecast_Original", "forecast_Adjusted")), names_to = "type", names_prefix = "forecast_", values_to = "forecast") %>%
    mutate(type = factor(type, levels = c("Original", "Adjusted")))

  forecast_best_plot[[i]] = ggplot(data = dat_long) +
    geom_point(aes(x = forecast, y = observed, shape = type, size = type)) +
    annotate(geom = "text", x = 2.2, y = 0.9, size = 8,
             label = bquote(paste(R^2, " = ", .(round(lm_r2, 3))))) +
    annotate(geom = "text", x = 2.2, y = 0.7, size = 8,
             label = bquote(paste(Ratio[min], " = ", .(round(ratio_min, 3))))) +
    annotate(geom = "text", x = 2.2, y = 0.5, size = 8,
             label = bquote(paste("Original MAE: ", .(round(mae_orig, 3))))) +
    annotate(geom = "text", x = 2.2, y = 0.3, size = 8,
             label = bquote(paste("Adjusted MAE: ", .(round(mae_adj, 3))))) +
    annotate(geom = "text", x = 2.2, y = 0.1, size = 8,
             label = bquote(paste("Observed = Forecast * ", .(round(lm_coef[2], 3)), " + ", .(round(lm_coef[1], 3))))) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_abline(intercept = lm_coef[1], slope = lm_coef[2], color = "red", linetype = "dotted") +
    scale_shape_manual(name = "Forecast type", values = c("Original" = 20, "Adjusted" = 1)) +
    scale_size_manual(name = "Forecast type", values = c("Original" = 2, "Adjusted" = 3)) +
    labs(title = paste0(yr_i, "-year target period"),
         x = bquote(paste("Forecasted carbon loss rate (MgC ", ha^-1, " ", yr^-1, ")")),
         y = bquote(paste("Observed carbon loss rate (MgC ", ha^-1, " ", yr^-1, ")"))) +
    scale_x_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 0.5)) +
    scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "cm"),
          plot.title = element_text(size = 28, hjust = 0.5),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          axis.ticks = element_blank(),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.key.size = unit(1.5, "cm"))
}
forecast_best_summ = forecast_best_summ %>%
  as.data.frame() %>%
  set_colnames(c("year", "r2", "mae_orig", "ratio_min", "mae_adj"))

#Figure 7. Forecasted vs observed values for the 5-year and 10-year target periods
plot_5_10 = (forecast_best_plot[[1]] + forecast_best_plot[[6]]) +
  plot_layout(axis_titles = "collect_y", guides = "collect") &
  theme(legend.position = "bottom")
ggsave(paste0(fig_path, "forecast_4_vs_obs.png"),
          width = 50, height = 30, unit = "cm")


#Figure 8. R2, ratio_min, original MAE and adjusted MAE over target periods
plot_r2 = ggplot(data = forecast_best_summ, aes(x = year, y = r2)) +
  geom_line() +
  labs(x = axis_label_target, y = expression(paste("Average ", R^2, sep = ""))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(size = 28, hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))

plot_ratiomin = ggplot(data = forecast_best_summ, aes(x = year, y = ratio_min)) +
  geom_line() +
  labs(x = axis_label_target, y = expression(Ratio[min])) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(size = 28, hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))

forecast_best_summ_long = forecast_best_summ %>%
  dplyr::select(year, mae_orig, mae_adj) %>%
  pivot_longer(mae_orig:mae_adj, names_to = "type", names_prefix = "mae_", values_to = "mae") %>%
  mutate(type = factor(type, levels = c("orig", "adj"), labels = c("Original", "Adjusted")))
plot_mae = ggplot(data = forecast_best_summ_long, aes(x = year, y = mae, group = type, color = type, linewidth = type)) +
  geom_line() +
  scale_color_manual(values = c("red", "black")) +
  scale_linewidth_manual(values = c(1, 2)) +
  labs(x = axis_label_target, y = "Mean absolute error (MAE)", color = "Forecast type", linewidth = "Forecast type") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(size = 28, hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20))

plot_all = (plot_r2 / plot_ratiomin / plot_mae) +
  plot_layout(axis_titles = "collect_x")
ggsave(paste0(fig_path, "forecast_5_vs_target_period.png"),
          width = 40, height = 40, unit = "cm")





# Supplementary: se GAM to look at how forecast r2 changes with forecasting parameters ----
forecast_gam = mgcv::gam(r2 ~ s(project_used, bs = "tp", k = 10) +
                              s(region_used, bs = "tp", k = 10) +
                              s(year, bs = "tp", k = 10), data = forecast_hyb_summ)
gam.check(forecast_gam)
summary(forecast_gam)
AIC(forecast_gam)


# -- Digression into discrepancy in GAM visualisation --

#variable "project_used": partial effects (line) and residuals (points) have the same pattern but are shifted
gam_plot = plot(forecast_gam, scale = 0, residuals = T, select = 1, cex = 5, shift = coef(forecast_gam)[1])

(gam_plot_gg = visreg(forecast_gam, "project_used", type = "conditional", partial = T, jitter = F, gg = T) +
    scale_x_continuous(breaks = c(-10, -5, -1)) +
    scale_y_continuous(limits = c(0.15, 0.6)) +
    labs(title = "Project carbon loss",
         x = "Interval start (year)",
         y = "Partial effect") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 22, hjust = 0.5),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)))


#variable "year": here the predicted values dips into the negative when "year" is low, even though all observed values are positive
gam_plot = plot(forecast_gam, scale = 0, residuals = T, select = 3, cex = 5, shift = coef(forecast_gam)[1])

(gam_plot_gg = visreg(forecast_gam, "year", type = "conditional", partial = T, jitter = F, gg = T) +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    labs(title = "Effect of forecasted interval",
         x = "Interval end (year)",
         y = "Partial effect") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 24, hjust = 0.5),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)))

ddd = forecast_summ %>%
    mutate(pred = forecast_gam$fitted.values,
           resid = forecast_gam$residuals,
           resid_from_plot = gam_plot[[1]]$p.resid,
           resid_from_visreg = gam_plot_gg$data$y,
           resid_diff = resid_from_visreg - resid_from_plot)
ggplot(data = ddd, aes(x = project_used, y = resid)) + geom_point()
ggplot(data = ddd, aes(x = project_used, y = resid_from_plot)) + geom_point()
ggplot(data = ddd, aes(x = project_used, y = resid_from_plot - resid)) +
    geom_line() +
    scale_y_continuous(limits = c(-0.05, 0.05), breaks = c(-0.05, seq(-0.04, 0.04, by = 0.02), 0.05))
ggplot(data = ddd, aes(x = project_used, y = resid_from_visreg)) + geom_point()
ggplot(data = ddd, aes(x = project_used, y = resid_from_visreg - resid_from_plot)) +
    geom_line()
ggplot(data = ddd, aes(x = project_used, y = resid_diff)) + geom_point()
plot(forecast_gam, scale = 0, select = 1, ylim = c(-0.05, 0.05))

pred_df = data.frame(project_used = -5.5, region_used = -5.5, year = seq(1, 10, len = 500))
pred_response = predict(forecast_gam, pred_df, type = "response")
pred_df = pred_df %>%
    mutate(response = pred_response)
ggplot(data = pred_df, aes(x = year, y = response)) +
    geom_line() +
    scale_y_continuous(limits = c(0.15, 0.6))

#residuals in gamObject (fitted GAM object): working residuals for the fitted model
#forecast_gam$fitted.values + forecast_gam$residuals = observed values
#p.resid in plot.gam output object: partial residuals (working residuals + partial effect)
#the shift argument in plot.gam doesn't change output of partial residuals (p.resid)
#visreg generates the predictions as from plot.gam and predict.gam(type = "response") and mean covariates
#the partial residuals in visreg also follow the same pattern as p.resid in plot.gam output object, but with a difference of 0.4138:
#my guess is that visreg includes the intercept coefficient 0.3215 but p.resid does not (even though it is included on the plot due to the shift argument),
#but removing the intercept coefficient leaves the difference of 0.09237 (which seems to be the gap shown on the graphs)
#where does this come from?

# -- End of digression --

# Plot partial effects of fitted GAM (ignoring the digression) ----
plot_proj = visreg(forecast_gam, "project_used", type = "conditional", partial = T, jitter = F, gg = T) +
    scale_x_continuous(breaks = c(-10, -5, -1)) +
    scale_y_continuous(limits = c(0.15, 0.6)) +
    labs(title = "Project carbon loss",
         x = "Interval start (year)",
         y = "Partial effect") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 22, hjust = 0.5),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
plot_reg = visreg(forecast_gam, "region_used", type = "conditional", partial = T, jitter = T, gg = T) +
    scale_x_continuous(breaks = c(-10, -5, -1)) +
    labs(title = "Regional carbon loss",
         x = "Interval start (year)",
         y = "Partial effect") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 22, hjust = 0.5),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
plot_year = visreg(forecast_gam, "year", type = "conditional", partial = T, jitter = T, gg = T) +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    labs(title = "Effect of forecasted interval",
         x = "Interval end (year)",
         y = "Partial effect") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 24, hjust = 0.5),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))

title_forecasting = ggplot() +
    ggtitle("Effect of forecasting interval chosen to estimate:") +
    theme_void() + theme(plot.title = element_text(size = 24, hjust = 0.5, margin = margin(t = 10)))

#use patchwork to combine plots
design = "AAAAAAAAAA
          BBBBBBBBBB
          ##CCCCC###"
plot_forecasting = (plot_proj + plot_reg) +
    plot_layout(axis_titles = "collect_y")
(plot_all = title_forecasting / plot_forecasting / plot_year +
    plot_layout(design = design, heights = c(0.02, 1, 1)))
ggsave(plot_all, filename = paste0(fig_path, "forecast_5_gam_summ.png"), width = 30, height = 30, unit = "cm")