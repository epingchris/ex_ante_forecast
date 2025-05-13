# 0. Setup ----
rm(list = ls())

#Load packages
library(tidyverse) #ggplot2, dplyr, and stringr used in plotPlacebo/plotBaseline.r: tibble to store labels with bquote()
library(magrittr) #pipe operators
#library(sf) #st_drop_geometry() used in GetCarbonLoss.r (runs on GDAL 3.10)
#library(arrow) #read_parquet()
#library(MatchIt) #matchit(), used in the customised function AssessBalance()
library(boot) #boot
library(scales) #trans_break
library(Metrics) #CalcError.r: rmse, mae
library(patchwork)

options(dplyr.summarise.inform = F) #remove dplyr summarise grouping message because it prints a lot

#Load pre-defined functions
source("BootOut.r")
source("CalcError.r")
source("plotPlacebo.r")
source("plotBaseline.r")

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

#Express year in relative terms to project start
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

# A. Generate forecasts from different intervals ----
forecast_df = data.frame(project = numeric(), year = numeric(), forecast = numeric(), observed = numeric())
for(i in -10:-1) {
  for(j in -10:-1) {
    for(k in seq_along(projects)) {
      project_k = projects[k]

      #simple forecasts using project or regional rates
      rate_project = ante_project_boot_rel %>%
        filter(year == i & project == project_k) %>%
        pull(mean)
      rate_region = ante_region_boot_rel %>%
        filter(year == j & project == project_k) %>%
        pull(mean)

      #hybrid forecasts using project and regional rates
      forerate_hybrid = ((rate_project ^ (9:0)) * (rate_region ^ (0:9))) ^ (1 / 10)
      observed = additionality_boot_rel %>%
        filter(year <= 10 & project == project_k) %>%
        pull(mean)
      if(length(observed) < 10) observed = c(observed, rep(NA, 10 - length(observed)))

      forecast_i = data.frame(project_used = i, region_used = j, project = project_k, year = 1:10,
                              forecast_proj = rate_project, forecast_reg = rate_region, forecast_hybr = forerate_hybrid,
                              observed = observed[1:10])
      forecast_df = rbind(forecast_df, forecast_i)
    }
    cat("Forecast from project rate in year", i, "and regional rate in year", j, "done\n")
  }
}

GoF = function(obs, pred) {
  sse = sum((obs - pred) ^ 2, na.rm = T)
  sst = sum((obs - mean(obs, na.rm = T)) ^ 2, na.rm = T)
  if(sst != 0) {
    return(1 - sse / sst)
  } else {
    return(NA)
  }
}

# Summarise across projects ----
forecast_summ = forecast_df %>%
  group_by(project_used, region_used, year) %>%
  summarise(r2_proj = summary(lm(observed ~ forecast_proj))$r.squared,
            r2_reg = summary(lm(observed ~ forecast_reg))$r.squared,
            r2_hybr = summary(lm(observed ~ forecast_hybr))$r.squared) %>%
  ungroup()

write.csv(forecast_df, paste0(out_path, "_forecast.csv"), row.names = F)
write.csv(forecast_summ, paste0(out_path, "_forecast_r2.csv"), row.names = F)

# Plot how forecast r2 changes with forecasting parameters ----
ggplot(data = forecast_summ) +
  geom_boxplot(aes(x = project_used, y = r2_hybr, group = project_used)) +
  scale_x_continuous(breaks = c(-10, -5, -1), labels = c(-10, -5, -1)) +
  labs(title = "Forecasting interval for project carbon loss", x = "Interval start (year)", y = expression(R^2)) +
  theme_bw() +
  theme(plot.title = element_text(size = 24),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid = element_blank())
ggsave(paste0(fig_path, "forecast_1_r2_vs_project_used.png"),
       width = 30, height = 20, unit = "cm")

ggplot(data = forecast_summ) +
  geom_boxplot(aes(group = region_used, y = r2_hybr)) +
  scale_x_continuous(breaks = c(-10, -5, -1), labels = c(-10, -5, -1)) +
  labs(title = "Forecasting interval for regional carbon loss", x = "Interval start (year)", y = expression(R^2)) +
  theme_bw() +
  theme(plot.title = element_text(size = 24),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid = element_blank())
ggsave(paste0(fig_path, "forecast_2_r2_vs_region_used.png"),
       width = 30, height = 20, unit = "cm")

ggplot(data = forecast_summ) +
  geom_boxplot(aes(group = year, y = r2_hybr)) +
  scale_x_continuous(breaks = c(1, 5, 10), labels = c(1, 5, 10)) +
  labs(title = "Forecasted interval", x = "Interval end (year)", y = expression(R^2)) +
  theme_bw() +
  theme(plot.title = element_text(size = 24),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid = element_blank())
ggsave(paste0(fig_path, "forecast_3_r2_vs_year.png"),
       width = 30, height = 20, unit = "cm")


# Use GAM to look at how forecast r2 changes with forecasting parameters ----
forecast_gam = mgcv::gam(r2_hybr ~ s(project_used, bs = "tp", k = 10) +
                                   s(region_used, bs = "tp", k = 10) +
                                   s(year, bs = "tp", k = 10), data = forecast_summ)
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
ggsave(plot_all, filename = paste0(fig_path, "forecast_4_gam_summ.png"), width = 30, height = 30, unit = "cm")


# ---- OLD RESULTS TO BE REMOVED ----

if(analysis_type == "control") {
  closs_placebo = read.csv(paste0(out_path, "_post_p_c_loss.csv"), header = T)
  closs_ante_regional = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
  closs_ante_project = read.csv(paste0(out_path, "_baseline_best_boot.csv"), header = T)
  closs_ante_historical = read.csv(paste0(out_path, "_baseline_lagged_new_boot.csv"), header = T)
  closs_post_counterfactual = read.csv(paste0(out_path, "_post_cf_c_loss.csv"), header = T)

  # Create plots
  out_regional = plotPlacebo(dat = rbind(closs_placebo, closs_ante_regional),
                             label_to_x = "p_c_loss", col = "#006CD1")
  out_project = plotPlacebo(dat = rbind(closs_placebo, closs_ante_project),
                            label_to_x = "p_c_loss", col = "#40B0A6")
  out_historical = plotPlacebo(dat = rbind(closs_placebo, closs_ante_historical),
                               label_to_x = "p_c_loss", col = "#CDAC60")
  out_post = plotPlacebo(dat = rbind(closs_placebo, closs_post_counterfactual),
                         label_to_x = "p_c_loss")

  p_regional = out_regional$p
  p_project = out_project$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_historical = out_historical$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_post = out_post$p

  # Create column labels
  col_ante = ggplot() +
    ggtitle(bquote("A." ~ italic("Ex ante") ~ "forecasts")) +
    theme_void() + theme(plot.title = element_text(size = 64, hjust = 0.5, margin = margin(t = 10)))
  col_post = ggplot() +
    ggtitle(bquote("B." ~ italic("Ex post") ~ "match")) +
    theme_void() + theme(plot.title = element_text(size = 64, hjust = 0.5, margin = margin(t = 10)))

  col_regional = ggplot() +
    ggtitle("Recent regional") +
    theme_void() + theme(plot.title = element_text(size = 60, hjust = 0.5, margin = margin(t = 10)))
  col_project = ggplot() +
    ggtitle("Recent project") +
    theme_void() + theme(plot.title = element_text(size = 60, hjust = 0.5, margin = margin(t = 10)))
  col_historical = ggplot() +
    ggtitle("Time-shifted historical match") +
    theme_void() + theme(plot.title = element_text(size = 60, hjust = 0.5, margin = margin(t = 10)))

# Figure 3a. ex ante forecasts with placebo areas
  plots = (p_regional + p_project + p_historical) +
    plot_layout(nrow = 1, guides = "collect", axis_titles = "collect")
  cols = (col_regional + col_project + col_historical) +
    plot_layout(nrow = 1)
  plot_complete =
    col_ante / cols / plots +
    plot_layout(nrow = 3, heights = c(0.02, 0.01, 1))
  ggsave(paste0(fig_path, "figure3a_placebo_ex_ante_new.png"), width = 48, height = 20, units = "in")

# Figure 3b. ex post estimation with placebo areas
  plot_complete =
    col_post / p_post +
    plot_layout(nrow = 2, heights = c(0.02, 1))
  ggsave(paste0(fig_path, "figure3b_placebo_ex_post.png"), width = 18, height = 20, units = "in")

# Figure 4. t-test result

  # load data
  out_defor = read.csv("/maps/jh2589/eping/project_rates.csv") %>%
    mutate(mae_region = abs(regional_exante_rate - k_rate),
           mae_project = abs(k_exante_rate - k_rate),
           mae_shifted = abs(s_exante_rate - k_rate),
           mae_expost = abs(s_expost_rate - k_rate),
           mape_region = mae_region / k_rate,
           mape_project = mae_project / k_rate,
           mape_shifted = mae_shifted / k_rate,
           mape_expost = mae_expost / k_rate)

  mape_region = mean(out_defor$mape_region)
  mape_project = mean(out_defor$mape_project)
  mape_shifted = mean(out_defor$mape_shifted)
  mape_expost = mean(out_defor$mape_expost)

  # Create a named vector for labels of different types
  label_types = c(
    "Recent regional" = "Recent regional",
    "Recent project" = "Recent project",
    "Time-shifted historical match" = "Time-shifted\nhistorical match",
    "Ex post match" = expression(italic("Ex post") ~ "\nmatch")
  )
  colors = c("#006CD1", "#40B0A6", "#CDAC60", "#C13C3C")

  t_out_list = list(t.test(out_defor$regional_exante_rate, out_defor$k_rate, paired = T),
                    t.test(out_defor$k_exante_rate, out_defor$k_rate, paired = T),
                    t.test(out_defor$s_exante_rate, out_defor$k_rate, paired = T),
                    t.test(out_defor$s_expost_rate, out_defor$k_rate, paired = T))
  t_out_df = lapply(t_out_list, function(x) {
    data.frame(estimate = as.numeric(x$estimate),
               ci_lower = as.numeric(x$conf.int[1]),
               ci_upper = as.numeric(x$conf.int[2]),
               pval = as.numeric(x$p.value))
    }) %>%
    list_rbind() %>%
    mutate(type = factor(names(label_types), levels = names(label_types)))

  p_t_out = ggplot(data = t_out_df, aes(x = type)) +
    geom_col(aes(y = estimate, fill = type)) +
    geom_hline(yintercept = 0, linewidth = 1, linetype = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.5) +
    geom_text(aes(y = pmax(0, ci_upper) + 0.05, label = ifelse(round(pval, 2) == 0, "< 0.01", round(pval, 2))), size = 5) +
    scale_fill_manual(values = colors, guide = NULL) +
    scale_x_discrete(labels = label_types) +
    labs(x = "Estimate type", y = "Estimate-to-observed difference") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.text.x = element_text(margin = margin(t = 30), hjust = 0.5, vjust = 0.5))
  ggsave(paste0(fig_path, "figure5_t_test_out.png"), width = 10, height = 10, units = "in")
  }


scale_color_manual(values = c(best = "#40B0A6", loose = "#006CD1", lagged = "#CDAC60")) +

# Figure in SI: ongoing projects
if(analysis_type == "ongoing") {
  closs_ante_regional = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
  closs_ante_project = read.csv(paste0(out_path, "_baseline_best_boot.csv"), header = T)
  closs_ante_historical = read.csv(paste0(out_path, "_baseline_lagged_boot.csv"), header = T)
  closs_post_counterfactual = read.csv(paste0(out_path, "_post_cf_c_loss.csv"), header = T)

  # Create plots
  out_regional = plotPlacebo(dat = rbind(closs_post_counterfactual, closs_ante_regional),
                             label_to_x = "cf_c_loss", col = "#006CD1")
  out_project = plotPlacebo(dat = rbind(closs_post_counterfactual, closs_ante_project),
                            label_to_x = "cf_c_loss", col = "#40B0A6")
  out_historical = plotPlacebo(dat = rbind(closs_post_counterfactual, closs_ante_historical),
                               label_to_x = "cf_c_loss", col = "#CDAC60")

  p_regional = out_regional$p
  p_project = out_project$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_historical = out_historical$p +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # Create column labels
  col_ante = ggplot() +
    ggtitle(bquote(italic("Ex ante") ~ "forecasts")) +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))

  col_regional = ggplot() +
    ggtitle("Recent regional") +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))
  col_project = ggplot() +
    ggtitle("Recent project") +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))
  col_historical = ggplot() +
    ggtitle("Time-shifted historical match") +
    theme_void() + theme(plot.title = element_text(size = 40, hjust = 0.5))

  plots = (p_regional + p_project + p_historical) +
    plot_layout(nrow = 1, guides = "collect", axis_titles = "collect")
  cols = (col_regional + col_project + col_historical) +
    plot_layout(nrow = 1)
  plot_complete =
    col_ante / cols / plots +
    plot_layout(nrow = 3, heights = c(0.02, 0.01, 1))
  ggsave(paste0(fig_path, "figure_s4_ongoing_projects_new.png"), width = 48, height = 16, units = "in")