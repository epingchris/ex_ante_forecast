# 0. Setup ----
rm(list = ls())

library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators
library(boot) #boot::boot

out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_ongoing") #where outputs are stored
fig_path = paste0("/maps/epr26/ex_ante_forecast_out/out_") #where outputs are stored
project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)
additionality_df = read.csv(paste0(out_path, "_additionality_estimates.csv"), header = T)

projects = project_var$project

#bootstrapped mean counterfactual C loss and mean additionality for ongoing projects
c_loss_ongoing_boot_df = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  area_i = project_var$area_ha[i]
  cf_loss_obs_i = additionality_df %>%
    filter(project == projects[i]) %>%
    dplyr::select(c_loss) %>%
    mutate(c_loss = c_loss / area_i)
  cf_loss_obs_boot_i = boot::boot(data = cf_loss_obs_i,
                                  statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                                  R = 1000)
  cf_loss_obs_boot_ci = boot::boot.ci(boot.out = cf_loss_obs_boot_i, type = "norm")
  cf_loss_obs_boot_df_i = data.frame(type = "cf_loss",
                                     mean = mean(as.vector(cf_loss_obs_boot_i$t)),
                                     ci_lower = cf_loss_obs_boot_ci$normal[2],
                                     ci_upper = cf_loss_obs_boot_ci$normal[3])

  add_i = additionality_df %>%
    filter(project == projects[i]) %>%
    dplyr::select(additionality) %>%
    mutate(additionality = additionality / area_i)
  add_boot_i = boot::boot(data = add_i,
                                 statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                                 R = 1000)
  add_boot_ci = boot::boot.ci(boot.out = add_boot_i, type = "norm")
  add_boot_df_i = data.frame(type = "additionality",
                             mean = mean(as.vector(add_boot_i$t)),
                             ci_lower = add_boot_ci$normal[2],
                             ci_upper = add_boot_ci$normal[3])
  b = Sys.time()
  cat(projects[i], ":", b - a, "\n")

  bootstrap_i = rbind(cf_loss_obs_boot_df_i, add_boot_df_i) %>%
    mutate(project = projects[i])
  return(bootstrap_i)
}) %>%
  do.call(rbind, .)

#Output: bootstrapped mean counterfactual C loss for ongoing projects
write.csv(c_loss_ongoing_boot_df, paste0(out_path, "_c_loss_ongoing_boot.csv"), row.names = F)
#c_loss_ongoing_boot_df = read.csv(paste0(out_path, "_c_loss_ongoing_boot.csv"), header = T)


#bootstrapped mean observed counterfactual C loss and additionality, combined with baseline
baseline_best_boot_df = read.csv("/maps/epr26/ex_ante_forecast_out/out_ongoing_baseline_best_boot.csv", header = T)
baseline_loose_boot_df = read.csv("/maps/epr26/ex_ante_forecast_out/out_ongoing_baseline_loose_boot.csv", header = T)
#baseline_offset_boot_df = read.csv("/maps/epr26/ex_ante_forecast_out/out_offset_baseline_best_boot.csv", header = T) %>%
#  mutate(type = "offset")

#c_loss_ongoing_stat = do.call(bind_rows, list(c_loss_ongoing_boot_df, baseline_best_boot_df, baseline_loose_boot_df, baseline_offset_boot_df))
c_loss_ongoing_stat = do.call(bind_rows, list(c_loss_ongoing_boot_df, baseline_best_boot_df, baseline_loose_boot_df))
c_loss_ongoing_mean = c_loss_ongoing_stat %>%
  dplyr::select(c("project", "type", "mean")) %>%
  pivot_wider(names_from = "type", values_from = "mean")


#Figure 5. show how baseline compares to counterfactual carbon loss
ggplot(data = c_loss_ongoing_mean, aes(x = best, y = cf_loss)) +
  geom_point() +
  geom_text(aes(label = project), hjust = -0.1, size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  labs(x = "Best-matched baseline carbon loss (MgC/ha/yr)",
       y = "Observed counterfactual carbon loss (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))
ggsave(paste0(fig_path, "figure5a_baseline_best_vs_c_loss_cf.png"), width = 2500, height = 2500, units = "px")

ggplot(data = c_loss_ongoing_mean, aes(x = loose, y = cf_loss)) +
  geom_point() +
  geom_text(aes(label = project), hjust = -0.1, size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
  labs(x = "Loose-matched baseline carbon loss (MgC/ha/yr)",
       y = "Observed counterfactual carbon loss (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))
ggsave(paste0(fig_path, "figure5b_baseline_loose_vs_c_loss_cf.png"), width = 2500, height = 2500, units = "px")

# ggplot(data = c_loss_ongoing_mean, aes(x = offset, y = cf_loss)) +
#   geom_point() +
#   geom_text(aes(label = project), hjust = -0.1, size = 4) +
#   geom_abline(intercept = 0, slope = 1, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 3) +
#   scale_x_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
#   scale_y_continuous(limits = c(0, 1.5), expand = c(0.01, 0.01)) +
#   labs(x = "Mean time-offsetted best-matched baseline carbon loss (MgC/ha/yr)",
#        y = "Mean observed counterfactual carbon loss (MgC/ha/yr)") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14))
# ggsave(paste0(fig_path, "figure5c_baseline_offset_vs_c_loss_cf.png"), width = 3500, height = 3500, units = "px")


# #Figure 6. show how baseline compares to additionality
# ggplot(data = c_loss_ongoing_mean, aes(x = best, y = additionality)) +
#   geom_point() +
#   geom_text(aes(label = project), hjust = -0.1, size = 4) +
#   geom_abline(intercept = 0, slope = 1, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 3) +
#   scale_x_continuous(limits = c(0, 1.1), expand = c(0.01, 0.01)) +
#   scale_y_continuous(limits = c(-0.2, 1.1), expand = c(0.01, 0.01)) +
#   labs(x = "Mean best-matched baseline carbon loss (MgC/ha/yr)",
#        y = "Mean additionality (MgC/ha/yr)") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14))
# ggsave(paste0(fig_path, "figure6a_baseline_best_vs_additionality.png"), width = 3500, height = 3500, units = "px")

# ggplot(data = c_loss_ongoing_mean, aes(x = loose, y = additionality)) +
#   geom_point() +
#   geom_text(aes(label = project), hjust = -0.1, size = 4) +
#   geom_abline(intercept = 0, slope = 1, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 3) +
#   scale_x_continuous(limits = c(0, 1.1), expand = c(0.01, 0.01)) +
#   scale_y_continuous(limits = c(-0.2, 1.1), expand = c(0.01, 0.01)) +
#   labs(x = "Mean loosely-matched baseline carbon loss (MgC/ha/yr)",
#        y = "Mean additionality (MgC/ha/yr)") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14))
# ggsave(paste0(fig_path, "figure6b_baseline_loose_vs_additionality.png"), width = 3500, height = 3500, units = "px")

# ggplot(data = c_loss_ongoing_mean, aes(x = offset, y = additionality)) +
#   geom_point() +
#   geom_text(aes(label = project), hjust = -0.1, size = 4) +
#   geom_abline(intercept = 0, slope = 1, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 3) +
#   scale_x_continuous(limits = c(0, 1.1), expand = c(0.01, 0.01)) +
#   scale_y_continuous(limits = c(-0.2, 1.1), expand = c(0.01, 0.01)) +
#   labs(x = "Mean time-offsetted best-matched baseline carbon loss (MgC/ha/yr)",
#        y = "Mean additionality (MgC/ha/yr)") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14))
# ggsave(paste0(fig_path, "figure6c_baseline_offset_vs_additionality.png"), width = 3500, height = 3500, units = "px")


#Figure 7. Comparing three types of baselines

ggplot(data = c_loss_ongoing_stat)