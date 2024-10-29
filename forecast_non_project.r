# 0. Setup ----
rm(list = ls())

library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators
library(boot) #boot::boot

BootDF = function(type, in_df, boot_n = 1000) {
  data.frame(type = type,
             val = as.vector(boot::boot(data = in_df,
                                        statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                                        R = boot_n)$t)) #R: number of bootstrapping repetitions
}

out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_control") #where outputs are stored
fig_path = paste0("/maps/epr26/ex_ante_forecast_out/out_") #where outputs are stored
project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)
additionality_df = read.csv(paste0(out_path, "_additionality_estimates.csv"), header = T)

projects = project_var$project

#bootstrap mean of observed C loss values for control areas (target and counterfactual pixels)
c_loss_control_boot_df = lapply(seq_along(projects), function(i) {
  a = Sys.time()
  area_i = project_var$area_ha[i]
  dat_i = additionality_df %>%
    filter(project == projects[i] & started == T) %>%
    mutate(t_loss = t_loss / area_i,
           c_loss = c_loss / area_i,
           additionality = additionality / area_i)

  p_loss_boot_i = boot::boot(data = dat_i %>% dplyr::select(t_loss),
                             statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                             R = 1000)
  p_loss_boot_ci = boot::boot.ci(boot.out = p_loss_boot_i, type = "norm")
  p_loss_boot_df_i = data.frame(type = "p_loss",
                                mean = mean(as.vector(p_loss_boot_i$t)),
                                ci_lower = p_loss_boot_ci$normal[2],
                                ci_upper = p_loss_boot_ci$normal[3])

  cf_loss_boot_i = boot::boot(data = dat_i %>% dplyr::select(c_loss),
                              statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                              R = 1000)
  cf_loss_boot_ci = boot::boot.ci(boot.out = cf_loss_boot_i, type = "norm")
  cf_loss_boot_df_i = data.frame(type = "cf_loss",
                                 mean = mean(as.vector(cf_loss_boot_i$t)),
                                 ci_lower = cf_loss_boot_ci$normal[2],
                                 ci_upper = cf_loss_boot_ci$normal[3])

  add_boot_i = boot::boot(data = dat_i %>% dplyr::select(additionality),
                          statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                          R = 1000)
  add_boot_ci = boot::boot.ci(boot.out = add_boot_i, type = "norm")
  add_boot_df_i = data.frame(type = "additionality",
                                 mean = mean(as.vector(add_boot_i$t)),
                                 ci_lower = add_boot_ci$normal[2],
                                 ci_upper = add_boot_ci$normal[3])

  add_boot_total_df_i = data.frame(type = "additionality_total",
                                   mean = add_boot_df_i$mean * area_i,
                                   ci_lower = add_boot_df_i$ci_lower * area_i,
                                   ci_upper = add_boot_df_i$ci_upper * area_i)

  b = Sys.time() # less than one second per run
  cat(projects[i], ":", b - a, "\n")

  bootstrap_i = rbind(p_loss_boot_df_i, cf_loss_boot_df_i, add_boot_df_i, add_boot_total_df_i) %>%
    mutate(project = projects[i])
  return(bootstrap_i)
}) %>%
  do.call(rbind, .)

#Output: bootstrapped observed carbon loss values for control areas (project and counterfactual pixels)
write.csv(c_loss_control_boot_df, paste0(out_path, "_c_loss_control_boot.csv"), row.names = F)
#c_loss_control_boot_df = read.csv(paste0(out_path, "_c_loss_obs_control_boot.csv"), header = T)

#Output: bootstrapped mean of C loss values of control area (target/counterfactual pixels), combined with baseline
c_loss_control_boot_mean = c_loss_control_boot_df %>%
  dplyr::select(type, mean, project) %>%
  pivot_wider(names_from = "type", values_from = "mean", id_expand = T)

continent_name = c(as = "Asia", af = "Africa", sa = "South America")
c_loss_control_boot_mean = c_loss_control_boot_mean %>%
  mutate(Continent = continent_name[str_sub(project, 1, 2)])

c_loss_control_boot_mean_merged = merge(c_loss_control_boot_mean, project_var[, c("project", "area_ha")], by = "project", all = T)

#Figure 4. show that there is no bias in additionality estimation
ggplot(data = c_loss_control_boot_mean, aes(x = p_loss, y = cf_loss)) +
  geom_point(aes(shape = Continent, color = Continent, fill = Continent), size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_shape_manual(values = c(Asia = 1, Africa = 3, `South America` = 18)) +
  scale_color_manual(values = c(Asia = "blue", Africa = "black", `South America` = "red")) +
  scale_fill_manual(values = c(Asia = NA, Africa = NA, `South America` = "red")) +
  scale_x_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) +
  labs(x = "Project carbon loss (MgC/ha/yr)",
       y = "Counterfactual carbon loss (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom")
ggsave(paste0(fig_path, "figure4_c_loss_p_vs_cf_post.png"), width = 2500, height = 2600, units = "px")


#Figure 4. show that there is no bias in additionality estimation
ggplot(data = c_loss_control_boot_mean, aes(x = cf_loss, y = additionality)) +
  geom_point(aes(shape = Continent, color = Continent, fill = Continent), size = 4) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(Asia = 1, Africa = 3, `South America` = 18)) +
  scale_color_manual(values = c(Asia = "blue", Africa = "black", `South America` = "red")) +
  scale_fill_manual(values = c(Asia = NA, Africa = NA, `South America` = "red")) +
  scale_x_continuous(limits = c(-0.2, 2.5), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(-0.7, 1), expand = c(0, 0)) +
  labs(x = "Counterfactual carbon loss (MgC/ha/yr)",
       y = "Additionality (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom")
ggsave(paste0(fig_path, "figure4a_add_vs_cf_post.png"), width = 2500, height = 2600, units = "px")

#Figure 4. show that there is no bias in additionality estimation
ggplot(data = c_loss_control_boot_mean, aes(x = cf_loss, y = additionality_total)) +
  geom_point(aes(shape = Continent, color = Continent, fill = Continent), size = 4) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(Asia = 1, Africa = 3, `South America` = 18)) +
  scale_color_manual(values = c(Asia = "blue", Africa = "black", `South America` = "red")) +
  scale_fill_manual(values = c(Asia = NA, Africa = NA, `South America` = "red")) +
  scale_x_continuous(limits = c(-0.2, 2.5), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(-45000, 35000), expand = c(0, 0)) +
  labs(x = "Counterfactual carbon loss (MgC/ha/yr)",
       y = "Total additionality (MgC/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom")
ggsave(paste0(fig_path, "figure4b_add_total_vs_cf_post.png"), width = 2500, height = 2600, units = "px")

#Figure 4. show that there is no bias in additionality estimation
ggplot(data = c_loss_control_boot_mean_merged, aes(x = area_ha, y = additionality)) +
  geom_point(aes(shape = Continent, color = Continent, fill = Continent), size = 4) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(Asia = 1, Africa = 3, `South America` = 18)) +
  scale_color_manual(values = c(Asia = "blue", Africa = "black", `South America` = "red")) +
  scale_fill_manual(values = c(Asia = NA, Africa = NA, `South America` = "red")) +
  scale_x_continuous(limits = c(0, 80000), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) +
  labs(x = "Project area (ha)",
       y = "Additionality (MgC/ha/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom")
ggsave(paste0(fig_path, "figure4c_add_vs_area.png"), width = 2500, height = 2600, units = "px")


#Figure 4. show that there is no bias in additionality estimation
ggplot(data = c_loss_control_boot_mean_merged, aes(x = area_ha, y = additionality_total)) +
  geom_point(aes(shape = Continent, color = Continent, fill = Continent), size = 4) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(Asia = 1, Africa = 3, `South America` = 18)) +
  scale_color_manual(values = c(Asia = "blue", Africa = "black", `South America` = "red")) +
  scale_fill_manual(values = c(Asia = NA, Africa = NA, `South America` = "red")) +
  scale_x_continuous(limits = c(0, 80000), expand = c(0, 0)) + #ensures no padding
  scale_y_continuous(limits = c(-45000, 35000), expand = c(0, 0)) +
  labs(x = "Project area (ha)",
       y = "Total additionality (MgC/yr)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom")
ggsave(paste0(fig_path, "figure4d_add_total_vs_area.png"), width = 2500, height = 2600, units = "px")
