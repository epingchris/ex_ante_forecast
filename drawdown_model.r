#multiple regression

rm(list = ls())

library(tidyverse)
library(corrplot)
library(magrittr)

path = paste0("/maps/epr26/tmf_pipe_out/")
fit_type = "normal"

#project-level variables:
#mean and sd of normal distribution used to fit additionality before/after project start
fit_before_param = readRDS(file.path(paste0(path, fit_type, "/drawdown_fit_before_param.rds")))
fit_after_param = readRDS(file.path(paste0(path, fit_type, "/drawdown_fit_after_param.rds")))

#project-level variables: area, ACD of undisturbed forest, country, ecoregion
#pair-level variables: min/median/max of elevation, slope, accessibility, cpc0/5/10_u, cpc0/5/10_d across the value in 100 pairs
#the value of each pair taken as the median of all pixels in each pair (control + treat)

project_var = readRDS(file.path(paste0(path, "project_var.rds")))

project_var = merge(fit_before_param, fit_after_param, by = "project") %>%
  merge(., project_var, by = "project") %>%
  rename(sd.before = sd.x, mean.after = mean.y, sd.after = sd.y) %>%
  dplyr::select(mean.after, sd.after, everything()) %>%
  mutate(mean.after.per.ha = mean.after / area_ha, sd.after.per.ha = sd.after / area_ha)



#1. Choose all independent variables
project_var_indep = project_var %>%
  dplyr::select(!c("project", "country", "mean.x")) %>%
  dplyr::select(sd.before:slope_max) %>%
  filter(complete.cases(.)) %>%
  as.matrix()

cor_var_indep = cor(project_var_indep)
testRes = cor.mtest(project_var_indep, conf.level = 0.95)
png("plot_corr_independent_var.png", width = 800, height = 1080)
corrplot(cor_var_indep, method = "circle", p.mat = testRes$p, sig.level = 0.05) %>%
  corrRect(c(seq(0, 27, by = 3) + 1, 30))
dev.off()

#2. Remove min and max statistics because they are always highly collinear with median statistics
project_var_indep = project_var %>%
  dplyr::select(!c("project", "country", "mean.x")) %>%
  dplyr::select(sd.before:area_ha, ends_with("median")) %>%
  filter(complete.cases(.)) %>%
  as.matrix()

cor_var_indep = cor(project_var_indep)
testRes = cor.mtest(project_var_indep, conf.level = 0.95)
png("plot_corr_independent_var_trimmed.png", width = 800, height = 1080)
corrplot(cor_var_indep, method = "number", p.mat = testRes$p, sig.level = 0.05, order = 'hclust', addrect = 5)
dev.off()
#cpc_u and cpc_d variables are highly negatively correlated
#cpc_u also positively correlated with higher inaccessbility
#elevation correlated with slope
#variability of carbon flux before project is positively correlated to project area

#3. Select ACD of undisturbed, area, accessibility, cpc0_d/u, elevation, slope
#both elevation and slope are meaningful so retained despite collinearity
project_var_indep = project_var %>%
  dplyr::select(acd_u, area_ha, accessibility_median, cpc0_d_median, cpc0_u_median,
                elevation_median, slope_median) %>%
  filter(complete.cases(.)) %>%
  as.matrix()

cor_var_indep = cor(project_var_indep)
testRes = cor.mtest(project_var_indep, conf.level = 0.95)
png("plot_corr_independent_var_final.png", width = 800, height = 1080)
corrplot(cor_var_indep, method = "number", p.mat = testRes$p, sig.level = 0.05)
dev.off()

(cor_var_res = cor(project_var$mean.after.per.ha, project_var$sd.after.per.ha))
(testRes = cor.mtest(project_var[, c("mean.after.per.ha", "sd.after.per.ha")], conf.level = 0.95))
#Mean.after.per.ha and sd.after.per.ha are positively correlated, but the correlation coefficient is not high (0.3547696)

project_var_tofit = project_var %>%
  mutate(mean.after.per.ha = mean.after / area_ha, sd.after.per.ha = sd.after / area_ha) %>%
  dplyr::select(mean.after.per.ha, sd.after.per.ha, acd_u, area_ha, accessibility_median,
                cpc0_d_median, cpc0_u_median, elevation_median, slope_median)

drawdown_lm = lm(cbind(mean.after.per.ha, sd.after.per.ha) ~ acd_u + ., data = project_var_tofit)
summary(drawdown_lm)
#positive effect of area, but it's to be expected because it measures project-level mean additionality; should be per ha
#higher inaccessibility, lower additionality?!
#higher cpc0_u (proportion of undisturbed forest cover at beginning of project), higher additionality