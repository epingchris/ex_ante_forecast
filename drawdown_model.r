#multiple regression

rm(list = ls())

library(tidyverse)
library(corrplot)
library(magrittr)

path = paste0("/maps/epr26/tmf_pipe_out/")

CorrSummary = function(x, method = "number", rect = NULL, plot_name, ...) {
  corRes = cor(x)
  testRes = cor.mtest(x, conf.level = 0.95)
  png(plot_name, width = 800, height = 1080, type = "cairo")
  corrplot(corRes, method = method, p.mat = testRes$p, sig.level = 0.05, ...) %>%
    {if(!is.null(rect)) corrRect(., rect) else .}
  dev.off()
}


# Load variables ----
#project-level variables:
#mean and sd of sampled distributions of c_loss/t_loss/additionality before/after project start
fit_param = readRDS(file.path(paste0(path, "fit_param.rds")))
loss_before_param = fit_param$before %>% dplyr::select(!starts_with("drawdown"))
loss_after_param = fit_param$after

#project-level variables: area, ACD of undisturbed forest, country, ecoregion
#pair-level variables: min/median/max of elevation, slope, accessibility, cpc0/5/10_u, cpc0/5/10_d across the value in 100 pairs
#the value of each pair taken as the median of all pixels in each pair (control + treat)
project_var = readRDS(file.path(paste0(path, "project_var.rds")))

var_all = merge(project_var, loss_before_param, by = "project") %>%
  merge(., loss_after_param, by = "project") %>%
  filter(complete.cases(.)) %>% #remove rows with NAs
  rename_with(~str_replace(.x, '\\.x', '\\.before')) %>%
  rename_with(~str_replace(.x, '\\.y', '\\.after'))


# Test collinearity of independent variables ----

#1. All variables
var_indep = var_all %>%
  dplyr::select(!country) %>% #remove categorical variables
  dplyr::select(acd_u:t_loss.sd.before) #remove categorical variables
cor_var_indep = CorrSummary(var_indep, method = "circle", rect = c(1, seq(3, 33, by = 3), 35),
                            plot_name = "plot_corr_independent_var.png")

#2. Remove min and max because collinear with median, and c_loss and loss mean because collinear with t_loss median
var_indep_trimmed = var_all %>%
  dplyr::select(!country) %>% #remove categorical variables
  dplyr::select(acd_u:area_ha, ends_with("_median"), t_loss.median.before:t_loss.sd.before)
cor_var_indep_trimmed = CorrSummary(var_indep_trimmed, method = "number", rect = NA,
                                    plot_name = "plot_corr_independent_var_trimmed.png",
                                    order = 'hclust', addrect = 5)
#cpc_u and cpc_d variables are highly negatively correlated, also highly collinear between different years
#cpc_u also positively correlated with higher inaccessbility
#elevation correlated with slope
#mean and sd of t_loss before project is positively correlated to project area
#normalize to per-area for project-level modelling, but not to be normalized for grid-level modelling (next-step)

#3. Select ACD of undisturbed, area, accessibility, cpc0_u, elevation, slope
#both elevation and slope are meaningful so retained despite collinearity
var_all %<>% mutate(t_loss.median.before.prop = t_loss.median.before / area_ha, t_loss.sd.before.prop = t_loss.sd.before / area_ha)

var_indep_final = var_all %>%
  dplyr::select(acd_u, area_ha, accessibility_median, cpc0_u_median,
                elevation_median, slope_median, t_loss.median.before.prop, t_loss.sd.before.prop)
cor_var_indep_final = CorrSummary(var_indep_final, method = "number", rect = NA, plot_name = "plot_corr_independent_var_final.png")


# Test collinearity of response variables ----

#1. All variables
var_res = var_all %>% dplyr::select(ends_with("after"), starts_with("drawdown"))
cor_var_res = CorrSummary(var_res, method = "number", rect = NA, plot_name = "plot_corr_response_var.png")
#all are highly collinear

#2. Use c_loss (baseline deforestation) and drawdown median and sd
#normalize to per-area for project-level modelling, but not to be normalized for grid-level modelling (next-step)
var_all %<>% mutate(c_loss.median.after.prop = c_loss.median.after / area_ha, c_loss.sd.after.prop = c_loss.sd.after / area_ha,
                    drawdown.median.prop = drawdown.median / area_ha, drawdown.sd.prop = drawdown.sd / area_ha)
var_res_final = var_all %>%
  dplyr::select(c_loss.median.after.prop, c_loss.sd.after.prop, drawdown.median.prop, drawdown.sd.prop)
cor_var_res_final = CorrSummary(var_res_final, method = "number", rect = NA, plot_name = "plot_corr_response_var_final.png")


# Run multivariate linear regression model ----
var_lm = var_all %>%
  dplyr::select(contains(colnames(var_res_final)), contains(colnames(var_indep_final)))

drawdown_lm = lm(cbind(c_loss.median.after.prop, c_loss.sd.after.prop, drawdown.median.prop, drawdown.sd.prop) ~ acd_u + ., data = var_lm)
summary(drawdown_lm)
#baseline deforestation after project
##median: inaccessibility(-), cpc_u(+), slope(- marginal), baseline deforestation before project (+)
##sd: inaccessibility(-), baseline deforestation before project (+)
#drawdown
#median: inaccessibility(-)
#sd: inaccessibility(-), cpc_u(+), baseline deforestation before project (+)
