#multiple regression

rm(list = ls())

library(tidyverse)
library(corrplot)
library(magrittr)
library(MASS)

path = paste0("/maps/epr26/tmf_pipe_out/")

CorrSummary = function(x, method = "number", rect = NULL, plot_name,
                       width = 640, height = 800, ...) {
  corRes = cor(x)
  testRes = cor.mtest(x, conf.level = 0.95)
  png(plot_name, width = width, height = height, type = "cairo")
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

#normalize to per-area for project-level modelling, but not to be normalized for grid-level modelling (next-step)
var_all %<>% mutate(t_loss.median.before.prop = t_loss.median.before / area_ha, t_loss.sd.before.prop = t_loss.sd.before / area_ha)

# Test collinearity of independent variables ----

#1. All variables
var_indep = var_all %>%
  dplyr::select(!country) %>% #remove categorical variables
  dplyr::select(acd_u:t_loss.sd.before, t_loss.median.before.prop:t_loss.sd.before.prop) #remove categorical variables
cor_var_indep = CorrSummary(var_indep, method = "circle", rect = c(1, seq(3, 33, by = 3), 35),
                            plot_name = "plot_corr_independent_var.png",
                            width = 800, height = 1080,
                            cl.cex = 1, tl.cex = 1)

#2. Remove min and max because collinear with median, and c_loss and loss mean because collinear with t_loss median
var_indep_trimmed = var_all %>%
  dplyr::select(!country) %>% #remove categorical variables
  dplyr::select(acd_u:area_ha, ends_with("_median"), t_loss.median.before.prop:t_loss.sd.before.prop)
cor_var_indep_trimmed = CorrSummary(var_indep_trimmed, method = "number", rect = NA,
                                    plot_name = "plot_corr_independent_var_trimmed.png",
                                    order = 'hclust', addrect = 5)
#2. Remove min and max because collinear with median, and c_loss and loss mean because collinear with t_loss median
#cpc_u and cpc_d variables are negatively correlated, also highly collinear between different years
#cpc_u also positively correlated with higher inaccessbility
#defor_5_0 is only correlated to cpc0_d/u, makes sense because it indicates more recent deforestation; it is also highly correlated to sd of t_loss
#defor_10_5 is collinear to all six CPC variables cpc10/5/0_d/u
#elevation correlated with slope

#3. Select ACD of undisturbed, area, accessibility, cpc0_d, defor_5_0, slope, median of t_loss before project

var_indep_final = var_all %>%
  dplyr::select(acd_u, area_ha, accessibility_median, cpc0_d_median, defor_5_0_median,
                slope_median, t_loss.median.before.prop)
cor_var_indep_final = CorrSummary(var_indep_final, method = "number", rect = NA, plot_name = "plot_corr_independent_var_final.png",
                                  width = 800, height = 1080,
                                  number.cex = 2, tl.cex = 2)


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
cor_var_res_final = CorrSummary(var_res_final, method = "number", rect = NA, plot_name = "plot_corr_response_var_final.png",
                       number.cex = 2, tl.cex = 2)


# Run multivariate linear regression model ----
var_lm = var_all %>%
  dplyr::select(contains(colnames(var_res_final)), contains(colnames(var_indep_final)))

drawdown_lm = lm(cbind(c_loss.median.after.prop, c_loss.sd.after.prop, drawdown.median.prop, drawdown.sd.prop) ~ acd_u + ., data = var_lm)
#baseline deforestation after project
##median: inaccessibility(-), cpc_u(+), slope(- marginal), baseline deforestation before project (+)
##sd: inaccessibility(-), baseline deforestation before project (+)
#drawdown
#median: inaccessibility(-)
#sd: inaccessibility(-), cpc_u(+), baseline deforestation before project (+)


# Select best model ----
#stepAIC does not work for multivariate linear model, so I did manual comparison using anova()
FindInsignif = function(lm_out) {
  summ = summary(lm_out)

  pval_list = lapply(summ, function(x) {
    coef = x$coefficients
    out_df = matrix(coef[, "Pr(>|t|)"], 1, nrow(coef)) %>% as.data.frame()
    colnames(out_df) = rownames(coef)
    return(out_df)
  })

  pval_df = do.call(rbind, pval_list)
  insignif_var = colnames(pval_df)[which(apply(pval_df, 2, function(x) !any(x < 0.05)))]
  return(insignif_var)
}

#starting from all the median statistics
var_indep_trimmed2 = var_indep_trimmed %>% dplyr::select(!starts_with(c("cpc10_", "cpc5_")))
var_lm = var_all %>% dplyr::select(contains(colnames(var_res_final)), contains(colnames(var_indep_trimmed2)))

drawdown_lm = lm(cbind(c_loss.median.after.prop, c_loss.sd.after.prop, drawdown.median.prop, drawdown.sd.prop) ~ acd_u + ., data = var_lm)
summary(drawdown_lm)
FindInsignif(drawdown_lm)

drawdown_lm2 = update(drawdown_lm, . ~ . - elevation_median - slope_median)
summary(drawdown_lm2)
anova(drawdown_lm2, drawdown_lm)
FindInsignif(drawdown_lm2)

drawdown_lm3 = update(drawdown_lm2, . ~ . - defor_10_5_median - defor_5_0_median)
summary(drawdown_lm3)
anova(drawdown_lm3, drawdown_lm2)
FindInsignif(drawdown_lm3)

drawdown_lm4 = update(drawdown_lm3, . ~ . - cpc0_d_median)
summary(drawdown_lm4)
anova(drawdown_lm4, drawdown_lm3)
FindInsignif(drawdown_lm4)

drawdown_lm5 = update(drawdown_lm4, . ~ . - acd_u - area_ha)
summary(drawdown_lm5)
anova(drawdown_lm5, drawdown_lm4)
FindInsignif(drawdown_lm5)
#all the insignificant variables have been removed


drawdown_lm7 = update(drawdown_lm6, . ~ . - accessibility_median)
anova(drawdown_lm7, drawdown_lm6)

drawdown_lm7 = update(drawdown_lm6, . ~ . - cpc0_u_median)
anova(drawdown_lm7, drawdown_lm6) #this one is only marginally significant

drawdown_lm7 = update(drawdown_lm6, . ~ . - t_loss.median.before.prop)
anova(drawdown_lm7, drawdown_lm6)

drawdown_lm7 = update(drawdown_lm6, . ~ . - t_loss.sd.before.prop)
anova(drawdown_lm7, drawdown_lm6)
