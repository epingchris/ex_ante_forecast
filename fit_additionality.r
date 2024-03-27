rm(list = ls())

# install.packages(c('arrow','configr', 'tidyverse', 'magrittr', 'sf', 'magrittr', 'MatchIt',
#                    'rnaturalearthdata', 'configr', 'terra', 'pbapply', 'cleangeo', 'doParallel',
#                    'foreach', 'readr', 'lwgeom', 'rnaturalearth', 'mclust', 'ggpubr', 'EnvStat'), depends = TRUE)

library(tidyverse)
library(magrittr)
library(doParallel)
library(foreach)
library(mclust)
library(ggpubr)
library(EnvStats)
library(MASS)

path = paste0('/maps/epr26/tmf_pipe_out/')

proj_list = readRDS(paste0(path, 'project_summaries.rds'))
proj_id_list = names(proj_list)

SampGMM = function(mclust_obj, n = 1000, onlyPositive = T) {
  if(is.null(mclust_obj)) return(NA)

  params = mclust_obj$param
  if(length(params$variance$sigmasq) == 1) params$variance$sigmasq[2] = params$variance$sigmasq[1]
  samp_vec = NULL
  while(length(samp_vec) < n) {
    distr_chosen = ifelse(runif(1) <= params$pro[1], 1, 2)
    val = rnorm(1, mean = params$mean[distr_chosen], sd = sqrt(params$variance$sigmasq[distr_chosen]))
    if(onlyPositive) {
      if(val > 0) samp_vec = c(samp_vec, val)
    } else {
      samp_vec = c(samp_vec, val)
    }
  }
  return(samp_vec)
}

fit_type = "normal" #GMM, normal

# Fit additionality distribution in the periods before and after project start ----
#drawdown_distr_list = lapply(seq_along(proj_list), function(i) { #used on Windows
drawdown_distr_list = mclapply(seq_along(proj_list), mc.cores = 30, function(i) {
  proj = proj_list[[i]]
  proj_name = names(proj_list[i])

  #obtain observed drawdown values
  drawdown_before = subset(proj, started == F)$additionality
  drawdown_after = subset(proj, started == T)$additionality

  #fit distributions to observed drawdown values with GMM
  fit_before = NULL
  fit_after = NULL
  if(fit_type == "GMM") {
    if(length(unique(drawdown_before)) > 1) fit_before = mclust::Mclust(drawdown_before, 2, verbose = F)
    if(length(unique(drawdown_after)) > 1) fit_after = mclust::Mclust(drawdown_after, 2, verbose = F)
  } else if(fit_type == "normal") {
    if(length(unique(drawdown_before)) > 1) fit_before = MASS::fitdistr(drawdown_before, "normal")
    if(length(unique(drawdown_after)) > 1) fit_after = MASS::fitdistr(drawdown_after, "normal")
  }

  #generate sampled drawdown values from fitted distributions
  if(fit_type == "GMM") {
    samp_before = SampGMM(fit_before, n = 1000, onlyPositive = F)
    samp_after = SampGMM(fit_after, n = 1000, onlyPositive = F)
  } else if(fit_type == "normal") {
    if(is.null(fit_before)) {
      samp_before = NA
    } else {
      samp_before = rnorm(1000, fit_before$estimate[1], fit_before$estimate[2])
    }
    if(is.null(fit_after)) {
      samp_after = NA
    } else {
      samp_after = rnorm(1000, fit_after$estimate[1], fit_after$estimate[2])
    }
  }

  #test goodness of fit
  gof_before = NULL
  gof_after = NULL
  if(sum(is.na(drawdown_before), is.na(samp_before)) == 0) gof_before = EnvStats::gofTest(drawdown_before, samp_before)
  if(sum(is.na(drawdown_after), is.na(samp_after)) == 0) gof_after = EnvStats::gofTest(drawdown_after, samp_after)

  gof_df = data.frame(project = proj_name,
                      period = c("Before", "After"),
                      stat = c(ifelse(is.null(gof_before), NA, gof_before$statistic),
                               ifelse(is.null(gof_after), NA, gof_after$statistic)),
                      pval = c(ifelse(is.null(gof_before), NA, gof_before$p.value),
                               ifelse(is.null(gof_after), NA, gof_after$p.value)))

  #plot observed and fitted density functions for both periods
  drawdown_summ = rbind(data.frame(val = drawdown_before, period = "Before", type = "Observed"),
                      data.frame(val = samp_before, period = "Before", type = "Fitted"),
                      data.frame(val = drawdown_after, period = "After", type = "Observed"),
                      data.frame(val = samp_after, period = "After", type = "Fitted")) %>%
                  mutate(project = proj_name)

  plot_fit = ggplot(data = drawdown_summ, aes(val)) +
    geom_density(aes(linetype = type, color = period)) +
    scale_linetype_manual(name = "Type", values = c(2, 1), labels = c("Fitted", "Observed")) +
    scale_color_manual(name = "Period", values = c("blue", "red"), labels = c("After", "Before")) +
    labs(title = paste("Project", proj_name),
         x = expression("Additionality (Mg CO"[2]*")")) +
    theme_bw() +
    theme(panel.grid = element_blank())

  return(list(fit_before = fit_before, fit_after = fit_after,
              gof_df = gof_df, drawdown_summ = drawdown_summ, plot_fit = plot_fit))
})

fit_before_list = lapply(drawdown_distr_list, function(x) x$fit_before) %>% `names<-`(proj_id_list)
fit_after_list = lapply(drawdown_distr_list, function(x) x$fit_after) %>% `names<-`(proj_id_list)
gof_df = lapply(drawdown_distr_list, function(x) x$gof_df) %>% do.call(rbind, .)
drawdown_df = lapply(drawdown_distr_list, function(x) x$drawdown_summ) %>% do.call(rbind, .)
plot_fit_list = lapply(drawdown_distr_list, function(x) x$plot_fit) %>% `names<-`(proj_id_list)

saveRDS(fit_before_list, file.path(paste0(path, fit_type, '/drawdown_fit_before.rds')))
saveRDS(fit_after_list, file.path(paste0(path, fit_type, '/drawdown_fit_after.rds')))
saveRDS(gof_df, file.path(paste0(path, fit_type, '/drawdown_gof.rds')))
saveRDS(drawdown_df, file.path(paste0(path, fit_type, '/drawdown_val.rds')))

if(fit_type == "normal") {
  fit_before_param = lapply(fit_before_list, function(x) x$estimate) %>% do.call(rbind, .)
  fit_after_param = lapply(fit_after_list, function(x) x$estimate) %>% do.call(rbind, .)

  saveRDS(fit_before_param, file.path(paste0(path, fit_type, '/drawdown_fit_before_param.rds')))
  saveRDS(fit_after_param, file.path(paste0(path, fit_type, '/drawdown_fit_after_param.rds')))
}

ggpubr::ggarrange(plotlist = plot_fit_list, ncol = 3, nrow = 3, common.legend = T, legend = "right") %>%
  ggpubr::ggexport(filename = paste0(path, fit_type, '/plot_drawdown_fit.pdf'), width = 18, height = 8.5, units = "in")
