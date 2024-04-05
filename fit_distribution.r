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


# Custom functions ----
SampGMM = function(mclust_obj, n = 1000, onlyPositive = T) {
  if(is.null(mclust_obj)) return(rep(NA, n))

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

SampNorm = function(fitdistr_obj, n = 1000) {
  if(is.null(fitdistr_obj)) return(rep(NA, n))
  samp = rnorm(n, fitdistr_obj$estimate[1], fitdistr_obj$estimate[2])
  return(samp)
}

gofWrapper = function(x, y){
  gof = NULL
  if(!all(is.na(x)) & !all(is.na(y))) gof = EnvStats::gofTest(x, y, na.action = na.omit)
  return(gof)
}

plot_drawdown = F
samp_n = 1000

# Fit additionality distribution in the periods before and after project start ----
#drawdown_distr_list = lapply(seq_along(proj_list), function(i) { #used on Windows
drawdown_distr_list = mclapply(seq_along(proj_list), mc.cores = 30, function(i) {
  proj = proj_list[[i]]
  proj_name = proj_id_list[i]

  #obtain observed values
  proj_before = proj %>% filter(started == F)
  c_loss_before = proj_before %>% pull(c_loss)
  t_loss_before = proj_before %>% pull(t_loss)
  drawdown_before = proj_before %>% pull(additionality)

  proj_after = proj %>% filter(started == T)
  c_loss_after = proj_after %>% pull(c_loss)
  t_loss_after = proj_after %>% pull(t_loss)
  drawdown_after = proj_after %>% pull(additionality)

  #fit observed carbon flux values with GMM, fit observed drawdown values to normal distribution
  fit_before = list(c_loss = NULL, t_loss = NULL, drawdown = NULL)
  if(length(unique(c_loss_before)) > 1) fit_before$c_loss = mclust::Mclust(c_loss_before, 2, verbose = F)
  if(length(unique(t_loss_before)) > 1) fit_before$t_loss = mclust::Mclust(t_loss_before, 2, verbose = F)
  if(length(unique(drawdown_before)) > 1) fit_before$drawdown = MASS::fitdistr(drawdown_before, "normal")

  fit_after = list(c_loss = NULL, t_loss = NULL, drawdown = NULL)
  if(length(unique(c_loss_after)) > 1) fit_after$c_loss = mclust::Mclust(c_loss_after, 2, verbose = F)
  if(length(unique(t_loss_after)) > 1) fit_after$t_loss = mclust::Mclust(t_loss_after, 2, verbose = F)
  if(length(unique(drawdown_after)) > 1) fit_after$drawdown = MASS::fitdistr(drawdown_after, "normal")

  #generate sampled drawdown values from fitted distributions
  samp_before = list(c_loss = SampGMM(fit_before$c_loss, n = samp_n, onlyPositive = F),
                     t_loss = SampGMM(fit_before$t_loss, n = samp_n, onlyPositive = F),
                     drawdown = SampNorm(fit_before$drawdown, n = samp_n))
  samp_after = list(c_loss = SampGMM(fit_after$c_loss, n = samp_n, onlyPositive = F),
                    t_loss = SampGMM(fit_after$t_loss, n = samp_n, onlyPositive = F),
                    drawdown = SampNorm(fit_after$drawdown, n = samp_n))

  #test goodness of fit
  gof_before = list(c_loss = gofWrapper(c_loss_before, samp_before$c_loss),
                    t_loss = gofWrapper(t_loss_before, samp_before$t_loss),
                    drawdown = gofWrapper(drawdown_before, samp_before$drawdown))
  gof_after = list(c_loss = gofWrapper(c_loss_after, samp_after$c_loss),
                    t_loss = gofWrapper(t_loss_after, samp_after$t_loss),
                    drawdown = gofWrapper(drawdown_after, samp_after$drawdown))

  gof_df = data.frame(project = proj_name,
                      period = rep(c("Before", "After"), each = 3),
                      var = rep(c("c_loss", "t_loss", "drawdown"), 2),
                      stat = c(sapply(gof_before, function(x) ifelse(is.null(x), NA, x$statistic)),
                               sapply(gof_after, function(x) ifelse(is.null(x), NA, x$statistic))),
                      pval = c(sapply(gof_before, function(x) ifelse(is.null(x), NA, x$p.value)),
                               sapply(gof_after, function(x) ifelse(is.null(x), NA, x$p.value))))

  #plot observed and fitted density functions for both periods
  drawdown_summ = NULL
  plot_fit = NULL
  if(plot_drawdown) {
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
  }

  return(list(fit_before = fit_before, fit_after = fit_after,
              samp_before = samp_before, samp_after = samp_after,
              gof_df = gof_df, drawdown_summ = drawdown_summ, plot_fit = plot_fit))
})

drawdown_distr_list %<>% `names<-`(proj_id_list)
saveRDS(drawdown_distr_list, file.path(paste0(path, 'drawdown_distr_list.rds')))

#project-level parameters:
#mean and sd of normal distribution used to fit additionality before/after project start
fit_param = lapply(seq_along(drawdown_distr_list), function(i) {
  proj_name = proj_id_list[i]
  x_before = drawdown_distr_list[[i]]$samp_before
  x_after = drawdown_distr_list[[i]]$samp_after

  df_before = lapply(seq_along(x_before), function(i) {
    x = x_before[[i]]
    nom = names(x_before[i])
    df = data.frame(mean(x), median(x), sd(x))
    colnames(df) = paste0(nom, ".", c("mean", "median", "sd"))
    return(df)
  }) %>%
    do.call(cbind, .) %>%
    mutate(project = proj_name)

  df_after = lapply(seq_along(x_after), function(i) {
    x = x_after[[i]]
    nom = names(x_after[i])
    df = data.frame(mean(x), median(x), sd(x))
    colnames(df) = paste0(nom, ".", c("mean", "median", "sd"))
    return(df)
  }) %>%
    do.call(cbind, .) %>%
    mutate(project = proj_name)

  return(list(before = df_before, after = df_after))
})

fit_param_before = lapply(fit_param, function(x) x$before) %>% do.call(rbind, .)
fit_param_after = lapply(fit_param, function(x) x$after) %>% do.call(rbind, .)
saveRDS(list(before = fit_param_before, after = fit_param_after), file.path(paste0(path, "fit_param.rds")))

if(plot_drawdown) {
  drawdown_df = lapply(drawdown_distr_list, function(x) x$drawdown_summ) %>% do.call(rbind, .)
  saveRDS(drawdown_df, file.path(paste0(path, fit_type, '/drawdown_val.rds')))

  plot_fit_list = lapply(drawdown_distr_list, function(x) x$plot_fit) %>% `names<-`(proj_id_list)
  ggpubr::ggarrange(plotlist = plot_fit_list, ncol = 3, nrow = 3, common.legend = T, legend = "right") %>%
    ggpubr::ggexport(filename = paste0(path, fit_type, '/plot_drawdown_fit.pdf'), width = 18, height = 8.5, units = "in")
}