FilterBaseline = function(analysis_type, proj_id) { #mclapply() does not work on Windows
  if(analysis_type == "old_source") {
    cat("Use new results from epr26 instead.\n")
    return(NULL)
  }
  a = Sys.time()

  t0 = switch(analysis_type,
              "full" = filter(proj_meta, ID == str_replace(proj_id, "a", ""))$t0,
              "grid" = filter(proj_meta, ID == "1201")$t0,
              "control" = 2011,
              "ac" = 2021)

  #column names to select from vicinity
  luc_t_10 = paste0("luc_", t0 - 10)
  luc_t0 = paste0("luc_", t0)

  #get biome variables and project vicinity
  file_prefix = paste0(project_dir, proj_id, "/",
                       switch(analysis_type,
                              "full" = "",
                              "grid" = "1201_",
                              "control" = "0000_",
                              "ac" = ""),
                       proj_id)

  #get ACD and ACD of undisturbed forest
  acd = read.csv(paste0(file_prefix, "carbon-density.csv"))
  for(class in 1:6) {
    if(class %in% acd$land.use.class == F) acd = rbind(acd, c(class, NA))
  }

  #load set K (sampled project pixels) and set M (potentially matched pixels)
  K = read_parquet(paste0(file_prefix, "k.parquet"))
  M = read_parquet(paste0(file_prefix, "matches.parquet"))

  #down-sample baseline to around 1/10 of the smallest set M size of the 13 projects (2559309)
  if(nrow(M) > 250000) M = M[sample(nrow(M), 250000), ]

  #obtain dependent variable of presence/absence of deforestation
  #only keep pixels who where undisturbed at t-10, because for pixels who aren't it becomes complicated to talk about their deforestation
  K = K %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0)) %>%
    dplyr::select(c("lat", "lng", all_of(var_vec), "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d")) %>%
    filter(luc10 == 1) %>%
    mutate(defor = ifelse(luc0 %in% c(2, 3, 4), 1, 0))
  M = M %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0)) %>%
    dplyr::select(c("lat", "lng", all_of(var_vec), "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d")) %>%
    filter(luc10 == 1) %>%
    mutate(defor = ifelse(luc0 %in% c(2, 3, 4), 1, 0))

  #perform logistic regression
  glm_k = glm(defor ~ slope + elevation + access, data = K, family = "binomial")
  glm_m = glm(defor ~ slope + elevation + access, data = M, family = "binomial")
#  summary(glm_k)
#  confint(glm_k)

  #retrieve model coefficient
  coef_k = summary(glm_k)$coefficients
  coef_m = summary(glm_m)$coefficients

  #predict deforestation probability using the logit function (predict() is far too slow)
  baseline = M %>%
    mutate(defor_prob_k = 1 / (1 + exp(-(coef_k[1, 1] + coef_k[2, 1] * slope + coef_k[3, 1] * elevation + coef_k[4, 1] * access))),
           defor_prob_m = 1 / (1 + exp(-(coef_m[1, 1] + coef_m[2, 1] * slope + coef_m[3, 1] * elevation + coef_m[4, 1] * access)))) %>%
    mutate(risk_k = ifelse(defor_prob_k < 0.01, "low", "high"),
           risk_m = ifelse(defor_prob_m < 0.01, "low", "high")) %>%
    mutate(risk_k = factor(risk_k, levels = c("low", "high")),
           risk_m = factor(risk_m, levels = c("low", "high")),
           acd10 = acd$carbon.density[match(luc10, acd$land.use.class)],
           acd0 = acd$carbon.density[match(luc0, acd$land.use.class)],
           c_loss = (acd10 - acd0) / 10)

  #apply same criteria on sampled project pixels
  project_filtered = K %>%
    mutate(defor_prob_k = 1 / (1 + exp(-(coef_k[1, 1] + coef_k[2, 1] * slope + coef_k[3, 1] * elevation + coef_k[4, 1] * access))),
           defor_prob_m = 1 / (1 + exp(-(coef_m[1, 1] + coef_m[2, 1] * slope + coef_m[3, 1] * elevation + coef_m[4, 1] * access)))) %>%
    mutate(risk_k = ifelse(defor_prob_k < 0.01, "low", "high"),
           risk_m = ifelse(defor_prob_m < 0.01, "low", "high")) %>%
    mutate(risk_k = factor(risk_k, levels = c("low", "high")),
           risk_m = factor(risk_m, levels = c("low", "high")))

  #plot difference in environmental variables between low vs high-risk pixels
  p_k = lapply(seq_along(var_vec), function(i) {
    ggplot(data = baseline, aes(x = risk_k, y = .data[[var_vec[i]]], group = risk_k)) +
      geom_boxplot() +
      scale_x_discrete(labels = c("Low", "High")) +
      labs(x = "Pixel deforestation probability", y = var_label[i], title = "Predicted by project pixels") +
      theme_classic() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 14))
  })
  p_m = lapply(seq_along(var_vec), function(i) {
    ggplot(data = baseline, aes(x = risk_m, y = .data[[var_vec[i]]], group = risk_m)) +
      geom_boxplot() +
      scale_x_discrete(labels = c("Low", "High")) +
      labs(x = "Pixel deforestation probability", y = var_label[i], title = "Predicted by baseline pixels") +
      theme_classic() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 14))
  })

  #summarise effect significance/sign and obtain ratio of pixels in baseline or project area that are predicted to be low-risk
  effects_k = case_when(
    coef_k[2:4, 4] < 0.05 & coef_k[2:4, 1] > 0 ~ "Pos.",
    coef_k[2:4, 4] < 0.05 & coef_k[2:4, 1] < 0 ~ "Neg.",
    coef_k[2:4, 4] >= 0.05 ~ "N.S.")
  effects_m = case_when(
    coef_m[2:4, 4] < 0.05 & coef_m[2:4, 1] > 0 ~ "Pos.",
    coef_m[2:4, 4] < 0.05 & coef_m[2:4, 1] < 0 ~ "Neg.",
    coef_m[2:4, 4] >= 0.05 ~ "N.S.")
  effects = data.frame(by_k = effects_k, by_m = effects_m)
  rownames(effects) = var_vec

  #summarise ratio of pixels in baseline or project area that are predicted to be low-risk
  low_risk_perc = data.frame(
    by_k = c(nrow(subset(baseline, risk_k == "low")) / nrow(baseline),
             nrow(subset(project_filtered, risk_k == "low")) / nrow(project_filtered)),
    by_m = c(nrow(subset(baseline, risk_m == "low")) / nrow(baseline),
             nrow(subset(project_filtered, risk_m == "low")) / nrow(project_filtered)))
  rownames(low_risk_perc) = c("baseline_low_risk_perc", "project_low_risk_perc")
  df_summary = rbind(effects, round(low_risk_perc * 100, 2))

  b = Sys.time()
  cat(proj_id, ":", b - a, "\n")
  return(list(glm_k = glm_k, glm_m = glm_m, plotlist = c(p_k, p_m),
              baseline = baseline, project_filtered = project_filtered, df_summary = df_summary))
}