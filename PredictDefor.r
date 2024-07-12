PredictDefor = function(proj_id, t0, acd, K, M, model_by = "M") { #mclapply() does not work on Windows
  a = Sys.time()

  #down-sample baseline to around 1/10 of the smallest set M size of the 13 projects (2559309)
  if(nrow(M) > 250000) M = M[sample(nrow(M), 250000), ]

  #obtain dependent variable of presence/absence of deforestation
  #only keep pixels who where undisturbed at t-10, because for pixels who aren't it becomes complicated to talk about their deforestation
  K = K %>%
    dplyr::select(c("lat", "lng", all_of(var_vec), "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d")) %>%
    filter(luc10 == 1) %>%
    mutate(defor = ifelse(luc0 %in% c(2, 3, 4), 1, 0))
  M = M %>%
    dplyr::select(c("lat", "lng", all_of(var_vec), "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d")) %>%
    filter(luc10 == 1) %>%
    mutate(defor = ifelse(luc0 %in% c(2, 3, 4), 1, 0))

  #perform logistic regression
  if(model_by == "M") {
    glm_out = glm(defor ~ slope + elevation + access, data = M, family = "binomial")
  } else {
    glm_out = glm(defor ~ slope + elevation + access, data = K, family = "binomial")
  }
#  summary(glm_out)
#  confint(glm_out)

  #retrieve model coefficient
  coef = summary(glm_out)$coefficients

  #predict deforestation probability using the logit function (predict() is far too slow)
  baseline = M %>%
    mutate(defor_prob = 1 / (1 + exp(-(coef[1, 1] + coef[2, 1] * slope + coef[3, 1] * elevation + coef[4, 1] * access)))) %>%
    mutate(risk = factor(ifelse(defor_prob < 0.01, "low", "high"), levels = c("low", "high"))) %>%
    mutate(acd10 = acd$carbon.density[match(luc10, acd$land.use.class)],
           acd0 = acd$carbon.density[match(luc0, acd$land.use.class)],
           c_loss = (acd10 - acd0) / 10)

  #apply same criteria on sampled project pixels
  project_defor_prob = K %>%
    mutate(defor_prob = 1 / (1 + exp(-(coef[1, 1] + coef[2, 1] * slope + coef[3, 1] * elevation + coef[4, 1] * access)))) %>%
    mutate(risk = factor(ifelse(defor_prob < 0.01, "low", "high"), levels = c("low", "high")))

  #summarise effect significance/sign and obtain ratio of pixels in baseline or project area that are predicted to be low-risk
  effects = case_when(
    coef[2:4, 4] < 0.05 & coef[2:4, 1] > 0 ~ "Pos.",
    coef[2:4, 4] < 0.05 & coef[2:4, 1] < 0 ~ "Neg.",
    coef[2:4, 4] >= 0.05 ~ "N.S.")
  names(effects) = var_vec

  #summarise ratio of pixels in baseline or project area that are predicted to be low-risk
  low_risk_ratio = c(nrow(subset(baseline, risk == "low")) / nrow(baseline),
                    nrow(subset(project_defor_prob, risk == "low")) / nrow(project_defor_prob))
  names(low_risk_ratio) = c("baseline", "project")

  #plot difference in environmental variables between low vs high-risk pixels
  p_var_by_risk = lapply(seq_along(var_vec), function(i) {
    label_position = max(baseline %>% pull(var_vec[i])) * 1.1
    label_text = paste(round(low_risk_ratio[1] * 100, 4), "% of all pixels")

    ggplot(data = baseline, aes(x = risk, y = .data[[var_vec[i]]], group = risk)) +
      geom_boxplot() +
      annotate("text", x = 1, y = label_position, label = label_text, size = 8) +
      scale_x_discrete(labels = c("Low", "High")) +
      labs(x = "Deforestation risk", y = var_label[i], title = proj_id) +
      theme_classic() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 14))
  })
  names(p_var_by_risk) = var_vec

  b = Sys.time()
  cat(proj_id, ":", b - a, "\n")
  return(list(glm_out = glm_out, plotlist = p_var_by_risk, effects = effects, low_risk_ratio = low_risk_ratio,
              baseline = baseline, project_defor_prob = project_defor_prob))
}