PredictDefor = function(var_type, coef, src) {
  df_pred = data.frame(slope = rep(mean(src$slope), 100),
                       elevation = rep(mean(src$elevation), 100),
                       access = rep(mean(src$access), 100))
  val = pull(src, all_of(var_type))
  df_pred[, var_type] = seq(min(val), max(val), len = 100)
  df_pred = df_pred %>%
    mutate(defor = 1 / (1 + exp(-(coef[1, 1] + coef[2, 1] * slope + coef[3, 1] * elevation + coef[4, 1] * access)))) #predict() is far too slow
  # roc(k$defor, k_pred2)

  effect = "N.S."
  exclude_ratio = NA
  thres = NA

  if(coef[var_type, 4] < 0.05) { #if effect is significant
    exclude_ratio = sum(df_pred$defor < 0.01) / 100 #how many pixels have deforestation probability < 1%
    effect = ifelse(coef[var_type, 1] > 0, "Pos.", "Neg.") #determine the effect labels based on coefficient signs

    if(exclude_ratio > 0 & exclude_ratio < 1) { #if proportion of pixels with <1% deforestation prob. is between 1 and 0
      if(effect == "Neg.") {
        thres = df_pred[head(which(df_pred$defor < 0.01), 1), var_type] #set maximum threshold if effect is negative
      } else {
        thres = df_pred[tail(which(df_pred$defor < 0.01), 1), var_type] #set minimum threshold if effect is positive
      }
    }
  }

  return(list(df_pred = df_pred, effect = effect, exclude_ratio = exclude_ratio, thres = thres))
}

ApplyThreshold = function(src, thres, effect_labels, variable) {
  if(is.na(thres[, variable])) {
    exclude = rep(F, nrow(src))
  } else {
    exclude = switch(effect_labels[, variable],
                           "Neg." = (pull(src, all_of(variable)) >= thres[, variable]),
                           "Pos." = (pull(src, all_of(variable)) <= thres[, variable]),
                           "N.S." = rep(F, nrow(src)))
  }
  return(exclude)
}

FilterVicinity = function(analysis_type, proj_id) { #mclapply() does not work on Windows
  if(analysis_type == "old_source") {
    cat("Use new results from epr26 instead.\n")
    return(NULL)
  }
  a = Sys.time()

  var_vec = c("slope", "elevation", "access")

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

  #sub-sample project vicinity down to around 1/10 of the smallest vicinity size of the 15 projects (2559309)
  K = read_parquet(paste0(file_prefix, "k.parquet"))
  M = read_parquet(paste0(file_prefix, "matches.parquet"))
  if(nrow(M) > 250000) M = M[sample(nrow(M), 250000), ]

  K = K %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0)) %>%
    dplyr::select(c("lat", "lng", "elevation", "slope", "access", "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d")) %>%
    filter(luc10 == 1) %>%
    mutate(defor = ifelse(luc0 %in% c(2, 3, 4), 1, 0))

  M = M %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0)) %>%
    dplyr::select(c("lat", "lng", "elevation", "slope", "access", "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d")) %>%
    filter(luc10 == 1) %>%
    mutate(defor = ifelse(luc0 %in% c(2, 3, 4), 1, 0))

  #perform logistic regression on pixels who where undisturbed at t-10
  glm_k = glm(defor ~ slope + elevation + access, data = K, family = "binomial")
  glm_m = glm(defor ~ slope + elevation + access, data = M, family = "binomial")

  #predict() is far too slow: use logit function to calculate predicted value directly
  coef_k = summary(glm_k)$coefficients
  coef_m = summary(glm_m)$coefficients
  df_pred_k = M %>%
    dplyr::select(all_of(var_vec)) %>%
    mutate(defor = 1 / (1 + exp(-(coef_k[1, 1] + coef_k[2, 1] * slope + coef_k[3, 1] * elevation + coef_k[4, 1] * access)))) %>%
    mutate(risk = ifelse(defor < 0.01, "low", "high")) %>%
    mutate(risk = factor(risk, levels = c("low", "high")))
  df_pred_m = M %>%
    dplyr::select(all_of(var_vec)) %>%
    mutate(defor = 1 / (1 + exp(-(coef_m[1, 1] + coef_m[2, 1] * slope + coef_m[3, 1] * elevation + coef_m[4, 1] * access)))) %>%
    mutate(risk = ifelse(defor < 0.01, "low", "high")) %>%
    mutate(risk = factor(risk, levels = c("low", "high")))

  var_label = c("Slope (dgree)", "Elevation (meter)", "Remoteness (minutes)")
  p_k = lapply(seq_along(var_vec), function(i) {
    ggplot(data = df_pred_k, aes(x = risk, y = .data[[var_vec[i]]], group = risk)) +
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
    ggplot(data = df_pred_m, aes(x = risk, y = .data[[var_vec[i]]], group = risk)) +
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

plot_k_m = cowplot::plot_grid(plotlist = c(p_k, p_m), byrow = F, nrow = 3, ncol = 2)

#  pred_k = lapply(var_vec, function(x) PredictDefor(x, coef = summary(glm_k)$coefficients, src = k))
#  pred_matches = lapply(var_vec, function(x) PredictDefor(x, coef = summary(glm_matches)$coefficients, src = matches))
#  summary(glm_k)
#  confint(glm_k)

  summ_filtering = lapply(seq_along(var_vec), function(i) {
    pred_out_k = pred_k[[i]]
    pred_out_matches = pred_matches[[i]]
    data.frame(project = proj_id,
               var = var_vec[i],
               k_effect = pred_out_k$effect,
               k_exclude_ratio = pred_out_k$exclude_ratio,
               k_thres = pred_out_k$thres,
               m_effect = pred_out_matches$effect,
               m_exclude_ratio = pred_out_matches$exclude_ratio,
               m_thres = pred_out_matches$thres)
  }) %>% do.call(rbind, .)

  plotlist = lapply(seq_along(var_vec), function(i) {
    thres_k = summ_filtering$k_thres[i]
    thres_matches = summ_filtering$m_thres[i]

    p_k = ggplot(data = pred_k[[i]]$df_pref, aes(x = .data[[x]], y = defor)) +
        geom_line() +
        geom_hline(yintercept = 0.01, linetype = 2) +
        labs(title = proj_id, x = x_label[i], y = "Deforestation probability") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))
    if(!is.na(thres_k)) p_k = p_k + geom_vline(xintercept = thres, color = "red")

    p_m = ggplot(data = pred_matches[[i]]$df_pref, aes(x = .data[[x]], y = defor)) +
        geom_line() +
        geom_hline(yintercept = 0.01, linetype = 2) +
        labs(title = proj_id, x = x_label[i], y = "Deforestation probability") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))
    if(!is.na(thres_k)) p_m = p_m + geom_vline(xintercept = thres, color = "red")

    return(list(p_k = p_k, p_m = p_m))
  })

  #filter matches to be baseline using model fitted with matches
  slope_exclude = ApplyThreshold(matches, thres = thres, effect_labels = effect_labels, variable = "slope")
  elevation_exclude = ApplyThreshold(matches, thres = thres, effect_labels = effect_labels, variable = "elevation")
  access_exclude = ApplyThreshold(matches, thres = thres, effect_labels = effect_labels, variable = "access")

  #filter matches to be baseline using model fitted with matches
  slope_exclude = ApplyThreshold(matches, thres = thres, effect_labels = effect_labels, variable = "slope")
  elevation_exclude = ApplyThreshold(matches, thres = thres, effect_labels = effect_labels, variable = "elevation")
  access_exclude = ApplyThreshold(matches, thres = thres, effect_labels = effect_labels, variable = "access")


  vicinity = matches %>%
    mutate(exclude = (slope_exclude | elevation_exclude | access_exclude),
           acd_t_10 = acd$carbon.density[match(luc10, acd$land.use.class)],
           acd_t0 = acd$carbon.density[match(luc0, acd$land.use.class)],
           c_loss = (acd_t_10 - acd_t0) / 10)

  vicinity_area = nrow(vicinity) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares
  vicinity_area_filtered = nrow(subset(vicinity, !exclude)) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares

  k_slope_exclude = ApplyThreshold(k, thres = thres, effect_labels = effect_labels, variable = "slope")
  k_elevation_exclude = ApplyThreshold(k, thres = thres, effect_labels = effect_labels, variable = "elevation")
  k_access_exclude = ApplyThreshold(k, thres = thres, effect_labels = effect_labels, variable = "access")

  k = k %>% mutate(exclude = (k_slope_exclude | k_elevation_exclude | k_access_exclude))
  k_area_ratio = nrow(subset(k, !exclude)) / nrow(k)

  b = Sys.time()
  cat(proj_id, ":", b - a, "\n")
  return(list(glm_logit = glm_logit, glm_k = glm_k, glm_mtaches = glm_mtaches, summ_filtering, plotlist = plotlist,
              vicinity_area = vicinity_area, vicinity_area_filtered = vicinity_area_filtered, vicinity = vicinity, k_area_ratio = k_area_ratio))
}