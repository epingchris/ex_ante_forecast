FilterVicinity = function(analysis_type, proj_id, by_source = "k") { #mclapply() does not work on Windows
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

  k = read_parquet(paste0(file_prefix, by_source, ".parquet")) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0)) %>%
    dplyr::select(c("elevation", "slope", "access", "luc10", "luc0")) %>%
    mutate(lucc = paste(luc10, luc0, sep = "_"),
           defor = ifelse(lucc %in% c("1_2", "1_3", "1_4"), 1, 0))

# %>% filter(.data[[luc_t_10]] == 1)

  #perform logistic regression on pixels who where undisturbed at t-10
  logit_k = glm(defor ~ slope + elevation + access, data = k %>% filter(luc10 == 1), family = "binomial")
#  summary(logit_k)
#  confint(logit_k)
  logit_k_coef = summary(logit_k)$coefficients

  #determine the effect labels based on p-values and coefficient signs
  coefficients = logit_k_coef[2:4, 1]
  pval = logit_k_coef[2:4, 4]
  effect_labels = case_when(
    pval < 0.05 & coefficients > 0 ~ "Pos.",
    pval < 0.05 & coefficients < 0 ~ "Neg.",
    pval >= 0.05 ~ "N.S."
  ) %>%
    matrix(., nrow = 1, ncol = 3) %>%
    as.data.frame()
  colnames(effect_labels) = var_vec

  out = lapply(var_vec, function(x) {
    x_label = switch(x,
                    "slope" = "Slope",
                    "elevation" = "Elevation (m)",
                    "access" = "Remoteness (minutes)")
    exclude_ratio = NA
    thres = NA

    p = ggplot(data = k, aes(x = .data[[x]], y = defor)) +
        geom_point(shape = "bullet", size = 1, color = "darkgray") +
        geom_hline(yintercept = 0.01, linetype = 2) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(title = proj_id, x = x_label, y = "Deforestation probability") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))

    if(pval[x] < 0.05) {
      k_pred = data.frame(slope = rep(mean(k$slope), 10000),
                          elevation = rep(mean(k$elevation), 10000),
                          access = rep(mean(k$access), 10000))
      k_pred[, x] = seq(min(k[, x]), max(k[, x]), len = 10000)
      k_pred$defor = predict(logit_k, newdata = k_pred, type = "response")

      # k_pred2 = predict(logit_k, newdata = k[, c("defor", "slope", "elevation", "access")],
      #                  type = "response")
      # roc(k$defor, k_pred2)

      exclude_ratio = sum(k_pred$defor < 0.01) / nrow(k_pred) #how many pixels have deforestation probability < 1%
      if(exclude_ratio > 0 & exclude_ratio < 1) {
        thres = switch(effect_labels[, x],
                        "Neg." = k_pred[head(which(k_pred$defor < 0.01), 1), x], #set maximum threshold if effect is negative
                        "Pos." = k_pred[tail(which(k_pred$defor < 0.01), 1), x]) #set minimum threshold if effect is positive
      }

      p = p + geom_line(data = k_pred, aes(x = .data[[x]], y = defor))
      if(!is.na(thres)) p = p + geom_vline(xintercept = thres, color = "red")
    }

    return(list(exclude_ratio = exclude_ratio, thres = thres, p = p))
  })

  exclude_ratio = sapply(out, function(x) x$exclude_ratio)
  thres = sapply(out, function(x) x$thres) %>%
    matrix(., nrow = 1, ncol = 3) %>%
    as.data.frame()
  colnames(thres) = var_vec
  plotlist = lapply(out, function(x) x$p)
  names(plotlist) = var_vec

  #filter matches to be vicinity
  matches = read_parquet(paste0(file_prefix, "matches.parquet")) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0))

  if(is.na(thres$slope)) {
    slope_exclude = rep(F, nrow(matches))
  } else {
    slope_exclude = switch(effect_labels$slope,
                           "Neg." = (matches$slope >= thres$slope),
                           "Pos." = (matches$slope <= thres$slope),
                           "N.S." = rep(F, nrow(matches)))
  }

  if(is.na(thres$elevation)) {
    elevation_exclude = rep(F, nrow(matches))
  } else {
    elevation_exclude = switch(effect_labels$elevation,
                               "Neg." = (matches$elevation >= thres$elevation),
                               "Pos." = (matches$elevation <= thres$elevation),
                               "N.S." = rep(F, nrow(matches)))
  }

  if(is.na(thres$access)) {
    access_exclude = rep(F, nrow(matches))
  } else {
    access_exclude = switch(effect_labels$access,
                            "Neg." = (matches$access >= thres$access),
                            "Pos." = (matches$access <= thres$access),
                            "N.S." = rep(F, nrow(matches)))
  }

  vicinity = matches %>%
    dplyr::select(c(var_vec, "luc10", "luc0")) %>%
    mutate(exclude = (slope_exclude | elevation_exclude | access_exclude),
           lucc = paste(luc10, luc0, sep = "_"),
           defor = ifelse(lucc %in% c("1_2", "1_3", "1_4"), 1, 0),
           acd_t_10 = acd$carbon.density[match(luc10, acd$land.use.class)],
           acd_t0 = acd$carbon.density[match(luc0, acd$land.use.class)],
           c_loss = (acd_t_10 - acd_t0) / 10)

  vicinity_area = nrow(vicinity) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares
  vicinity_area_filtered = nrow(subset(vicinity, !exclude)) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares

  #sub-sample project vicinity down to around 1/10 of the smallest vicinity size of the 15 projects (2559309)
  if(nrow(vicinity) > 250000) vicinity = vicinity[sample(nrow(vicinity), 250000), ]

  b = Sys.time()
  b - a
  return(list(k = k, logit_k = logit_k, coefficients = coefficients, pval = pval, effect_labels = effect_labels,
              exclude_ratio = exclude_ratio, thres = thres, plotlist = plotlist,
              vicinity_area = vicinity_area, vicinity_area_filtered = vicinity_area_filtered, vicinity = vicinity))
}

#1 / (1 + exp(-(logit_k_coef[1, 1] + logit_k_coef[2, 1] * 242.6669 + logit_k_coef[3, 1] * 6.0885 + logit_k_coef[4, 1] * 22)))
