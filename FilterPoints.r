FilterPoints = function(analysis_type, proj_id) { #mclapply() does not work on Windows
  if(analysis_type == "old_source") {
    cat("Use new results from epr26 instead.\n")
    return(NULL)
  }

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

  k = read_parquet(paste0(file_prefix, "k.parquet")) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0)) %>%
    dplyr::select(c("elevation", "slope", "access", luc10, luc0)) %>%
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
  effect_labels = ifelse(pval < 0.05,
                         ifelse(coefficients > 0, "Pos.", "Neg."),
                         "N.S.") %>%
    matrix(., nrow = 1, ncol = 3) %>%
    as.data.frame()
  colnames(effect_labels) = c("slope", "elevation", "access")

  out = lapply(c("slope", "elevation", "access"), function(x) {
    x_label = switch(x,
                    "slope" = "Slope",
                    "elevation" = "Elevation (m)",
                    "access" = "Remoteness (minutes)")
    exclude_ratio = NA
    thres = NA

    p = ggplot(data = k, aes_string(x = x, y = "defor")) +
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

      p = p + geom_line(data = k_pred, aes_string(x = x, y = "defor"))
      if(!is.na(thres)) p = p + geom_vline(xintercept = thres, color = "red")
    }

    return(list(exclude_ratio = exclude_ratio, thres = thres, p = p))
  })

  exclude_ratio = sapply(out, function(x) x$exclude_ratio)
  thres = sapply(out, function(x) x$thres)
  plotlist = lapply(out, function(x) x$p)
  names(plotlist) = c("slope", "elevation", "access")

  return(list(logit_k = logit_k, coefficients = coefficients, pval = pval, effect_labels = effect_labels,
              exclude_ratio = exclude_ratio, thres = thres, plotlist = plotlist))
}

#1 / (1 + exp(-(logit_k_coef[1, 1] + logit_k_coef[2, 1] * 242.6669 + logit_k_coef[3, 1] * 6.0885 + logit_k_coef[4, 1] * 22)))
