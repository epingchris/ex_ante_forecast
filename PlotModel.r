PlotModel = function(yr, type, model = "naive") {
  #select data to use
  yr_excl = switch(as.character(yr),
    "5" = "10",
    "10" = "5")
  figtitle = switch(as.character(yr),
    "5" = "Five-year prediction",
    "10" = "Ten-year prediction")
  if(type == "cf") {
    model_df_selected = model_df_scaled %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & starts_with(c("forecast", "closs_obs_cf"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = closs_obs_cf) %>%
      mutate(forecast = forecast * 100, observed = observed * 100) #turn into percentage
    x_lab = paste("Predicted annual counterfactual carbon loss (%)")
    model = "naive"
    max_val = 6
    break_val = 0:6
    text_x = 4.5
    text_y = c(2, 1.5, 1)
  } else if(type == "p") {
    model_df_selected = model_df_scaled %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & !starts_with(c("project", "closs_obs_cf", "add_rate", "add_obs"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = closs_obs_p) %>%
      mutate(forecast = forecast * 100, observed = observed * 100) #turn into percentage
    x_lab = paste("Predicted annual project carbon loss (%)")
    max_val = 6
    break_val = 0:6
    text_x = 4.5
    text_y = c(2, 1.5, 1)
  } else if(type == "add") {
    model_df_selected = model_df_scaled %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & !starts_with(c("project", "closs_obs", "add_rate"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = add_obs) %>%
      mutate(forecast = forecast * 100) #turn into percentage only for carbon loss forecast
    x_lab = bquote(paste("Predicted annual carbon credit production (MgC ", ha^-1, " ", yr^-1, ")"))
    max_val = 1.75
    break_val = seq(0, 1.75, 0.25)
    text_x = 1.25
    text_y = c(0.5, 0.375, 0.25)
  } else if(type == "add_rate") {
    model_df_selected = model_df_scaled %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & !starts_with(c("project", "closs_obs_", "add_obs"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = add_rate) %>%
      mutate(forecast = forecast * 100, observed = observed * 100) #turn into percentage
    x_lab = bquote(paste("Predicted difference in carbon loss rate (%)"))
    max_val = 6
    break_val = 0:6
    text_x = 4.5
    text_y = c(2, 1.5, 1)
  }

  #run linear model
  if(model == "full") {
    forecast_lm = lm(observed ~ ., data = model_df_selected) #full model using best forecast and socio-environmental variables
  } else if(model == "naive") {
    forecast_lm = lm(observed ~ forecast, data = model_df_selected) #naive model using only best forecast
  } else if(model == "sel") {
    forecast_lm_full = lm(observed ~ ., data = model_df_selected) #full model predicting project counterfactual carbon loss
    forecast_lm = stepAIC(forecast_lm_full,
                          scope = list(upper = ~ ., lower = ~ forecast),
                          direction = "backward", trace = 1) #backward selection
  }

  #print diagnostic plots
  par(mfrow = c(2, 2))
  png(paste0(fig_path, "figure_diagnostic_", yr, "_", type, "_", model, ".png"), width = 600, height = 600)
  plot(forecast_lm)
  dev.off()

  #Calculate predictions and predictive performance
  pred_df = data.frame(pred = predict(forecast_lm),
                       observed = forecast_lm$model$observed)
  R2 = GOF(pred_df$pred, pred_df$observed) #goodness-of-fit (R2 over 1:1 line)
  mape = MAPE(pred_df$pred, pred_df$observed) #mean absolute percentage error (MAPE)
  mpb = MPB(pred_df$pred, pred_df$observed) #mean percentage bias (MPB)

  #Plot model (observed vs predicted project carbon loss)
  plot_model = ggplot(data = pred_df) +
    geom_point(aes(x = pred, y = observed), size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    annotate(geom = "text", x = text_x, y = text_y[1], size = 10,
             label = bquote(paste("MAPE: ", .(round(mape, 3)), "%"))) +
    annotate(geom = "text", x = text_x, y = text_y[2], size = 10,
             label = bquote(paste("MPB: ", .(round(mpb, 3)), "%"))) +
    annotate(geom = "text", x = text_x, y = text_y[3], size = 10,
             label = bquote(paste("Goodness-of-fit: ", .(round(R2, 3))))) +
    labs(title = figtitle,
         x = x_lab,
         y = "") +
    scale_x_continuous(limits = c(0, max_val), breaks = break_val) +
    scale_y_continuous(limits = c(0, max_val), breaks = break_val) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "cm"),
          plot.title = element_text(size = 32, hjust = 0.5),
          axis.title.x = element_text(size = 28),
          axis.text.x = element_text(size = 24),
          axis.title.y = element_text(size = 28),
          axis.text.y = element_text(size = 24),
          axis.ticks = element_blank(),
          axis.line = element_line(color = "black"))
  
  return(list(model = forecast_lm, pred = pred_df, R2 = R2, mape = mape, mpb = mpb, plot = plot_model))
}