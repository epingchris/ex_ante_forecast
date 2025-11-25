PlotModel = function(input, yr, type, model) {
  #select data to use
  yr_excl = switch(as.character(yr),
    "5" = "10",
    "10" = "5")
  figtitle = switch(as.character(yr),
    "5" = "Five-year prediction",
    "10" = "Ten-year prediction")

  lim_val = c(-1, 5)
  break_val = seq(-1, 5)
  text_x = -0.75
  text_y = c(4, 4.5, 5)
  text_size = 12
  x_lab = "Predicted rate (%)"
  y_lab = "Observed rate (%)"

  if(type == "cf") {
    model_df_selected = input %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & starts_with(c("forecast", "closs_obs_cf"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = closs_obs_cf) %>%
      mutate(observed = observed * 100) #turn into percentage
    model = "naive"
    title_type = "A. Counterfactual\ncarbon loss rate"
    point_col = "darkred"
  } else if(type == "p") {
    model_df_selected = input %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & !starts_with(c("project", "closs_obs_cf", "add_", "credit"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = closs_obs_p) %>%
      mutate(observed = observed * 100) #scale forecasts and turn observed into percentage
    title_type = "B. Project\ncarbon loss rate"
    point_col = "blue"
  } else if(type == "reduc") {
    model_df_selected = input %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & !starts_with(c("project", "closs_obs_", "add_obs", "credit"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = add_rate) %>%
      mutate(observed = observed * 100) #turn into percentage
    title_type = "C. Reduction of\ncarbon loss rate"
    point_col = "darkgreen"
  } else if(type == "credit") {
    model_df_selected = input %>%
      dplyr::select(!ends_with(as.character(yr_excl)) & !starts_with(c("project", "closs_obs_", "add_"))) %>%
      rename_with(~ gsub("_[0-9]+$", "", .x)) %>%
      rename(observed = credit)
    title_type = "\nD. Carbon credit production"
    lim_val = c(-2, 8)
    break_val = seq(-2, 8, 2)
    text_x = -1.5
    text_y = c(6, 7, 8)
    x_lab = expression(paste("Predicted credits (MgC", " ", ha^-1, " ", yr^-1, ")", sep = " "))
    y_lab = expression(paste("Observed credits (MgC", " ", ha^-1, " ", yr^-1, ")", sep = " "))
    point_col = "darkgoldenrod3"
  }


  #run linear model
  if(length(model) > 1) {
    model_text = "custom"
    formula = reformulate(model, response = "observed")
    forecast_lm = lm(formula, data = model_df_selected) #full model using best forecast and socio-environmental variables
  } else {
    model_text = model
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
  }

  #Calculate predictions and predictive performance
  pred_df = data.frame(pred = predict(forecast_lm),
                       observed = forecast_lm$model$observed)
  R2 = GOF(pred_df$pred, pred_df$observed) #goodness-of-fit (R2 over 1:1 line)
  mae = MAE(pred_df$pred, pred_df$observed) #mean absolute  error (MAE)
  rmse = RMSE(pred_df$pred, pred_df$observed) #root mean squared error (RMSE)

  #Plot model (observed vs predicted project carbon loss)
  plot_model = ggplot(data = pred_df) +
    geom_point(aes(x = pred, y = observed), size = 5, col = point_col) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    annotate(geom = "text", x = text_x, y = text_y[1], size = text_size,
             label = paste("MAE:", round(mae, 2)), hjust = 0) +
    annotate(geom = "text", x = text_x, y = text_y[2], size = text_size,
             label = paste("RMSE:", round(rmse, 2)), hjust = 0) +
    annotate(geom = "text", x = text_x, y = text_y[3], size = text_size,
             label = paste("Goodness-of-fit:", round(R2, 3)), hjust = 0) +
    labs(title = title_type,
         x = x_lab,
         y = y_lab) +
    scale_x_continuous(limits = lim_val, breaks = break_val) +
    scale_y_continuous(limits = lim_val, breaks = break_val) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "cm"),
          plot.title = element_text(size = 38, hjust = 0.5),
          axis.title = element_text(size = 36),
          axis.text = element_text(size = 36),
          axis.ticks = element_blank(),
          axis.line = element_line(color = "black"))
  
  return(list(model = forecast_lm, pred = pred_df, R2 = R2, mae = mae, rmse = rmse, plot = plot_model))
}