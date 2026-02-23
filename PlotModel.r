PlotModel = function(input, obs_var, model) {

  #select data to use
  envir_var = c("area_ha", "prj_slope", "prj_remote", "gdppc_mean", "gdppc_rate", "wgicc_mean")
  #slope and elevation collinear: remove elevation
  #initial carbon density and forecasted counterfactual C loss collinear: remove initial carbon density

  model_df_selected = input %>%
    dplyr::select(starts_with(c(envir_var, obs_var, "forecast"))) %>%
    rename(observed = all_of(obs_var)) %>%
    mutate(forecast = forecast * 100) #convert to percentage
  if (obs_var != "credit") {
    model_df_selected = model_df_selected %>%
      mutate(observed = observed * 100) #convert to percentage
  }

  #set graphical parameters
  lim_x = c(0, 1.2)
  break_x = seq(0, 1.2, by = 0.2)
  text_x = 1.2
  axtitle_x = expression(atop("Forecasted counterfactual", "carbon loss rate (%)"))

  lim_y = switch(obs_var,
                 "closs_obs_cf" = c(-1, 5),
                 "closs_obs_p" = c(-1, 5),
                 "add_rate" = c(-1, 5),
                 "credit" = c(-2, 8))
  break_y = switch(obs_var,
                   "closs_obs_cf" = seq(-1, 5),
                   "closs_obs_p" = seq(-1, 5),
                   "add_rate" = seq(-1, 5),
                   "credit" = seq(-2, 8))
  text_y = switch(obs_var,
                  "closs_obs_cf" = c(0.2, -0.2, -0.6, -1),
                  "closs_obs_p" = c(0.2, -0.2, -0.6, -1),
                  "add_rate" = c(0.2, -0.2, -0.6, -1),
                  "credit" = c(0.1, -0.6, -1.3, -2))
  axtitle_y = switch(obs_var,
                     "closs_obs_cf" = expression(atop("Observed counterfactual", "carbon loss (%)")),
                     "closs_obs_p" = expression(atop("Observed project", "carbon loss (%)")),
                     "add_rate" = expression(atop("Observed counterfactual - project", "carbon loss (%)")),
                     "credit" = expression(atop("Observed carbon credit production", paste("(MgC", " ", ha^-1, " ", yr^-1, ")", sep = " "))))

  label_pos_y = switch(obs_var,
                       "closs_obs_cf" = 4.75,
                       "closs_obs_p" = 4.75,
                       "add_rate" = 4.75,
                       "credit" = 7.5)
  label_text = switch(obs_var,
                      "closs_obs_cf" = "A",
                      "closs_obs_p" = "B",
                      "add_rate" = "C",
                      "credit" = "D")

  point_col = switch(obs_var,
                     "closs_obs_cf" = "darkred",
                     "closs_obs_p" = "blue",
                     "add_rate" = "darkgreen",
                     "credit" = "darkgoldenrod3")
  text_size = 10

  #run linear model
  if(model == "full") {
    #full model using best forecast and socio-environmental variables
    forecast_lm = lm(observed ~ ., data = model_df_selected) 
  } else if(model == "naive") {
    #simple model using only best forecast
    forecast_lm = lm(observed ~ forecast, data = model_df_selected)
  } else if(model == "sel") {
    #backward model selection from full model
    forecast_lm_full = lm(observed ~ ., data = model_df_selected) 
    forecast_lm = stepAIC(forecast_lm_full,
                          scope = list(upper = ~ ., lower = ~ forecast),
                          direction = "backward", trace = 1)
  } else {
    formula = reformulate(model, response = "observed")
    forecast_lm = lm(formula, data = model_df_selected) #full model using best forecast and socio-environmental variables
  }

  #extract model coefficients
  intercept_val = round(coef(forecast_lm)[1], 2)
  slope_val = round(coef(forecast_lm)[2], 2)
  formula_text = paste0("Observed = ", intercept_val, " + Forecast * ", slope_val)

  #calculate predictions and predictive performance
  pred_df = data.frame(pred = predict(forecast_lm),
                       observed = forecast_lm$model$observed)
  R2 = GOF(pred_df$pred, pred_df$observed) #goodness-of-fit (R2 over 1:1 line)
  mae = MAE(pred_df$pred, pred_df$observed) #mean absolute  error (MAE)
  rmse = RMSE(pred_df$pred, pred_df$observed) #root mean squared error (RMSE)

  #plot model (observed vs predicted project carbon loss)
  plot_model = ggplot(data = model_df_selected) +
    geom_point(aes(x = forecast, y = observed), size = 5, col = point_col) +
    geom_abline(intercept = intercept_val, slope = slope_val, linewidth = 1, col = point_col) +
    {if (obs_var != "credit") geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2)} +
    geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray", linewidth = 2) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray", linewidth = 2) +
    annotate(geom = "text", x = 0.05, y = label_pos_y, size = text_size * 2, fontface = "bold",
             label = label_text) +
    annotate(geom = "text", x = text_x, y = text_y[1], size = text_size,
             label = paste("MAE:", round(mae, 2)), hjust = 1) +
    annotate(geom = "text", x = text_x, y = text_y[2], size = text_size,
             label = paste("RMSE:", round(rmse, 2)), hjust = 1) +
    annotate(geom = "text", x = text_x, y = text_y[3], size = text_size,
             label = paste("Goodness-of-fit:", round(R2, 3)), hjust = 1) +
    annotate(geom = "text", x = text_x, y = text_y[4], size = text_size,
             label = formula_text, hjust = 1) +
    labs(x = axtitle_x, y = axtitle_y) +
    scale_x_continuous(limits = lim_x, breaks = break_x) +
    scale_y_continuous(limits = lim_y, breaks = break_y) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          plot.margin = margin(t = 20, b = 20, l = 50),
          plot.title = element_text(size = 38, hjust = 0.5),
          axis.title.x = element_text(size = 34, margin = margin(t = 0, b = 15)),
          axis.title.y = element_text(size = 34, vjust = 0.5),
          axis.text.x = element_text(size = 30, margin = margin(t = 10, b = 20)),
          axis.text.y = element_text(size = 30), margin = margin(l = 10, r = 20),
          axis.ticks = element_blank(),
          axis.line = element_line(color = "black"))

  return(list(model = forecast_lm, pred = pred_df, R2 = R2, mae = mae, rmse = rmse, plot = plot_model))
}