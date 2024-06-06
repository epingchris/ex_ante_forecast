CalcExAnte = function(proj_id, vicinity, obs_val) {
  c = Sys.time()

  #initialise plots
  baseline_c_loss_filtered_val = NA

  #calculate pre-project vicinity baseline deforestation and C loss
  meanCLoss = function(dat, ind) mean(dat[ind, ]$c_loss, na.rm = T)
  boot_n = 1000
  baseline_c_loss = boot::boot(vicinity, meanCLoss, R = boot_n) #around 700 seconds
  baseline_c_loss_val = as.vector(baseline_c_loss$t)

#  vicinity_filtered = subset(vicinity, !exclude)
#  vicinity_filtered = subset(vicinity, !(slope_exclude | elevation_exclude | access_exclude))
  vicinity_filtered = subset(vicinity, !exclude)
  if(nrow(vicinity_filtered) == 0) {
    return(list(baseline_c_loss = baseline_c_loss, baseline_c_loss_filtered = NULL,
                ses = NULL, plot_df = NULL, risk_overcrediting = NULL, risk_reversal = NULL))
  }
  baseline_c_loss_filtered = boot::boot(vicinity_filtered, meanCLoss, R = boot_n) #around 700 seconds
  baseline_c_loss_filtered_val = as.vector(baseline_c_loss_filtered$t)
  ses = (mean(baseline_c_loss_filtered_val) - mean(baseline_c_loss_val)) / sd(baseline_c_loss_val)

  #compile plotting data frame
  plot_df = data.frame(Value = c(baseline_c_loss_val, baseline_c_loss_filtered_val),
                       Type = c(rep("base_c_loss", boot_n),
                                rep("base_c_loss_filt", boot_n))) %>%
    rbind(., data.frame(Value = obs_val$c_loss,
                        Type = rep("obs_c_loss", nrow(obs_val))))

  forecast_median = median(baseline_c_loss_filtered_val)
  n_val = nrow(obs_val)
  scenario = c(1, 0.75, 0.5, 0.25)
  risk_overcrediting = as.data.frame(matrix(NA, 4, 3))
  for(i in 1:4) {
    risk_overcrediting[i, ] = c(scenario[i],
                                forecast_median * scenario[i],
                                sum(obs_val$additionality < forecast_median * scenario[i]) / n_val)
  }
  colnames(risk_overcrediting) = c("scenario", "forecast", "risk")
  risk_reversal = sum(obs_val$additionality < 0) / n_val  #risk of reversal

  d = Sys.time()
  cat(proj_id, ":", d - c, "\n")
  return(list(baseline_c_loss = baseline_c_loss, baseline_c_loss_filtered = baseline_c_loss_filtered,
              ses = ses, plot_df = plot_df, risk_overcrediting = risk_overcrediting, risk_reversal = risk_reversal))
}
