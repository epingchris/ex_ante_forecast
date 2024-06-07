CalcExAnte = function(proj_id, vicinity, obs_val) {
  c = Sys.time()

  #calculate pre-project vicinity baseline deforestation and C loss
  meanCLoss = function(dat, ind) mean(dat[ind, ]$c_loss, na.rm = T)
  boot_n = 1000
  baseline_c_loss = boot::boot(vicinity, meanCLoss, R = boot_n) #around 700 seconds

  vicinity_filtered = subset(vicinity, !exclude)
  if(nrow(vicinity_filtered) == 0) {
    return(list(baseline_c_loss = baseline_c_loss, baseline_c_loss_filtered = NULL,
                ses = NULL, plot_df = NULL, risk_overcrediting = NULL, risk_reversal = NULL))
  }
  baseline_c_loss_filtered = boot::boot(vicinity_filtered, meanCLoss, R = boot_n) #around 700 seconds

  baseline_val = data.frame(base_c_loss = as.vector(baseline_c_loss$t),
                            base_c_loss_filt = as.vector(baseline_c_loss_filtered$t))
  ses = (mean(baseline_val$base_c_loss_filt) - mean(baseline_val$base_c_loss)) / sd(baseline_val$base_c_loss)

  forecast_median = median(baseline_val$base_c_loss_filt)
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

  #compile plotting data frame
  baseline_long = pivot_longer(baseline_val, everything(), names_to = "Type", values_to = "Value")
  obs_long = obs_val %>%
    rename(obs_c_loss = c_loss, obs_t_loss = t_loss) %>%
    pivot_longer(everything(), names_to = "Type", values_to = "Value")
  plot_df = rbind(baseline_long, obs_long)

  d = Sys.time()
  cat(proj_id, ":", d - c, "\n")
  return(list(baseline_c_loss = baseline_c_loss, baseline_c_loss_filtered = baseline_c_loss_filtered,
              ses = ses, plot_df = plot_df, risk_overcrediting = risk_overcrediting, risk_reversal = risk_reversal))
}
