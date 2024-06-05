CalcExAnte = function(proj_id, vicinity, obs_val) {
  c = Sys.time()

  #initialise plots
  p_filter = NULL
  p_c_loss = NULL
  baseline_c_loss_filtered = NULL
  baseline_c_loss_filtered_val = NA

  #calculate pre-project vicinity baseline deforestation and C loss
  meanCLoss = function(dat, ind) mean(dat[ind, ]$c_loss, na.rm = T)
  boot_n = 1000
  baseline_c_loss = boot::boot(vicinity, meanCLoss, R = boot_n) #around 700 seconds
  baseline_c_loss_val = as.vector(baseline_c_loss$t)

#  vicinity_filtered = subset(vicinity, !exclude)
  vicinity_filtered = subset(vicinity, !(slope_exclude | elevation_exclude | access_exclude))
  if(nrow(vicinity_filtered) == 0) {
    return(list(baseline_c_loss = baseline_c_loss, baseline_c_loss_filtered = NULL,
                plot_df = NULL, p_filter = NULL, p_c_loss = NULL,
                risk_overcrediting = NULL, risk_reversal = NULL))
  }
  baseline_c_loss_filtered = boot::boot(vicinity_filtered, meanCLoss, R = boot_n) #around 700 seconds
  baseline_c_loss_filtered_val = as.vector(baseline_c_loss_filtered$t)

  #compile plotting data frame
  plot_df = data.frame(Value = c(baseline_c_loss_val, baseline_c_loss_filtered_val),
                       Type = c(rep("base_c_loss", boot_n),
                                rep("base_c_loss_filt", boot_n)))

  if(!is.null(obs_val)) {
    plot_df = plot_df %>%
      rbind(., data.frame(Value = obs_val$c_loss,
                          Type = rep("obs_c_loss", nrow(obs_val))))
  }

  risk_overcrediting = sum(obs_val$c_loss < median(baseline_c_loss_filtered_val)) / nrow(obs_val)
  risk_reversal = sum(obs_val$c_loss < 0) / nrow(obs_val)  #risk of reversal

  #make plot title
  proj_id_bare = proj_id
#  proj_id_bare = proj_id %>% str_remove("a")
  plot_title = ifelse(analysis_type == "control", paste("Control", proj_id_bare), proj_id_bare)

  #make plots comparing before/after filter
  p_filter = ggplot(data = filter(plot_df, Type != "obs_c_loss"), aes(x = Type, y = Value)) +
      geom_boxplot(aes(color = Type)) +
#      ggpubr::stat_compare_means(method = "t.test", aes(label = ..p.signif..),
#                                 label.x = 1.5, size = 8) +
      scale_x_discrete(labels = c("Unfiltered", "Filtered")) +
      scale_color_manual(values = c("red", "blue"),
                         labels = c("Unfiltered", "Filtered")) +
      labs(x = "", y = "Annual carbon loss (Mg/ha)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))

  #make plots comparing before/after filter
  if(!is.null(obs_val)) {
    p_c_loss = ggplot(data = filter(plot_df, Type != "base_c_loss"), aes(x = Type, y = Value)) +
      geom_boxplot(aes(color = Type)) +
      scale_x_discrete(labels = c("Baseline", "Observed")) +
      scale_color_manual(values = c("red", "black"),
                         labels = c("Baseline", "Observed")) +
      scale_linetype_manual(values = c(2, 1),
                            labels = c("Baseline", "Observed")) +
      labs(x = "", y = "Annual carbon loss (Mg/ha)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  }

  d = Sys.time()
  cat(proj_id, ":", d - c, "\n")
  return(list(baseline_c_loss = baseline_c_loss, baseline_c_loss_filtered = baseline_c_loss_filtered,
              plot_df = plot_df, p_filter = p_filter, p_c_loss = p_c_loss,
              risk_overcrediting = risk_overcrediting, risk_reversal = risk_reversal))
}
