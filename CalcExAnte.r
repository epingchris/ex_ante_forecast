CalcExAnte = function(proj_id, area_ha, obs_val, path, forecast = F) {
  c = Sys.time()

  #load ACD and vicinitiy (potentially matched points) data
  acd = read.csv(paste0(path, "carbon-density.csv"))
    matches = read_parquet(paste0(path, "matches.parquet")) %>%
      slice_sample(n = 2500000) #sample down to around the smallest vicinity size of the 15 projects (2559309)

  #compile data of pre-project vicinity and during-project observed values
  vicinity_area = nrow(matches) * 900 / 10000 #convert from numbers of 30x30m2 pixels to hectares

  #find carbon loss associated with LUC change undisturbed -> deforested (MgC / ha)
  #exit if not available for classes 1 and 3
  acd_1 = acd %>% filter(land.use.class == 1) %>% pull(carbon.density)
  acd_3 = acd %>% filter(land.use.class == 3) %>% pull(carbon.density)
  acd_change = acd_1 - acd_3
  if(length(acd_change) == 0) return(list(vicinity_area = vicinity_area, plot_df = NULL, over_df = NULL,
                                          p0 = NULL, p1 = NULL, p2 = NULL, p_perc = NULL, p_overcredit = NULL))

  #calculate pre-project vicinity baseline C loss
  plot_df = matches %>%
    #convert from change in proportion of 30x30m2 pixels to change in MgC/ha
    mutate(carbon_loss_5_0 = acd_change * (cpc5_u - cpc0_u) / 5,
           carbon_loss_10_5 = acd_change * (cpc10_u - cpc5_u) / 5,
           carbon_loss_10_0 = acd_change * (cpc10_u - cpc0_u) / 10) %>%
    #since previous data shows that the three time intervals are almost the same, only use one (t-10 to t0)
    dplyr::select(carbon_loss_5_0:carbon_loss_10_0) %>%
    pivot_longer(cols = everything(), values_to = "Value", names_to = "Period") %>%
    mutate(Type = "baseline_c_loss")

  if(!forecast) {
    #add during-project observed values
    plot_df = plot_df %>% 
    rbind(., data.frame(Value = obs_val$c_loss,
                        Period = "after",
                        Type = "obs_c_loss")) %>% #Observed counterfactual carbon loss (MgC/ha)
    rbind(., data.frame(Value = obs_val$additionality,
                        Period = "after",
                        Type = "obs_add")) #Observed additionality (MgC/ha)

    #obtain different quantiles of baseline carbon loss and calculate probability of overcrediting
    len_pr = length(pr_vec)
    carbon_loss_q = quantile(plot_df %>% filter(Period == "carbon_loss_10_0") %>% pull(Value), pr_vec)
    p_over = sapply(carbon_loss_q, function(x) length(which(obs_val$additionality < x)) / nrow(obs_val))
    forecast_df = data.frame(proj = proj_id,
                             Pr = pr_vec,
                             Variable = rep(c("Percentile", "Overcrediting risk"), each = length(pr_vec)),
                             Value = c(carbon_loss_q, p_over))

    pr_max_zero = pr_vec[tail(which(carbon_loss_q == 0), 1)]
  }

  #make plo title
  proj_id_bare = proj_id %>% str_remove("a")
  plot_title = ifelse(analysis_type == "control", paste("Control", proj_id_bare), proj_id_bare)

  #make plots
  p0 = ggplot(data = plot_df %>% filter(Period != "after"), aes(Value, after_stat(density))) +
    geom_freqpoly(aes(color = Period)) +
    scale_color_manual(values = c("red", "black", "blue"),
                       labels = c(expression(t [-10]*" to "*t [0]),
                                  expression(t [-10]*" to "*t [-5]),
                                  expression(t [-5]*" to "*t [0]))) +
    xlab("Annual carbon flux (Mg/ha)") +
    ggtitle(plot_title) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.direction = "vertical")

  if(!forecast) {
    p_legend = ggplot(data = plot_df %>% filter(str_detect(Period, "10_0") | str_detect(Type, "obs")), aes(Value, after_stat(density))) +
      geom_freqpoly(aes(color = Type, linetype = Type)) +
      scale_color_manual(values = c("red", "black", "blue"),
                        labels = c("Pre-project carbon loss",
                                   "During-project additionality",
                                   "During-project counterfactual carbon loss")) +
      scale_linetype_manual(values = c(2, 1, 1),
                            labels = c("Pre-project carbon loss",
                                       "During-project additionality",
                                       "During-project counterfactual carbon loss")) +
      theme(legend.key = element_rect(fill = "white"),
            legend.position = "bottom",
            legend.direction = "vertical")
    p_legend_grob = ggpubr::get_legend(p_legend)

    p1 = ggplot(data = plot_df %>% filter(Period == "carbon_loss_10_0" | Type == "obs_c_loss"), aes(Value, after_stat(density))) +
      geom_freqpoly(aes(color = Type, linetype = Type)) +
      scale_color_manual(values = c("red", "black", "blue"),
                        labels = c("Pre-project carbon loss",
                                   "During-project counterfactual carbon loss")) +
      scale_linetype_manual(values = c(2, 1),
                            labels = c("Pre-project carbon loss",
                                       "During-project additionality",
                                       "During-project counterfactual carbon loss")) +
      xlab("Annual carbon flux (Mg/ha)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "none")

    p2 = ggplot(data = plot_df %>% filter(Period == "carbon_loss_10_0" | Type == "obs_add"), aes(Value, after_stat(density))) +
      geom_freqpoly(aes(color = Type, linetype = Type)) +
      scale_color_manual(values = c("red", "black"),
                        labels = c("Pre-project carbon loss",
                                    "During-project additionality")) +
      scale_linetype_manual(values = c(2, 1),
                            labels = c("Pre-project carbon loss",
                                      "During-project additionality")) +
      xlab("Annual carbon flux (Mg/ha)") +
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "none")

    p_perc = ggplot(data = forecast_df %>% filter(Variable == "Percentile")) +
      geom_line(aes(x = Pr, y = Value), color = "blue") +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_vline(xintercept = pr_max_zero, linetype = 2) +
      scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
      scale_y_continuous(limits = c(0, 5)) +
      labs(title = plot_title, x = "Percentage", y = "C loss percentile") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))

    p_overcredit = ggplot(data = forecast_df %>% filter(Variable == "Overcrediting risk")) +
      geom_line(aes(x = Pr, y = Value), color = "red") +
      geom_hline(yintercept = c(0.05, 0.1), linetype = 2) +
      geom_vline(xintercept = pr_max_zero, linetype = 2) +
      scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "", x = "Percentage", y = "Overcrediting risk") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  }

  d = Sys.time()
  cat(proj_id, ":", d - c, "\n")
  if(!forecast) {
    return(list(vicinity_area = vicinity_area, plot_df = plot_df, forecast_df = forecast_df,
                p0 = p0, p1 = p1, p2 = p2, p_legend_grob = p_legend_grob, p_perc = p_perc, p_overcredit = p_overcredit))
  } else {
    return(list(vicinity_area = vicinity_area, plot_df = plot_df, p0 = p0))
  }
}