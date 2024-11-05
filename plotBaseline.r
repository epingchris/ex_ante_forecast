plotBaseline = function(dat, baseline_used, use_log10) {

    dat_long = dat %>%
        pivot_longer(best:offset, names_to = "baseline_type", values_to = "baseline")

    plot_partial = baseline_used %in% c("best", "loose", "offset")
    if(plot_partial) {
        dat_long = subset(dat_long, baseline_type == baseline_used)
    }

    fig_num = switch(analysis_type,
        "control" = "5",
        "ongoing" = "6")
    fig_ind = switch(baseline_used,
        "all" = "",
        "best" = "a",
        "loose" = "b",
        "offset" = "c")
    fig_baseline_type = switch(baseline_used,
        "all" = "baseline",
        "best" = "baseline_best",
        "loose" = "baseline_loose",
        "offset" = "baseline_offset")
    fig_name = paste0("figure", fig_num, fig_ind, "_", analysis_type, "_", fig_baseline_type, "_vs_cf_c_loss")

    fig_color = switch(baseline_used,
        "all" = "black",
        "best" = "blue",
        "loose" = "red",
        "offset" = "purple")
    fig_title = switch(baseline_used,
        "all" = "",
        "best" = "Best-matched",
        "loose" = "Loosely-matched",
        "offset" = "Time-lagged")

    if(use_log10) {
        dat_long = dat_long %>%
            filter(cf_c_loss > 0 & baseline > 0) %>%
            mutate(cf_c_loss_log10 = log10(cf_c_loss),
                   baseline_log10 = log10(baseline))
        lm_out = lm(cf_c_loss_log10 ~ baseline_log10, data = dat_long)
        hjust_val = 1.3
        if(analysis_type == "control") {
            x_limit = c(10 ^ (-4), 10 ^ 0.4)
            y_limit = c(10 ^ (-2.2), 10 ^ 0.4)
        } else if(analysis_type == "ongoing") {
            x_limit = c(10 ^ (-4), 10 ^ 0.2)
            y_limit = c(10 ^ (-1.8), 10 ^ 0.2)
        }
    } else {
        lm_out = lm(cf_c_loss ~ baseline, data = dat_long)
        hjust_val = -0.5
        if(analysis_type == "control") {
            x_limit = c(-0.02, 2.36)
            y_limit = c(-0.02, 2.36)
        } else if(analysis_type == "ongoing") {
            x_limit = c(0, 1.45)
            y_limit = c(0, 1.45)
        }
    }

    if(plot_partial) {
        dat = dat %>%
            mutate(text_pos = .data[[baseline_used]])
    } else {
        dat = dat %>%
            mutate(text_pos = c_loss_min)
    }

    p = ggplot(data = dat_long, aes(x = baseline, y = cf_c_loss)) +
        geom_segment(data = dat, aes(x = c_loss_min, xend = c_loss_max, y = cf_c_loss), linetype = ifelse(plot_partial, 0, 3))
    if(analysis_type == "control") {
        p = p +
            geom_point(aes(shape = Continent, color = baseline_type, fill = baseline_type), size = 3) +
            scale_shape_manual(values = c(Asia = 1, Africa = 3, `South America` = 18)) +
            labs(shape = "Continent")
    } else {
        p = p +
            geom_point(aes(shape = baseline_type, color = baseline_type, fill = baseline_type), size = 3) +
            geom_text(data = dat, aes(x = text_pos, y = cf_c_loss, label = project),
                    hjust = hjust_val, size = 7) +
            scale_shape_manual(values = c(best = 16, loose = 18, offset = 17),
                            labels = c("Best-matched", "Loosely-matched", "Time-lagged")) +
            labs(shape = "Baseline type")
    }
    p = p +
        geom_abline(intercept = 0, slope = 1, linetype = 3) +
        geom_abline(intercept = lm_out$coefficients[1], slope = lm_out$coefficients[2], linetype = ifelse(plot_partial, 2, 0), color = fig_color) +
        scale_color_manual(values = c(best = "blue", loose = "red", offset = "purple"),
                           labels = c("Best-matched", "Loosely-matched", "Time-lagged"), guide = ifelse(plot_partial, "none", "legend")) +
        scale_fill_manual(values = c(best = "blue", loose = "red", offset = "purple"),
                          labels = c("Best-matched", "Loosely-matched", "Time-lagged"), guide = ifelse(plot_partial, "none", "legend")) +
        scale_x_continuous(limits = x_limit, expand = c(0.01, 0.01),
                           transform = ifelse(use_log10, "log10", "identity"),
                           breaks = if (use_log10) trans_breaks("log10", function(x) 10 ^ x) else waiver(),
                           labels = if (use_log10) trans_format("log10", math_format(10 ^ .x)) else waiver()) +
        scale_y_continuous(limits = y_limit, expand = c(0.01, 0.01),
                           transform = ifelse(use_log10, "log10", "identity"),
                           breaks = if (use_log10) trans_breaks("log10", function(x) 10 ^ x) else waiver(),
                           labels = if (use_log10) trans_format("log10", math_format(10 ^ .x)) else waiver()) +
        labs(title = fig_title,
             x = "Baseline carbon loss (MgC/ha/yr)",
             y = "Observed counterfactual carbon loss (MgC/ha/yr)",
             color = "Baseline type",
             fill = "Baseline type") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title = element_text(size = 22),
              axis.text = element_text(size = 20),
              legend.title = element_text(size = 22),
              legend.text = element_text(size = 20),
              legend.position = "bottom",
              legend.box = "vertical",
              plot.title = element_text(size = 24, hjust = 0.5))

    if(plot_partial & analysis_type == "ongoing") {
        p = p +
            theme(legend.position = "none")
    }

    file_path = paste0(fig_path, fig_name, ifelse(use_log10, "_log10", ""), ".png")
#    cat(file_path, "\n")
    ggsave(file_path, plot = p, width = 2500, height = 5000, units = "px")
    return(p)
}
