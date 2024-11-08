plotBaseline = function(dat, baseline_used, use_log10, size_vec) {
    dat_y_mean = dat %>%
        filter(type == "cf_c_loss") %>%
        dplyr::select(project, y_mean = mean)

    dat_ci_x = dat %>%
        filter(type == baseline_used) %>%
        dplyr::select(project, ci_lower, ci_upper, baseline_type = type) %>%
        left_join(dat_y_mean, by = "project")

    dat_x_mean = dat %>%
        filter(type == baseline_used) %>%
        dplyr::select(project, x_mean = mean, baseline_type = type)

    dat_ci_y = dat %>%
        filter(type == "cf_c_loss") %>%
        dplyr::select(project, ci_lower, ci_upper) %>%
        left_join(dat_x_mean, by = "project")

    dat_wide = dat %>%
        dplyr::select(-c(ci_lower, ci_upper)) %>%
        pivot_wider(names_from = "type", values_from = "mean", id_expand = T) %>%
        mutate(c_loss_min = pmin(best, loose, offset, na.rm = T),
               c_loss_max = pmax(best, loose, offset, na.rm = T))

    dat_long = dat_wide %>%
        pivot_longer(best:offset, names_to = "baseline_type", values_to = "baseline")

    plot_partial = baseline_used %in% c("best", "loose", "offset")
    if(plot_partial) {
        dat_long = subset(dat_long, baseline_type == baseline_used)
        dat_wide = dat_wide %>%
            mutate(text_pos = .data[[baseline_used]])
    } else {
        dat_wide = dat_wide %>%
            mutate(text_pos = c_loss_min)
    }

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
            y_limit = c(10 ^ (-1.6), 10 ^ 0.3)
        }
    } else {
        lm_out = lm(cf_c_loss ~ baseline, data = dat_long)
        hjust_val = -0.5
        if(analysis_type == "control") {
            x_limit = c(-0.02, 2.36)
            y_limit = c(-0.02, 2.36)
        } else if(analysis_type == "ongoing") {
            x_limit = c(-0.2, 2.55)
            y_limit = c(-0.2, 2.55)
        }
    }

    p = ggplot(data = dat_long, aes(x = baseline, y = cf_c_loss)) +
        geom_segment(data = dat_wide, aes(x = c_loss_min, xend = c_loss_max, y = cf_c_loss), linetype = ifelse(plot_partial, 0, 3)) +
        geom_segment(data = dat_ci_x, aes(x = ci_lower, xend = ci_upper, y = y_mean, color = baseline_type)) +
        geom_segment(data = dat_ci_y, aes(x = x_mean, y = ci_lower, yend = ci_upper, color = baseline_type))
    if(analysis_type == "control") {
        p = p +
            geom_point(aes(shape = Continent, color = baseline_type, fill = baseline_type), size = 3) +
            scale_shape_manual(values = c(Asia = 1, Africa = 5, `South America` = 18)) +
            labs(shape = "Continent")
    } else {
        p = p +
            geom_point(aes(shape = baseline_type, color = baseline_type, fill = baseline_type), size = 3) +
            geom_text(data = dat_wide, aes(x = text_pos, y = cf_c_loss, label = project),
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
              legend.position = "bottom",
              legend.box = "vertical",
              axis.title = element_text(size = size_vec[2]),
              axis.text = element_text(size = size_vec[3]),
              legend.title = element_text(size = size_vec[2]),
              legend.text = element_text(size = size_vec[3]),
              plot.title = element_text(size = size_vec[1], hjust = 0.5))

    if(plot_partial & analysis_type == "ongoing") {
        p = p +
            theme(legend.position = "none")
    }

    return(p)
}
