plotBaseline = function(dat_wide = c_loss_ongoing_summ_wide, type, use_log10) {

    dat_long = dat_wide %>%
        pivot_longer(best:offset, names_to = "baseline_type", values_to = "baseline")

    plot_partial = type %in% c("best", "loose", "offset")
    if(plot_partial) {
        dat_long = subset(dat_long, baseline_type == type)
    }

    fig_name = case_when(
        type == "best" ~ "figure6a_ongoing_baseline_best_vs_c_loss_cf",
        type == "loose" ~ "figure6b_ongoing_baseline_loose_vs_c_loss_cf",
        type == "offset" ~ "figure6c_ongoing_baseline_offset_vs_c_loss_cf",
        .default = "figure6_ongoing_baseline_vs_c_loss_cf")

    fig_color = case_when(
        type == "best" ~ "blue",
        type == "loose" ~ "red",
        type == "offset" ~ "purple",
        .default = "black")

    hjust_val = ifelse(use_log10, 1.3, -0.5)

    if(use_log10) {
        dat_long = dat_long %>%
            mutate(cf_c_loss_log10 = log10(cf_c_loss),
                   baseline_log10 = log10(baseline))
        lm_out = lm(cf_c_loss_log10 ~ baseline_log10, data = dat_long)
        #loose: skewed by 1201 - Gola; offset: skewed by 1541
    } else {
        lm_out = lm(cf_c_loss ~ baseline, data = dat_long)
        #loose: skewed by 1201 - Gola; offset: skewed by 1541
    }

    p = ggplot(data = dat_long, aes(x = baseline, y = cf_c_loss))
    if(plot_partial) {
        p = p +
            geom_point(aes(shape = baseline_type, color = baseline_type, fill = baseline_type), size = 4) +
            geom_text(data = dat_wide, aes(x = .data[[type]], y = cf_c_loss, label = project),
                      hjust = hjust_val, size = 5)
    } else {
        p = p +
            geom_segment(data = dat_wide, aes(x = c_loss_min, xend = c_loss_max, y = cf_c_loss), linetype = 3) +
            geom_point(aes(shape = baseline_type, color = baseline_type, fill = baseline_type), size = 4) +
            geom_text(data = dat_wide, aes(x = c_loss_min, y = cf_c_loss, label = project),
                      hjust = hjust_val, size = 5)
    }
    p = p +
        geom_abline(intercept = 0, slope = 1, linetype = 3) +
        scale_shape_manual(values = c(best = 16, loose = 18, offset = 17),
                        labels = c("Best-matched", "Loosely-matched", "Time-lagged")) +
        scale_color_manual(values = c(best = "blue", loose = "red", offset = "purple"),
                        labels = c("Best-matched", "Loosely-matched", "Time-lagged")) +
        scale_fill_manual(values = c(best = "blue", loose = "red", offset = "purple"),
                        labels = c("Best-matched", "Loosely-matched", "Time-lagged")) +
        labs(x = "Baseline carbon loss (MgC/ha/yr)",
                y = "Observed counterfactual carbon loss (MgC/ha/yr)",
                shape = "Baseline type",
                color = "Baseline type",
                fill = "Baseline type") +
        theme_bw() +
        theme(panel.grid = element_blank(),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 16))

    if(plot_partial) {
        p = p +
            geom_abline(intercept = lm_out$coefficients[1], slope = lm_out$coefficients[2], linetype = 2, color = fig_color) +
            theme(legend.position = "none")
    } else {
        p = p +
            theme(legend.title = element_text(size = 16),
                  legend.text = element_text(size = 14),
                  legend.position = "bottom")
    }

    if(use_log10) {
        p = p +
            scale_x_continuous(limits = c(10 ^ (-4), 10 ^ 0.2), expand = c(0.01, 0.01),
                               transform = transform_log10(),
                               breaks = trans_breaks("log10", function(x) 10 ^ x),
                               labels = trans_format("log10", math_format(10 ^ .x))) +
            scale_y_continuous(limits = c(10 ^ (-1.8), 10 ^ 0.2), expand = c(0.01, 0.01),
                               transform = transform_log10(),
                               breaks = trans_breaks("log10", function(x) 10 ^ x),
                               labels = trans_format("log10", math_format(10 ^ .x)))
    } else {
        p = p +
            scale_x_continuous(limits = c(0, 1.45), expand = c(0.01, 0.01)) +
            scale_y_continuous(limits = c(0, 1.45), expand = c(0.01, 0.01))
    }
    file_path = paste0(fig_path, fig_name, ifelse(use_log10, "_log10", ""), ".png")
    cat(file_path, "\n")
    ggsave(file_path, plot = p, width = 3000, height = 3500, units = "px")

        # scale_x_continuous(limits = c(exp(-8), exp(0.4)), expand = c(0.01, 0.01),
        #                 transform = transform_log(),
        #                 breaks = trans_breaks("log", function(x) exp(x)),
        #                 labels = trans_format("log", math_format(e ^ .x))) +
        # scale_y_continuous(limits = c(exp(-8), exp(0.4)), expand = c(0.01, 0.01),
        #                 transform = transform_log(),
        #                 breaks = trans_breaks("log", function(x) exp(x)),
        #                 labels = trans_format("log", math_format(e ^ .x))) +
        # scale_x_continuous(limits = c(2 ^ (-12), 2 ^ 0.6), expand = c(0.01, 0.01),
        #                 transform = transform_log2(),
        #                 breaks = trans_breaks("log2", function(x) 2 ^ x),
        #                 labels = trans_format("log2", math_format(2 ^ .x))) +
        # scale_y_continuous(limits = c(2 ^ (-6), 2 ^ 0.6), expand = c(0.01, 0.01),
        #                 transform = transform_log2(),
        #                 breaks = trans_breaks("log2", function(x) 2 ^ x),
        #                 labels = trans_format("log2", math_format(2 ^ .x))) +
}
