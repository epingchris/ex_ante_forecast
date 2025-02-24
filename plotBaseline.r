plotBaseline = function(dat, baseline_used, metrics, corr = F, add_label = F) {
    summ_wide_mean = dat %>%
        dplyr::select(type, mean, project) %>%
        pivot_wider(names_from = "type", values_from = "mean") %>%
        arrange(project)

    summ_wide_ci_lower = dat %>%
        dplyr::select(type, ci_lower, project) %>%
        pivot_wider(names_from = "type", values_from = "ci_lower") %>%
        arrange(project)

    summ_wide_ci_upper = dat %>%
        dplyr::select(type, ci_upper, project) %>%
        pivot_wider(names_from = "type", values_from = "ci_upper") %>%
        arrange(project)

    summ_ci_x = data.frame(project = summ_wide_mean$project,
                           baseline_type = baseline_used,
                           y = summ_wide_mean$cf_c_loss,
                           x0 = summ_wide_ci_lower[, baseline_used],
                           x1 = summ_wide_ci_upper[, baseline_used])
    colnames(summ_ci_x) = c("project", "baseline_type", "y", "x0", "x1")

    summ_ci_y = data.frame(project = summ_wide_mean$project,
                           baseline_type = baseline_used,
                           x = summ_wide_mean[, baseline_used],
                           y0 = summ_wide_ci_lower$cf_c_loss,
                           y1 = summ_wide_ci_upper$cf_c_loss)
    colnames(summ_ci_y) = c("project", "baseline_type", "x", "y0", "y1")

    dat_wide = dat %>%
        dplyr::select(type, mean, project, code) %>%
        pivot_wider(names_from = "type", values_from = "mean") %>%
        mutate(c_loss_min = pmin(best, loose, lagged, na.rm = T),
               c_loss_max = pmax(best, loose, lagged, na.rm = T))

    dat_long = dat_wide %>%
        pivot_longer(best:lagged, names_to = "baseline_type", values_to = "baseline")

    dat_long = subset(dat_long, baseline_type == baseline_used)
    dat_wide = dat_wide %>%
        mutate(text_pos = .data[[baseline_used]])

    lm_out = lm(cf_c_loss ~ baseline - 1, data = dat_long)
    lm_slope = summary(lm_out)$coefficients[1]
    x_limit = c(-0.2, 2.55)
    y_limit = c(-0.2, 2.55)


    if(corr) {
        note = tibble(
            x = 2.5,
            y = c(0, -0.2),
            text = list(bquote("MAE: " * .(metrics[1])),
                        bquote("Correction factor: " * .(metrics[2]) * " [" * .(metrics[3]) * ", " * .(metrics[4]) * "]")
            )
        )
    } else {
        note = tibble(
            x = 2.5,
            y = c(0, -0.2),
            text = list(bquote("MAE: " * .(metrics[1])), NULL)
        )
    }

    p = ggplot(data = dat_long, aes(x = baseline, y = cf_c_loss))
    if(add_label) {
        p = p +
            geom_text(aes(label = code), hjust = -0.5, size = 5, parse = T)
    }
    p = p +
        geom_segment(data = summ_ci_x, aes(x = x0, xend = x1, y = y, color = baseline_type), linewidth = 1.2) +
        geom_segment(data = summ_ci_y, aes(x = x, y = y0, yend = y1, color = baseline_type), linewidth = 1.2) +
        geom_abline(aes(intercept = 0, slope = lm_slope, color = baseline_type), linetype = 1, linewidth = 0) +
        geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 1.5) +
        geom_point(aes(shape = baseline_type, color = baseline_type, fill = baseline_type), size = 5) +
        geom_text(data = note, aes(x = x, y = y, label = text), hjust = 1, vjust = 0, size = 9, parse = T) +
        scale_shape_manual(values = c(best = 16, loose = 18, lagged = 17)) +
        scale_color_manual(values = c(best = "#40B0A6", loose = "#006CD1", lagged = "#CDAC60")) +
        scale_fill_manual(values = c(best = "#40B0A6", loose = "#006CD1", lagged = "#CDAC60")) +
        scale_x_continuous(limits = x_limit, expand = c(0.01, 0.01)) +
        scale_y_continuous(limits = y_limit, expand = c(0.01, 0.01)) +
        labs(x = bquote("Forecast carbon loss in baseline (MgC" ~ ha^"-1" ~ yr^"-1" * ")"),
             y = bquote(italic("ex post") ~ "carbon loss in baseline (MgC" ~ ha^"-1" ~ yr^"-1" * ")")) +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", fill = NA),
              panel.grid = element_blank(),
              plot.title = element_blank(),
              axis.title = element_text(size = 32),
              axis.title.x = element_text(margin = margin(t = 30)),
              axis.title.y = element_text(margin = margin(r = 30)),
              axis.text = element_text(size = 28),
              axis.text.x = element_text(margin = margin(t = 15)),
              axis.text.y = element_text(margin = margin(r = 15)),
              axis.ticks = element_line(linewidth = 2),
              axis.ticks.length = unit(.5, "cm"),
              legend.position = "none")

    return(p)
}
