plotPlacebo = function(dat, label_to_x, col = "black") {
    dat = dat %>%
        mutate(type = replace(type, type != label_to_x, "estimate")) %>%
        mutate(type = replace(type, type == label_to_x, "measure"))

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
                           y = summ_wide_mean$estimate,
                           x0 = summ_wide_ci_lower$measure,
                           x1 = summ_wide_ci_upper$measure)

    summ_ci_y = data.frame(project = summ_wide_mean$project,
                           x = summ_wide_mean$measure,
                           y0 = summ_wide_ci_lower$estimate,
                           y1 = summ_wide_ci_upper$estimate)

    lm_out = lm(estimate ~ measure - 1, data = summ_wide_mean)
    slope = summary(lm_out)$coefficients[1, 1] %>% round(., 2)
    slope_se = summary(lm_out)$coefficients[1, 2]
    slope_df = summary(lm_out)$df[2]
    slope_t = qt(0.975, slope_df)
    slope_ci = c(slope - slope_se * slope_t, slope + slope_se * slope_t) %>% round(., 2)

    r2 = summary(lm_out)$r.squared %>% round(., 3)
    note_df = tibble(
        x = 0,
        y = c(2.4, 2.2),
        text = list(bquote("Slope: " * .(slope) * " [" * .(slope_ci[1]) * ", " * .(slope_ci[2]) * "]"),
                    bquote(R^2 * ": " * .(r2)))
    )

    t_out = t.test(summ_wide_mean$estimate, summ_wide_mean$measure, paired = T)

    p = ggplot(data = summ_wide_mean, aes(x = measure, y = estimate)) +
        geom_point(size = 6, color = col) +
        geom_segment(data = summ_ci_x, aes(x = x0, xend = x1, y = y), color = col, linewidth = 1.2) +
        geom_segment(data = summ_ci_y, aes(x = x, y = y0, yend = y1), color = col, linewidth = 1.2) +
        geom_text(data = note_df, aes(x = x, y = y, label = text), hjust = 0, size = 17, parse = T) +
        geom_abline(intercept = 0, slope = 1, linetype = 3, linewidth = 1.5) +
        scale_x_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) + #ensures no padding
        scale_y_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) +
        labs(x = bquote("Observed carbon loss (MgC" ~ ha^"-1" ~ yr^"-1" * ")"),
             y = bquote("Estimated carbon loss (MgC" ~ ha^"-1" ~ yr^"-1" * ")")) +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", fill = NA),
              panel.grid = element_blank(),
              axis.title = element_text(size = 48),
              axis.title.x = element_text(margin = margin(t = 20)),
              axis.title.y = element_text(margin = margin(r = 20)),
              axis.text = element_text(size = 44),
              axis.text.x = element_text(margin = margin(t = 10)),
              axis.text.y = element_text(margin = margin(r = 10)),
              axis.ticks = element_line(linewidth = 2),
              axis.ticks.length = unit(.5, "cm"),
              legend.key.spacing = unit(1, "cm"),
              legend.title = element_text(size = 32, margin = margin(t = 10)),
              legend.text = element_text(size = 32, margin = margin(t = 10)))

    return(list(p = p, lm = lm_out, t = t_out))
}
