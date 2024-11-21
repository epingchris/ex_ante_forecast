plotPlacebo = function(dat, period_used) {
    dat_used = dat %>%
        filter(period == period_used)

    summ_wide_mean = dat_used %>%
        dplyr::select(type, mean, project, Continent) %>%
        pivot_wider(names_from = "type", values_from = "mean") %>%
        arrange(project)

    summ_wide_ci_lower = dat_used %>%
        dplyr::select(type, ci_lower, project, Continent) %>%
        pivot_wider(names_from = "type", values_from = "ci_lower") %>%
        arrange(project)

    summ_wide_ci_upper = dat_used %>%
        dplyr::select(type, ci_upper, project, Continent) %>%
        pivot_wider(names_from = "type", values_from = "ci_upper") %>%
        arrange(project)

    summ_ci_x = data.frame(project = summ_wide_mean$project,
                           Continent = summ_wide_mean$Continent,
                           y = summ_wide_mean$cf_c_loss,
                           x0 = summ_wide_ci_lower$p_c_loss,
                           x1 = summ_wide_ci_upper$p_c_loss)

    summ_ci_y = data.frame(project = summ_wide_mean$project,
                           Continent = summ_wide_mean$Continent,
                           x = summ_wide_mean$p_c_loss,
                           y0 = summ_wide_ci_lower$cf_c_loss,
                           y1 = summ_wide_ci_upper$cf_c_loss)

    lm_control = lm(p_c_loss ~ cf_c_loss, data = summ_wide_mean)
    r2 = summary(lm_control)$r.squared %>% round(., 3)
    rmse = rmse(na.omit(summ_wide_mean)$p_c_loss, na.omit(summ_wide_mean)$cf_c_loss) %>% round(., 3)
    mae = mae(na.omit(summ_wide_mean)$p_c_loss, na.omit(summ_wide_mean)$cf_c_loss) %>% round(., 3)
    note_df = tibble(
        x = 2.5,
        y = c(0.4, 0.2, 0.0),
        text = list(bquote(R^2 ~ ": " ~ .(r2)),
                    bquote("RMSE: " ~ .(rmse)),
                    bquote("MAE: " ~ .(mae)))
    )

    title_text = ifelse(period_used == "pre", "a. Pre-start period", "b. Post-start period")

    p = ggplot(data = summ_wide_mean, aes(x = p_c_loss, y = cf_c_loss)) +
        geom_point(aes(shape = Continent, color = Continent, fill = Continent), size = 6) +
        geom_segment(data = summ_ci_x, aes(x = x0, xend = x1, y = y, color = Continent), linewidth = 1.2) +
        geom_segment(data = summ_ci_y, aes(x = x, y = y0, yend = y1, color = Continent), linewidth = 1.2) +
        geom_text(data = note_df, aes(x = x, y = y, label = text), hjust = 1, vjust = 0, size = 15, parse = T) +
        geom_abline(intercept = 0, slope = 1, linetype = 3) +
        scale_shape_manual(values = c(Asia = 16, Africa = 17, `South America` = 18)) +
        scale_color_manual(values = c(Asia = "#006CD1", Africa = "#40B0A6", `South America` = "#CDAC60")) +
        scale_fill_manual(values = c(Asia = "#006CD1", Africa = "#40B0A6", `South America` = "#CDAC60")) +
        scale_x_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) + #ensures no padding
        scale_y_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) +
        labs(title = title_text,
             x = "Placebo carbon loss (MgC/ha/yr)",
             y = "Counterfactual carbon loss (MgC/ha/yr)",
             color = "Continent  ",
             fill = "Continent  ",
             shape = "Continent  ") +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", fill = NA),
              panel.grid = element_blank(),
              plot.title = element_text(size = 48, hjust = 0.5, margin = margin(b = 10)),
              axis.title = element_text(size = 40),
              axis.title.x = element_text(margin = margin(t = 20)),
              axis.title.y = element_text(margin = margin(r = 20)),
              axis.text = element_text(size = 36),
              axis.text.x = element_text(margin = margin(t = 10)),
              axis.text.y = element_text(margin = margin(r = 10)),
              axis.ticks = element_line(linewidth = 2),
              axis.ticks.length = unit(.5, "cm"),
              legend.key.spacing = unit(1, "cm"),
              legend.title = element_text(size = 32, margin = margin(t = 10)),
              legend.text = element_text(size = 32, margin = margin(t = 10)))

    return(p)
}
