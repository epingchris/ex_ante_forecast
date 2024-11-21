CalcError = function(dat) {
    dat_wide = dat %>%
        dplyr::select(type, mean, project) %>%
        pivot_wider(names_from = "type", values_from = "mean", id_expand = T)

    dat_long = dat_wide %>%
        pivot_longer(best:lagged, names_to = "baseline_type", values_to = "baseline")

    dat_long_best = subset(dat_long, baseline_type == "best") %>% na.omit()
    dat_long_loose = subset(dat_long, baseline_type == "loose") %>% na.omit()
    dat_long_lagged = subset(dat_long, baseline_type == "lagged") %>% na.omit()

    error_df = matrix(c(rmse(dat_long_best$cf_c_loss, dat_long_best$baseline),
                        rmse(dat_long_loose$cf_c_loss, dat_long_loose$baseline),
                        rmse(dat_long_lagged$cf_c_loss, dat_long_lagged$baseline),
                        mae(dat_long_best$cf_c_loss, dat_long_best$baseline),
                        mae(dat_long_loose$cf_c_loss, dat_long_loose$baseline),
                        mae(dat_long_lagged$cf_c_loss, dat_long_lagged$baseline)),
                      nrow = 3, ncol = 2) %>%
        as.data.frame()
    colnames(error_df) = c("RMSE", "MAE")
    rownames(error_df) = c("Best-matched", "Loosely-matched", "Time-lagged")
    return(error_df)
}