BootOut = function(in_df, column, from, to, boot_n = 1000) {
  n_intrvs = to - from + 1

  boot_summ = data.frame(year = from:to,
                         mean = rep(NA, n_intrvs),
                         ci_lower = rep(NA, n_intrvs),
                         ci_upper = rep(NA, n_intrvs))
  for(i in seq_len(n_intrvs)) {
    data_i = in_df %>%
      filter(year == from + i - 1) %>%
      dplyr::select(all_of(column))

    boot_out = boot::boot(data = data_i,
                          statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #bootstrapped mean
                          R = boot_n)

    boot_mean = mean(boot_out$t)
    if(length(unique(boot_out$t)) > 1) {
      boot_ci = boot::boot.ci(boot.out = boot_out, type = "perc")$percent[4:5] #percentile intervals
    } else {
      cat("Warning: all values of t are equal for year", from + i - 1, "\n")
      boot_ci = c(boot_mean, boot_mean)
    }

    boot_summ$mean[i] = boot_mean
    boot_summ$ci_lower[i] = boot_ci[1]
    boot_summ$ci_upper[i] = boot_ci[2]
  }

  return(boot_summ)
}