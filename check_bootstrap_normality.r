#need the variable "projects"

norm_pval = matrix(NA, nrow = length(projects), ncol = 6)

BootOut = function(type, in_df, boot_n = 1000) {
  boot_out = boot::boot(data = in_df,
                        statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                        R = boot_n)
  return(boot_out$t)
}

for(i in seq_along(projects)) {
    area_i = project_var$area_ha[i]
    observed = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T) %>%
      filter(started == T) %>%
      mutate(t_loss = t_loss / area_i,
             c_loss = c_loss / area_i,
             additionality = additionality / area_i)
    baseline_best = read.csv(paste0(out_path, "_baseline_best_", projects[i], ".csv"), header = T)
    baseline_loose = read.csv(paste0(out_path, "_baseline_loose_", projects[i], ".csv"), header = T)
    baseline_offset = read.csv(paste0(out_path, "_baseline_offset_", projects[i], ".csv"), header = T)

    bootsum_out = BootOut(type = "cf_c_loss", in_df = dplyr::select(observed, c_loss))
    norm_pval[i, 1] = shapiro.test(bootsum_out)$p.value

    bootsum_out = BootOut(type = "p_c_loss", in_df = dplyr::select(observed, t_loss))
    norm_pval[i, 2] = shapiro.test(bootsum_out)$p.value

    bootsum_out = BootOut(type = "additionality", in_df = dplyr::select(observed, additionality))
    norm_pval[i, 3] = shapiro.test(bootsum_out)$p.value

    bootsum_out = BootOut(type = "best", in_df = dplyr::select(baseline_best, c_loss))
    norm_pval[i, 4] = shapiro.test(bootsum_out)$p.value

    bootsum_out = BootOut(type = "loose", in_df = dplyr::select(baseline_loose, c_loss))
    norm_pval[i, 5] = shapiro.test(bootsum_out)$p.value

    if(nrow(baseline_offset) > 0) {
        bootsum_out = BootOut(type = "offset", in_df = dplyr::select(baseline_offset, c_loss))
        norm_pval[i, 6] = shapiro.test(bootsum_out)$p.value
    } else {
      norm_pval[i, 6] = NA
    }
    cat("Project", i, "\n")
}
