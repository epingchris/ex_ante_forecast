cpc_rename <- function(x, t0) {
  x %<>%
    as.data.frame() %>%
    dplyr::select(starts_with('cpc'))
  mycolnames <- colnames(x)
  years <- mycolnames %>%
    stringr::str_extract('[:digit:]+') %>% as.numeric()
  suffix <- mycolnames %>%
    stringr::str_extract('_u|_d') %>%
    stringr::str_replace('u', '1') %>%
    stringr::str_replace('d', '3')
  newnames <- paste0('JRC', t0-years, suffix)
  colnames(x) <- newnames
  return(x)
}

tmfemi_reformat <- function(df, t0) {
  df %<>%
    st_as_sf(coords = c("lng", "lat")) %>%
    dplyr::rename(accessibility = access) %>%
    dplyr::rename_with(~ gsub("luc_", "JRC", .x, fixed = TRUE))

  other <- df %>% dplyr::select(-starts_with('cpc'))

  cpcs <- df %>%
    dplyr::select(starts_with('cpc')) %>%
    cpc_rename(t0 = t0)

  df <- cbind(df, cpcs)
  # JRC2002_1 = cpc10_u,
  # JRC2007_1 = cpc5_u,
  # JRC2012_1 = cpc0_u,
  # JRC2002_3 = cpc10_d,
  # JRC2007_3 = cpc5_d,
  # JRC2012_3 = cpc0_d) %>%

  if(any(colnames(df) %>% stringr::str_detect('ecoregion')))
    df %<>% dplyr::rename(biome = ecoregion)
  return(df)
}

#functions for multi-page plots
SaveMultiPagePlot = function(plotlist, suffix, n = 4, width = 4000, height = 4000) {
  if(length(plotlist) > n ^ 2) {
    plot_all = plotlist %>%
      ggpubr::ggarrange(plotlist = ., ncol = n, nrow = n, common.legend = T, legend = "bottom")
    lapply(seq_along(plot_all), function(i) {
      ggsave(paste0(out_path, "_", suffix, "_", i, ".png"), plot_all[[i]],
             width = width, height = height, units = "px", bg = "white")
    })
  } else {
    plot_nrow = ceiling(length(plotlist) / n)
    plot_all = plotlist %>%
      ggpubr::ggarrange(plotlist = ., ncol = n, nrow = plot_nrow, common.legend = T, legend = "bottom")
    height_adjusted = height * plot_nrow / n
    ggsave(paste0(out_path, "_", suffix, ".png"), plot_all,
           width = width, height = height_adjusted, units = "px", bg = "white")
  }
}

cat("These functions require ggplot2, dplyr and stars packages")

#function to turn empty object as NA
fillNA = function(x) if (length(x) == 0) return(NA) else return(x)

#function to find LUC transitions from undisturbed to disturbed
findLUCC = function(x) sum(x %in% c("1_2", "1_3", "1_4")) / length(x)
