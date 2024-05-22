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
SaveMultiPagePlot = function(plots, suffix, width = 5000, height = 3000) {
  if(class(plots)[1] == "list") {
    lapply(seq_along(plots), function(i) {
      ggsave(paste0(out_path, "_", suffix, "_", i, ".png"), plots[[i]],
             width = width, height = height, units = "px", bg = "white")
    })
  } else {
    ggsave(paste0(out_path, "_", suffix, ".png"), plots,
           width = width, height = height, units = "px", bg = "white")
  }
}

cat("These functions require ggplot2, dplyr and stars packages")

#function to turn empty object as NA
fillNA = function(x) if (length(x) == 0) return(NA) else return(x)