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

  return(df)
}

simulate_area_series<-function(pts_matched,
                               class_prefix, t0, match_years, match_classes,
                               exp_n, area, verbose = TRUE){
  match_assess<-summary(assess_balance(pts_matched, class_prefix, t0, match_years, match_classes))
  balance_test<-all(abs(match_assess$sum.matched[, 'Std. Mean Diff.'])<=0.2)

  if(verbose){
    print(match_assess$sum.matched)
    print(balance_test)
  }
  # if(balance_test){
    # Make an adjustment of the area represented based on the number of points not matched
  match_sample_adj<-(nrow(pts_matched)/2)/exp_n
  area <- area * match_sample_adj
  # print(area)

  area_series<-make_area_series(pts_matched %>% na.omit(), area_ha = area, class_prefix = class_prefix)
  out<-list(area_series, match_assess$sum.matched)
  names(out) <- c('series', 'balance')

  return(out)
  # }
  # else
  #   return(NULL)
}

make_area_series<-function(pts, area_ha, class_prefix){
  pts<-pts[str_detect(colnames(pts), paste(class_prefix, '[:digit:]{4}$|treatment', sep = ''))]

  pts %>%
    # select(-matches(paste('_[1-6]$', sep = ''))) %>% # exclude proportional cover values
    as.data.frame() %>%
    # select(-id, -X2010_agb, -(accessibility:transition), treatment) %>% 
    group_by(treatment) %>%
    # Add up the number of points in each habitat class within each treatment:
    summarise(across(where(is.numeric), list(`1` = ~ sum(.x ==1),
                                             `2` = ~ sum(.x ==2),
                                             `3` = ~ sum(.x ==3),
                                             `4` = ~ sum(.x ==4),
                                             `5` = ~ sum(.x ==5),
                                             `6` = ~ sum(.x ==6) # we do not have an AGB value for class 4 because
                                             # it is so rare
    )
    )) %>%
    # Convert to wide format:
    pivot_longer(cols = contains(class_prefix),
                 names_to = c('year', 'class'),
                 names_prefix = class_prefix,
                 names_sep = '_',
                 values_to = 'n') %>%
    # Convert from character to numeric
    mutate(across(year, as.numeric)) %>%
    # Add on the agb data for each habitat class by joining with the hab_class_agb object
    # left_join(class_agb %>% select(class, agb), by = 'class') %>%
    # Group by treatment and year and then generate stats for these groups:
    group_by(treatment, year) %>%
    mutate(n_total = sum(n),                   # total number of points in sample
           class_prop = n / n_total,           # the proportion of points in each class
           class_area = class_prop * area_ha,  # the area of the project this represents
           # class_agb = class_area * agb,       # the above ground biomass
           # class_co2e = class_agb * cf_c* cf_co2e
           )
}

assess_balance<-function(pts, class_prefix, t0, match_years, match_classes){
  fmla<-make_match_formula(prefix = class_prefix,
                           t0 = t0,
                           match_years = match_years,
                           match_classes = match_classes,
                           suffix = '')

  matchit(
    fmla,
    method = 'nearest',
    distance = 'mahalanobis',
    ratio = 1,
    order = 'smallest',
    replace = FALSE,
    discard = 'none',
    data = pts %>%
      as.data.frame %>%
      mutate(treatment = ifelse(treatment == "treatment", 1, 0))
  )
}

make_match_formula<-function(prefix,
                             t0,
                             match_years,
                             match_classes,
                             suffix){

  # generate the match variables:
  match_var_grid<-expand.grid(prefix = prefix,
                              years = t0 + match_years, # the years to match on
                              # match_years should be zero 
                              # or negative
                              classes = match_classes,
                              suffix = "")

  match_vars<-apply(match_var_grid, 1, function(x){
    paste(x['prefix'], x['years'], '_', x['classes'], x['suffix'], sep='')
    # %>% str_replace('(^_|_$)', '')
  }) %>%
    c("accessibility", "elevation", "slope")

  # the match formula
  fmla <- as.formula(paste("treatment ~ ", paste(match_vars, collapse= "+")))

  return(fmla)
}

#functions for multi-page plots
# SaveMultiPagePlot = function(plotlist, suffix, n = 4, width = 4000, height = 4000) {
#   if(length(plotlist) > n ^ 2) {
#     plot_all = plotlist %>%
#       ggpubr::ggarrange(plotlist = ., ncol = n, nrow = n, common.legend = T, legend = "bottom")
#     lapply(seq_along(plot_all), function(i) {
#       ggsave(paste0(out_path, "_", suffix, "_", i, ".png"), plot_all[[i]],
#              width = width, height = height, units = "px", bg = "white")
#     })
#   } else {
#     plot_nrow = ceiling(length(plotlist) / n)
#     plot_all = plotlist %>%
#       ggpubr::ggarrange(plotlist = ., ncol = n, nrow = plot_nrow, common.legend = T, legend = "bottom")
#     height_adjusted = height * plot_nrow / n
#     ggsave(paste0(out_path, "_", suffix, ".png"), plot_all,
#            width = width, height = height_adjusted, units = "px", bg = "white")
#   }
# }

cat("These functions require tidyverse, MatchIt, and sf packages\n")