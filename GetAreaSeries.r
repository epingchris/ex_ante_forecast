#function to generate yearly time series of the proportion and area of each land class based on the sample
GetAreaSeries = function(pixels, t0, area_ha) {
  area_series = pixels %>%
    #add up the number of pixels per year and land class within each treatment
    summarise(across(starts_with("JRC"), list(`1` = ~ sum(.x == 1),
                                              `2` = ~ sum(.x == 2),
                                              `3` = ~ sum(.x == 3),
                                              `4` = ~ sum(.x == 4),
                                              `5` = ~ sum(.x == 5),
                                              `6` = ~ sum(.x == 6)))) %>%
    pivot_longer(starts_with("JRC"),
                 names_to = c("year", "class"),
                 names_prefix = "JRC",
                 names_sep = "_",
                 values_to = "count") %>%
    mutate(year = as.numeric(year)) %>%
    group_by(year) %>%
    mutate(class_prop = count / nrow(pixels), #proportion of pixels in each class
           class_area = class_prop * area_ha) #area (ha) of the project this sample represents

  return(area_series)
}