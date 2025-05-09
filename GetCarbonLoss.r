#function to generate time series of annual change in project-level average carbon density
GetCarbonLoss = function(pixels, t0, area_ha, area_adj_ratio, cdens, adjustArea = T, pair) {
  #make an adjustment of the area represented based on the proportion of unmatched pixels
  if(adjustArea) area_ha = area_ha * area_adj_ratio

  #calculate yearly time series of the proportion and area of each land class
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

  #calculate annual time series of average carbon density of the project (MgC/ha)
  carbon_series = area_series %>%
    merge(., cdens, by.x = "class", by.y = "land.use.class", all.x = T) %>%
    mutate(carbon_density = class_prop * carbon.density) %>% #MgC/ha
    group_by(year) %>%
    summarise(carbon_density = sum(carbon_density, na.rm = T)) %>%
    ungroup() %>%
    mutate(pair = pair) %>%
    sf::st_drop_geometry()

  return(carbon_series)
}