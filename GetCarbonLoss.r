GetCarbonLoss = function(pixels, t0, area_ha, area_adj_ratio, cdens, adjustArea = T, pair) {
  #make an adjustment of the area represented based on the proportion of unmatched pixels
  if(adjustArea) area_ha = area_ha * area_adj_ratio

  #calculate yearly time series of the proportion and area of each land class
  area_series = GetAreaSeries(pixels, t0 = t0, area_ha = area_ha)

  #calculate annual time series of average carbon density of the project (MgC/ha)
  carbon_series = area_series %>%
    merge(., cdens, by.x = "class", by.y = "land.use.class", all.x = T) %>%
    mutate(carbon_density = class_prop * carbon.density) %>% #MgC/ha
    group_by(year) %>%
    summarise(carbon_density = sum(carbon_density, na.rm = T)) %>%
    ungroup() %>%
    mutate(pair = pair)

  #calculate annual fluxes of carbon density
  carbon_loss = data.frame(year = carbon_series$year[-1],
                           c_loss = -diff(carbon_series$carbon_density),
                           pair = pair)

  return(list(carbon_series = carbon_series, carbon_loss = carbon_loss))
}