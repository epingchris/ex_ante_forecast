#Compare between E-Ping and Jody's output
rm(list = ls())

library(tidyverse) #ggplot2, dplyr, stringr, plotPlacebo/plotBaseline.r: tibble to store labels with bquote()
library(magrittr) #pipe operators

path = paste0("/maps/epr26/ex_ante_forecast_out/") #where outputs are stored
observed_add_eping = read.csv(paste0(path, "out_ongoing_observed_add.csv"), header = T)
observed_add_jody = read.csv(paste0(path, "merged_additionality_by_jody.csv"), header = T)
year_max = observed_add_jody[nrow(observed_add_jody), ]$year

observed_add_jody_mean = observed_add_jody %>%
  filter(year == year_max) %>%
  rename_with(function(x) gsub("X", "", gsub("_additionality", "", x))) %>%
  pivot_longer(cols = !year, names_to = "ID", values_to = "cumul_add") %>%
  mutate(ID = as.numeric(ID)) %>%
  left_join(x = ., y = project_var %>% dplyr::select("ID", "t0", "area_ha"), by = "ID") %>%
  filter(!is.na(t0)) %>%
  mutate(additionality_jody = cumul_add / ((year - t0) * area_ha))

observed_add_compare = observed_add_eping %>%
  left_join(x = ., y = project_var %>% dplyr::select("ID", "t0", "area_ha"), by = join_by(project == ID)) %>%
  left_join(x = ., y = observed_add_jody_mean %>% dplyr::select("ID", "additionality_jody"), by = join_by(project == ID))
write.csv(observed_add_boot_compare, paste0(path, "additionality_compare.csv"), row.names = F)
