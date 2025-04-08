#create project_var_basic.csv containing only basic variables
#create condensed additionality distribution data files
#to send to Ofir for his project
rm(list = ls())

library(tidyverse) #ggplot2, dplyr, stringr, plotPlacebo/plotBaseline.r: tibble to store labels with bquote()
library(magrittr) #pipe operators

path = "/maps/epr26/ex_ante_forecast_out/"
project_var = read.csv(paste0(path, "out_ongoing_project_var.csv"), header = T)
write.table(project_var %>% dplyr::select(project, t0, country, area_ha),
            paste0(path, "ofir_project_var_basic.csv"), sep = ",", row.names = F)

projects = project_var$ID

additionality_distribution = lapply(seq_along(projects), function(i) {
    add = read.csv(paste0(path, "out_ongoing_additionality_", projects[i], ".csv"), header = T) %>%
        filter(started) %>%
        dplyr::select(year, additionality, pair) %>%
        mutate(project = projects[i])
}) %>%
    list_rbind()
write.csv(additionality_distribution, paste0(path, "ofir_additionality_distribution.csv"), row.names = F)