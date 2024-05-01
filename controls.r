rm(list = ls())

library(tidyverse)
library(sf)
library(magrittr)

#stored in "tmf-data-grid/0000"

controls = st_read("/maps/epr26/tmf-data-grid/0000/controls.geojson")

controls$id = stringr::str_sub(controls$id, -2)
controls$id_2 = seq_len(nrow(controls))

st_write(controls, "/maps/epr26/tmf-data-grid/0000/controls_new.geojson", driver = "GeoJSON")

#st_relate(controls, pattern = "****0****") find overlapping controls

for (i in seq_len(nrow(controls))) {
  st_write(controls[i, ], paste0("/maps/epr26/tmf-data-grid/0000/0000_", i, ".geojson"), driver = "GeoJSON")
}