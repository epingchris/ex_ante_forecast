rm(list = ls())

install.packages(c("arrow","configr", "tidyverse", "magrittr", "sf", "magrittr", "MatchIt",
                   "rnaturalearthdata", "configr", "terra", "pbapply", "cleangeo", "doParallel",
                   "foreach", "readr", "lwgeom", "rnaturalearth", "stars"), depends = TRUE)

install.packages("languageserver")
library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators
library(arrow) #arrow::read_parquet

loose = read_parquet("/maps/epr26/tmf_pipe_out_luc_t/1201/matches.parquet")

acd_i = read.csv("/maps/epr26/tmf_pipe_out_luc_t/1201/carbon-density.csv")

block = read_parquet("/maps/epr26/tmf_pipe_out_luc_t/1201/block_baseline.parquet") %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
    as.data.frame() %>%
    filter(elevation >= 0 & access >= 0 & slope >= 0) %>%
    dplyr::select(-starts_with("luc_")) %>%
    dplyr::select(-starts_with("cpc"))
if(nrow(block) > 250000) block = block[sample(nrow(block), 250000), ]

baseline_i = block %>%
    mutate(acd10 = acd_i$carbon.density[match(luc10, acd_i$land.use.class)],
           acd0 = acd_i$carbon.density[match(luc0, acd_i$land.use.class)],
           c_loss = (acd10 - acd0) / 10)


library(terra)
block_tif = terra::raster("/maps/epr26/tmf_pipe_out/1047/block_baseline.tif")