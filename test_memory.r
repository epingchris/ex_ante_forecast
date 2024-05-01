rm(list = ls())
library(tidyverse)
library(magrittr)
library(pryr)

proj = "562"

acd = read.csv(paste0('/maps/pf341/results/live-pipeline/', proj, '-carbon-density.csv'))
acd_change = ifelse(sum(is.na(acd$carbon.density[c(1, 3)])) == 0,
                    acd$carbon.density[1] - acd$carbon.density[3], NA)

defor = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "matches.parquet"))
dim(defor) #[1] 68302317       41
print(object.size(defor), units = "MB", standard = "SI") #22130 MB
object_size(defor) #273.22 MB

defor = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "matches.parquet")) %>%
    dplyr::select(cpc0_u:cpc10_d)
dim(defor) #[1] 68302317        6
print(object.size(defor), units = "MB", standard = "SI") #3278.5 MB
object_size(defor) #2.52 kB

defor = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "matches.parquet")) %>%
    dplyr::select(cpc0_u:cpc10_d) %>%
    mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10)
dim(defor) #[1] 68302317        9
print(object.size(defor), units = "MB", standard = "SI") #4917.8 MB
object_size(defor) #1.64 GB

defor = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "matches.parquet")) %>%
    dplyr::select(cpc0_u:cpc10_d) %>%
    mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
    mutate(acd_defor_5_0 = defor_5_0 * acd_change, acd_defor_10_5 = defor_10_5 * acd_change, acd_defor_10_0 = defor_10_0 * acd_change)
dim(defor) #[1] 68302317       12
print(object.size(defor), units = "MB", standard = "SI") #6557 MB
object_size(defor) #3.28 GB
