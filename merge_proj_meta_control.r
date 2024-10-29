library(tidyverse)
library(magrittr)

asn_info = read.csv("asian_tropics_controls.csv")
sa_info = read.csv("neotropics_controls.csv")
af_info = read.csv("afrotropics_controls.csv")
proj_info = do.call(bind_rows, list(asn_info, sa_info, af_info)) %>%
  mutate(ID = proj_name, COUNTRY = name, t0 = 2011) %>%
  filter(ID != "") %>%
  dplyr::select(ID, COUNTRY, t0)
write.csv(proj_info, "proj_meta_control.csv")