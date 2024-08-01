#Code to find the subset of projects that are waiting to be run

library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators

#projects that are already done
project_status = read.csv("/maps/epr26/tmf_pipe_out/project_status.csv", header = T)
done_id = subset(project_status, done)$project

#projects without complete ACD: from past and present results
acd_dir = "/maps/pf341/results/live-pipeline/"
acd_id = list.files(acd_dir) %>%
  str_subset("carbon-density") %>%
  str_replace("-carbon-density.csv", "")
acd_list = lapply(acd_id, function(x) {
  acd = read.csv(paste0(acd_dir, x, "-carbon-density.csv"))
})
incomplete_acd_id = acd_id[which(sapply(acd_list, nrow) < 4)] %>%
  union(subset(project_status, !full_acd)$project %>% str_replace("a", ""))

#select from available shapefiles
shp_id = list.files("/maps/epr26/tmf-data/projects") %>%
  str_replace(".geojson", "") %>%
  str_subset("ac", negate = T) %>%
  str_subset("as", negate = T) %>%
  str_replace("a", "")
shp_id = shp_id[!shp_id %in% c("0000", "9999")] #test polygons
shp_id = shp_id[!shp_id %in% union(incomplete_acd_id, done_id)] #incomplete ACD data, done
shp_id = shp_id[!shp_id %in% c("1359", "902", "3347", "3335", "3114")] #withdrawn and on hold and under development
#shp_id = shp_id[!shp_id %in% c("1566", "1067", "958", "1133")] #problematic results: why did I say that?

#find t0
proj_candidate = read.csv("proj_meta.csv") %>%
  filter(ID %in% shp_id) %>%
  filter(t0 >= 2007 & t0 <= 2017)

proj_done = read.csv("proj_meta.csv") %>%
  filter(ID %in% done_id)
