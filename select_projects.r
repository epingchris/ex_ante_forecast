acd_dir = "/maps/pf341/results/live-pipeline/"
acd_id = list.files(acd_dir, full = TRUE) %>%
  str_subset("carbon-density") %>%
  basename() %>%
  str_replace("-carbon-density.csv", "")

project_dir = "/maps/pf341/results/2024-january-pipeline"
project_id_full = list.files(project_dir, full = TRUE) %>%
  str_subset("pairs") %>%
  basename() %>%
  str_replace("_pairs", "")
project_id = project_id_full[project_id_full %in% acd_id]

acd_list = lapply(project_id, function(x) {
  acd = read.csv(paste0(acd_dir, x, "-carbon-density.csv"))
})
names(acd_list) = project_id
acd_list[which(sapply(acd_list, nrow) < 4)]
acd_exclude = project_id[which(sapply(acd_list, nrow) < 4)]

shp_id = list.files("/maps/epr26/tmf-data/projects", full = TRUE) %>%
  basename() %>%
  str_replace(".geojson", "")
shp_id = shp_id[!shp_id %in% c("0000", "9999")] #test polygons
shp_id = shp_id[!shp_id %in% c("1566", "1067", "958", "1133")] #problematic results
shp_id = shp_id[!shp_id %in% acd_exclude] #incomplete ACD data
shp_id = shp_id[!shp_id %in% c("1359", "902", "3347", "3335", "3114")] #withdrawn and on hold and under development

proj_meta = read.csv("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv") #for t0
proj_meta %>% filter(ID %in% as.numeric(shp_id))
