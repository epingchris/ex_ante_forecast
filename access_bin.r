rm(list = ls())

library(tidyverse)
library(magrittr)
library(stars)
library(arrow)
library(parallel)

#load functions
source("./functions.r") #cpc_rename, tmfemi_reformat

#load the list of projects to be run
projects = c("562", "674", "944", "1047", "1112", "1113", "1201", "1311", "1396", "1399")

#load metadata
proj_meta = read.csv("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv") #for t0
proj_var = readRDS("/maps/epr26/tmf_pipe_out/project_var.rds") #for area
proj_estimate = readRDS("/maps/epr26/tmf_pipe_out/project_estimates.rds") #for during-project additionality

access_bin_list = lapply(seq_along(projects), function(i) {
  c = Sys.time()
  proj = projects[i]
  t0 = proj_meta %>% filter(ID == proj) %>% pull(t0)
  area_ha = proj_var %>% filter(project == proj) %>% pull(area_ha)

  acd = read.csv(paste0('/maps/pf341/results/live-pipeline/', proj, '-carbon-density.csv'))
  acd_change = ifelse(sum(is.na(acd$carbon.density[c(1, 3)])) == 0,
                      acd$carbon.density[1] - acd$carbon.density[3], NA)

  #obtain vicinity baseline carbon loss rate (undisturbed -> deforested LUC change)
  defor = read_parquet(paste0("/maps/epr26/tmf_pipe_out/", proj, "/", proj, "matches.parquet")) %>%
    dplyr::select(access, cpc0_u:cpc10_d) %>%
    mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10) %>%
    mutate(acd_defor_5_0 = defor_5_0 * acd_change, acd_defor_10_5 = defor_10_5 * acd_change, acd_defor_10_0 = defor_10_0 * acd_change) %>%
    mutate(access_bins = cut(access, breaks = seq(0, 2400, by = 60)))

   #get the proportion of pixels with zero deforestation in each accessibility bin
  access_bin = data.frame(interval = 0:39,
                          dens = hist(defor$access / 60, breaks = 0:40)$density,
                          prop_zero_defor = tapply(defor$acd_defor_5_0, defor$access_bins,
                                                   function(x) length(which(x == 0)) / length(x)),
                          proj = proj)

  max_access = max(defor$access, na.rm = T)
  thres_vec = seq(0.5, 0.9, by = 0.05)
  for(i in seq_along(thres_vec)) {
    defor_filt = defor %>% filter(access < max_access * thres_vec[i])
  }

  d = Sys.time()
  cat(proj, ":", d - c, "\n")
  return(access_bin)
# return(list(counts = counts, prop_zero = prop_zero))
})

access_bin_df = do.call(rbind, access_bin_list) %>%
  mutate(proj = factor(proj, levels = projects))

ggplot(data = access_bin_df, aes(x = interval)) +
  geom_line(aes(y = dens, color = "prop")) +
  geom_line(aes(y = prop_zero_defor, color = "prop_zero")) + 
  scale_color_manual(values = c("black", "red"),
                     labels = c("Proportion of total pixels",
                                "Proportion of pixels with zero baseline C loss in each hourly bin")) +
  facet_wrap(vars(proj), ncol = 3) +
  xlab("Inaccessibility (hours)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
ggsave("access_bin.png", width = 2500, height = 3000, units = "px")
