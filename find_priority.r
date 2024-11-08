rm(list = ls())

library(tidyverse)
library(magrittr)

#Load basic info about projects from Tom's directory: for country and t0 input (E-Ping's current workflow)
proj_meta = read.csv(paste0("/home/tws36/4c_evaluations/data/project_metadata/proj_meta.csv"))

#block of code to determine which projects to run next
exclude_id = c("612", "1340", "2363", "3057", "2974", "2975", "2976")
#612, 1340 done under 612a, 1340a
#2363 withdrawn
#3057, 2974-2976 are mangrove

done_id = list.files("/maps/epr26/tmf_pipe_out/") %>% #full = T and basename() negates one another
    str_subset("\\.", negate = T) %>%
    str_subset("\\_", negate = T) %>%
    str_subset("ac", negate = T) #already done

#projects with shapefiles, aren't excluded or already done
candidates_id = list.files("/maps/epr26/tmf-data/projects/") %>%
  gsub("\\.geojson", "", .) %>%
  str_subset("ac", negate = T) %>%
  setdiff(exclude_id) %>%
  setdiff(done_id)

#countries with signficant amount of JRC-TMF pixels
country_tmf = unique(read.csv("country_tmf.csv")$SOVEREIGNT)
exclude_country = c("China", "Uruguay", "United States", "Argentina", "Pakistan", "Chile", "Zimbabwe", "Australia", "Mozambique", "India", "Zambia", "Benin")
#China, Uruguay, US, Argentina and Pakistan explicitly excluded
#Pakistan, Chile and Zimbabwe de facto excluded
#Australia only in northeast; Mozambique only in north; India only in southwest; Zambia only in northwest; Benin only in south
proj_candidate = proj_meta %>%
  filter(ID %in% candidates_id) %>%
  filter(COUNTRY %in% country_tmf) %>%
  filter(COUNTRY %in% exclude_country == F) %>%
  filter(t0 >= 2009 & t0 < 2018) # projects with t0 in [2009, 2017], in order to have available luc columns

#table(subset(proj_meta, ID %in% projects)$t0)
proj_priority = filter(proj_candidate, t0 %in% c(2012, 2014, 2015, 2016, 2017)) #priortise running projects in these years