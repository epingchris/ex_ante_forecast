# set out a table of parameters you want to change
# read the config
# write the config with the parameters
# Run the script

rm(list = ls())

# install.packages(c('arrow','configr', 'tidyverse', 'magrittr', 'sf', 'magrittr', 'MatchIt',
#                    'rnaturalearthdata', 'configr', 'terra', 'pbapply', 'cleangeo', 'doParallel',
#                    'foreach', 'readr', 'lwgeom', 'rnaturalearth', 'mclust', 'ggpubr', 'EnvStat'), depends = TRUE)

library(tidyverse)
library(configr)
library(magrittr)
library(readr)
library(arrow)
library(sf)
library(MatchIt)
library(rnaturalearthdata)
library(terra)
library(pbapply)
library(cleangeo)
library(doParallel)
library(foreach)
library(readr)
library(lwgeom)
library(countrycode)
library(mclust)
library(ggpubr)

in_path = '/home/tws36/4c_evaluations/'

# The list of projects to be run in this evaluation:
proj_meta <- read.csv(paste0(in_path, 'data/project_metadata/proj_meta.csv'))
#proj_to_eval <- read.table('./data/project_metadata/proj_to_eval.txt') %>% unlist() %>% as.numeric()
#projects_agb <- read_csv('./data/GEDI/project_agb.csv')

# For testing on local
# data_suffix<- '230313'
# proj_id <- 1201
# site <- 'Gola'

cpc_rename <- function(x, t0) {
  x %<>%
    as.data.frame() %>%
    dplyr::select(starts_with('cpc'))
  mycolnames <- colnames(x)
  years <- mycolnames %>%
    stringr::str_extract('[:digit:]+') %>% as.numeric()
  suffix <- mycolnames %>%
    stringr::str_extract('_u|_d') %>%
    stringr::str_replace('u', '1') %>%
    stringr::str_replace('d', '3')
  newnames <- paste0('JRC', t0-years, suffix)
  colnames(x) <- newnames
  return(x)
}

tmfemi_reformat <- function(df, t0) {
  df %<>%
    st_as_sf(coords = c("lng", "lat")) %>%
    dplyr::rename(accessibility = access) %>%
    dplyr::rename_with(~ gsub("luc_", "JRC", .x, fixed = TRUE))

  other <- df %>% dplyr::select(-starts_with('cpc'))

  cpcs <- df %>%
    dplyr::select(starts_with('cpc')) %>%
    cpc_rename(t0 = t0)

  df <- cbind(df, cpcs)
  # JRC2002_1 = cpc10_u,
  # JRC2007_1 = cpc5_u,
  # JRC2012_1 = cpc0_u,
  # JRC2002_3 = cpc10_d,
  # JRC2007_3 = cpc5_d,
  # JRC2012_3 = cpc0_d) %>%

  if(any(colnames(df) %>% stringr::str_detect('ecoregion')))
    df %<>% dplyr::rename(biome = ecoregion)
  return(df)
}

config <- read.config(paste0(in_path, 'config/fixed_config_sherwood.ini'))
config$USERPARAMS$data_path <- '/maps/pf341/tom'
#write.config(config, './config/fixed_config_tmp.ini') #error: permission denied

# Load user-defined functions that Tom wrote
sapply(list.files(paste0(in_path, 'R'), full.names = TRUE, pattern = '.R$'), source)

# Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = FALSE)

# theme_set(theme_minimal())
# theme_replace(panel.grid.minor = element_line(colour = "red"))

source(paste0(in_path, 'R/scripts/setup_and_reusable/load_config.R'))
# source('./R/scripts/0.2_load_project_details.R')

# match_years<-c(0, -5, -10)

# # Setup AGB values:
# # Gola PDD
# co2<-c(606, 237, 127, 0, 0, 21)
# # Setup AGB values:
# class_agb <- mk_class_agb(agb = NULL, acd = NULL, co2 = co2)
# class_agb <-adjust_stock(class_agb, bgb_adj = 0.2, deadwood_adj = 0.11)

# # GEDI
# acd<-c(130.68634598999023, 24.467604306411744, 12.470903750324249, 17.916133327388764, 0, 10.001206178855895)
# class_agb <- mk_class_agb(agb = NULL, acd = acd, co2 = NULL)
# # class_agb <-adjust_stock(class_agb, bgb_adj = 0.2, deadwood_adj = 0.11)

# project_area_ha<-702200042.54 / 10000

# # Patrick's candidates

# control<-read_parquet('../../../maps/pf341/tom-candidates/controls.parquet') # all potential matches
# treat<-read_parquet('../../../maps/pf341/tom-candidates/treatment.parquet') # treatment subset


class_prefix <- 'JRC'
match_years <- c(0, -5, -10)
match_classes <- c(1, 3)

# Find projects with a non-empty carbon_density.csv
acd_dir = '/maps/pf341/results/live-pipeline/'
acd_paths <- list.files(acd_dir, full = TRUE) %>%
  str_subset('carbon-density') %>%
  sapply(., function(x) ifelse(nrow(read.csv(x)) == 0, NA, x)) %>%
  na.omit() %>%
  as.vector()

acd_proj_id <- basename(acd_paths) %>%
  str_replace('-carbon-density', '') %>%
  str_replace('.csv', '')

# Find projects with matched pair data
# Select those with carbon density data, unselect problematic projects (1566, 1067, 958, 1133)
exclude_id = c(1566, 1067, 958, 1133)
pair_dir <- '/maps/pf341/results/2024-january-pipeline/'
# '/maps/pf341/tom-add-paper'
project_paths <- list.files(pair_dir, full = TRUE) %>%
  str_subset('pairs') %>%
  sapply(., function(x) {
    proj_id <- basename(x) %>% str_replace('_pairs', '')
    ifelse(proj_id %in% acd_proj_id & proj_id %in% exclude_id == F, x, NA)
  }) %>%
  na.omit() %>%
  as.vector()


# Obtain annual carbon loss and additionality values ----
#i <- 1 #used to test just one project
#proj_list <- lapply(1:10, function(i) { #used to loop just ten projects
#proj_list <- lapply(seq_along(project_paths), function(i) { #used on Windows
proj_list <- mclapply(seq_along(project_paths), mc.cores = 30, function(i) {
  myproject_path <- project_paths[i]
  proj_id <- basename(myproject_path) %>% str_replace('_pairs', '')

  # Extract project start date:
  myproj <- proj_meta %>%
    filter(ID == proj_id)

  t0 <- myproj$t0
  eval_end <- myproj$eval_end
  eval_end <- ifelse(is.na(eval_end), 2021, eval_end)

  site <- paste('VCS', proj_id, sep = '_')
  if(!str_detect(supplier_path, '.shp')) {
    aoi_path <- file.path(supplier_path, site, 'GIS', 'aoi.shp')
    if(!file.exists(aoi_path)) aoi_path <- file.path(supplier_path, site, paste(site, '.shp', sep = ''))
  }
  aoi_project <- read_sf(aoi_path) %>% # NEED TO CONSISTENTLY SET AOI NAME
    st_make_valid() %>%
    st_union()

  # Transform the projection:
  aoi_project <- aoi_project %>% st_transform(4326)

  # Find the area of the region:
  project_area_ha <- st_area_ha(aoi_project)

  # Find paths to match and unmatached points:
  pair_paths <- list.files(myproject_path, full = TRUE)
  matchless_ind <- pair_paths %>% str_detect('matchless')
  matchless_paths <- pair_paths[matchless_ind]
  matched_paths <- pair_paths[!matchless_ind]

  # Read and analyse pairs:
  #j <- 4 #to test just one pair
  project_estimates <- lapply(seq_along(matched_paths), function(j) {
    pairs <- read_parquet(matched_paths[j])
    unmatched_pairs <- read_parquet(matchless_paths[j])

    control <- pairs %>%
      select(starts_with('s_')) %>%
      rename_with(~str_replace(.x, 's_', '')) %>%
      mutate(treatment = 'control') %>%
      tmfemi_reformat(t0 = t0)

    treat <- pairs %>%
      select(starts_with('k_')) %>%
      rename_with(~str_replace(.x, 'k_', '')) %>%
      mutate(treatment = 'treatment') %>%
      tmfemi_reformat(t0 = t0)

    exp_n_pairs <- nrow(treat) + nrow(unmatched_pairs)

    pts_matched <- rbind(treat, control)

    # m.out<-assess_balance(pts_matched, class_prefix = class_prefix, t0 = t0,
    #                       match_years = match_years, match_classes = match_classes)
    # summary(m.out, standardize = TRUE)

    control_series <- simulate_area_series(pts_matched,
                                           class_prefix, t0 = t0, match_years, match_classes,
                                           exp_n_pairs, project_area_ha,
                                           verbose = FALSE)
    # control_series_post_pairs$series %>%
    # left_join(class_agb %>% select(class, agb), by = 'class') %>%
    # mutate(class_agb = class_area * agb, # the above ground biomass
    #         class_co2e = class_agb * cf_c * cf_co2e
    #         ) %>%
    #         additionality_total_series()

    # Extract ACD per LUC
    acd = read.csv(paste0(acd_dir, proj_id, '-carbon-density.csv'))

    y <- control_series$series %>%
      #left_join(class_agb %>% select(class, agb), by = 'class') %>%
      #mutate(class_agb = class_area * agb,
      #       class_co2e = class_agb * cf_c * cf_co2e) %>% #the above ground biomass
      #filter(class == 1) %>% #I do not want to just use class 1
      #filter(year %in% c(t0, eval_end)) #I do not want to filter to only t0 and eval_end
      merge(., acd, by.x = "class", by.y = "land.use.class", all.x = T) %>%
      mutate(carbon_content = class_area * carbon.density) %>%
      group_by(treatment, year) %>%
      summarise(carbon_content = sum(carbon_content, na.rm = T)) %>%
      ungroup()

    year <- y %>% filter(treatment == 'control') %>% pull(year)
    yc <- y %>% filter(treatment == 'control') %>% pull(carbon_content)
    yt <- y %>% filter(treatment == 'treatment')  %>% pull(carbon_content)

    out_df <- data.frame(pair = j, year = year[-1], c_loss = -diff(yc), t_loss = -diff(yt)) %>%
      mutate(additionality = c_loss - t_loss)

    #out_df <- data.frame(
    #  ID = proj_id,
    #  Start = t0,
    #  End = eval_end,
    #  control_undisturbed = yc[2],
    #  project_undisturbed = yt[2],
    #  avoided_disturbance_ha = undisturbed_additionality,
    #  avoided_disturbance_ha_yr = undisturbed_additionality / (eval_end - t0)
    #)

    return(out_df)
  })


  project_estimates <- do.call(rbind, project_estimates) %>%
    mutate(started = ifelse(year > t0, T, F))

  return(project_estimates)
})

proj_id_list = basename(project_paths) %>% str_replace('_pairs', '')
names(proj_list) = proj_id_list

out_path = paste0('/maps/epr26/tmf_pipe_out/')

saveRDS(proj_list, file.path(paste0(out_path, 'project_summaries.rds')))
