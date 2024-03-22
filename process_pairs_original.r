# set out a table of parameters you want to change
# read the config
# write the config with the parameters
# Run the script

rm(list= ls())

# install.packages(c('arrow','configr', 'tidyverse', 'magrittr', 'sf', 'magrittr', 'MatchIt', 'rnaturalearthdata', 'configr', 'terra', 'pbapply', 'cleangeo', 'doParallel', 'foreach', 'readr', 'lwgeom', 'rnaturalearth'), depends = TRUE)

library(tidyverse)
library(configr)
library(magrittr)
library(readr)
library(arrow)
library(sf)
library(MatchIt)
library(rnaturalearthdata)
library(configr)
library(terra)
library(pbapply)
library(cleangeo)
library(doParallel)
library(foreach)
library(readr)
library(lwgeom)
library(countrycode)


setwd("/home/tws36/4c_evaluations")

proj_meta<-read.csv('./data/project_metadata/proj_meta.csv') 

# The list of projects to be run in this evaluation:
proj_meta<-read.csv('./data/project_metadata/proj_meta.csv') 
proj_to_eval<-read.table('./data/project_metadata/proj_to_eval.txt') %>%
  unlist() %>% as.numeric()
projects_agb<-read_csv('./data/GEDI/project_agb.csv')

# For testing on local
# data_suffix<- '230313'
# proj_id <- 1201
# site <- 'Gola'

cpc_rename<-function(x, t0){
    x %<>% 
        as.data.frame() %>%
        select(starts_with('cpc'))
    mycolnames<-colnames(x)
    years<-mycolnames %>% str_extract('[:digit:]+') %>% as.numeric()
    suffix<-mycolnames %>% str_extract('_u|_d') %>%
        str_replace('u', '1') %>%
        str_replace('d', '3')
    newnames<-paste0('JRC', t0-years, suffix)
    colnames(x)<-newnames
    return(x)
}

tmfemi_reformat<-function(df, t0){
    df %<>%
        st_as_sf(coords = c("lng","lat")) %>%
        rename(accessibility = access) %>%
        rename_with(~ gsub("luc_", "JRC", .x, fixed = TRUE))
    
    other<-df %>% select(-starts_with('cpc'))

    cpcs<-df %>%
        select(starts_with('cpc')) %>%
        cpc_rename(t0 = t0)
    
    df<-cbind(df, cpcs)
        # JRC2002_1 = cpc10_u,
        # JRC2007_1 = cpc5_u,
        # JRC2012_1 = cpc0_u,
        # JRC2002_3 = cpc10_d,
        # JRC2007_3 = cpc5_d,
        # JRC2012_3 = cpc0_d) %>%

  if(any(colnames(df) %>% str_detect('ecoregion')))
    df %<>% rename(biome = ecoregion)
  return(df)

}

config<-read.config('./config/fixed_config_sherwood.ini')
config$USERPARAMS$data_path <- '/maps/pf341/tom'
#write.config(config, './config/fixed_config_tmp.ini')

sapply(list.files('./R', full.names = TRUE, pattern = '.R$'), source)

# Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = FALSE)

# theme_set(theme_minimal())
# theme_replace(panel.grid.minor = element_line(colour = "red"))

source('./R/scripts/setup_and_reusable/load_config.R')
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
match_classes <- c(1,3)


readdir<-'/maps/pf341/results/2024-january-pipeline/'
# '/maps/pf341/tom-add-paper'
project_paths<-list.files(readdir, full = TRUE) %>% 
    str_subset('pairs')

i<-3
mylist<-mclapply(1:length(project_paths), mc.cores = 30, function(i){
#mylist<-mclapply(1:10, mc.cores = 30, function(i){
   
    myproject_path<-project_paths[i]
    proj_id<-basename(myproject_path) %>% str_replace('_pairs', '')

        # extract project start date: 
    myproj<-proj_meta %>%
        filter(ID == proj_id)
    
    print(myproj)

    t0 <- myproj$t0
    eval_end <- myproj$eval_end
    eval_end<-ifelse(is.na(eval_end), 2021, eval_end)

    site<-paste('VCS', proj_id, sep = '_')
    if(!str_detect(supplier_path, '.shp')){
    aoi_path<-file.path(supplier_path, site, 'GIS', 'aoi.shp')
    if(!file.exists(aoi_path))
        aoi_path<-file.path(supplier_path, site, paste(site, '.shp', sep = ''))
    }
    aoi_project<-read_sf(aoi_path) %>% # NEED TO CONSISTENTLY SET AOI NAME
    st_make_valid() %>%
    st_union()

    aoi_project<-aoi_project %>% st_transform(4326) 
    # Find the area of the region:
    project_area_ha <- st_area_ha(aoi_project)


    # find paths to match and unmatached points:
    pair_paths<-list.files(myproject_path, full = TRUE)
    matchless_ind<-pair_paths %>% str_detect('matchless')
    matchless_paths<-pair_paths[matchless_ind]
    matched_paths<-pair_paths[!matchless_ind]

    # Read and analyse pairs:
    #my_matched_path<-matched_paths[1]
    #j<-4
    project_estimates<-lapply(1:length(matched_paths), function(j){

        pairs<-read_parquet(matched_paths[j]) 
        unmatched_pairs<-read_parquet(matchless_paths[j]) 

        control<-pairs %>%
        select(starts_with('s_')) %>%
        rename_with(~str_replace(.x, 's_', '')) %>%
        mutate(treatment = 'control') %>%
        tmfemi_reformat(t0=t0)

        treat<-pairs %>%
        select(starts_with('k_')) %>%
        rename_with(~str_replace(.x, 'k_', '')) %>%
        mutate(treatment = 'treatment') %>%
        tmfemi_reformat(t0=t0)

        exp_n_pairs<-nrow(treat) + nrow(unmatched_pairs)

        pts_matched<-rbind(treat,
            control
            )

        # m.out<-assess_balance(pts_matched, class_prefix = class_prefix, t0 = t0, match_years = match_years, match_classes = match_classes)
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

        y <- control_series$series %>%
            #   left_join(class_agb %>% select(class, agb), by = 'class') %>%
            #   mutate(class_agb = class_area * agb,
            #          # the above ground biomass
            #          class_co2e = class_agb * cf_c * cf_co2e) %>%
            filter(class ==1) %>%
            filter(year %in% c(t0, eval_end))
            
        yc <- y %>% filter(treatment == 'control') %>% pull(class_area)
        yt <- y %>% filter(treatment == 'treatment')  %>% pull(class_area)

        c_undisturbed_loss <- yc[1] - yc[2]
        t_undisturbed_loss <- yt[1] - yt[2]

        undisturbed_additionality <- c_undisturbed_loss - t_undisturbed_loss

        out_df<-data.frame(
            ID = proj_id,
            Start = t0, 
            End = eval_end,
            control_undisturbed = yc[2],
            project_undisturbed = yt[2],
            avoided_disturbance_ha = undisturbed_additionality,
            avoided_disturbance_ha_yr = undisturbed_additionality / (eval_end - t0)
                    )

        return(out_df)
    })
    project_estimates<-do.call(rbind, project_estimates)
})

# summarise_all(mylist[[1]], mean)

# item<-project_estimates

project_summaries<-lapply(mylist[-c(3, 33)], function(item){
    ID<-item$ID[1]
    item_means <- item %>% 
        select(-ID) %>%
        summarise_all(mean) %>%
        mutate(ID = ID) %>%
        relocate(ID)
    item_sd <-item %>%
        select(control_undisturbed, 
            project_undisturbed, 
            avoided_disturbance_ha, 
            avoided_disturbance_ha_yr) %>%
            summarise_all(list(sd = sd))
    item<-bind_cols(item_means, item_sd)
    return(item)
}) %>%
    bind_rows() %>%
    arrange(as.numeric(as.character(ID)))

project_summaries

project_summaries %>%
    ggplot(aes(x = avoided_disturbance_ha_yr)) +
    geom_histogram() +
    geom_vline(xintercept = 0)

write_parquet(project_summaries, 
    file.path(readdir, 'project_summaries.parquet'))

# this is where the files are: 
list.files("../../../maps/4C/cf_candidate/country/counterfactuals")