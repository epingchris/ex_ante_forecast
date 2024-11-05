# 0. Setup ----
rm(list = ls())

#Load packages
# install.packages(c("arrow","configr", "tidyverse", "magrittr", "sf", "magrittr", "MatchIt",
#                    "rnaturalearthdata", "configr", "terra", "pbapply", "cleangeo", "doParallel",
#                    "foreach", "readr", "lwgeom", "rnaturalearth", "stars", "Metrics"), depends = TRUE)

library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators
library(units) #units::set_units
library(sf) #sf::st_area
library(arrow) #arrow::read_parquet
library(MatchIt) #MatchIt::matchit
library(boot) #boot::boot
library(Metrics) #rmse, mae
library(scales) #scales::trans_break
#library(pryr) #pryr::object_size
#library(ggpubr) #ggpubr::ggarrange
#library(cowplot)
#library(parallel) #parallel::mclapply

#Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = F)

#Load pre-defined functions
source("functions.r") #cpc_rename, tmfemi_reformat, simulate_area_series, make_area_series, assess_balance, make_match_formula
source("AdditionalityPair.r")
source("plotBaseline.r")

FindFiles = function(dir, pattern, full = F, negate = F) {
  file_matched = list.files(dir, full = full) %>% str_subset(pattern, negate = negate)
  if(length(file_matched) == 0) {
    return(NA)
  } else {
    return(file_matched)
  }
}

BootSumm = function(type, in_df, boot_n = 1000) {
  boot_out = boot::boot(data = in_df,
                        statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                        R = boot_n)
  boot_ci = boot::boot.ci(boot.out = boot_out, type = "norm")
  boot_summ = data.frame(type = type,
                         mean = mean(as.vector(boot_out$t)),
                         ci_lower = boot_ci$normal[2],
                         ci_upper = boot_ci$normal[3])
  return(boot_summ)
}


#Define input variables needed to read TMF implementation output and other data

# It requires the following input variables to read TMF implementation output and other data.
# All variables are vectors containing one value for each project to be analysed:

#1. projects: an index of all projects to be analysed
# This should correspond to the filenames of the shapefiles and to the -p argument in the implementation code
# It is usually be the ongoing projects' VCS ID or customised (e.g. prefixed series of integers)

#2. pair_dirs: absolute paths of the directories containing all matched pair pixel sets (typically "/pairs/xxx.parquet" and  "/pairs/xxx_matchless.parquet")
# The directory should containing pairs of parquet files with the same file name, with and without the "_matchless" suffix.
# This is used to calculate estimated observed additionality.

#3. k_paths: absolute paths of the set K (typically "k.parquet")
#4. m_paths: absolute paths of the set M (typically "matches.parquet")
# Both should be in parquet format, containing the following columns:
# "lat", "lng" (degrees), "slope" (degrees), "elevation" (metres), "access" (remoteness, in minutes), "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d" (from 0 to 1), "luc_[t-10]" to "luc_2021" (categorical, 1-6, based on the JRC-TMF dataset, Vancutsem et al. 2021), "ecoregion" (categorical, based on the RESOLVE dataset, Dinerstein et al. 2017)

#5. acd_paths: absolute paths of the carbon density per LUC (typically "carbon-density.csv")
# This should be an csv file (although the script could be modiified in the future to support txt format) containing columns "land.use.class" (categorical, 1-6) and "carbon.density" (MgC/ha) for all six LUCs, although the script checks and fill missing LUC with NAs

#6. polygon_paths: absolute paths of the shapefile of the project extent
# This should be a geojson file containing valid geometries in WGS84 (EPSG: 4326), although the script checks for both conditions.
# This is currently only used to calculate project area (ha), but could be useful for other purposes in the future.

#7. country: country of the project
#8. t0: year of start of the project (real or hypothetical)
#9. OPTIONAL: proj_ID: ID used to identify country t0 from the proj_info_all data frame
# This is only used to remove the trailing "a" that I added to the filename of some shapefiles that needed fixing
# If unspecified, the projects variable will be used
#10. out_path: absolute paths of the directory where outputs are to be saved; include file prefix if desired


# 0a. E-Ping's workflow to obtain input variables ----

#Define analysis type
analysis_type = "control"
#"ongoing": ongoing REDD+ projects (best-matched and loosely-matched baselines)
#"control": non-project polygons
#"ac": Amazonian Collective polygons
ofir = F

polygon_dir = "/maps/epr26/tmf-data/projects/" #where polygons are stored
if(analysis_type == "control") { #where the results for the offsetted baselines are stored
    offset_dir = "/maps/epr26/tmf_pipe_out_offset_new/"
} else {
    offset_dir = "/maps/epr26/tmf_pipe_out_offset/"
}
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored
fig_path = paste0("/maps/epr26/ex_ante_forecast_out/out_") #where figures are stored
projects_to_exclude = c("674", "934", "2502", "1408") #which projects to exclude manually

#Based on analysis type, find project names and store in vector "projects", and project basic info and store in dataframe "proj_info"
if(analysis_type == "ongoing") {
  project_dir = "/maps/epr26/tmf_pipe_out/" #define where implementation code results are stored

  #Load basic information (csv file copied from Tom's directory)
  proj_info = read.csv("proj_meta.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #Find all project IDs
  exclude_strings = c("slopes", "elevation", "srtm", "ac", "as", "\\.", "\\_", "0000", "9999")
  projects = map(exclude_strings, function(x) FindFiles(project_dir, x, negate = T)) %>%
    reduce(intersect)

} else if(analysis_type == "control") {
  project_dir = "/maps/epr26/tmf_pipe_out_luc_t/" #define where implementation code results are stored

  #Load basic information
  proj_info = read.csv("proj_meta_control.csv") %>%
    dplyr::select(ID, COUNTRY, t0)

  #Find all project IDs
  projects = map(c("asn", "af", "sa"), function(x) FindFiles(project_dir, x)) %>%
    reduce(union) %>%
    str_subset("\\.", negate = T)

} else if(analysis_type == "ac") {
  project_dir = "/maps/epr26/tmf_pipe_out/" #define where implementation code results are stored

  #Load basic information
  proj_info = data.frame(ID = paste0("ac", sprintf("%02d", c(1:6))), COUNTRY = "Brazil", t0 = 2021)

  #Amazonian Collective polygons
  projects = paste0("ac", sprintf("%02d", c(1:6)))

}

done_vec = sapply(projects, function(x) FindFiles(paste0(project_dir, x, "/pairs"), ".parquet") %>% length() == 200)

#only keep projects with complete ACD values for LUC 1, 2, 3, and 4
full_acd_vec = sapply(projects, function(x) {
  acd_path = FindFiles(paste0(project_dir, x), "carbon-density", full = T)
  if(!is.na(acd_path)) acd = read.csv(acd_path)
  acd = acd[, 1:2]
  colnames(acd) = c("land.use.class", "carbon.density")
  return(Reduce("&", 1:4 %in% acd$land.use.class))
})

projects_df = data.frame(project = projects, done = done_vec, full_acd = full_acd_vec)
projects_df$to_exclude = projects_df$project %in% projects_to_exclude
write.csv(projects_df, paste0(out_path, "_project_status.csv"), row.names = F)
projects = subset(projects_df, done & full_acd & !to_exclude)$project

#Produce input variables needed for the analysis
pair_dirs = paste0(project_dir, projects, "/pairs/")
k_paths = rep(NA, length(projects))
m_paths = rep(NA, length(projects))
acd_paths = rep(NA, length(projects))
for(i in seq_along(projects)) {
  project_out_dir = paste0(project_dir, projects[i])
  k_paths[i] = FindFiles(project_out_dir, "k.parquet", full = T)
  m_paths[i] = FindFiles(project_out_dir, "matches.parquet", full = T)
  acd_paths[i] = FindFiles(project_out_dir, "carbon-density.csv", full = T)
}

pair_dirs_offset = paste0(offset_dir, projects, "/pairs/")
k_paths_offset = rep(NA, length(projects))
m_paths_offset = rep(NA, length(projects))
for(i in seq_along(projects)) {
  project_out_dir_offset = paste0(offset_dir, projects[i])
  k_paths_offset[i] = FindFiles(project_out_dir_offset, "k.parquet", full = T)
  m_paths_offset[i] = FindFiles(project_out_dir_offset, "matches.parquet", full = T)
}

polygon_paths = paste0(polygon_dir, projects, ".geojson")
proj_ID = gsub("(?<=[0-9])a$", "", projects, perl = T)
project_var = proj_info[match(proj_ID, proj_info$ID), ]
country = project_var$COUNTRY
t0_vec = project_var$t0


# 0b. User-defined input variables ----

#projects = NULL
#pair_dirs = NULL
#k_paths = NULL
#m_paths = NULL
#acd_paths = NULL
#polygon_paths = NULL
#country = NULL
#t0 = NULL
#proj_ID = NULL
#out_path = NULL


# A. Read data ----

#vector containing area (ha) of every project
area_ha_vec = sapply(seq_along(projects), function(i) {
  area_ha_i = st_read(polygon_paths[i]) %>%
    st_make_valid() %>%
    st_union() %>%
    st_transform(4326) %>%
    st_area() %>% #area in m^2
    set_units(ha) #convert into hectares
  return(area_ha_i)
})
project_var$area_ha = area_ha_vec

#list containing data frame of ACD (MgC/ha) per LUC of every project
acd_list = lapply(seq_along(projects), function(i) {
  acd_i = read.csv(acd_paths[i])
  acd_i = acd_i[, 1:2]
  colnames(acd_i) = c("land.use.class", "carbon.density")
  for(class in 1:6) {
    if(class %in% acd_i$land.use.class == F) acd_i = rbind(acd_i, c(class, NA))
  }
  return(acd_i)
})
project_var$acd_undisturbed = sapply(acd_list, function(x) filter(x, land.use.class == 1)$carbon.density)

#list containing set K of every project
setK = lapply(seq_along(projects), function(i) {
  luc_t_10 = paste0("luc_", t0_vec[i] - 10)
  luc_t0 = paste0("luc_", t0_vec[i])
  read_parquet(k_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
    as.data.frame()
})

#list containing set M of every project
setM = lapply(seq_along(projects), function(i) {
  luc_t_10 = paste0("luc_", t0_vec[i] - 10)
  luc_t0 = paste0("luc_", t0_vec[i])
  read_parquet(m_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
    as.data.frame()
})

#list containing set K of every project: for offsetted baseline
setK_offset = lapply(seq_along(projects), function(i) {
  if(is.na(k_paths_offset[i])) return(NULL)
  luc_t_10 = paste0("luc_", t0_vec[i] - 10)
  luc_t0 = paste0("luc_", t0_vec[i])
  read_parquet(k_paths[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
    as.data.frame()
})

#list containing set M of every project: for offsetted baseline
setM_offset = lapply(seq_along(projects), function(i) {
  if(is.na(m_paths_offset[i])) return(NULL)
  luc_t_10 = paste0("luc_", t0_vec[i] - 20) #offsetted baseline is matched to ten years
  luc_t0 = paste0("luc_", t0_vec[i] - 10)
  read_parquet(m_paths_offset[i]) %>%
    rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
    as.data.frame()
})


#Output: project-level variables
write.csv(project_var, paste0(out_path, "_project_var.csv"), row.names = F)
#project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)

# B. Get additionality and baseline ----
for(i in seq_along(projects)) {

    t0 = t0_vec[i]
    area_ha = area_ha_vec[i]
    acd = acd_list[[i]]

    pair_dir = pair_dirs[i]
    k = dplyr::select(setK[[i]], c("lat", "lng", "k_ecoregion"))
    matches = dplyr::select(setM[[i]], c("lat", "lng", "s_ecoregion"))
    a = Sys.time()
    pairs_best = AdditionalityPair(pair_dir = pair_dir, t0 = t0, area_ha = area_ha, acd = acd, k = k, matches = matches, offset = F)
    b = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- pairs_best :", b - a, "\n")

    additionality = lapply(pairs_best, function(x) x$out_df) %>%
        do.call(dplyr::bind_rows, .) %>%
        mutate(started = ifelse(year > t0, T, F)) %>%
        mutate(project = projects[i])

    baseline_best = lapply(pairs_best, function(x) {
        filter(x$out_df, year <= t0 & year > t0 - 10) %>%
        dplyr::select(c_loss) %>%
        mutate(c_loss = c_loss / area_ha)
    }) %>%
        do.call(rbind, .) %>%
        mutate(project = projects[i])

    write.csv(additionality, paste0(out_path, "_additionality_", projects[i], ".csv"), row.names = F)
    write.csv(baseline_best, paste0(out_path, "_baseline_best_", projects[i], ".csv"), row.names = F)

    matches = setM[[i]] %>%
      dplyr::select(-starts_with("luc_")) %>%
      dplyr::select(-starts_with("cpc")) %>%
      as.data.frame()
    if(nrow(matches) > 250000) matches = matches[sample(nrow(matches), 250000), ]

    baseline_loose = matches %>%
        mutate(acd10 = acd$carbon.density[match(luc10, acd$land.use.class)],
               acd0 = acd$carbon.density[match(luc0, acd$land.use.class)],
               c_loss = (acd10 - acd0) / 10) %>%
        dplyr::select(c_loss) %>%
        mutate(project = projects[i])

    write.csv(baseline_loose, paste0(out_path, "_baseline_loose_", projects[i], ".csv"), row.names = F)

    baseline_offset = data.frame(c_loss = numeric())
    if(!(is.null(setK_offset[[i]]) | is.null(setM_offset[[i]]))) {
        pair_dir_offset = pair_dirs_offset[i]
        k_offset = dplyr::select(setK_offset[[i]], c("lat", "lng", "k_ecoregion"))
        matches_offset = dplyr::select(setM_offset[[i]], c("lat", "lng", "s_ecoregion"))

        a = Sys.time()
        pairs_offset = AdditionalityPair(pair_dir = pair_dir_offset, t0 = t0, area_ha = area_ha, acd = acd,
                                             k = k_offset, matches = matches_offset, offset = T)
        b = Sys.time()
        cat("Project", i, "/", length(projects), "-", projects[i], "- pairs_offset :", b - a, "\n")

        baseline_offset = lapply(pairs_offset, function(x) {
            filter(x$out_df, year <= 0 & year > -10) %>%
            dplyr::select(c_loss) %>%
            mutate(c_loss = c_loss / area_ha)
        }) %>%
            do.call(rbind, .) %>%
            mutate(project = projects[i])
    }

    write.csv(baseline_offset, paste0(out_path, "_baseline_offset_", projects[i], ".csv"), row.names = F)
}


# C. Bootstrap baselines ----
observed_cf_c_loss_boot_list = vector("list", length(projects))
observed_p_c_loss_boot_list = vector("list", length(projects))
baseline_best_boot_list = vector("list", length(projects))
baseline_loose_boot_list = vector("list", length(projects))
baseline_offset_boot_list = vector("list", length(projects))

for(i in seq_along(projects)) {
    area_i = project_var$area_ha[i]
    additionality = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T) %>%
        mutate(t_loss = t_loss / area_i,
               c_loss = c_loss / area_i)
    baseline_best = read.csv(paste0(out_path, "_baseline_best_", projects[i], ".csv"), header = T)
    baseline_loose = read.csv(paste0(out_path, "_baseline_loose_", projects[i], ".csv"), header = T)
    baseline_offset = read.csv(paste0(out_path, "_baseline_offset_", projects[i], ".csv"), header = T)

    time_a = Sys.time()
    observed_cf_c_loss_boot_list[[i]] = BootSumm(type = "cf_c_loss", in_df = dplyr::select(additionality, c_loss)) %>%
        mutate(project = projects[i])
    time_b = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- cf_c_loss_boot :", time_b - time_a, "\n")

    observed_p_c_loss_boot_list[[i]] = BootSumm(type = "p_c_loss", in_df = dplyr::select(additionality, t_loss)) %>%
        mutate(project = projects[i])
    time_c = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- p_c_loss_boot :", time_c - time_b, "\n")

    baseline_best_boot_list[[i]] = BootSumm(type = "best", in_df = dplyr::select(baseline_best, c_loss)) %>%
        mutate(project = projects[i])
    time_d = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_best_boot :", time_d - time_c, "\n")

    baseline_loose_boot_list[[i]] = BootSumm(type = "loose", in_df = dplyr::select(baseline_loose, c_loss)) %>%
        mutate(project = projects[i])
    time_e = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_loose_boot :", time_e - time_d, "\n")

    if(nrow(baseline_offset) > 0) {
        baseline_offset_boot_list[[i]] = BootSumm(type = "offset", in_df = dplyr::select(baseline_offset, c_loss)) %>%
            mutate(project = projects[i])
    } else {
        baseline_offset_boot_list[[i]] = data.frame(type = character(), mean = numeric(), ci_lower = numeric(), ci_upper = numeric(), project = character())
    }
    time_f = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_offset_boot :", time_f - time_e, "\n")
}

observed_cf_c_loss_boot_df = do.call(rbind, observed_cf_c_loss_boot_list)
observed_p_c_loss_boot_df = do.call(rbind, observed_p_c_loss_boot_list)
baseline_best_boot_df = do.call(rbind, baseline_best_boot_list)
baseline_loose_boot_df = do.call(rbind, baseline_loose_boot_list)
baseline_offset_boot_df = do.call(rbind, baseline_offset_boot_list)

write.csv(observed_cf_c_loss_boot_df, paste0(out_path, "_observed_cf_c_loss.csv"), row.names = F)
write.csv(observed_p_c_loss_boot_df, paste0(out_path, "_observed_p_c_loss.csv"), row.names = F)
write.csv(baseline_best_boot_df, paste0(out_path, "_baseline_best_boot.csv"), row.names = F)
write.csv(baseline_loose_boot_df, paste0(out_path, "_baseline_loose_boot.csv"), row.names = F)
write.csv(baseline_offset_boot_df, paste0(out_path, "_baseline_offset_boot.csv"), row.names = F)

observed_cf_c_loss_boot_df = read.csv(paste0(out_path, "_observed_cf_c_loss.csv"), header = T)
observed_p_c_loss_boot_df = read.csv(paste0(out_path, "_observed_p_c_loss.csv"), header = T)
baseline_best_boot_df = read.csv(paste0(out_path, "_baseline_best_boot.csv"), header = T)
baseline_loose_boot_df = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
baseline_offset_boot_df = read.csv(paste0(out_path, "_baseline_offset_boot.csv"), header = T)


# D. Generate results ----

if(analysis_type == "control") {
    continent_name = c(as = "Asia", af = "Africa", sa = "South America")
    c_loss_control_summ = do.call(bind_rows, list(observed_cf_c_loss_boot_df,
                                                       observed_p_c_loss_boot_df,
                                                       baseline_best_boot_df,
                                                       baseline_loose_boot_df,
                                                       baseline_offset_boot_df)) %>%
        dplyr::select(type, mean, project) %>%
        pivot_wider(names_from = "type", values_from = "mean", id_expand = T) %>%
        mutate(Continent = continent_name[str_sub(project, 1, 2)],
               c_loss_min = pmin(best, loose, offset, na.rm = T),
               c_loss_max = pmax(best, loose, offset, na.rm = T))

    #Figure 4. show that there is no bias in our counterfactuals
    ggplot(data = c_loss_control_summ, aes(x = p_c_loss, y = cf_c_loss)) +
        geom_point(aes(shape = Continent, color = Continent, fill = Continent), size = 4) +
        geom_abline(intercept = 0, slope = 1, linetype = 2) +
        scale_shape_manual(values = c(Asia = 1, Africa = 3, `South America` = 18)) +
        scale_color_manual(values = c(Asia = "blue", Africa = "black", `South America` = "red")) +
        scale_fill_manual(values = c(Asia = NA, Africa = NA, `South America` = "red")) +
        scale_x_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) + #ensures no padding
        scale_y_continuous(limits = c(-0.2, 2.6), expand = c(0, 0)) +
        labs(x = "Observed project carbon loss (MgC/ha/yr)",
             y = "Observed counterfactual carbon loss (MgC/ha/yr)") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title = element_text(size = 18),
              axis.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14),
              legend.position = "bottom")
    ggsave(paste0(fig_path, "figure4_c_loss_p_vs_cf_post.png"), width = 2500, height = 2600, units = "px")

    #calculate RMSE and MAE
    p_c_loss_val = na.omit(c_loss_control_summ)$p_c_loss
    cf_c_loss_val = na.omit(c_loss_control_summ)$cf_c_loss

    Metrics::rmse(p_c_loss_val, cf_c_loss_val)
    Metrics::mae(p_c_loss_val, cf_c_loss_val)


    #Figure 5. show how baseline compares to counterfactual carbon loss in non-project areas
    plotBaseline(dat = c_loss_control_summ, baseline_used = "all", use_log10 = T)
    p1 = plotBaseline(dat = c_loss_control_summ, baseline_used = "best", use_log10 = T)
    p2 = plotBaseline(dat = c_loss_control_summ, baseline_used = "loose", use_log10 = T)
    p3 = plotBaseline(dat = c_loss_control_summ, baseline_used = "offset", use_log10 = T)

    #calculate RMSE and MAE of each baseline compared to observed counterfactual C loss
    p_c_loss_val = na.omit(c_loss_control_summ)$p_c_loss
    cf_c_loss_val = na.omit(c_loss_control_summ)$cf_c_loss
    best_val = na.omit(c_loss_control_summ)$best
    loose_val = na.omit(c_loss_control_summ)$loose
    offset_val = na.omit(c_loss_control_summ)$offset

    error_df = data.frame(error_type = rep(c("rmse", "mae"), each = 3),
                          baseline_type = rep(c("best", "loose", "offset"), 2),
                          val = c(rmse(cf_c_loss_val, best_val), rmse(cf_c_loss_val, loose_val), rmse(cf_c_loss_val, offset_val),
                                  mae(cf_c_loss_val, best_val), mae(cf_c_loss_val, loose_val), mae(cf_c_loss_val, offset_val)))
    write.csv(error_df, paste0(out_path, "_error.csv"), row.names = F)
}

if(analysis_type == "ongoing") {
    c_loss_ongoing_summ = do.call(bind_rows, list(observed_cf_c_loss_boot_df,
                                                  baseline_best_boot_df,
                                                  baseline_loose_boot_df,
                                                  baseline_offset_boot_df)) %>%
        dplyr::select(type, mean, project) %>%
        pivot_wider(names_from = "type", values_from = "mean", id_expand = T) %>%
        mutate(c_loss_min = pmin(best, loose, offset, na.rm = T),
               c_loss_max = pmax(best, loose, offset, na.rm = T))

    #Figure 6. show how baseline compares to counterfactual carbon loss in ongoing projects
    plotBaseline(dat = c_loss_ongoing_summ, baseline_used = "all", use_log10 = T)
    plotBaseline(dat = c_loss_ongoing_summ, baseline_used = "best", use_log10 = T)
    plotBaseline(dat = c_loss_ongoing_summ, baseline_used = "loose", use_log10 = T)
    plotBaseline(dat = c_loss_ongoing_summ, baseline_used = "offset", use_log10 = T)

    #calculate RMSE and MAE of each baseline compared to observed counterfactual C loss
    cf_c_loss_val = na.omit(c_loss_ongoing_summ)$cf_c_loss
    best_val = na.omit(c_loss_ongoing_summ)$best
    loose_val = na.omit(c_loss_ongoing_summ)$loose
    offset_val = na.omit(c_loss_ongoing_summ)$offset

    error_df = data.frame(error_type = rep(c("rmse", "mae"), each = 3),
                          baseline_type = rep(c("best", "loose", "offset"), 2),
                          val = c(rmse(cf_c_loss_val, best_val), rmse(cf_c_loss_val, loose_val), rmse(cf_c_loss_val, offset_val),
                                  mae(cf_c_loss_val, best_val), mae(cf_c_loss_val, loose_val), mae(cf_c_loss_val, offset_val)))
    write.csv(error_df, paste0(out_path, "_error.csv"), row.names = F)
}


#OPTIONAL output: only basic variables, additionality distribution data to send to Ofir
if(ofir) {
    write.table(project_var %>% dplyr::select(project, t0, country, area_ha),
                paste0(out_path, "_project_var_basic.csv"), sep = ",", row.names = F)

    additionality_distribution = lapply(seq_along(projects), function(i) {
        additionality = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T)

        additionality %>%
            filter(started) %>%
            dplyr::select(year, additionality, pair) %>%
            mutate(project = projects[i])
    }) %>%
        do.call(rbind, .)
    write.csv(additionality_distribution, paste0("/maps/epr26/tmf_pipe_out/additionality_distribution.csv"), row.names = F)
}