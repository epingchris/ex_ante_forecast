# 0. Setup ----
rm(list = ls())

#Load packages
# install.packages(c("arrow","configr", "tidyverse", "magrittr", "sf", "magrittr", "MatchIt",
#                    "rnaturalearthdata", "configr", "terra", "pbapply", "cleangeo", "doParallel",
#                    "foreach", "readr", "lwgeom", "rnaturalearth", "stars", "Metrics", "patchwork"), depends = TRUE)

library(tidyverse) #ggplot2, dplyr, stringr
library(magrittr) #pipe operators
library(units) #units::set_units
library(sf) #sf::st_area
library(arrow) #arrow::read_parquet
library(MatchIt) #MatchIt::matchit
library(boot) #boot::boot
library(Metrics) #rmse, mae
library(tibble) #tibble to store labels with bquote()
library(scales) #scales::trans_break
library(patchwork)
#library(pryr) #pryr::object_size
#library(parallel) #parallel::mclapply

#Remove dplyr summarise grouping message because it prints a lot
options(dplyr.summarise.inform = F)

#Load pre-defined functions
source("functions.r") #cpc_rename, tmfemi_reformat, simulate_area_series, make_area_series, assess_balance, make_match_formula
source("AdditionalityPair.r")
source("plotPlacebo.r")
source("plotBaseline.r")

FindFiles = function(dir, pattern, full = F, negate = F) {
  file_matched = list.files(dir, full = full) %>% str_subset(pattern, negate = negate)
  if(length(file_matched) == 0) {
    return(NA)
  } else {
    return(file_matched)
  }
}

BootOut = function(type, in_df, boot_n = 1000, summ = T) {
  boot_out = boot::boot(data = in_df,
                        statistic = function(dat, ind) mean(dat[ind, ], na.rm = T), #function for bootstrapped mean
                        R = boot_n)
  if(summ) {
    boot_ci = boot::boot.ci(boot.out = boot_out, type = "perc")
    boot_summ = data.frame(type = type,
                           mean = mean(boot_out$t),
                           ci_lower = boot_ci$percent[4],
                           ci_upper = boot_ci$percent[5])
    return(boot_summ)
  } else {
    return(boot_out$t)
  }
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
#9. OPTIONAL: proj_ID: this is only used to remove the trailing "a" that I added to the filename of some shapefiles that needed fixing,
# So that I can retrieve country and t0 information from the the proj_info data frame
#10. out_path: absolute paths of the directory where outputs are to be saved; include file prefix if desired


# A. Read input (E-Ping's workflow) ----

#Define analysis type
analysis_type = "ongoing"
#"control": non-project polygons
#"ongoing": ongoing REDD+ projects (best-matched and loosely-matched baselines)
#"ac": Amazonian Collective polygons
ofir = F

polygon_dir = "/maps/epr26/tmf-data/projects/" #where polygons are stored
lagged_dir = "/maps/epr26/tmf_pipe_out_lagged/" #where the results for the lagged baselines are stored
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored
fig_path = paste0("/maps/epr26/ex_ante_forecast_out/out_") #where figures are stored
projects_to_exclude = c("674", "934", "2502", "1408", "sa11", "1122", "1133") #which projects to exclude manually; 1122 not done yet

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

pair_dirs_lagged = paste0(lagged_dir, projects, "/pairs/")
k_paths_lagged = rep(NA, length(projects))
m_paths_lagged = rep(NA, length(projects))
for(i in seq_along(projects)) {
  project_out_dir_lagged = paste0(lagged_dir, projects[i])
  k_paths_lagged[i] = FindFiles(project_out_dir_lagged, "k.parquet", full = T)
  m_paths_lagged[i] = FindFiles(project_out_dir_lagged, "matches.parquet", full = T)
}

polygon_paths = paste0(polygon_dir, projects, ".geojson")
proj_ID = gsub("(?<=[0-9])a$", "", projects, perl = T)
project_var = proj_info[match(proj_ID, proj_info$ID), ]
country = project_var$COUNTRY
t0_vec = project_var$t0

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

#Output: project-level variables
write.csv(project_var, paste0(out_path, "_project_var.csv"), row.names = F)
#project_var = read.csv(paste0(out_path, "_project_var.csv"), header = T)


# A2. Read input (user-defined) ----

#projects = NULL
#pair_dirs = NULL
#k_paths = NULL
#m_paths = NULL
#pair_dirs_lagged = NULL
#k_paths_lagged = NULL
#m_paths_lagged = NULL
#polygon_paths = NULL
#country = NULL
#t0_vec = NULL
#area_ha_vec = NULL
#acd_list = NULL
#out_path = NULL
#fig_path = NULL

# B. Get additionality and baseline ----
for(i in seq_along(projects)) {

    t0 = t0_vec[i]
    luc_t_20 = paste0("luc_", t0_vec[i] - 20)
    luc_t_10 = paste0("luc_", t0_vec[i] - 10)
    luc_t0 = paste0("luc_", t0_vec[i])

    area_ha = area_ha_vec[i]
    acd = acd_list[[i]]
    pair_dir = pair_dirs[i]

    k = read_parquet(k_paths[i]) %>%
      rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
      as.data.frame() %>%
      dplyr::select(lat, lng, k_ecoregion)

    setM = read_parquet(m_paths[i]) %>%
      rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), s_ecoregion = ecoregion) %>%
      as.data.frame()

    matches = setM %>%
      dplyr::select(lat, lng, s_ecoregion)

    a = Sys.time()
    pairs_best = AdditionalityPair(pair_dir = pair_dir, t0 = t0, area_ha = area_ha, acd = acd,
                                   k = k, matches = matches, lagged = F)
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

    setM = setM %>%
      dplyr::select(-starts_with("luc_")) %>%
      dplyr::select(-starts_with("cpc"))
    if(nrow(setM) > 250000) setM = setM[sample(nrow(setM), 250000), ]

    baseline_loose = setM %>%
        mutate(acd10 = acd$carbon.density[match(luc10, acd$land.use.class)],
               acd0 = acd$carbon.density[match(luc0, acd$land.use.class)],
               c_loss = (acd10 - acd0) / 10) %>%
        dplyr::select(c_loss) %>%
        mutate(project = projects[i])

    write.csv(baseline_loose, paste0(out_path, "_baseline_loose_", projects[i], ".csv"), row.names = F)

    baseline_lagged = data.frame(c_loss = numeric())
    if(!(is.na(k_paths_lagged[i]) | is.na(m_paths_lagged[i]))) {
        pair_dir_lagged = pair_dirs_lagged[i]
        k_lagged = read_parquet(k_paths_lagged[i]) %>%
          rename(luc10 = all_of(luc_t_10), luc0 = all_of(luc_t0), k_ecoregion = ecoregion) %>%
          as.data.frame() %>%
          dplyr::select(lat, lng, k_ecoregion)
        matches_lagged = read_parquet(m_paths_lagged[i]) %>%
          rename(luc10 = all_of(luc_t_20), luc0 = all_of(luc_t_10), s_ecoregion = ecoregion) %>%
          as.data.frame() %>%
          dplyr::select(lat, lng, s_ecoregion)

        a = Sys.time()
        pairs_lagged = AdditionalityPair(pair_dir = pair_dir_lagged, t0 = t0, area_ha = area_ha, acd = acd,
                                         k = k_lagged, matches = matches_lagged, lagged = T)
        b = Sys.time()
        cat("Project", i, "/", length(projects), "-", projects[i], "- pairs_lagged :", b - a, "\n")

        baseline_lagged = lapply(pairs_lagged, function(x) {
            filter(x$out_df, year <= 0 & year > -10) %>%
            dplyr::select(c_loss) %>%
            mutate(c_loss = c_loss / area_ha)
        }) %>%
            do.call(rbind, .) %>%
            mutate(project = projects[i])
    }

    write.csv(baseline_lagged, paste0(out_path, "_baseline_lagged_", projects[i], ".csv"), row.names = F)
}


# C. Bootstrap baselines ----
pre_cf_c_loss_boot_list = vector("list", length(projects))
pre_p_c_loss_boot_list = vector("list", length(projects))
post_cf_c_loss_boot_list = vector("list", length(projects))
post_p_c_loss_boot_list = vector("list", length(projects))
observed_add_boot_list = vector("list", length(projects))
baseline_best_boot_list = vector("list", length(projects))
baseline_loose_boot_list = vector("list", length(projects))
baseline_lagged_boot_list = vector("list", length(projects))

for(i in seq_along(projects)) {
    area_i = project_var$area_ha[i]
    observed = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T) %>%
      mutate(t_loss = t_loss / area_i,
             c_loss = c_loss / area_i,
             additionality = additionality / area_i)
    obs_pre = observed %>% filter(started == F)
    obs_post = observed %>% filter(started == T)

    baseline_best = read.csv(paste0(out_path, "_baseline_best_", projects[i], ".csv"), header = T)
    baseline_loose = read.csv(paste0(out_path, "_baseline_loose_", projects[i], ".csv"), header = T)
    baseline_lagged = read.csv(paste0(out_path, "_baseline_lagged_", projects[i], ".csv"), header = T)

    time_a = Sys.time()
    pre_cf_c_loss_boot_list[[i]] = BootOut(type = "cf_c_loss", in_df = dplyr::select(obs_pre, c_loss)) %>%
      mutate(project = projects[i])
    time_b = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- pre_cf_c_loss_boot :", time_b - time_a, "\n")

    pre_p_c_loss_boot_list[[i]] = BootOut(type = "p_c_loss", in_df = dplyr::select(obs_pre, t_loss)) %>%
      mutate(project = projects[i])
    time_c = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- pre_p_c_loss_boot :", time_c - time_b, "\n")

    post_cf_c_loss_boot_list[[i]] = BootOut(type = "cf_c_loss", in_df = dplyr::select(obs_post, c_loss)) %>%
      mutate(project = projects[i])
    time_d = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- post_cf_c_loss_boot :", time_d - time_c, "\n")

    post_p_c_loss_boot_list[[i]] = BootOut(type = "p_c_loss", in_df = dplyr::select(obs_post, t_loss)) %>%
      mutate(project = projects[i])
    time_e = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- post_p_c_loss_boot :", time_e - time_d, "\n")

    observed_add_boot_list[[i]] = BootOut(type = "additionality", in_df = dplyr::select(observed, additionality)) %>%
      mutate(project = projects[i])
    time_f = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- add_boot :", time_f - time_e, "\n")

    baseline_best_boot_list[[i]] = BootOut(type = "best", in_df = dplyr::select(baseline_best, c_loss)) %>%
        mutate(project = projects[i])
    time_g = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_best_boot :", time_g - time_f, "\n")

    baseline_loose_boot_list[[i]] = BootOut(type = "loose", in_df = dplyr::select(baseline_loose, c_loss)) %>%
        mutate(project = projects[i])
    time_h = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_loose_boot :", time_h - time_g, "\n")

    if(nrow(baseline_lagged) > 0) {
      baseline_lagged_boot_list[[i]] = BootOut(type = "lagged", in_df = dplyr::select(baseline_lagged, c_loss)) %>%
        mutate(project = projects[i])
    } else {
      baseline_lagged_boot_list[[i]] = data.frame(type = character(), mean = numeric(), ci_lower = numeric(), ci_upper = numeric(), project = character())
    }
    time_i = Sys.time()
    cat("Project", i, "/", length(projects), "-", projects[i], "- baseline_lagged_boot :", time_i - time_h, "\n")
}

pre_cf_c_loss_boot_df = do.call(rbind, pre_cf_c_loss_boot_list)
pre_p_c_loss_boot_df = do.call(rbind, pre_p_c_loss_boot_list)
post_cf_c_loss_boot_df = do.call(rbind, post_cf_c_loss_boot_list)
post_p_c_loss_boot_df = do.call(rbind, post_p_c_loss_boot_list)
observed_add_boot_df = do.call(rbind, observed_add_boot_list)
baseline_best_boot_df = do.call(rbind, baseline_best_boot_list)
baseline_loose_boot_df = do.call(rbind, baseline_loose_boot_list)
baseline_lagged_boot_df = do.call(rbind, baseline_lagged_boot_list)

append_result = T
write.table(pre_cf_c_loss_boot_df, paste0(out_path, "_pre_cf_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_pre_cf_c_loss.csv")), row.names = F, append = append_result)
write.table(pre_p_c_loss_boot_df, paste0(out_path, "_pre_p_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_pre_p_c_loss.csv")), row.names = F, append = append_result)
write.table(post_cf_c_loss_boot_df, paste0(out_path, "_post_cf_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_post_cf_c_loss.csv")), row.names = F, append = append_result)
write.table(post_p_c_loss_boot_df, paste0(out_path, "_post_p_c_loss.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_post_p_c_loss.csv")), row.names = F, append = append_result)
write.table(observed_add_boot_df, paste0(out_path, "_observed_add.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_observed_add.csv")), row.names = F, append = append_result)
write.table(baseline_best_boot_df, paste0(out_path, "_baseline_best_boot.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_baseline_best_boot.csv")), row.names = F, append = append_result)
write.table(baseline_loose_boot_df, paste0(out_path, "_baseline_loose_boot.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_baseline_loose_boot.csv")), row.names = F, append = append_result)
write.table(baseline_lagged_boot_df, paste0(out_path, "_baseline_lagged_boot.csv"), sep = ",",
            col.names = !file.exists(paste0(out_path, "_baseline_lagged_boot.csv")), row.names = F, append = append_result)


# D. Generate results ----

# Figure 4. show that there is no bias in our counterfactuals using placebo areas
if(analysis_type == "control") {
  continent_name = c(as = "Asia", af = "Africa", sa = "South America")
  pre_cf_c_loss_boot_df = read.csv(paste0(out_path, "_pre_cf_c_loss.csv"), header = T)
  pre_p_c_loss_boot_df = read.csv(paste0(out_path, "_pre_p_c_loss.csv"), header = T)
  post_cf_c_loss_boot_df = read.csv(paste0(out_path, "_post_cf_c_loss.csv"), header = T)
  post_p_c_loss_boot_df = read.csv(paste0(out_path, "_post_p_c_loss.csv"), header = T)

  summ = do.call(bind_rows, list(pre_cf_c_loss_boot_df %>% mutate(period = "pre"),
                                 pre_p_c_loss_boot_df %>% mutate(period = "pre"),
                                 post_cf_c_loss_boot_df %>% mutate(period = "post"),
                                 post_p_c_loss_boot_df %>% mutate(period = "post"))) %>%
    mutate(Continent = continent_name[str_sub(project, 1, 2)])


  p1 = plotPlacebo(dat = summ, period_used = "pre")
  p2 = plotPlacebo(dat = summ, period_used = "post")
  (p1 + p2) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom")
  ggsave(paste0(fig_path, "figure4_placebo_c_loss.png"), width = 7000, height = 3600, units = "px")

  #linear regression analyses
  summ_wide_mean = summ %>%
    filter(period == "post") %>%
    dplyr::select(type, mean, project, Continent) %>%
    pivot_wider(names_from = "type", values_from = "mean")

  lm_control = lm(p_c_loss - cf_c_loss ~ cf_c_loss, data = summ_wide_mean)
  lm_control_as = lm(p_c_loss - cf_c_loss ~ cf_c_loss, data = subset(summ_wide_mean, Continent == "Asia"))
  lm_control_af = lm(p_c_loss - cf_c_loss ~ cf_c_loss, data = subset(summ_wide_mean, Continent == "Africa"))
  lm_control_sa = lm(p_c_loss - cf_c_loss ~ cf_c_loss, data = subset(summ_wide_mean, Continent == "South America"))
  summary(lm_control)
  summary(lm_control_as)
  summary(lm_control_af)
  summary(lm_control_sa)
}


# Figure 5. show how baseline compares to counterfactual carbon loss in ongoing projects
if(analysis_type == "ongoing") {
  post_cf_c_loss_boot_df = read.csv(paste0(out_path, "_post_cf_c_loss.csv"), header = T)
  baseline_best_boot_df = read.csv(paste0(out_path, "_baseline_best_boot.csv"), header = T)
  baseline_loose_boot_df = read.csv(paste0(out_path, "_baseline_loose_boot.csv"), header = T)
  baseline_lagged_boot_df = read.csv(paste0(out_path, "_baseline_lagged_boot.csv"), header = T)

  summ = do.call(bind_rows, list(post_cf_c_loss_boot_df,
                                 baseline_best_boot_df,
                                 baseline_loose_boot_df,
                                 baseline_lagged_boot_df))

  p1 = plotBaseline(dat = summ, baseline_used = "best", note_x = 1.5, note_y = c(0.7, 0.5))
  p2 = plotBaseline(dat = summ, baseline_used = "loose", note_x = 1.5, note_y = c(0.7, 0.5))
  p3 = plotBaseline(dat = summ, baseline_used = "lagged", note_x = 1.5, note_y = c(0.7, 0.5))
  (p1 + p2 + p3) +
    plot_layout(axes = "collect", axis_titles = "collect")
  ggsave(paste0(fig_path, "figure5_ongoing_baseline_vs_cf_c_loss.png"), width = 7500, height = 4000, units = "px")

  #calculate RMSE and MAE
  error_df = CalcError(summ)
  write.csv(error_df, paste0(out_path, "_error.csv"))

  #fit LM while forcing through the origin (0, 0)
  summ_wide = summ %>%
    dplyr::select(type, mean, project) %>%
    pivot_wider(names_from = "type", values_from = "mean", id_expand = T)

  lm_best = lm(cf_c_loss ~ best - 1, data = summ_wide)
  lm_loose = lm(cf_c_loss ~ loose - 1, data = summ_wide)
  lm_lagged = lm(cf_c_loss ~ lagged - 1, data = summ_wide)
  adj_val = c(summary(lm_best)$coefficients[1],
              summary(lm_loose)$coefficients[1],
              summary(lm_lagged)$coefficients[1])

  summ_adj = do.call(bind_rows, list(post_cf_c_loss_boot_df,
                                     baseline_best_boot_df %>% mutate(across(mean:ci_upper, function(x) x * adj_val[1])),
                                     baseline_loose_boot_df %>% mutate(across(mean:ci_upper, function(x) x * adj_val[2])),
                                     baseline_lagged_boot_df %>% mutate(across(mean:ci_upper, function(x) x * adj_val[3]))))

  p4 = plotBaseline(dat = summ_adj, baseline_used = "best", note_x = 1.5, note_y = c(0.4, 0.2))
  p5 = plotBaseline(dat = summ_adj, baseline_used = "loose", note_x = 1.5, note_y = c(0.4, 0.2))
  p6 = plotBaseline(dat = summ_adj, baseline_used = "lagged", note_x = 1.5, note_y = c(0.4, 0.2))
  (p4 + p5 + p6) +
    plot_layout(axes = "collect", axis_titles = "collect")
  ggsave(paste0(fig_path, "figure6_ongoing_baseline_vs_cf_c_loss_corrected.png"), width = 7500, height = 4000, units = "px")

  #calculate RMSE and MAE
  error_adj_df = CalcError(summ_adj)
  write.csv(error_adj_df, paste0(out_path, "_error_adj.csv"))


  #Figure 7. compare additionality with baselines (best or time-lagged) to evaluate "project effectiveness"
  effectiveness_list = vector("list", length(projects))

  for(i in seq_along(projects)) {
      area_i = project_var$area_ha[i]
      baseline_best = read.csv(paste0(out_path, "_baseline_best_", projects[i], ".csv"), header = T)
      observed = read.csv(paste0(out_path, "_additionality_", projects[i], ".csv"), header = T) %>%
        filter(started == T) %>%
        mutate(t_loss = t_loss / area_i,
               c_loss = c_loss / area_i,
               additionality = additionality / area_i)

      baseline_best_boot = BootOut(type = "best", in_df = dplyr::select(baseline_best, c_loss), summ = F)
      observed_add_boot = BootOut(type = "additionality", in_df = dplyr::select(observed, additionality), summ = F)
      effectiveness = observed_add_boot / baseline_best_boot

      effectiveness_list[[i]] = data.frame(project = projects[i],
                                           best = mean(baseline_best_boot, na.rm = T),
                                           best_lower = quantile(baseline_best_boot, probs = 0.025, na.rm = T),
                                           best_upper = quantile(baseline_best_boot, probs = 0.975, na.rm = T),
                                           add = mean(observed_add_boot, na.rm = T),
                                           add_lower = quantile(observed_add_boot, probs = 0.025, na.rm = T),
                                           add_upper = quantile(observed_add_boot, probs = 0.975, na.rm = T),
                                           eff = mean(effectiveness, na.rm = T),
                                           eff_lower = quantile(effectiveness, probs = 0.025, na.rm = T),
                                           eff_upper = quantile(effectiveness, probs = 0.975, na.rm = T))

      cat("Project", i, "/", length(projects), "-", projects[i], "- effectiveness boostrapped\n")
  }

  effectiveness_df = effectiveness_list %>%
    do.call(rbind, .)
  write.csv(effectiveness_df, paste0(out_path, "_effectiveness.csv"), row.names = F)

  ggplot(data = effectiveness_df, aes(x = add, y = eff)) +
    geom_segment(aes(y = eff_lower, yend = eff_upper), color = "blue") +
    geom_segment(aes(x = add_lower, xend = add_upper), color = "blue") +
    geom_point(shape = 16, color = "blue", fill = "blue", size = 3) +
    geom_text(aes(label = project), vjust = -0.5, hjust = -0.1, size = 10) +
    geom_hline(yintercept = c(0.5, 1, 1.5), linetype = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-0.3, 0.8), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-6, 26), breaks = c(-5, 0, 0.5, 1, 1.5, 2, 5, seq(10, 25, by = 5)),
                       expand = c(0, 0),
                       transform = scales::pseudo_log_trans(base = 10, sigma = 0.1)) +
    labs(x = "Observed additionality (MgC/ha/yr)",
        y = "Project effectiveness") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 40),
          axis.text = element_text(size = 36),
          legend.position = "none")
  ggsave(paste0(fig_path, "figure7_effectiveness_vs_add_transformed.png"), width = 8000, height = 8000, units = "px")

  ggplot(data = subset(effectiveness_df, eff > 0), aes(x = add, y = eff)) +
    geom_segment(aes(y = eff_lower, yend = eff_upper), color = "blue") +
    geom_segment(aes(x = add_lower, xend = add_upper), color = "blue") +
    geom_point(shape = 16, color = "blue", fill = "blue", size = 3) +
    geom_text(aes(label = project), vjust = -0.5, hjust = -0.1, size = 10) +
    geom_hline(yintercept = c(0.5, 1, 1.5), linetype = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.2, 26), breaks = c(0.1, 0.2, 0.5, 1, 1.5, 2, 3, 5, 10, 25),
                       transform = scales::transform_log10()) +
    labs(x = "Observed additionality (MgC/ha/yr)",
         y = "Project effectiveness") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 40),
          axis.text = element_text(size = 36),
          legend.position = "none")
  ggsave(paste0(fig_path, "figure7_effectiveness_vs_add_excluded.png"), width = 8000, height = 8000, units = "px")

  ggplot(data = subset(effectiveness_df, eff > 0), aes(x = best, y = eff)) +
    geom_segment(aes(y = eff_lower, yend = eff_upper), color = "blue") +
    geom_segment(aes(x = best_lower, xend = best_upper), color = "blue") +
    geom_point(shape = 16, color = "blue", fill = "blue", size = 3) +
    geom_text(aes(label = project), vjust = -0.5, hjust = -0.1, size = 10) +
    geom_hline(yintercept = c(0.5, 1, 1.5), linetype = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.2, 25), breaks = c(0.1, 0.2, 0.5, 1, 1.5, 2, 3, 5, 10), expand = c(0, 0),
                       transform = scales::transform_log10()) +
    labs(x = "Predicted maximum additionality (MgC/ha/yr)",
         y = "Project effectiveness") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 40),
          axis.text = element_text(size = 36),
          legend.position = "none")
  ggsave(paste0(fig_path, "figure7_effectiveness_vs_baseline_excluded.png"), width = 8000, height = 8000, units = "px")
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


#Compare with Jody's output
analysis_type = "ongoing"
out_path = paste0("/maps/epr26/ex_ante_forecast_out/out_", analysis_type) #where outputs are stored
observed_add_boot_df = read.csv(paste0(out_path, "_observed_add.csv"), header = T)

observed_add_jody = read.csv("/maps/epr26/ex_ante_forecast_out/merged_additionality_by_jody.csv", header = T)
year_max = observed_add_jody[nrow(observed_add_jody), ]$year

observed_add_jody_mean = observed_add_jody %>%
  filter(year == year_max) %>%
  rename_with(function(x) gsub("X", "", gsub("_additionality", "", x))) %>%
  pivot_longer(cols = !year, names_to = "ID", values_to = "cumul_add") %>%
  mutate(ID = as.numeric(ID)) %>%
  left_join(x = ., y = project_var %>% dplyr::select("ID", "t0", "area_ha"), by = "ID") %>%
  filter(!is.na(t0)) %>%
  mutate(additionality_jody = cumul_add / ((year - t0) * area_ha))

observed_add_boot_compare = observed_add_boot_df %>%
  left_join(x = ., y = project_var %>% dplyr::select("ID", "t0", "area_ha"), by = join_by(project == ID)) %>%
  left_join(x = ., y = observed_add_jody_mean %>% dplyr::select("ID", "additionality_jody"), by = join_by(project == ID))
write.csv(observed_add_boot_compare, paste0("/maps/epr26/ex_ante_forecast_out/additionality_compare.csv"), row.names = F)
