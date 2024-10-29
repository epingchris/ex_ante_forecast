#retrieve time series of observed carbon loss and additionality per matched pair
AdditionalityPair = function(matched_path, matchless_path, k, matches, t0, area_ha, acd, pair_id, offset = F) {
    c = Sys.time()

    pairs = read_parquet(matched_path) %>%
      dplyr::select(-dplyr::any_of(c("k_x", "k_y", "s_x", "s_y")))

    #add ecoregion information to each matched pixel if it is not already there
    if("k_ecoregion" %in% colnames(pairs) == F) pairs = dplyr::left_join(pairs, k, by = join_by(k_lat == lat, k_lng == lng))
    if("s_ecoregion" %in% colnames(pairs) == F) pairs = dplyr::left_join(pairs, matches, by = join_by(s_lat == lat, s_lng == lng))

    unmatched_pairs = read_parquet(matchless_path)

    control = pairs %>%
        dplyr::select(starts_with("s_")) %>%
            rename_with(~str_replace(.x, "s_", "")) %>%
            mutate(treatment = "control") %>%
            tmfemi_reformat(t0 = t0)

    treat = pairs %>%
        dplyr::select(starts_with("k_")) %>%
            rename_with(~str_replace(.x, "k_", "")) %>%
            mutate(treatment = "treatment") %>%
            tmfemi_reformat(t0 = t0)

    # For when offset = T: remove columns "JRC2021" "JRC2022" which contain only NA's
    if(offset) {
        control = control[!str_detect(colnames(control), "JRC[:digit:]{4}$|JRC[1-9]$|JRC10$")]
        treat = treat[!str_detect(colnames(treat), "JRC[:digit:]{4}$|JRC[1-9]$|JRC10$")]
    }

    exp_n_pairs = nrow(treat) + nrow(unmatched_pairs)
    pts_matched = rbind(treat, control) %>%
        mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5,
               defor_10_5 = (cpc10_u - cpc5_u) / 5,
               defor_10_0 = (cpc10_u - cpc0_u) / 10,
               pair = pair_id)

    #calculate annual proportion of each land use class
    class_prefix = "JRC"
    match_years = c(0, -5, -10)
    match_classes = c(1, 3)
    luc_series = simulate_area_series(pts_matched, class_prefix = class_prefix, t0 = t0,
                                      match_years = match_years, match_classes = match_classes,
                                      exp_n = exp_n_pairs, area = area_ha, verbose = F, offset = offset)

    #calculate annual carbon stock (MgC)
    carbon_series = luc_series$series %>%
        merge(., acd, by.x = "class", by.y = "land.use.class", all.x = T) %>%
        mutate(carbon_content = class_area * carbon.density) %>% #MgC
        group_by(treatment, year) %>%
        summarise(carbon_content = sum(carbon_content, na.rm = T)) %>%
        ungroup()

    #calculate annual carbon flux and additionality (CO2e)
    carbon_wide = pivot_wider(carbon_series, names_from = "treatment", values_from = "carbon_content")
    out_df = data.frame(year = carbon_wide$year[-1],
                        c_loss = -diff(carbon_wide$control),
                        t_loss = -diff(carbon_wide$treatment)) %>%
        mutate(additionality = c_loss - t_loss, pair = pair_id)
    d = Sys.time()
    cat(pair_id, ":", d - c, "\n")
    return(list(out_df = out_df, pts_matched = pts_matched, exp_n_pairs = exp_n_pairs))
}
