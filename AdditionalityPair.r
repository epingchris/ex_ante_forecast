#retrieve time series of observed carbon loss and additionality per matched pair
AdditionalityPair = function(matched_path, matchless_path, k, matches, t0, area_ha, acd, pair_id) {
    c = Sys.time()

    #add ecoregion information to each matched pixel
    pairs = read_parquet(matched_path) %>%
    dplyr::left_join(., k, by = join_by(k_lat == lat, k_lng == lng)) %>%
    dplyr::left_join(., matches, by = join_by(s_lat == lat, s_lng == lng))
    #mutate(s_id = 1:n(), k_id = 1:n())

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

    exp_n_pairs = nrow(treat) + nrow(unmatched_pairs)
    pts_matched = rbind(treat, control) %>%
      mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5,
             defor_10_5 = (cpc10_u - cpc5_u) / 5,
             defor_10_0 = (cpc10_u - cpc0_u) / 10,
             pair = pair_id)

    # Pair-level independent variables: median of all pixels in each pair (control + treat), then min/median/max across 100 pairs
    # elevation, slope, accessibility, cpc0/5/10_u, cpc0/5/10_d, defor_5_0 = cpc5_u - cpc0_u, defor_10_5 = cpc10_u - cpc5_u
    pair_var = pts_matched %>%
        reframe(elevation = median(elevation),
                slope = median(slope),
                accessibility = median(accessibility),
                cpc0_u = median(cpc0_u),
                cpc0_d = median(cpc0_d),
                cpc5_u = median(cpc5_u),
                cpc5_d = median(cpc5_d),
                cpc10_u = median(cpc10_u),
                cpc10_d = median(cpc10_d),
                defor_5_0 = median(defor_5_0),
                defor_10_5 = median(defor_10_5),
                defor_10_0 = median(defor_10_0)) %>%
        pivot_longer(cols = elevation:defor_10_0, names_to = "var", values_to = "val") %>%
        mutate(pair = pair_id)

    #calculate annual proportion of each land use class
    class_prefix = "JRC"
    match_years = c(0, -5, -10)
    match_classes = c(1, 3)
    luc_series = simulate_area_series(pts_matched, class_prefix = class_prefix, t0 = t0,
                                      match_years = match_years, match_classes = match_classes,
                                      exp_n_pairs, area_ha, verbose = F)

    #project pixel and matched pixel's pre-project and during-project LUC change
    lucc = luc_series$series %>%
        group_by(treatment, class) %>%
        reframe(year = year[-1],
                prop_df = -diff(class_prop)) %>%
        ungroup() %>%
        mutate(started = (year > t0)) %>%
        filter(class == "1") %>%
        dplyr::select(c("treatment", "prop_df", "started"))

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
    return(list(pair_var = pair_var, out_df = out_df, lucc = lucc, pts_matched = pts_matched, exp_n_pairs = exp_n_pairs))
}
