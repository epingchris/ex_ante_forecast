RetrievePoints = function(matched_path, matchless_path, k, matches, t0) {
    k = k %>% dplyr::select(c("lat", "lng", "ecoregion"))
    matches = matches %>% dplyr::select(c("lat", "lng", "ecoregion"))

    pairs = read_parquet(matched_path) %>%
        dplyr::left_join(., k, by = join_by(k_lat == lat, k_lng == lng)) %>%
        rename(k_ecoregion = ecoregion) %>%
        dplyr::left_join(., matches, by = join_by(s_lat == lat, s_lng == lng)) %>%
        rename(s_ecoregion = ecoregion) %>%
        mutate(s_id = 1:n(), k_id = 1:n())

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

    pts_matched = rbind(treat, control)

    return(list(exp_n_pairs = exp_n_pairs, pts_matched = pts_matched))
}

ProcessPairs = function(matched_path, matchless_path, k, matches, t0, area_ha, acd, pair_id) { #loop through all sampled pairs
    c = Sys.time()

    points = RetrievePoints(matched_path, matchless_path, k, matches, t0)
    exp_n_pairs = points$exp_n_pairs
    pts_matched = points$pts_matched %>%
        mutate(defor_5_0 = (cpc5_u - cpc0_u) / 5, defor_10_5 = (cpc10_u - cpc5_u) / 5, defor_10_0 = (cpc10_u - cpc0_u) / 10)


    #project pixel and matched pixel's pre-project CPC change
    cpcc = list(treatment = pts_matched %>% filter(treatment == "treatment") %>% pull(defor_10_0),
                control = pts_matched %>% filter(treatment == "control") %>% pull(defor_10_0))

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

    pair_biome = pts_matched$biome

    #calculate annual proportion of each land use class
    luc_series = simulate_area_series(pts_matched,
                                      class_prefix, t0 = t0, match_years, match_classes,
                                      exp_n_pairs, area_ha,
                                      verbose = F)

    #project pixel and matched pixel's pre-project and during-project LUC change
    luc_diff = luc_series$series %>%
        group_by(treatment, class) %>%
        reframe(year = year[-1],
                prop_df = -diff(class_prop)) %>%
        ungroup() %>%
        mutate(started = (year > t0))

    lucc = list(treatment_pre = luc_diff %>% filter(treatment == "treatment" & !started & class == "1") %>% pull(prop_df),
                control_pre = luc_diff %>% filter(treatment == "control" & !started & class == "1") %>% pull(prop_df),
                treatment_during = luc_diff %>% filter(treatment == "treatment" & started & class == "1") %>% pull(prop_df),
                control_during = luc_diff %>% filter(treatment == "control" & started & class == "1") %>% pull(prop_df))

    #calculate annual carbon stock (Mg)
    carbon_series = luc_series$series %>%
        merge(., acd, by.x = "class", by.y = "land.use.class", all.x = T) %>%
        mutate(carbon_content = class_area * carbon.density) %>% #Mg
        group_by(treatment, year) %>%
        summarise(carbon_content = sum(carbon_content, na.rm = T)) %>%
        ungroup()

    #calculate annual carbon flux and additionality
    carbon_wide = pivot_wider(carbon_series, names_from = "treatment", values_from = "carbon_content")
    out_df = data.frame(year = carbon_wide$year[-1],
                        c_loss = -diff(carbon_wide$control),
                        t_loss = -diff(carbon_wide$treatment)) %>%
        mutate(additionality = c_loss - t_loss, pair = pair_id)
    d = Sys.time()
    cat(pair_id, ":", d - c, "\n")
    return(list(pair_var = pair_var, out_df = out_df, cpcc = cpcc, lucc = lucc, biome = pair_biome))
}