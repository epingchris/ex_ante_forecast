ReformatPixels = function(in_df, prefix, t0, treatment, pair) {
  in_df = in_df %>%
    as.data.frame()

  if(nchar(prefix) > 0) {
    in_df = in_df %>%
      dplyr::select(starts_with(prefix)) %>%
      rename_with(~str_replace(.x, prefix, ""))
  }

  luc_df = in_df %>%
    dplyr::select(!starts_with("cpc")) %>%
    st_as_sf(coords = c("lng", "lat")) %>%
    dplyr::rename(remoteness = access) %>%
    dplyr::rename_with(~ gsub("luc_", "JRC", .x, fixed = T)) %>%
    mutate(treatment = treatment, pair = pair)

  hfc_df = in_df %>%
    dplyr::select(starts_with("cpc"))
  colnames_hfc = colnames(hfc_df)
  years_bfr = colnames_hfc %>%
    stringr::str_extract("[:digit:]+") %>%
    as.numeric()
  suffixes = colnames_hfc %>%
    stringr::str_extract("_u|_d") %>%
    stringr::str_replace("u", "1") %>%
    stringr::str_replace("d", "3")
  colnames(hfc_df) = paste0("HFC", t0 - years_bfr, suffixes)

  out_df = cbind(luc_df, hfc_df)

  return(out_df)
}