#function to assess matching balance
AssessBalance = function(pixels, t0) {
  match_years = c(0, -5, -10)
  match_classes = c(1, 3)
  prefix = "HFC"

  #generate the match variables
  match_var_grid = expand.grid(prefix = prefix,
                               years = t0 + match_years,
                               #the years to match on match_years should be zero or negative
                               classes = match_classes,
                               suffix = "",
                               stringsAsFactors = F)

  match_vars = apply(match_var_grid, 1, function(x) {
      paste0(x["prefix"], x["years"], "_", x["classes"], x["suffix"])
    }) %>%
    c("remoteness", "elevation", "slope") %>%
    gsub("\\s+", "", .) %>%
    gsub("\\-", "\\.", .)

  #generate match formula
  fmla = as.formula(paste("treatment ~ ", paste(match_vars, collapse = "+")))

  #match
  match_out = matchit(
    fmla,
    method = "nearest",
    distance = "mahalanobis",
    ratio = 1,
    order = "smallest",
    replace = FALSE,
    discard = "none",
    data = pixels %>%
      mutate(treatment = ifelse(treatment == "project", 1, 0))
  )
  match_summary = summary(match_out)
  balance_test = all(abs(match_summary$sum.matched[, "Std. Mean Diff."]) <= 0.2)

  out = list(out = match_out, summary = match_summary, balance_test = balance_test)
  return(out)
}