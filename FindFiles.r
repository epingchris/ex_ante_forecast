#Wrapper function to search files or folders based on inclusion or exclusion keywords
FindFiles = function(dir, include = NULL, exclude = NULL, full = F) {
  files = list.files(dir, full = full)

  if (!is.null(include)) {
    include_pattern = paste(include, collapse = "|")
    files = files %>% str_subset(include_pattern)
  }

  if (!is.null(exclude)) {
    exclude_pattern = paste(exclude, collapse = "|")
    files = files %>% str_subset(exclude_pattern, negate = T)
  }

  if(length(files) == 0) {
    return(NA)
  } else {
    return(files)
  }
}