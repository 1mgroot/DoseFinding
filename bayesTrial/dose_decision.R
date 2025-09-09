get_admissible_set <- function(posterior_summaries, config) {
  # Placeholder
  return(1:length(config$dose_levels))
}

adaptive_randomization <- function(admissible_set, posterior_summaries) {
  # Placeholder
  n_admissible <- length(admissible_set)
  if (n_admissible == 0) {
    return(numeric(0))
  } else {
    return(rep(1/n_admissible, n_admissible))
  }
}

select_final_od <- function(admissible_set, posterior_summaries, config) {
  # Placeholder
  if (length(admissible_set) > 0) {
    return(admissible_set[1])
  } else {
    return(NA)
  }
}