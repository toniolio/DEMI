# Minimal package installer for behavior-only scripts (00–03)
# Run: source("scripts/install_deps.R")

pkgs <- c(
  # core tidyverse
  "tidyverse", "readr", "dplyr", "tidyr", "purrr", "tibble", "stringr",
  # modeling
  "brms", "bayestestR",
  # optional (present in headers / legacy)
  "lme4", "lmerTest",
  # preprocessing utilities seen in your scripts
  "TSEntropies",
  # Stan interface used in 02_imagery_clean.R
  "cmdstanr"
)

# install missing
is_missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(is_missing)) {
  message("Installing: ", paste(is_missing, collapse=", "))
  install.packages(is_missing)
}

# cmdstanr from mc-stan repo if needed
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}

message("All set. You may still need to run cmdstanr::install_cmdstan() once.")
