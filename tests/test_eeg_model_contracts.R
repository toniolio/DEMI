# Focused R contracts for DEMI EEG Stage 19.
#
# These tests validate exact formulas, random-effect structure, Student-t
# family, named fixed-effect mapping, T/A/B contrast algebra by two routes,
# deterministic synthetic generation, exact prior classes, synthetic-input
# outcome rejection, direct prior-predictive finiteness, and atomic directory
# publication. They use balanced synthetic fixtures only and never read
# accepted EEG outcomes or fit an accepted EEG model.
#
# Run from the repository root with:
# `Rscript --no-init-file tests/test_eeg_model_contracts.R`.

options(stringsAsFactors = FALSE)

repo_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo_root, "analysis/eeg_models/model_contracts.R"), local = FALSE)
stopifnot(requireNamespace("brms", quietly = TRUE))
stopifnot(requireNamespace("yaml", quietly = TRUE))

model_config <- yaml::read_yaml(file.path(repo_root, "analysis/eeg_models/model_config_v1.yaml"))
prior_config <- yaml::read_yaml(file.path(repo_root, "analysis/eeg_models/prior_config_v1.yaml"))

rating <- generate_rating_synthetic(model_config$computation$synthetic_rating_seed)
beta <- generate_beta_synthetic(model_config$computation$synthetic_beta_seed)
stopifnot(identical(rating$data, generate_rating_synthetic(model_config$computation$synthetic_rating_seed)$data))
stopifnot(identical(beta$data, generate_beta_synthetic(model_config$computation$synthetic_beta_seed)$data))

for (model_id in c("H1", "H2", "H3")) {
  fixture <- if (model_id %in% c("H1", "H2")) rating$data else beta$data
  matrix <- validate_formula_contract(model_id, fixture, model_config$models[[model_id]]$formula)
  stopifnot(identical(colnames(matrix), expected_fixed_terms(model_id)))
  family <- model_family()
  stopifnot(family$family == "student", family$link == "identity")
  prior <- make_model_priors(model_id, prior_config)
  stopifnot(identical(as.character(prior$class), c("Intercept", "b", "sd", "sigma", "nu")))
  validated <- brms::validate_prior(
    prior, formula = model_formulas(model_id)$full,
    data = fixture, family = family
  )
  stopifnot(inherits(validated, "brmsprior"))
  standata <- brms::make_standata(
    formula = model_formulas(model_id)$full, data = fixture,
    family = family, prior = prior, sample_prior = "yes", backend = "cmdstanr"
  )
  stopifnot(is.list(standata), length(standata) > 0L)
}

rating_algebra <- validate_estimand_algebra("H1", rating$fixed)
alpha_algebra <- validate_estimand_algebra("H2", rating$fixed)
beta_algebra <- validate_estimand_algebra("H3", beta$fixed)
stopifnot(max(rating_algebra$absolute_difference) <= 1e-12)
stopifnot(max(alpha_algebra$absolute_difference) <= 1e-12)
stopifnot(max(beta_algebra$absolute_difference) <= 1e-12)
stopifnot(identical(rating_algebra$estimand_id, c("T1", "T2", "T3")))
stopifnot(identical(alpha_algebra$estimand_id, c("A1", "A2", "A3")))
stopifnot(identical(beta_algebra$estimand_id, paste0("B", 1:6)))
stopifnot(max(abs(rating_algebra$direct_algebra - c(1.0, 1.6, 0.6))) <= 1e-12)
stopifnot(max(abs(beta_algebra$direct_algebra - c(1.0, 1.8, 0.8, 0.4, 0.8, 0.4))) <= 1e-12)

registry <- estimand_contrast_registry()
primary <- unique(registry$estimand_id[registry$classification == "primary"])
supplementary <- unique(registry$estimand_id[registry$classification == "supplementary"])
stopifnot(identical(primary, c("T1", "T2", "T3", "A1", "A2", "A3", "B1", "B2", "B3")))
stopifnot(identical(supplementary, c("B4", "B5", "B6")))

unsafe <- rating$data
unsafe$value_db <- 0
stopifnot(inherits(try(assert_synthetic_fit_input_safe(unsafe), silent = TRUE), "try-error"))
unsafe <- rating$data
unsafe$copied_value_log_power <- 0
stopifnot(inherits(try(assert_synthetic_fit_input_safe(unsafe), silent = TRUE), "try-error"))
stopifnot(isTRUE(assert_synthetic_fit_input_safe(rating$data)))

prior_small <- simulate_prior_predictive(
  "H1", rating$data[names(rating$data) != "value_synthetic"], prior_config,
  model_config$engineering_criteria, draws = 30, seed = 9919, figure_sample = 1000
)
stopifnot(prior_small$summary$finite_proportion == 1)
stopifnot(all(is.finite(prior_small$prediction_sample)))

temporary_parent <- tempfile("demi-stage19-atomic-")
dir.create(temporary_parent)
target <- file.path(temporary_parent, "model_validation_v1")
staging <- make_staging_directory(target)
writeLines("complete", file.path(staging, "proof.txt"))
publish_staging_directory(staging, target)
stopifnot(readLines(file.path(target, "proof.txt")) == "complete")

cat("PASS focused R EEG model contracts\n")
