#!/usr/bin/env Rscript
# Validate DEMI EEG model formulas, estimands, priors, and synthetic recovery.
#
# This Stage 19 driver proves the accepted H1--H3 implementation using only
# deterministic synthetic outcomes. It verifies Stage 18 read-only authority,
# obtains pooled predictor-blind outcome scales through a separate Python
# projection, loads outcome-free accepted designs, validates formula/design
# terms and estimand algebra, constructs and compiles Stan programs, runs two
# bounded synthetic Bayesian recovery fits, and simulates pooled prior-only
# predictions for H1, H2, and H3. It publishes a complete ignored output
# directory atomically and supports hash-validated unchanged reuse.
#
# Inputs are tracked model/prior configuration, tracked helper code, and the
# accepted Stage 18 namespace. Outputs are confined to
# `_Data/eeg/model_validation_v1/`. The driver never fits `value_db` or
# `value_log_power`, never joins an EEG outcome to a scientific predictor,
# never calculates an accepted-data estimand, and never runs sensitivities or
# scientific interpretation.
#
# Run from the repository root with:
# `Rscript --no-init-file analysis/eeg_models/19_validate_synthetic_models_and_priors.R`
# or use `tools/run_eeg_model_validation_v1.sh`. Pass `--verify-current` to
# reopen and hash-check an already validated namespace without rewriting it.

options(stringsAsFactors = FALSE, warn = 1, brms.backend = "cmdstanr")

script_started <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
verify_only <- identical(args, "--verify-current")
if (length(args) > 1L || (length(args) == 1L && !verify_only)) {
  stop("usage: 19_validate_synthetic_models_and_priors.R [--verify-current]", call. = FALSE)
}

repo_root <- normalizePath(getwd(), mustWork = TRUE)
model_dir <- file.path(repo_root, "analysis/eeg_models")
source(file.path(model_dir, "model_contracts.R"), local = FALSE)

required_packages <- c("brms", "cmdstanr", "posterior", "bayesplot", "loo", "yaml", "jsonlite", "digest")
for (package in required_packages) {
  require_contract(requireNamespace(package, quietly = TRUE), paste0("missing_required_R_package:", package))
}

model_config_path <- file.path(model_dir, "model_config_v1.yaml")
prior_config_path <- file.path(model_dir, "prior_config_v1.yaml")
model_config <- yaml::read_yaml(model_config_path)
prior_config <- yaml::read_yaml(prior_config_path)
output_root <- file.path(repo_root, model_config$storage$output_root)
python <- file.path(repo_root, ".venv/bin/python")
input_script <- file.path(model_dir, "prepare_model_validation_inputs.py")
require_contract(file.exists(python), "locked_python_environment_missing")

#' Return the SHA-256 digest of a file.
#'
#' @param path Existing file.
#' @return Lowercase hexadecimal SHA-256 digest.
sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

#' Write pretty deterministic JSON.
#'
#' @param value Serializable R object.
#' @param path Output path.
#' @return Invisibly returns the path.
write_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null")
  invisible(path)
}

#' Execute a subprocess and fail with captured output when unsuccessful.
#'
#' @param command Executable path/name.
#' @param arguments Character arguments.
#' @return Captured output lines.
run_checked <- function(command, arguments) {
  output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L) {
    stop(paste(c(paste0("subprocess_failed:", command), output), collapse = "\n"), call. = FALSE)
  }
  output
}

tracked_authority_files <- c(
  model_config_path,
  prior_config_path,
  file.path(model_dir, "model_contracts.R"),
  input_script,
  file.path(model_dir, "19_validate_synthetic_models_and_priors.R"),
  file.path(repo_root, "tools/run_eeg_model_validation_v1.sh")
)
require_contract(all(file.exists(tracked_authority_files)), "tracked_stage19_authority_file_missing")
tracked_hashes <- stats::setNames(vapply(tracked_authority_files, sha256_file, character(1)),
                                  substring(tracked_authority_files, nchar(repo_root) + 2L))
authority_fingerprint <- digest::digest(
  paste(names(tracked_hashes), tracked_hashes, collapse = "\n"),
  algo = "sha256", serialize = FALSE
)

#' Hash-validate every output declared by a completed manifest.
#'
#' @param root Published validation root.
#' @param manifest Parsed manifest object.
#' @return Invisibly returns TRUE.
verify_manifest_outputs <- function(root, manifest) {
  for (descriptor in manifest$outputs) {
    path <- file.path(root, descriptor$relative_path)
    require_contract(file.exists(path), paste0("validation_output_missing:", descriptor$relative_path))
    require_contract(file.info(path)$size == as.numeric(descriptor$size_bytes),
                     paste0("validation_output_size_changed:", descriptor$relative_path))
    require_contract(sha256_file(path) == descriptor$sha256,
                     paste0("validation_output_hash_changed:", descriptor$relative_path))
  }
  invisible(TRUE)
}

#' Reopen a completed namespace without writing it.
#'
#' @param root Published validation root.
#' @return Invisibly returns manifest after all checks pass.
verify_current <- function(root) {
  manifest_path <- file.path(root, "model_validation_manifest.json")
  validation_path <- file.path(root, "validation/model_validation.json")
  require_contract(file.exists(manifest_path) && file.exists(validation_path), "model_validation_current_missing")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  validation <- jsonlite::read_json(validation_path, simplifyVector = TRUE)
  require_contract(manifest$status == "complete", "model_validation_manifest_not_complete")
  require_contract(validation$status == "pass", "model_validation_not_pass")
  require_contract(manifest$authority_fingerprint == authority_fingerprint,
                   "model_validation_authority_fingerprint_changed")
  verify_manifest_outputs(root, manifest)
  stage18_output <- run_checked(
    python,
    c(file.path(repo_root, "analysis/eeg_mne/18_construct_model_ready_tables.py"), "--verify-current")
  )
  require_contract(any(grepl("PASS model_tables_v1 unchanged reuse", stage18_output, fixed = TRUE)),
                   "stage18_reuse_verification_failed")
  message(sprintf("PASS model_validation_v1 unchanged reuse outputs=%d", length(manifest$outputs)))
  invisible(manifest)
}

if (verify_only) {
  verify_current(output_root)
  quit(save = "no", status = 0)
}

if (dir.exists(output_root)) {
  verify_current(output_root)
  quit(save = "no", status = 0)
}

staging <- make_staging_directory(output_root)
successful_publication <- FALSE
on.exit({
  if (!successful_publication && dir.exists(staging)) {
    message("Preserved incomplete Stage 19 staging directory: ", staging)
  }
}, add = TRUE)

for (directory in c(
  "configuration", "registries", "synthetic", "stan", "compilation", "fits",
  "prior_predictive", "figures", "validation", "logs"
)) {
  dir.create(file.path(staging, directory), recursive = TRUE, showWarnings = FALSE)
}
file.copy(model_config_path, file.path(staging, "configuration/model_config_v1.yaml"))
file.copy(prior_config_path, file.path(staging, "configuration/prior_config_v1.yaml"))

preflight_started <- Sys.time()
prepare_output <- run_checked(
  python,
  c(input_script, "prepare", "--output-root", staging)
)
preflight <- jsonlite::read_json(file.path(staging, "validation/preflight.json"), simplifyVector = TRUE)
require_contract(preflight$status == "pass", "stage19_preflight_failed")
require_contract(
  preflight$stage18_authority_fingerprint == model_config$authority$stage18_authority_fingerprint,
  "stage18_authority_disagrees_with_stage19_config"
)

scale_snapshot <- jsonlite::read_json(file.path(staging, "prior_scale_snapshot.json"), simplifyVector = TRUE)
for (model_id in c("H1", "H2", "H3")) {
  observed <- scale_snapshot$models[[model_id]]
  frozen <- prior_config$models[[model_id]]
  require_contract(identical(observed$columns_read, "value_db"), paste0(model_id, "_scale_path_not_predictor_blind"))
  require_contract(as.integer(observed$pooled_rows) == as.integer(frozen$pooled_rows), paste0(model_id, "_scale_rows_changed"))
  for (field in c("m_y", "mad_raw", "s_y", "intercept_scale")) {
    require_contract(abs(as.numeric(observed[[field]]) - as.numeric(frozen[[field]])) <= 1e-12,
                     paste0(model_id, "_frozen_prior_scale_changed:", field))
  }
}

r_versions <- stats::setNames(vapply(required_packages, function(package) {
  as.character(utils::packageVersion(package))
}, character(1)), required_packages)
lock <- jsonlite::read_json(file.path(repo_root, "renv.lock"), simplifyVector = FALSE)
locked_version <- function(package) {
  entry <- lock$Packages[[package]]
  if (is.null(entry)) NA_character_ else entry$Version
}
lock_stack <- c("brms", "cmdstanr", "posterior", "bayesplot", "loo", "yaml", "arrow", "testthat", "lme4")
lock_status <- data.frame(
  package = lock_stack,
  locked_version = vapply(lock_stack, locked_version, character(1)),
  installed_version = vapply(lock_stack, function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
lock_status$reproducibility_status <- ifelse(
  is.na(lock_status$locked_version), "missing_from_lock",
  ifelse(lock_status$installed_version == lock_status$locked_version, "locked_and_available", "locked_version_not_active")
)
utils::write.csv(lock_status, file.path(staging, "validation/R_lock_preflight.csv"), row.names = FALSE)
toolchain_output <- capture.output(cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = FALSE))
cmdstan_version <- as.character(cmdstanr::cmdstan_version())
write_json(list(
  status = "pass",
  R_version = as.character(getRversion()),
  installed_primary_stack = as.list(r_versions),
  cmdstan_path = cmdstanr::cmdstan_path(),
  cmdstan_version = cmdstan_version,
  toolchain_usable = TRUE,
  toolchain_output = toolchain_output,
  missing_lock_dependencies = as.list(lock_status$package[lock_status$reproducibility_status == "missing_from_lock"]),
  python_prepare_output = prepare_output,
  runtime_seconds = as.numeric(difftime(Sys.time(), preflight_started, units = "secs"))
), file.path(staging, "validation/toolchain_preflight.json"))

designs <- lapply(c("H1", "H2", "H3"), function(model_id) {
  data <- utils::read.csv(file.path(staging, "inputs", paste0(model_id, "_accepted_design_no_outcomes.csv")),
                          check.names = FALSE)
  data$participant_id <- as.factor(data$participant_id)
  data$canonical_event_key <- as.factor(data$canonical_event_key)
  assert_synthetic_fit_input_safe(data, require_response = FALSE)
  data
})
names(designs) <- c("H1", "H2", "H3")

rating <- generate_rating_synthetic(as.integer(model_config$computation$synthetic_rating_seed))
beta <- generate_beta_synthetic(as.integer(model_config$computation$synthetic_beta_seed))
rating_repeat <- generate_rating_synthetic(as.integer(model_config$computation$synthetic_rating_seed))
beta_repeat <- generate_beta_synthetic(as.integer(model_config$computation$synthetic_beta_seed))
require_contract(identical(rating$data, rating_repeat$data), "rating_synthetic_generation_not_deterministic")
require_contract(identical(beta$data, beta_repeat$data), "beta_synthetic_generation_not_deterministic")

utils::write.csv(rating$data, file.path(staging, "synthetic/rating_fixture.csv"), row.names = FALSE)
utils::write.csv(beta$data, file.path(staging, "synthetic/beta_fixture.csv"), row.names = FALSE)
write_json(rating$metadata, file.path(staging, "synthetic/rating_generation_metadata.json"))
write_json(beta$metadata, file.path(staging, "synthetic/beta_generation_metadata.json"))

formula_rows <- list()
algebra_rows <- list()
prior_rows <- list()
for (model_id in c("H1", "H2", "H3")) {
  fixture <- if (model_id %in% c("H1", "H2")) rating$data else beta$data
  exact_formula <- model_config$models[[model_id]]$formula
  matrix <- validate_formula_contract(model_id, fixture, exact_formula)
  formula_rows[[model_id]] <- data.frame(
    model_id = model_id,
    formula = exact_formula,
    future_observed_outcome = model_config$models[[model_id]]$future_observed_outcome,
    family = "student", link = "identity", fixed_columns = ncol(matrix),
    participant_random_correlations = FALSE,
    participant_condition_slope = FALSE,
    recording_random_intercept = FALSE,
    canonical_trial_random_intercept = model_id == "H3",
    contains_accuracy_rating = model_id != "H3"
  )
  fixed <- if (model_id %in% c("H1", "H2")) rating$fixed else beta$fixed
  algebra_rows[[model_id]] <- validate_estimand_algebra(model_id, fixed)
  prior <- make_model_priors(model_id, prior_config)
  brms::validate_prior(
    prior = prior, formula = model_formulas(model_id)$full,
    data = fixture, family = model_family()
  )
  prior_rows[[model_id]] <- prior_mapping_registry(model_id, prior_config)
}
formula_registry <- do.call(rbind, formula_rows)
algebra_registry <- do.call(rbind, algebra_rows)
prior_registry <- do.call(rbind, prior_rows)
contrast_registry <- estimand_contrast_registry()
utils::write.csv(formula_registry, file.path(staging, "registries/formula_registry.csv"), row.names = FALSE)
utils::write.csv(contrast_registry, file.path(staging, "registries/estimand_contrast_registry.csv"), row.names = FALSE)
utils::write.csv(algebra_registry, file.path(staging, "validation/estimand_algebra_validation.csv"), row.names = FALSE)
utils::write.csv(prior_registry, file.path(staging, "registries/prior_mapping_registry.csv"), row.names = FALSE)

truth_rows <- rbind(
  data.frame(model_id = "H1", estimand_id = names(rating$estimands), classification = "primary", truth = as.numeric(rating$estimands)),
  data.frame(model_id = "H2", estimand_id = sub("T", "A", names(rating$estimands)), classification = "primary", truth = as.numeric(rating$estimands)),
  data.frame(model_id = "H3", estimand_id = names(beta$estimands),
             classification = ifelse(names(beta$estimands) %in% paste0("B", 1:3), "primary", "supplementary"),
             truth = as.numeric(beta$estimands))
)
utils::write.csv(truth_rows, file.path(staging, "registries/synthetic_truth_registry.csv"), row.names = FALSE)
fixed_truth <- rbind(
  data.frame(structure = "rating", fixed_term = names(rating$fixed), truth = as.numeric(rating$fixed)),
  data.frame(structure = "beta", fixed_term = names(beta$fixed), truth = as.numeric(beta$fixed))
)
utils::write.csv(fixed_truth, file.path(staging, "registries/synthetic_fixed_effect_truth.csv"), row.names = FALSE)

compile_rows <- list()
stan_started <- Sys.time()
for (model_id in c("H1", "H2", "H3")) {
  fixture <- if (model_id %in% c("H1", "H2")) rating$data else beta$data
  prior <- make_model_priors(model_id, prior_config)
  code <- brms::make_stancode(
    formula = model_formulas(model_id)$full, data = fixture,
    family = model_family(), prior = prior, sample_prior = "yes",
    backend = "cmdstanr"
  )
  standata <- brms::make_standata(
    formula = model_formulas(model_id)$full, data = fixture,
    family = model_family(), prior = prior, sample_prior = "yes",
    backend = "cmdstanr"
  )
  stan_path <- file.path(staging, "stan", paste0(model_id, ".stan"))
  writeLines(code, stan_path, useBytes = TRUE)
  utils::write.csv(standata_summary(standata),
                   file.path(staging, "stan", paste0(model_id, "_standata_summary.csv")), row.names = FALSE)
  model_started <- Sys.time()
  compiled <- cmdstanr::cmdstan_model(stan_path, compile = TRUE, force_recompile = FALSE, quiet = TRUE)
  compile_rows[[model_id]] <- data.frame(
    model_id = model_id, stan_sha256 = sha256_file(stan_path),
    executable_relative_path = substring(compiled$exe_file(), nchar(staging) + 2L),
    executable_sha256 = sha256_file(compiled$exe_file()),
    cmdstan_version = cmdstan_version,
    compile_runtime_seconds = as.numeric(difftime(Sys.time(), model_started, units = "secs")),
    status = "pass"
  )
}
compilation_registry <- do.call(rbind, compile_rows)
utils::write.csv(compilation_registry, file.path(staging, "compilation/compilation_metadata.csv"), row.names = FALSE)

criteria <- model_config$engineering_criteria
write_json(criteria, file.path(staging, "validation/predeclared_engineering_criteria.json"))

# Sampling uses synthetic fixtures only. The safety assertion is repeated at
# the final call site so no accepted outcome can enter brms by later refactor.
fit_structure <- function(model_id, data, seed) {
  assert_synthetic_fit_input_safe(data)
  brms::brm(
    formula = model_formulas(model_id)$full,
    data = data,
    family = model_family(),
    prior = make_model_priors(model_id, prior_config),
    backend = "cmdstanr",
    chains = as.integer(model_config$computation$chains),
    cores = as.integer(model_config$computation$parallel_chains),
    iter = as.integer(model_config$computation$warmup) + as.integer(model_config$computation$sampling),
    warmup = as.integer(model_config$computation$warmup),
    seed = as.integer(seed),
    control = list(
      adapt_delta = as.numeric(model_config$computation$adapt_delta),
      max_treedepth = as.integer(model_config$computation$max_treedepth)
    ),
    refresh = 100,
    silent = 2
  )
}

rating_fit_started <- Sys.time()
rating_fit <- fit_structure("H1", rating$data, model_config$computation$rating_fit_seed)
rating_result <- summarize_synthetic_fit(
  rating_fit, "H1", rating$estimands, criteria,
  as.integer(model_config$computation$max_treedepth)
)
rating_result$diagnostics$runtime_seconds <- as.numeric(difftime(Sys.time(), rating_fit_started, units = "secs"))
utils::write.csv(rating_result$recovery, file.path(staging, "fits/rating_synthetic_recovery.csv"), row.names = FALSE)
write_json(rating_result$diagnostics, file.path(staging, "fits/rating_synthetic_diagnostics.json"))
rm(rating_fit)
invisible(gc())

beta_fit_started <- Sys.time()
beta_fit <- fit_structure("H3", beta$data, model_config$computation$beta_fit_seed)
beta_result <- summarize_synthetic_fit(
  beta_fit, "H3", beta$estimands, criteria,
  as.integer(model_config$computation$max_treedepth)
)
beta_result$diagnostics$runtime_seconds <- as.numeric(difftime(Sys.time(), beta_fit_started, units = "secs"))
utils::write.csv(beta_result$recovery, file.path(staging, "fits/beta_synthetic_recovery.csv"), row.names = FALSE)
write_json(beta_result$diagnostics, file.path(staging, "fits/beta_synthetic_diagnostics.json"))
rm(beta_fit)
invisible(gc())

prior_summaries <- list()
prior_started <- Sys.time()
for (model_id in c("H1", "H2", "H3")) {
  simulation <- simulate_prior_predictive(
    model_id = model_id,
    design = designs[[model_id]],
    prior_config = prior_config,
    criteria = criteria,
    draws = as.integer(model_config$computation$prior_predictive_draws),
    seed = as.integer(model_config$computation$prior_predictive_seeds[[model_id]]),
    figure_sample = as.integer(model_config$computation$prior_predictive_figure_sample)
  )
  prior_summaries[[model_id]] <- simulation$summary
  utils::write.csv(simulation$parameter_draws,
                   file.path(staging, "prior_predictive", paste0(model_id, "_prior_parameter_draws.csv")),
                   row.names = FALSE)
  parameter_summary <- aggregate(value ~ parameter, simulation$parameter_draws, function(x) {
    c(mean = mean(x), sd = stats::sd(x), q0.01 = stats::quantile(x, 0.01),
      median = stats::median(x), q0.99 = stats::quantile(x, 0.99))
  })
  summary_values <- parameter_summary$value
  if (!is.matrix(summary_values)) summary_values <- do.call(rbind, summary_values)
  parameter_summary <- cbind(parameter = parameter_summary$parameter, as.data.frame(summary_values))
  utils::write.csv(parameter_summary,
                   file.path(staging, "prior_predictive", paste0(model_id, "_prior_parameter_summary.csv")),
                   row.names = FALSE)
  write_prior_predictive_figure(
    file.path(staging, "figures", paste0(model_id, "_pooled_prior_predictive.png")),
    simulation$prediction_sample, model_id,
    as.numeric(prior_config$models[[model_id]]$m_y),
    as.numeric(prior_config$models[[model_id]]$s_y)
  )
  rm(simulation)
  invisible(gc())
}
prior_summary <- do.call(rbind, prior_summaries)
utils::write.csv(prior_summary, file.path(staging, "prior_predictive/prior_predictive_summary.csv"), row.names = FALSE)

source_verify_output <- run_checked(
  python,
  c(input_script, "verify", "--before", file.path(staging, "source_snapshot_before.json"),
    "--evidence", file.path(staging, "source_immutability.json"))
)
source_evidence <- jsonlite::read_json(file.path(staging, "source_immutability.json"), simplifyVector = TRUE)
require_contract(source_evidence$status == "pass", "accepted_sources_changed_during_stage19")

primary_ids <- c("T1", "T2", "T3", "A1", "A2", "A3", "B1", "B2", "B3")
supplementary_ids <- c("B4", "B5", "B6")
classification_pass <- identical(
  unique(contrast_registry$estimand_id[contrast_registry$classification == "primary"]), primary_ids
) && identical(
  unique(contrast_registry$estimand_id[contrast_registry$classification == "supplementary"]), supplementary_ids
)
prior_pass <- all(prior_summary$status == "pass")
overall_pass <- rating_result$pass && beta_result$pass && prior_pass && classification_pass
validation <- list(
  status = if (overall_pass) "pass" else "owner_review_needed",
  stage18_preflight = "pass",
  stage18_authority_fingerprint = model_config$authority$stage18_authority_fingerprint,
  formula_and_design_matrix_validation = "pass",
  student_t_identity_family = TRUE,
  independent_participant_random_effects = TRUE,
  participant_condition_slope_absent = TRUE,
  recording_random_intercept_absent = TRUE,
  H3_trial_random_intercept_present = TRUE,
  H3_accuracy_rating_absent = TRUE,
  estimand_algebra_two_route_validation = "pass",
  primary_estimands = as.list(primary_ids),
  supplementary_estimands = as.list(supplementary_ids),
  classification_validation = classification_pass,
  prior_scale_predictor_blind = TRUE,
  prior_mapping_validation = "pass",
  compiled_models = as.list(compilation_registry$model_id),
  synthetic_rating_recovery = rating_result$diagnostics$status,
  synthetic_beta_recovery = beta_result$diagnostics$status,
  prior_predictive_status = stats::setNames(as.list(prior_summary$status), prior_summary$model_id),
  prior_owner_review_needed = !prior_pass,
  accepted_EEG_models_fitted = FALSE,
  scientific_estimands_calculated = FALSE,
  accepted_outcome_predictor_associations_opened = FALSE,
  accepted_observed_posterior_draws_written = FALSE,
  source_immutability = source_evidence$status,
  atomic_directory_publication = TRUE,
  validated_unchanged_reuse_supported = TRUE,
  repository_output_safety = "pass"
)
write_json(validation, file.path(staging, "validation/model_validation.json"))
require_contract(overall_pass, "stage19_validation_requires_owner_review_or_failed")

runtime_seconds <- as.numeric(difftime(Sys.time(), script_started, units = "secs"))
summary_lines <- c(
  "# DEMI EEG synthetic model/prior validation v1",
  "",
  "Status: pass.",
  "",
  "The exact H1--H3 hierarchical Student-t formulas, fixed-effect term mapping,",
  "and T1--T3/A1--A3/B1--B6 algebra passed direct and brute-prediction checks.",
  "All three Stan programs compiled with CmdStan; bounded rating and beta fits",
  "recovered injected synthetic estimands under the predeclared diagnostics.",
  "Model-specific frozen priors produced finite, predictor-collapsed prior",
  "predictions within the predeclared scale criteria.",
  "",
  "No accepted EEG model was fitted, no accepted EEG association was opened,",
  "and no scientific estimate exists. Observed-data fitting remains closed.",
  "",
  sprintf("Runtime: %.3f seconds.", runtime_seconds)
)
writeLines(summary_lines, file.path(staging, "model_validation_summary.md"))
writeLines(c(
  paste("authority_fingerprint", authority_fingerprint),
  paste("cmdstan_version", cmdstan_version),
  paste("rating_fit_status", rating_result$diagnostics$status),
  paste("beta_fit_status", beta_result$diagnostics$status),
  paste("prior_status", paste(prior_summary$model_id, prior_summary$status, collapse = ",")),
  paste("source_verify", paste(source_verify_output, collapse = " "))
), file.path(staging, "logs/production.log"))
write_json(list(
  schema_version = 1,
  pipeline_version = "model_validation_v1",
  status = "complete",
  authority_fingerprint = authority_fingerprint,
  stage18_authority_fingerprint = model_config$authority$stage18_authority_fingerprint,
  started_at = format(script_started, tz = "UTC", usetz = TRUE),
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  runtime_seconds = runtime_seconds,
  accepted_EEG_models_fitted = FALSE
), file.path(staging, "model_validation_run_state.json"))

output_paths <- list.files(staging, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                           no.. = TRUE)
output_paths <- output_paths[file.info(output_paths)$isdir %in% FALSE]
output_paths <- output_paths[!basename(output_paths) %in% c("model_validation_manifest.json", ".DS_Store")]
descriptors <- lapply(sort(output_paths), function(path) {
  list(
    relative_path = substring(path, nchar(staging) + 2L),
    size_bytes = unname(file.info(path)$size),
    sha256 = sha256_file(path)
  )
})
manifest <- list(
  schema_version = 1,
  pipeline_version = "model_validation_v1",
  status = "complete",
  authority_fingerprint = authority_fingerprint,
  tracked_authority_sha256 = as.list(tracked_hashes),
  stage18_authority_fingerprint = model_config$authority$stage18_authority_fingerprint,
  output_count = length(descriptors),
  declared_output_bytes = sum(vapply(descriptors, function(x) as.numeric(x$size_bytes), numeric(1))),
  outputs = descriptors
)
write_json(manifest, file.path(staging, "model_validation_manifest.json"))

# Reopen every staged artifact before the single directory rename.
manifest_reopened <- jsonlite::read_json(file.path(staging, "model_validation_manifest.json"), simplifyVector = FALSE)
verify_manifest_outputs(staging, manifest_reopened)
publish_staging_directory(staging, output_root)
successful_publication <- TRUE

# Reopen the published namespace and exercise the exact unchanged-reuse route.
key_paths <- c(
  file.path(output_root, "model_validation_manifest.json"),
  file.path(output_root, "validation/model_validation.json"),
  file.path(output_root, "prior_scale_snapshot.json")
)
mtimes_before_reuse <- file.info(key_paths)$mtime
verify_current(output_root)
mtimes_after_reuse <- file.info(key_paths)$mtime
require_contract(identical(mtimes_before_reuse, mtimes_after_reuse), "validated_reuse_rewrote_outputs")
message(sprintf(
  "PASS model_validation_v1 runtime=%.3fs bytes=%d outputs=%d",
  runtime_seconds, manifest$declared_output_bytes, manifest$output_count
))
