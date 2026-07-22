# DEMI EEG model formula, estimand, prior, and synthetic-validation helpers.
#
# This module implements the Stage 19 engineering contract without fitting any
# accepted EEG response. It owns the exact H1--H3 brms formulas, Student-t
# family, scale-aware prior mappings, fixed-effect contrast algebra, balanced
# synthetic generators, direct prior-predictive simulation, sampling
# diagnostics, and atomic-output helpers. Inputs to sampling must contain a
# generated `value_synthetic` column and must not contain accepted EEG outcome
# columns. The module reads no Parquet and writes nothing on import.
#
# Expected callers are `19_validate_synthetic_models_and_priors.R` and the
# focused R contract tests. Generated outputs are confined to the ignored
# `_Data/eeg/model_validation_v1/` namespace by the driver.

#' Stop with a stable contract label when a condition is false.
#'
#' @param condition Logical scalar.
#' @param label Stable failure label.
#' @return Invisibly returns TRUE; otherwise stops.
require_contract <- function(condition, label) {
  if (!isTRUE(condition)) stop(label, call. = FALSE)
  invisible(TRUE)
}

#' Return exact full and fixed formulas for one model.
#'
#' @param model_id One of H1, H2, or H3.
#' @return A list containing `full` and `fixed` formula objects.
model_formulas <- function(model_id) {
  rating_full <- paste0(
    "value_synthetic ~ condition_c * accuracy_rating_within + ",
    "accuracy_rating_between + familiarity_c + ",
    "(1 + accuracy_rating_within + familiarity_c || participant_id)"
  )
  rating_fixed <- "~ condition_c * accuracy_rating_within + accuracy_rating_between + familiarity_c"
  beta_full <- paste0(
    "value_synthetic ~ condition_c * phase_c * familiarity_c + ",
    "(1 + phase_c * familiarity_c || participant_id) + ",
    "(1 | canonical_event_key)"
  )
  beta_fixed <- "~ condition_c * phase_c * familiarity_c"
  require_contract(model_id %in% c("H1", "H2", "H3"), "unknown_model_id")
  if (model_id %in% c("H1", "H2")) {
    return(list(full = stats::as.formula(rating_full), fixed = stats::as.formula(rating_fixed)))
  }
  list(full = stats::as.formula(beta_full), fixed = stats::as.formula(beta_fixed))
}

#' Return the accepted Student-t family with identity link.
#'
#' @return A brms family object.
model_family <- function() {
  brms::student(link = "identity")
}

#' Return exact expected fixed-effect term identities.
#'
#' @param model_id One of H1, H2, or H3.
#' @return Character vector in R model-matrix order.
expected_fixed_terms <- function(model_id) {
  if (model_id %in% c("H1", "H2")) {
    return(c(
      "(Intercept)", "condition_c", "accuracy_rating_within",
      "accuracy_rating_between", "familiarity_c",
      "condition_c:accuracy_rating_within"
    ))
  }
  require_contract(model_id == "H3", "unknown_model_id")
  c(
    "(Intercept)", "condition_c", "phase_c", "familiarity_c",
    "condition_c:phase_c", "condition_c:familiarity_c",
    "phase_c:familiarity_c", "condition_c:phase_c:familiarity_c"
  )
}

#' Validate formula text, fixed terms, and prohibited structures.
#'
#' @param model_id One of H1, H2, or H3.
#' @param data A safe design frame containing needed columns.
#' @param expected_formula Exact formula string from tracked configuration.
#' @return The fixed-effect model matrix.
validate_formula_contract <- function(model_id, data, expected_formula) {
  formulas <- model_formulas(model_id)
  observed_formula <- paste(deparse(formulas$full, width.cutoff = 500L), collapse = "")
  require_contract(identical(observed_formula, expected_formula), paste0(model_id, "_formula_changed"))
  matrix <- stats::model.matrix(formulas$fixed, data = data)
  require_contract(
    identical(colnames(matrix), expected_fixed_terms(model_id)),
    paste0(model_id, "_fixed_term_mapping_changed")
  )
  require_contract(grepl("\\|\\| participant_id", observed_formula), paste0(model_id, "_participant_effects_not_independent"))
  require_contract(!grepl("recording", observed_formula, ignore.case = TRUE), paste0(model_id, "_recording_random_effect_present"))
  require_contract(!grepl("condition_c \\|\\|", observed_formula), paste0(model_id, "_participant_condition_slope_present"))
  if (model_id == "H3") {
    require_contract(!grepl("accuracy_rating", observed_formula), "H3_accuracy_rating_term_present")
    require_contract(grepl("\\(1 \\| canonical_event_key\\)", observed_formula), "H3_trial_intercept_missing")
  } else {
    require_contract(!grepl("canonical_event_key", observed_formula), paste0(model_id, "_trial_intercept_present"))
  }
  matrix
}

#' Construct exact model-specific brms priors.
#'
#' @param model_id One of H1, H2, or H3.
#' @param prior_config Parsed tracked prior configuration.
#' @return A brms prior object with class-level mappings.
make_model_priors <- function(model_id, prior_config) {
  constants <- prior_config$models[[model_id]]
  require_contract(!is.null(constants), "prior_model_missing")
  m_y <- as.numeric(constants$m_y)
  s_y <- as.numeric(constants$s_y)
  c(
    brms::set_prior(sprintf("student_t(3, %.17g, %.17g)", m_y, 2.5 * s_y), class = "Intercept"),
    brms::set_prior(sprintf("normal(0, %.17g)", s_y), class = "b"),
    brms::set_prior(sprintf("student_t(3, 0, %.17g)", s_y), class = "sd"),
    brms::set_prior(sprintf("student_t(3, 0, %.17g)", s_y), class = "sigma"),
    brms::set_prior("gamma(2, 0.1)", class = "nu")
  )
}

#' Expand prior classes to the parameters they govern.
#'
#' @param model_id One of H1, H2, or H3.
#' @param prior_config Parsed tracked prior configuration.
#' @return Data frame recording exact distribution-to-parameter mappings.
prior_mapping_registry <- function(model_id, prior_config) {
  constants <- prior_config$models[[model_id]]
  s_y <- as.numeric(constants$s_y)
  rows <- list(data.frame(
    model_id = model_id, parameter = "(Intercept)", brms_class = "Intercept",
    group = NA_character_, distribution = "student_t", arg1 = 3,
    arg2 = as.numeric(constants$m_y), arg3 = 2.5 * s_y,
    lower_bound = NA_real_, directional = FALSE
  ))
  fixed <- setdiff(expected_fixed_terms(model_id), "(Intercept)")
  rows[[length(rows) + 1L]] <- data.frame(
    model_id = model_id, parameter = fixed, brms_class = "b",
    group = NA_character_, distribution = "normal", arg1 = 0,
    arg2 = s_y, arg3 = NA_real_, lower_bound = NA_real_, directional = FALSE
  )
  participant_terms <- if (model_id %in% c("H1", "H2")) {
    c("(Intercept)", "accuracy_rating_within", "familiarity_c")
  } else {
    c("(Intercept)", "phase_c", "familiarity_c", "phase_c:familiarity_c")
  }
  rows[[length(rows) + 1L]] <- data.frame(
    model_id = model_id, parameter = participant_terms, brms_class = "sd",
    group = "participant_id", distribution = "student_t", arg1 = 3,
    arg2 = 0, arg3 = s_y, lower_bound = 0, directional = FALSE
  )
  if (model_id == "H3") {
    rows[[length(rows) + 1L]] <- data.frame(
      model_id = model_id, parameter = "(Intercept)", brms_class = "sd",
      group = "canonical_event_key", distribution = "student_t", arg1 = 3,
      arg2 = 0, arg3 = s_y, lower_bound = 0, directional = FALSE
    )
  }
  rows[[length(rows) + 1L]] <- data.frame(
    model_id = model_id, parameter = "sigma", brms_class = "sigma",
    group = NA_character_, distribution = "student_t", arg1 = 3,
    arg2 = 0, arg3 = s_y, lower_bound = 0, directional = FALSE
  )
  rows[[length(rows) + 1L]] <- data.frame(
    model_id = model_id, parameter = "nu", brms_class = "nu",
    group = NA_character_, distribution = "gamma", arg1 = 2,
    arg2 = 0.1, arg3 = NA_real_, lower_bound = 0, directional = FALSE
  )
  do.call(rbind, rows)
}

#' Return the long fixed-effect contrast registry for all 12 estimands.
#'
#' @return Data frame with one row per nonzero fixed-effect weight.
estimand_contrast_registry <- function() {
  make_rating <- function(model_id, prefix) {
    ids <- paste0(prefix, 1:3)
    terms <- c(
      "accuracy_rating_within", "condition_c:accuracy_rating_within",
      "accuracy_rating_within", "condition_c:accuracy_rating_within",
      "condition_c:accuracy_rating_within"
    )
    data.frame(
      model_id = model_id,
      estimand_id = c(ids[1], ids[1], ids[2], ids[2], ids[3]),
      classification = "primary",
      fixed_term = terms,
      weight = c(1, -0.5, 1, 0.5, 1)
    )
  }
  beta <- data.frame(
    model_id = "H3",
    estimand_id = c("B1", "B1", "B2", "B2", "B3", "B4", "B4", "B5", "B5", "B6"),
    classification = c(rep("primary", 5), rep("supplementary", 5)),
    fixed_term = c(
      "phase_c", "condition_c:phase_c", "phase_c", "condition_c:phase_c",
      "condition_c:phase_c", "phase_c:familiarity_c",
      "condition_c:phase_c:familiarity_c", "phase_c:familiarity_c",
      "condition_c:phase_c:familiarity_c", "condition_c:phase_c:familiarity_c"
    ),
    weight = c(1, -0.5, 1, 0.5, 1, 1, -0.5, 1, 0.5, 1)
  )
  rbind(make_rating("H1", "T"), make_rating("H2", "A"), beta)
}

#' Calculate named estimands from a named fixed-effect vector.
#'
#' @param model_id One of H1, H2, or H3.
#' @param fixed Named numeric vector using explicit model-matrix term names.
#' @return Named numeric estimand vector.
estimands_from_fixed <- function(model_id, fixed) {
  require_contract(!is.null(names(fixed)), "fixed_effect_vector_must_be_named")
  require_contract(
    identical(names(fixed), expected_fixed_terms(model_id)),
    paste0(model_id, "_fixed_effect_identity_mismatch")
  )
  registry <- estimand_contrast_registry()
  registry <- registry[registry$model_id == model_id, , drop = FALSE]
  ids <- unique(registry$estimand_id)
  stats::setNames(vapply(ids, function(id) {
    rows <- registry[registry$estimand_id == id, , drop = FALSE]
    sum(fixed[rows$fixed_term] * rows$weight)
  }, numeric(1)), ids)
}

#' Compute one fixed prediction for a numeric coding combination.
#'
#' @param model_id One of H1, H2, or H3.
#' @param fixed Named fixed-effect vector.
#' @param newdata Single-row data frame.
#' @return Numeric fixed-effect prediction.
fixed_prediction <- function(model_id, fixed, newdata) {
  matrix <- stats::model.matrix(model_formulas(model_id)$fixed, newdata)
  matrix <- matrix[, names(fixed), drop = FALSE]
  as.numeric(matrix %*% fixed)
}

#' Calculate estimands by brute-force fixed predictions.
#'
#' @param model_id One of H1, H2, or H3.
#' @param fixed Named fixed-effect vector.
#' @return Named numeric estimand vector.
estimands_from_predictions <- function(model_id, fixed) {
  if (model_id %in% c("H1", "H2")) {
    slope <- function(condition) {
      high <- data.frame(condition_c = condition, accuracy_rating_within = 0.5,
                         accuracy_rating_between = 0, familiarity_c = 0)
      low <- high
      low$accuracy_rating_within <- -0.5
      fixed_prediction(model_id, fixed, high) - fixed_prediction(model_id, fixed, low)
    }
    overt <- slope(-0.5)
    imagery <- slope(0.5)
    ids <- if (model_id == "H1") c("T1", "T2", "T3") else c("A1", "A2", "A3")
    return(stats::setNames(c(overt, imagery, imagery - overt), ids))
  }
  phase_contrast <- function(condition, familiarity) {
    high <- data.frame(condition_c = condition, phase_c = 0.5, familiarity_c = familiarity)
    low <- high
    low$phase_c <- -0.5
    fixed_prediction("H3", fixed, high) - fixed_prediction("H3", fixed, low)
  }
  overt_random <- phase_contrast(-0.5, -0.5)
  overt_repeated <- phase_contrast(-0.5, 0.5)
  imagery_random <- phase_contrast(0.5, -0.5)
  imagery_repeated <- phase_contrast(0.5, 0.5)
  overt <- mean(c(overt_random, overt_repeated))
  imagery <- mean(c(imagery_random, imagery_repeated))
  stats::setNames(c(
    overt,
    imagery,
    imagery - overt,
    overt_repeated - overt_random,
    imagery_repeated - imagery_random,
    (imagery_repeated - imagery_random) - (overt_repeated - overt_random)
  ), paste0("B", 1:6))
}

#' Require direct contrast algebra and brute predictions to agree.
#'
#' @param model_id One of H1, H2, or H3.
#' @param fixed Named fixed-effect vector.
#' @param tolerance Absolute numerical tolerance.
#' @return Data frame showing both routes and their difference.
validate_estimand_algebra <- function(model_id, fixed, tolerance = 1e-12) {
  direct <- estimands_from_fixed(model_id, fixed)
  brute <- estimands_from_predictions(model_id, fixed)
  require_contract(identical(names(direct), names(brute)), "estimand_route_names_differ")
  difference <- direct - brute
  require_contract(max(abs(difference)) <= tolerance, paste0(model_id, "_estimand_routes_disagree"))
  data.frame(
    model_id = model_id, estimand_id = names(direct),
    direct_algebra = as.numeric(direct), brute_force_prediction = as.numeric(brute),
    absolute_difference = abs(as.numeric(difference)), tolerance = tolerance
  )
}

#' Return accepted EEG outcome-like columns found in a frame.
#'
#' @param columns Character column names.
#' @return Character vector of prohibited columns.
observed_eeg_columns <- function(columns) {
  lower <- tolower(columns)
  columns[
    columns %in% c("value_db", "value_log_power") |
      grepl("value_db|value_log_power", lower) |
      grepl("^(observed_eeg|accepted_eeg)", lower)
  ]
}

#' Enforce the synthetic fitting boundary.
#'
#' @param data Proposed fit frame.
#' @param require_response Whether `value_synthetic` must exist and be finite.
#' @return Invisibly returns TRUE; otherwise stops before fitting.
assert_synthetic_fit_input_safe <- function(data, require_response = TRUE) {
  forbidden <- observed_eeg_columns(names(data))
  require_contract(length(forbidden) == 0L, paste0("synthetic_fit_input_contains_observed_eeg:", paste(forbidden, collapse = ",")))
  if (require_response) {
    require_contract("value_synthetic" %in% names(data), "synthetic_response_missing")
    require_contract(all(is.finite(data$value_synthetic)), "synthetic_response_nonfinite")
  }
  invisible(TRUE)
}

#' Generate the deterministic balanced rating-model fixture.
#'
#' @param seed Deterministic synthetic-generation seed.
#' @return List containing data, fixed truth, estimands, and generation metadata.
generate_rating_synthetic <- function(seed) {
  set.seed(seed)
  participants <- sprintf("P%02d", seq_len(24))
  data <- expand.grid(trial_index = seq_len(40), participant_id = participants,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  participant_number <- match(data$participant_id, participants)
  data$condition_c <- ifelse(participant_number <= 12, -0.5, 0.5)
  data$familiarity_c <- ifelse(data$trial_index %% 2L == 0L, 0.5, -0.5)
  rating_pattern <- rep(c(-2, -1, 0, 1, 2), each = 8)
  data$accuracy_rating_within <- rating_pattern[data$trial_index]
  between_pattern <- rep(c(-1, -0.5, 0, 0.5, 1, 0), length.out = 24)
  data$accuracy_rating_between <- between_pattern[participant_number]
  data$canonical_event_key <- sprintf("rating:%s:%02d", data$participant_id, data$trial_index)
  fixed <- c(
    "(Intercept)" = -1.0,
    "condition_c" = 0.5,
    "accuracy_rating_within" = 1.2,
    "accuracy_rating_between" = 0.3,
    "familiarity_c" = -0.2,
    "condition_c:accuracy_rating_within" = 0.8
  )
  random_sd <- c("(Intercept)" = 0.4, "accuracy_rating_within" = 0.25, "familiarity_c" = 0.2)
  random <- sapply(random_sd, function(sd) stats::rnorm(length(participants), 0, sd))
  rownames(random) <- participants
  fixed_matrix <- stats::model.matrix(model_formulas("H1")$fixed, data)
  random_matrix <- cbind(1, data$accuracy_rating_within, data$familiarity_c)
  colnames(random_matrix) <- names(random_sd)
  participant_rows <- match(data$participant_id, participants)
  random_value <- rowSums(random_matrix * random[participant_rows, , drop = FALSE])
  residual_sigma <- 0.5
  residual_nu <- 7
  data$value_synthetic <- as.numeric(fixed_matrix %*% fixed) + random_value +
    residual_sigma * stats::rt(nrow(data), df = residual_nu)
  assert_synthetic_fit_input_safe(data)
  list(
    data = data,
    fixed = fixed,
    estimands = estimands_from_fixed("H1", fixed),
    metadata = list(seed = seed, participants = 24L, rows = nrow(data),
                    residual_sigma = residual_sigma, residual_nu = residual_nu,
                    participant_random_sd = random_sd)
  )
}

#' Generate the deterministic balanced beta phase-model fixture.
#'
#' @param seed Deterministic synthetic-generation seed.
#' @return List containing data, fixed truth, estimands, and generation metadata.
generate_beta_synthetic <- function(seed) {
  set.seed(seed)
  participants <- sprintf("P%02d", seq_len(24))
  trials <- expand.grid(trial_index = seq_len(20), participant_id = participants,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  participant_number <- match(trials$participant_id, participants)
  trials$condition_c <- ifelse(participant_number <= 12, -0.5, 0.5)
  trials$familiarity_c <- ifelse(trials$trial_index %% 2L == 0L, 0.5, -0.5)
  trials$canonical_event_key <- sprintf("beta:%s:%02d", trials$participant_id, trials$trial_index)
  data <- trials[rep(seq_len(nrow(trials)), each = 2L), ]
  data$phase_c <- rep(c(-0.5, 0.5), times = nrow(trials))
  fixed <- c(
    "(Intercept)" = -2.0,
    "condition_c" = 0.4,
    "phase_c" = 1.4,
    "familiarity_c" = 0.2,
    "condition_c:phase_c" = 0.8,
    "condition_c:familiarity_c" = -0.1,
    "phase_c:familiarity_c" = 0.6,
    "condition_c:phase_c:familiarity_c" = 0.4
  )
  participant_sd <- c("(Intercept)" = 0.35, "phase_c" = 0.25,
                      "familiarity_c" = 0.2, "phase_c:familiarity_c" = 0.15)
  random <- sapply(participant_sd, function(sd) stats::rnorm(length(participants), 0, sd))
  rownames(random) <- participants
  z <- cbind(1, data$phase_c, data$familiarity_c, data$phase_c * data$familiarity_c)
  colnames(z) <- names(participant_sd)
  participant_rows <- match(data$participant_id, participants)
  participant_value <- rowSums(z * random[participant_rows, , drop = FALSE])
  trial_sd <- 0.3
  trial_random <- stats::setNames(stats::rnorm(nrow(trials), 0, trial_sd), trials$canonical_event_key)
  fixed_matrix <- stats::model.matrix(model_formulas("H3")$fixed, data)
  residual_sigma <- 0.45
  residual_nu <- 7
  data$value_synthetic <- as.numeric(fixed_matrix %*% fixed) + participant_value +
    trial_random[data$canonical_event_key] + residual_sigma * stats::rt(nrow(data), df = residual_nu)
  assert_synthetic_fit_input_safe(data)
  list(
    data = data,
    fixed = fixed,
    estimands = estimands_from_fixed("H3", fixed),
    metadata = list(seed = seed, participants = 24L, trials = nrow(trials), rows = nrow(data),
                    residual_sigma = residual_sigma, residual_nu = residual_nu,
                    participant_random_sd = participant_sd, trial_random_sd = trial_sd)
  )
}

#' Convert a brms fixed term to its posterior variable name.
#'
#' @param term Explicit model-matrix term.
#' @return brms posterior variable name.
brms_fixed_variable <- function(term) {
  if (term == "(Intercept)") "b_Intercept" else paste0("b_", term)
}

#' Derive estimand draws through the production contrast registry.
#'
#' @param fit A fitted brms object using only synthetic outcomes.
#' @param model_id H1 or H3 for the bounded recovery fits.
#' @return posterior draws_df containing only named estimands.
posterior_estimand_draws <- function(fit, model_id) {
  draws <- posterior::as_draws_df(fit)
  expected <- expected_fixed_terms(model_id)
  variables <- vapply(expected, brms_fixed_variable, character(1))
  require_contract(all(variables %in% names(draws)), paste0(model_id, "_posterior_fixed_terms_missing"))
  fixed_draws <- stats::setNames(lapply(variables, function(variable) draws[[variable]]), expected)
  registry <- estimand_contrast_registry()
  registry <- registry[registry$model_id == model_id, , drop = FALSE]
  result <- data.frame(.chain = draws$.chain, .iteration = draws$.iteration, .draw = draws$.draw)
  for (id in unique(registry$estimand_id)) {
    rows <- registry[registry$estimand_id == id, , drop = FALSE]
    value <- rep(0, nrow(result))
    for (index in seq_len(nrow(rows))) {
      value <- value + fixed_draws[[rows$fixed_term[index]]] * rows$weight[index]
    }
    result[[id]] <- value
  }
  posterior::as_draws_df(result)
}

#' Summarize sampler diagnostics and synthetic estimand recovery.
#'
#' @param fit Fitted synthetic brms object.
#' @param model_id H1 or H3.
#' @param truth Named injected estimand values.
#' @param criteria Parsed engineering criteria.
#' @param max_treedepth Configured sampler maximum.
#' @return List with recovery table, diagnostic counts, and pass flag.
summarize_synthetic_fit <- function(fit, model_id, truth, criteria, max_treedepth) {
  estimand_draws <- posterior_estimand_draws(fit, model_id)
  ids <- names(truth)
  estimand_draws <- posterior::subset_draws(estimand_draws, variable = ids)
  summary <- posterior::summarise_draws(
    estimand_draws,
    mean = mean,
    sd = stats::sd,
    q2.5 = ~stats::quantile(.x, 0.025),
    q97.5 = ~stats::quantile(.x, 0.975),
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  )
  summary <- as.data.frame(summary)
  summary_names <- names(summary)
  summary_names[summary_names == "variable"] <- "estimand_id"
  summary_names[summary_names == "2.5%"] <- "q2.5"
  summary_names[summary_names == "97.5%"] <- "q97.5"
  summary_names[summary_names == "posterior::rhat"] <- "rhat"
  summary_names[summary_names == "posterior::ess_bulk"] <- "ess_bulk"
  summary_names[summary_names == "posterior::ess_tail"] <- "ess_tail"
  names(summary) <- summary_names
  require_contract(
    all(c("estimand_id", "q2.5", "q97.5", "rhat", "ess_bulk", "ess_tail") %in% names(summary)),
    "posterior_summary_columns_changed"
  )
  summary$truth <- as.numeric(truth[summary$estimand_id])
  summary$truth_in_interval <- summary$truth >= summary$q2.5 & summary$truth <= summary$q97.5
  summary$absolute_mean_error <- abs(summary$mean - summary$truth)
  summary$mean_close_relative_to_uncertainty <- summary$absolute_mean_error <=
    as.numeric(criteria$posterior_mean_maximum_standard_errors_from_truth) * summary$sd
  require_contract(identical(fit$backend, "cmdstanr"), "synthetic_fit_backend_not_cmdstanr")
  sampler <- brms::nuts_params(fit)
  divergences <- sum(sampler$Value[sampler$Parameter == "divergent__"])
  treedepth_saturations <- sum(
    sampler$Value[sampler$Parameter == "treedepth__"] >= max_treedepth
  )
  pass <- divergences <= as.integer(criteria$maximum_divergences) &&
    treedepth_saturations <= as.integer(criteria$maximum_treedepth_saturations) &&
    all(summary$rhat <= as.numeric(criteria$maximum_rhat)) &&
    all(summary$ess_bulk >= as.numeric(criteria$minimum_bulk_ess)) &&
    all(summary$ess_tail >= as.numeric(criteria$minimum_tail_ess)) &&
    all(summary$truth_in_interval) && all(summary$mean_close_relative_to_uncertainty)
  list(
    recovery = summary,
    diagnostics = list(
      model_id = model_id, backend = fit$backend, divergences = divergences,
      max_treedepth_saturations = treedepth_saturations,
      maximum_rhat = max(summary$rhat), minimum_bulk_ess = min(summary$ess_bulk),
      minimum_tail_ess = min(summary$ess_tail), status = if (pass) "pass" else "fail"
    ),
    pass = pass
  )
}

#' Run direct prior-only predictive simulation over an accepted safe design.
#'
#' @param model_id One of H1, H2, or H3.
#' @param design Safe accepted predictor/group/key design with no outcomes.
#' @param prior_config Parsed tracked prior constants.
#' @param criteria Parsed engineering criteria.
#' @param draws Number of prior draws.
#' @param seed Deterministic prior-predictive seed.
#' @param figure_sample Maximum pooled predictions retained for plotting.
#' @return List with predictor-collapsed summaries, parameter draws, and sample.
simulate_prior_predictive <- function(model_id, design, prior_config, criteria,
                                      draws, seed, figure_sample) {
  assert_synthetic_fit_input_safe(design, require_response = FALSE)
  set.seed(seed)
  constants <- prior_config$models[[model_id]]
  m_y <- as.numeric(constants$m_y)
  s_y <- as.numeric(constants$s_y)
  x <- stats::model.matrix(model_formulas(model_id)$fixed, design)
  require_contract(identical(colnames(x), expected_fixed_terms(model_id)), "prior_fixed_matrix_terms_changed")
  participant <- as.integer(factor(design$participant_id, levels = unique(design$participant_id)))
  if (model_id %in% c("H1", "H2")) {
    z <- cbind(1, design$accuracy_rating_within, design$familiarity_c)
    random_names <- c("Intercept", "accuracy_rating_within", "familiarity_c")
  } else {
    z <- cbind(1, design$phase_c, design$familiarity_c, design$phase_c * design$familiarity_c)
    random_names <- c("Intercept", "phase_c", "familiarity_c", "phase_c:familiarity_c")
    trial <- as.integer(factor(design$canonical_event_key, levels = unique(design$canonical_event_key)))
  }
  total <- nrow(design) * draws
  pooled <- numeric(total)
  parameters <- vector("list", draws)
  cursor <- 1L
  for (draw in seq_len(draws)) {
    fixed <- c(m_y + 2.5 * s_y * stats::rt(1, 3), stats::rnorm(ncol(x) - 1L, 0, s_y))
    names(fixed) <- colnames(x)
    participant_sd <- abs(s_y * stats::rt(ncol(z), 3))
    participant_effects <- matrix(
      stats::rnorm(length(unique(participant)) * ncol(z)),
      nrow = length(unique(participant)), ncol = ncol(z)
    )
    participant_effects <- sweep(participant_effects, 2, participant_sd, `*`)
    mu <- as.numeric(x %*% fixed) + rowSums(z * participant_effects[participant, , drop = FALSE])
    trial_sd <- NA_real_
    if (model_id == "H3") {
      trial_sd <- abs(s_y * stats::rt(1, 3))
      trial_effects <- stats::rnorm(length(unique(trial)), 0, trial_sd)
      mu <- mu + trial_effects[trial]
    }
    sigma <- abs(s_y * stats::rt(1, 3))
    nu <- stats::rgamma(1, shape = 2, rate = 0.1)
    predicted <- mu + sigma * stats::rt(nrow(design), df = nu)
    destination <- cursor:(cursor + nrow(design) - 1L)
    pooled[destination] <- predicted
    cursor <- cursor + nrow(design)
    parameter <- c(fixed, stats::setNames(participant_sd, paste0("sd_participant_", random_names)),
                   sigma = sigma, nu = nu)
    if (model_id == "H3") parameter <- c(parameter, sd_canonical_event_key = trial_sd)
    parameters[[draw]] <- data.frame(draw = draw, parameter = names(parameter), value = as.numeric(parameter))
  }
  parameter_draws <- do.call(rbind, parameters)
  finite <- is.finite(pooled)
  standardized <- abs((pooled[finite] - m_y) / s_y)
  quantiles <- stats::quantile(pooled[finite], c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999), names = FALSE)
  central90_width_sy <- (quantiles[5] - quantiles[3]) / s_y
  proportion_20 <- mean(standardized > 20)
  proportion_50 <- mean(standardized > 50)
  sd_values <- parameter_draws$value[grepl("^sd_", parameter_draws$parameter)] / s_y
  nu_values <- parameter_draws$value[parameter_draws$parameter == "nu"]
  pass <- mean(finite) >= as.numeric(criteria$prior_predictive_minimum_finite_proportion) &&
    central90_width_sy >= as.numeric(criteria$prior_predictive_minimum_central90_width_sy) &&
    proportion_20 <= as.numeric(criteria$prior_predictive_maximum_proportion_beyond_20_sy) &&
    proportion_50 <= as.numeric(criteria$prior_predictive_maximum_proportion_beyond_50_sy) &&
    stats::median(sd_values) >= as.numeric(criteria$prior_sd_median_minimum_sy) &&
    stats::median(sd_values) <= as.numeric(criteria$prior_sd_median_maximum_sy) &&
    mean(nu_values < 10) >= as.numeric(criteria$nu_minimum_probability_below_10) &&
    mean(nu_values > 30) >= as.numeric(criteria$nu_minimum_probability_above_30)
  sample_size <- min(length(pooled), as.integer(figure_sample))
  sample_index <- if (sample_size < length(pooled)) sample.int(length(pooled), sample_size) else seq_along(pooled)
  list(
    summary = data.frame(
      model_id = model_id, draws = draws, predictions = length(pooled),
      finite_proportion = mean(finite), q0.001 = quantiles[1], q0.01 = quantiles[2],
      q0.05 = quantiles[3], median = quantiles[4], q0.95 = quantiles[5],
      q0.99 = quantiles[6], q0.999 = quantiles[7],
      central90_width_sy = central90_width_sy,
      proportion_beyond_20_sy = proportion_20,
      proportion_beyond_50_sy = proportion_50,
      median_group_sd_sy = stats::median(sd_values),
      probability_nu_below_10 = mean(nu_values < 10),
      probability_nu_above_30 = mean(nu_values > 30),
      status = if (pass) "pass" else "owner_review_needed"
    ),
    parameter_draws = parameter_draws,
    prediction_sample = pooled[sample_index],
    pass = pass
  )
}

#' Write a predictor-collapsed prior-predictive density figure.
#'
#' @param path PNG output path.
#' @param values Pooled sampled prior predictions.
#' @param model_id Model label.
#' @param m_y Frozen pooled outcome median.
#' @param s_y Frozen pooled robust scale.
#' @return Invisibly returns the path after closing the device.
write_prior_predictive_figure <- function(path, values, model_id, m_y, s_y) {
  grDevices::png(path, width = 1200, height = 800, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  standardized <- (values - m_y) / s_y
  clipped <- standardized[is.finite(standardized) & abs(standardized) <= 50]
  graphics::hist(clipped, breaks = 160, probability = TRUE, col = "#5B8FF9",
                 border = NA, main = paste(model_id, "pooled prior-predictive scale"),
                 xlab = "Prior prediction relative to pooled median (s_y units)")
  graphics::abline(v = 0, lwd = 2, col = "#222222")
  graphics::abline(v = c(-20, 20), lty = 2, col = "#D62728")
  invisible(path)
}

#' Return a compact structural summary for generated brms standata.
#'
#' @param standata Named list returned by `brms::make_standata`.
#' @return Data frame of element classes, dimensions, and lengths.
standata_summary <- function(standata) {
  data.frame(
    element = names(standata),
    class = vapply(standata, function(x) paste(class(x), collapse = "/"), character(1)),
    dimensions = vapply(standata, function(x) paste(dim(x), collapse = "x"), character(1)),
    length = vapply(standata, length, integer(1)),
    stringsAsFactors = FALSE
  )
}

#' Create an unused sibling staging directory for atomic publication.
#'
#' @param target Final validation directory.
#' @return Newly created staging path.
make_staging_directory <- function(target) {
  parent <- dirname(target)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- file.path(parent, paste0(".", basename(target), ".staging-", Sys.getpid()))
  require_contract(!file.exists(staging), "model_validation_staging_already_exists")
  dir.create(staging, recursive = FALSE)
  staging
}

#' Atomically publish a complete staging directory.
#'
#' @param staging Complete validated staging path.
#' @param target Final production path, which must not exist.
#' @return Invisibly returns the final path.
publish_staging_directory <- function(staging, target) {
  require_contract(dir.exists(staging), "model_validation_staging_missing")
  require_contract(!file.exists(target), "model_validation_target_exists")
  require_contract(file.rename(staging, target), "model_validation_atomic_rename_failed")
  invisible(target)
}
