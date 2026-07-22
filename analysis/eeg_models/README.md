# EEG model contract validation

This directory implements the predeclared H1--H3 hierarchical Student-t
formulas, fixed-effect estimand algebra, model-specific priors, synthetic
recovery checks, and predictor-collapsed prior-predictive validation. It does
not fit accepted EEG outcomes or calculate a scientific estimate.

The data boundary is deliberate. `prepare_model_validation_inputs.py` has one
predictor-blind route that projects only pooled `value_db` for prior scales and
a separate allowlisted route that projects design/group/key columns while
excluding observed EEG outcomes. Every R fit input must contain
`value_synthetic` and is rejected if an accepted outcome column is present.

Run the focused contracts with:

```sh
Rscript --no-init-file tests/test_eeg_model_contracts.R
PATH="$(pwd)/.venv/bin:$PATH" python3 -m pytest tests/test_eeg_model_input_safety.py
```

Run or reopen the ignored validation namespace with:

```sh
tools/run_eeg_model_validation_v1.sh
tools/run_eeg_model_validation_v1.sh --verify-current
```

The R lock records the primary `brms`/`cmdstanr`/`posterior`/`bayesplot`/`loo`
stack. It does not currently lock `arrow`, `testthat`, or fallback `lme4`;
therefore Parquet projection uses the repository's locked Python environment,
focused R contracts use base R assertions, and the missing lock entries are
reported as reproducibility dependencies rather than silently added here.

