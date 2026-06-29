# DEMI EEG Reanalysis Repository Audit

Date: 2026-06-29

## Purpose

This audit maps the current DEMI repository before implementing the planned Python/MNE-first EEG reanalysis. It is grounded in the existing repository conventions and the public EEG reanalysis plan in `docs/eeg_reanalysis_plan.md`.

The main conclusion is that the repository already has a coherent historical R-first pipeline under `_Scripts/`, but the active EEG path is visually mixed with behaviour scripts and historical EEG/GAM work. The behavioural analysis should remain frozen. The prior EEG/GAM analysis should be preserved as reference, but not left as the obvious active path once the MNE reanalysis begins.

## Concise Repository Map

```text
.
|-- _Scripts/
|   |-- 00_import.R
|   |-- 01_preprocessing.R
|   |-- 02_imagery_clean.R
|   |-- 03_behav_analysis.R
|   |-- 04_import_eeg.R
|   |-- 05_eeg_prep.R
|   |-- 06_eeg_plots.R
|   |-- 07_gam_prep.R
|   |-- 08_run_gams.R
|   |-- 09_maxT_inference.R
|   |-- 10_posthoc_group_maps.R
|   |-- 11_results_tables.R
|   |-- 12_main_figure_maxT.R
|   |-- 13a_posthoc_raw_difference_maps.R
|   |-- 13b_posthoc_raw_condition_maps.R
|   |-- _functions/
|   |-- _misc/
|   `-- _stan/
|-- _Data/
|   `-- eeg/BESA-81.csv
|-- docs/
|   `-- eeg_reanalysis_plan.md
|-- external/
|   `-- DEMI_EEG_Pipeline
|-- legacy/
|   `-- dissertation/
|-- media/
|   |-- behaviour/
|   `-- eegplots/
|-- renv/
|-- .gitignore
|-- README.md
|-- renv.lock
`-- requirements-eeg.txt
```

## File and Folder Classification

### Behaviour Analysis

- `_Scripts/00_import.R`: imports task text files and figure/tracing ZIP files from `_Data/task/` and `_Data/figure/`.
- `_Scripts/01_preprocessing.R`: computes figure complexity, cleans tracing responses, computes behavioural error metrics, and writes behavioural/intermediate RDS files.
- `_Scripts/02_imagery_clean.R`: flags imagery trials with implausible durations using a Stan mixture model.
- `_Scripts/03_behav_analysis.R`: published behavioural analyses and figures. This should remain frozen unless a technical issue blocks EEG trial metadata export.
- `_Scripts/_functions/complexity.R`, `_Scripts/_functions/filters.R`, `_Scripts/_functions/physical.R`, `_Scripts/_functions/visualization.R`: helper functions for behavioural import, tracing cleanup, accuracy/error metrics, and filter plots.
- `_Scripts/_stan/outlier_mixture.stan`: model used by the imagery cleanup script.

### Current / Old EEG Analysis

- `_Scripts/04_import_eeg.R`: imports externally preprocessed EDF files, joins behavioural timing data, updates events, epochs EEG/EMG, and caches `eeg_imported.Rds`.
- `_Scripts/05_eeg_prep.R`: EMG checks, wavelet time-frequency decomposition, dB normalization, participant shards, and band/epoch/sensor aggregation.
- `_Scripts/06_eeg_plots.R`: participant and group time-frequency QC plots.
- `_Scripts/07_gam_prep.R`: prepares band/epoch/sensor data for GAMs.
- `_Scripts/08_run_gams.R`: fits the hierarchical GAM-style scalp model.
- `_Scripts/09_maxT_inference.R`: global max-|T| inference for planned hypotheses.
- `_Scripts/10_posthoc_group_maps.R`: descriptive per-group post-hoc maps.
- `_Scripts/11_results_tables.R`: summary tables from maxT/post-hoc outputs.
- `_Scripts/12_main_figure_maxT.R`: main maxT figure.
- `_Scripts/13a_posthoc_raw_difference_maps.R`: raw subject-level difference maps.
- `_Scripts/13b_posthoc_raw_condition_maps.R`: raw subject-level condition maps.
- `_Scripts/_functions/eeg.R`: wavelet helpers and the event-update logic that must be preserved or carefully ported.
- `_Scripts/_functions/emg.R`: EMG filtering/rectification helpers.
- `_Scripts/_misc/sensor_latlong_chan_map.csv`: sensor coordinates for scalp plotting/GAMs.
- `_Data/eeg/BESA-81.csv`: tracked EEG coordinate/channel map.
- `external/DEMI_EEG_Pipeline`: git submodule pointer to the older external EEG preprocessing pipeline.

### Legacy Dissertation Material

- `legacy/dissertation/`: archived dissertation-era code and preprocessing. Its README marks it as archived/reference-only.
- `legacy/dissertation/_Scripts/`: older R analysis scripts with similar numbering.
- `legacy/dissertation/_Preprocessing/`: older Python preprocessing pipeline, including EDF-to-BIDS conversion and preprocessing scripts.
- `legacy/dissertation/_Data/`: placeholder/local-data conventions plus tracked `BESA-81.csv`.

### Generated Artifacts

These are mostly ignored by `.gitignore` and should not be treated as source:

- `_Scripts/_rds/`: RDS caches, model objects, behavioural intermediates, EEG shards, GAM-ready data.
- `_Plots/`: behavioural, EEG QC, post-hoc, and figure outputs.
- `_Tables/`: generated CSV result tables.
- `filters/`: behavioural tracing filter diagnostic plots.
- `_Scripts/_plots/`: ignored preprocessing/QC artifacts.
- `_Scripts/_stan/*` build artifacts except tracked `.stan` sources.
- Raw/preprocessed EEG formats such as `.edf`, `.fif`, `.set`, `.bdf`, `.vhdr`, `.vmrk`, `.eeg`.

### Docs and Media

- `README.md`: public repository overview. It currently describes `_Scripts/04+` as the EEG paper path and frames the EEG analysis as a GAM-based WIP.
- `docs/eeg_reanalysis_plan.md`: public plan for the new raw-to-MNE reanalysis.
- `media/behaviour/`: curated behavioural figures used by the README.
- `media/eegplots/`: curated EEG result summary image used by the README. This should eventually be checked for consistency with the new active analysis status.

### Local Data Expectations

- `_Data/` is local-only by default.
- `_Data/eeg/BESA-81.csv` is the only tracked file currently under `_Data/`.
- Current task data are expected locally under `_Data/task/`.
- Current figure/tracing data are expected locally under `_Data/figure/`.
- Old EEG import expects externally preprocessed EEG files at:
  - `_Data/eeg/prep_info.csv`
  - `_Data/eeg/edfs/sub-XXX_eeg_prepped.edf`

## Current Behaviour Pipeline and Trial Metadata Sources

The behavioural pipeline order is:

1. `_Scripts/00_import.R`
   - Reads `_Data/task/*.txt`.
   - Reads `_Data/figure/**/*.zip`.
   - Parses trial-level task data, figure vertices, figure segments, animation frames, and physical tracing samples.

2. `_Scripts/01_preprocessing.R`
   - Sources `00_import.R`.
   - Computes stimulus path length, sinuosity, curvature, entropy-style shape measures.
   - Cleans physical tracing samples with multiple filters: no-shape, failed-end trimming, touchscreen glitches, false starts, hand noise, incomplete traces, excessive gaps, and edge-hit flags.
   - Computes tracing duration/path metrics and multiple error metrics, including DTW-based error.
   - Joins task, stimulus, and tracing summaries.
   - Writes:
     - `_Scripts/_rds/bdat.rds`
     - `_Scripts/_rds/tracings.rds`
     - `_Scripts/_rds/figtrace.rds`
     - `_Scripts/_rds/trace_filter_info.rds`

3. `_Scripts/02_imagery_clean.R`
   - Reads `_Scripts/_rds/bdat.rds`.
   - Flags implausibly fast/unreasonable imagery trials.
   - Writes `_Scripts/_rds/bad_imagery_trials.rds`.

4. `_Scripts/03_behav_analysis.R`
   - Reads `bdat.rds` and `bad_imagery_trials.rds`.
   - Excludes pilot and bad trials.
   - Creates analysis variables used downstream: `group`, `condition`, `rep`, `complexity`, `error`, corrected `mt`, `skill`, z-scored predictors, and `trial`.
   - Writes `_Scripts/_rds/bdat2.rds`.
   - Fits/saves published behavioural models and behavioural figures.

For the new EEG pipeline, `bdat2.rds` is the closest existing source for stable trial metadata. A dedicated export should be added later, likely from a small wrapper that reads `bdat2.rds` and writes `derivatives/behaviour/trial_metadata_for_eeg.csv` without changing the behavioural analysis scripts. Candidate fields already available or derivable include participant, group, condition, block/trial numbers, trial index, repetition/random status, stimulus/movement duration, complexity, accuracy rating, physical error, expected/derived performance fields, handedness, and inclusion flags.

## Current EEG/GAM Pipeline Order, Inputs, and Outputs

The current EEG/GAM pipeline is R-first and starts after a separate preprocessing pipeline has already produced EDF+ files.

1. `_Scripts/04_import_eeg.R`
   - Inputs:
     - `_Scripts/_rds/bdat2.rds`
     - `_Scripts/_rds/tracings.rds`
     - `_Scripts/_rds/trace_filter_info.rds`
     - `_Data/eeg/prep_info.csv`
     - `_Data/eeg/edfs/sub-XXX_eeg_prepped.edf`
   - Uses `eeguana::read_edf()`.
   - Drops practice/unusable events, adds `real_trace_start` and `real_trace_end`, separates EEG from EOG/EMG, drops bad channels, epochs baseline/tracing/post-tracing EEG, and epochs EMG for imagery participants.
   - Output:
     - `_Scripts/_rds/eeg_imported.Rds`

2. `_Scripts/05_eeg_prep.R`
   - Inputs:
     - `_Scripts/_rds/eeg_imported.Rds`
     - `_Scripts/_rds/bdat2.rds`
     - `_Data/eeg/BESA-81.csv`
   - Performs EMG activity summaries and optional imagery trial exclusion.
   - Runs wavelet transforms for baseline, tracing, and post-tracing epochs over `wt_frequencies <- 1:40`.
   - Uses either pre-trial baseline dB normalization or an epoch-local baseline depending on `_Scripts/_settings.R`.
   - Renames tracing/post-trace epochs to `during`/`after`.
   - Outputs:
     - `_Scripts/_rds/participants/*_eeg_processed.rds`
     - `_Scripts/_rds/participants_agg/*_agg.rds`
     - `_Scripts/_rds/dat_bands_notime.rds`
     - optionally `_Scripts/_rds/all_dat.rds`

3. `_Scripts/06_eeg_plots.R`
   - Inputs:
     - `_Scripts/_rds/participants/*_eeg_processed.rds`
     - `_Scripts/_rds/bdat2.rds`
     - `_Data/eeg/BESA-81.csv`
   - Outputs:
     - `_Scripts/_rds/participants_mean/*_mean.rds`
     - `_Plots/eeg_timefreq/...`

4. `_Scripts/07_gam_prep.R`
   - Input:
     - `_Scripts/_rds/dat_bands_notime.rds`
     - `_Scripts/_rds/bdat2.rds` for handedness lookup
   - Drops reference channels, scales `accuracy_rating`, mirrors left-handed coordinates, normalizes factors, derives block if needed.
   - Outputs:
     - `_Scripts/_rds/dat_gam.rds`
     - `_Plots/pmbr_C3_beta.png`

5. `_Scripts/08_run_gams.R`
   - Input:
     - `_Scripts/_rds/dat_gam.rds`
   - Bins accuracy at the 20th/80th percentiles and fits the GAM model.
   - Outputs:
     - `_Scripts/_rds/accuracy_bin_cutpoints.rds`
     - `_Scripts/_rds/gam_full.rds`
     - `_Scripts/_rds/gam_full_summary.rds`
     - `_Scripts/_rds/gam_full_sp.rds`

6. `_Scripts/09_maxT_inference.R`
   - Inputs:
     - `_Scripts/_rds/gam_full.rds`
     - `_Scripts/_rds/dat_gam.rds`
   - Tests H1-H3 using a global max-|T| threshold across hypotheses and sensors.
   - Outputs:
     - `_Tables/maxT_results_h1.csv`
     - `_Tables/maxT_results_h2.csv`
     - `_Tables/maxT_results_h3.csv`
     - `_Plots/maxT_topo_h1.pdf`
     - `_Plots/maxT_topo_h2.pdf`
     - `_Plots/maxT_topo_h3.pdf`

7. `_Scripts/10_posthoc_group_maps.R`
   - Inputs:
     - `_Scripts/_rds/gam_full.rds`
     - `_Scripts/_rds/dat_gam.rds`
     - `_Scripts/_misc/sensor_latlong_chan_map.csv` if present
   - Outputs:
     - `_Tables/posthoc_*.csv`
     - `_Plots/posthoc_*.pdf`

8. `_Scripts/11_results_tables.R`
   - Inputs:
     - `_Tables/maxT_results_*.csv`
     - `_Tables/posthoc_*.csv`
   - Outputs:
     - `_Tables/main_hypotheses_summary.csv`
     - `_Tables/posthoc_maps_summary.csv`

9. `_Scripts/12_main_figure_maxT.R`
   - Inputs:
     - `_Tables/maxT_results_h1.csv`
     - `_Tables/maxT_results_h2.csv`
     - `_Tables/maxT_results_h3.csv`
     - `_Scripts/_misc/sensor_latlong_chan_map.csv`
   - Outputs:
     - `_Plots/figure_main_maxT.pdf`
     - `_Plots/figure_main_maxT.png`

10. `_Scripts/13a_posthoc_raw_difference_maps.R`
    - Inputs:
      - `_Scripts/_rds/dat_gam.rds`
      - `_Scripts/_rds/accuracy_bin_cutpoints.rds`
      - `_Scripts/_misc/sensor_latlong_chan_map.csv`
    - Outputs:
      - `_Tables/posthoc_raw_sensor_summaries.csv`
      - `_Plots/figure_posthoc_raw_maps.pdf`
      - `_Plots/figure_posthoc_raw_maps.png`

11. `_Scripts/13b_posthoc_raw_condition_maps.R`
    - Inputs:
      - `_Scripts/_rds/dat_gam.rds`
      - `_Scripts/_rds/accuracy_bin_cutpoints.rds`
      - `_Scripts/_misc/sensor_latlong_chan_map.csv`
    - Outputs:
      - `_Tables/posthoc_raw_condition_sensor_summaries.csv`
      - `_Plots/figure_posthoc_raw_conditions.pdf`
      - `_Plots/figure_posthoc_raw_conditions.png`

## Event-Alignment Logic to Preserve or Carefully Port

The critical event-alignment logic lives in `_Scripts/04_import_eeg.R` and `_Scripts/_functions/eeg.R`.

Important behaviours to preserve:

- Trial numbering is reconstructed from EEG triggers using `stim_on` and a red-on bad-start rule.
- Practice trials are identified by long figure durations and removed.
- Trials missing a valid `stim_on`/start pattern are removed.
- EEG trial numbers are reindexed after practice/bad starts are dropped.
- Behavioural trial count is `trial + (block - 1) * 20`.
- Physical trial timing is corrected using behavioural/tracing data:
  - `trace_onset` is the first raw tracing sample time.
  - `trace_end` is the max raw tracing sample time.
  - false-start `start_shift` from `trace_filter_info$false_start` is added when present.
  - physical `real_start` is based on task inter-trial timing, raw trace onset, and false-start shift.
  - physical `real_end` is `real_start + mt`.
  - imagery `real_start` is treated as 0 and imagery `real_end` as `mt`.
- New events are created:
  - `real_trace_start`
  - `real_trace_end`
- For physical trials, new event times are offsets from the original `trace_start` trigger.
- For imagery trials, the start remains the original start and the end uses the original end.
- The code validates alignment:
  - flags estimated trace-end mismatch over 150 ms;
  - warns when EEG figure duration differs from task `stimulus_mt` by more than 250 ms.
- Trials without `real_trace_end` are filtered from the final event table.
- Updated event times are converted from milliseconds to samples using the recording sampling rate.

For MNE, this should become an explicit, tested event-construction module. It should produce event-count and event-duration QC tables before any large preprocessing run.

## Ignored and Local Data Conventions

The root `.gitignore` is already conservative for local data:

- Ignores all of `_Data/`.
- Re-allows `_Data/eeg/`.
- Ignores everything inside `_Data/eeg/`.
- Re-allows only `_Data/eeg/BESA-81.csv`.
- Explicitly ignores:
  - `_Data/eeg/prep_info.csv`
  - `_Data/eeg/edfs/`
  - `_Data/eeg/edfs/**`
- Ignores common EEG formats anywhere:
  - `.edf`, `.EDF`, `.eeg`, `.EEG`, `.vhdr`, `.vmrk`, `.fif`, `.set`, `.bdf`, `.BDF`
- Ignores generated model/data artifacts:
  - `.rds`, `*_rds/`, `.RData`
- Ignores `_Plots/`, `_Tables/`, filter diagnostics, and most preprocessing/QC artifacts.

Recommendations for the MNE reanalysis:

- Keep raw EEG local-only under `_Data/eeg/`.
- Prefer a BIDS-like raw layout if the raw recordings can support it, for example `_Data/eeg/raw_bids/`.
- Keep source raw recordings and generated derivatives separate.
- Do not place new source data or derived `.fif`/EDF files under tracked paths.
- Use a public, tracked `derivatives/README.md` or `docs` note only if it describes expected local paths without committing data.
- Use ignored derivative paths for generated files, for example:
  - `derivatives/behaviour/trial_metadata_for_eeg.csv`
  - `derivatives/eeg/manifest/`
  - `derivatives/eeg/preprocessed/`
  - `derivatives/eeg/epochs/`
  - `derivatives/eeg/time_frequency/`
  - `derivatives/eeg/features/`
  - `derivatives/eeg/qc/`
  - `derivatives/eeg/stats/`
  - `derivatives/eeg/figures/`
- Before creating those paths, decide whether `derivatives/` should be fully ignored with tracked `.gitkeep`/README exceptions, or whether only heavyweight binary outputs should be ignored.
- Add explicit ignore rules for MNE reports and derivative tables once filenames are chosen.

## Archive / Retain / Active-Until-Replaced Assessment

### Likely Archive as Prior EEG/GAM Analysis

These should likely move together later under a clearly labelled archive path such as `archive/eeg_gam_analysis/`, with a small README explaining that they are superseded:

- `_Scripts/04_import_eeg.R`
- `_Scripts/05_eeg_prep.R`
- `_Scripts/06_eeg_plots.R`
- `_Scripts/07_gam_prep.R`
- `_Scripts/08_run_gams.R`
- `_Scripts/09_maxT_inference.R`
- `_Scripts/10_posthoc_group_maps.R`
- `_Scripts/11_results_tables.R`
- `_Scripts/12_main_figure_maxT.R`
- `_Scripts/13a_posthoc_raw_difference_maps.R`
- `_Scripts/13b_posthoc_raw_condition_maps.R`
- `_Scripts/_functions/eeg.R`
- `_Scripts/_functions/emg.R`
- `_Scripts/_misc/sensor_latlong_chan_map.csv`
- old generated `_Scripts/_rds/`, `_Plots/`, and `_Tables/` products if they are retained locally.

Archive timing should wait until the MNE pipeline has replacements for event construction, epoching, feature extraction, planned inference, and figures.

### Retain as Frozen Behaviour Foundation

These should remain in place unless owner explicitly decides to reorganize the published behavioural analysis:

- `_Scripts/00_import.R`
- `_Scripts/01_preprocessing.R`
- `_Scripts/02_imagery_clean.R`
- `_Scripts/03_behav_analysis.R`
- `_Scripts/_functions/complexity.R`
- `_Scripts/_functions/filters.R`
- `_Scripts/_functions/physical.R`
- `_Scripts/_functions/visualization.R`
- `_Scripts/_stan/outlier_mixture.stan`

If the MNE pipeline needs behavioural metadata, add a new export/wrapper rather than editing these scripts first.

### Treat as Active Only Until Replaced

- `_Scripts/04_import_eeg.R` and `_Scripts/_functions/eeg.R`: active only as the source of event-alignment logic until a tested MNE implementation exists.
- `_Scripts/05_eeg_prep.R`: active only as a reference for wavelet settings, baseline handling, bands, epochs, and variable joins.
- `_Scripts/07_gam_prep.R` through `_Scripts/13b_posthoc_raw_condition_maps.R`: active only as historical/reference analysis until the new primary statistical path is implemented.
- `external/DEMI_EEG_Pipeline`: retain as provenance for the previous preprocessed EDF path. It should not be the active preprocessing implementation for the new raw-MNE reanalysis unless explicitly revived.

## Proposed Next-Step Repo Structure

The current repo uses `_Scripts/` and `_Data/`, so the next structure should evolve from those conventions rather than replacing them immediately.

Recommended near-term structure:

```text
_Scripts/
  00_import.R
  01_preprocessing.R
  02_imagery_clean.R
  03_behav_analysis.R
  04_export_trial_metadata_for_eeg.R

analysis/
  eeg_mne/
    00_make_manifest.py
    01_import_raw.py
    02_build_events.py
    03_preprocess.py
    04_epoch.py
    05_qc_reports.py
    06_time_frequency.py
    07_extract_features.py
    08_primary_models.py
    09_cluster_support.py
    10_figures.py

docs/
  eeg_reanalysis_plan.md
  eeg_reanalysis_repo_audit.md
  eeg_methods_decisions.md

derivatives/
  behaviour/
    trial_metadata_for_eeg.csv
  eeg/
    manifest/
    events/
    preprocessed/
    epochs/
    time_frequency/
    features/
    qc/
    stats/
    figures/

archive/
  eeg_gam_analysis/
```

Notes:

- `_Scripts/04_export_trial_metadata_for_eeg.R` would be the bridge from the frozen R behavioural pipeline to the Python/MNE EEG pipeline.
- `analysis/eeg_mne/02_build_events.py` should be the port of the current R event-alignment logic.
- `derivatives/eeg/events/` should store auditable event tables before epoching.
- `archive/eeg_gam_analysis/` should be created only when archival movement is explicitly approved.
- The old EEG scripts should not be moved until the MNE path has enough functionality that future readers are not left without a runnable/reference EEG workflow.

## Open Questions and Blockers

1. What exact raw EEG file format and local path should the MNE pipeline consume?
2. Are the raw recordings already organized as BIDS, partially BIDS, or loose EDF files?
3. Should `_Data/eeg/raw_bids/` be the local raw path, or should another local convention be used?
4. Should `derivatives/` be ignored wholesale, or should lightweight CSV/QC summaries be tracked selectively?
5. Is CSD primary, sensitivity-only, or historical-only for the new analysis?
6. Which montage/channel coordinate source should be authoritative in MNE: `BESA-81.csv`, `_Scripts/_misc/sensor_latlong_chan_map.csv`, recording metadata, or a standard montage?
7. What are the final onset and offset event definitions for movement and imagery epochs?
8. Should false-start correction and physical-tracing onset logic be ported exactly before any methodological changes?
9. What baseline window should be primary for time-frequency features?
10. What ROI channel lists define frontal theta, posterior alpha, and sensorimotor beta?
11. Should accuracy be modeled continuously, binned, or both in the primary model?
12. Should the primary model be Bayesian hierarchical, frequentist mixed-effects, subject-level contrasts, or a layered approach?
13. Which parts of the prior GAM output should be retained as supplemental/reference figures or tables?
14. How should the README be updated later so the old GAM analysis no longer appears to be the active EEG pipeline?
15. Should `external/DEMI_EEG_Pipeline` remain as a submodule after the MNE pipeline is established, or be archived as historical provenance?
