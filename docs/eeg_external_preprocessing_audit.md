# External EEG Preprocessing Pipeline Audit

Date: 2026-06-29

## Purpose

This audit documents the prior external EEG preprocessing pipeline used before the current DEMI R EEG scripts. The external pipeline is included as the initialized git submodule `external/DEMI_EEG_Pipeline` at commit `b5ebb0230bc5cc204b821306b8ad0c6f8152f16a`.

The goal is not to revive that pipeline as the active analysis path. The DEMI EEG analysis is moving to a new Python/MNE-first reanalysis from raw EEG recordings. This document identifies what the old pipeline handled, what the downstream R scripts handled, and what should be preserved, reconsidered, or replaced.

## Repository Scope Reviewed

Primary files reviewed:

- `docs/eeg_reanalysis_plan.md`
- `docs/eeg_reanalysis_repo_audit.md`
- `external/DEMI_EEG_Pipeline/README.md`
- `external/DEMI_EEG_Pipeline/eeg_pipeline.py`
- `external/DEMI_EEG_Pipeline/edf2bids.py`
- `external/DEMI_EEG_Pipeline/Pipfile`
- `external/DEMI_EEG_Pipeline/utils/save_edf.py`
- `external/DEMI_EEG_Pipeline/utils/EDF.py`
- `_Scripts/04_import_eeg.R`
- `_Scripts/05_eeg_prep.R`
- `_Scripts/_functions/eeg.R`
- `_Scripts/_functions/emg.R`
- `_Scripts/_settings.R`

## What The External Pipeline Appears To Do

The old external pipeline has two layers: an optional EDF-to-BIDS conversion script and the main preprocessing script.

### Optional raw EDF to BIDS conversion

`external/DEMI_EEG_Pipeline/edf2bids.py` appears to:

1. Read loose EDF files from an `edfs/` directory inside the external pipeline environment.
2. Select files whose names match `demi_(\d+)[\s\.]`, padding subject IDs to three digits.
3. Ignore split files by design because of the filename pattern.
4. Read each EDF with MNE.
5. Assign channel types for non-EEG channels:
   - `HEO`, `VEO` as EOG;
   - `EMG-L`, `EMG-A` as EMG.
6. Round sampling-rate and filter metadata and set line frequency to 60 Hz.
7. Convert numeric trigger annotations or stim-channel events to named annotations:
   - `28` -> `stim_on`
   - `30` -> `red_on`
   - `44` -> `trace_start`
   - `46` -> `trace_end`
   - `60` -> `accuracy_submit`
   - `62` -> `vividness_submit`
   - `100` -> `file start`
8. Write a BIDS-like dataset under `bids/` using `mne_bids.write_raw_bids()`.
9. Update the EEG JSON sidecar with project-specific metadata, some of which is explicitly incomplete or marked for confirmation.

This conversion script does not perform signal cleaning. It organizes raw recordings and makes event/channel metadata easier for the preprocessing script to consume.

### Main preprocessing script

`external/DEMI_EEG_Pipeline/eeg_pipeline.py` appears to process each BIDS subject in this order:

1. Discover subjects under `bids/sub-*`.
2. Skip subjects already listed in `output/prep_info.csv`.
3. Read the raw BIDS EEG recording with `mne_bids.read_raw_bids()`.
4. Mark a recording as "complete" only if it has at least 600 annotations. Incomplete recordings are recorded in `prep_info.csv` but not preprocessed.
5. Apply the `standard_1005` montage after MNE EEGBCI-style channel-name standardization.
6. Load the data, and for recordings containing a `file start` annotation, crop at the first such onset to remove duplicated data. A code comment says this was needed for `sub-005`.
7. Split off EOG and EMG channels before PREP-style EEG-only cleaning.
8. Repair event annotations:
   - remove very close duplicate triggers;
   - advance labels for triggers judged to have the wrong label;
   - retain only the expected event names.
9. Save pre-cleaning PSD and channel plots for complete recordings.
10. Remove line noise using a CleanLine-like approach:
   - detrend with PyPREP `removeTrend`;
   - notch at 60 Hz harmonics up to Nyquist using MNE `notch_filter(..., method="spectrum_fit")`;
   - combine the cleaned trend with the original residual.
11. Run PyPREP robust re-referencing manually via `pyprep.reference.Reference`:
   - reference channels are all EEG channels;
   - RANSAC is enabled;
   - bad channels are interpolated when `interpolate_bads = True`.
12. Record original bad, interpolated, and still-noisy channels in `prep_info.csv`.
13. Flag sessions with more than 25% interpolated EEG channels as too noisy and route their final EDF output to a bad-output folder.
14. Re-append the EOG and EMG channels after EEG re-referencing/interpolation.
15. Filter the full prepped data from 1 to 50 Hz using MNE FIR filtering.
16. Fit ICA with 20 components using the Picard method, decimated by 5.
17. Detect ocular components with `ica.find_bads_eog()`.
18. If no EOG-related ICA components are found, save a partially processed FIF file under `output/eeg/bad/ica_err/` and skip EDF export for that subject.
19. Save ICA diagnostic plots.
20. Apply ICA exclusion to remove blink-related components.
21. If `perform_csd = True`, compute current source density after dropping still-noisy channels.
22. Save final preprocessed data as EDF+ using the local writer in `utils/save_edf.py`.
23. Append one row per processed subject to `output/prep_info.csv`.

## Expected Inputs And Produced Outputs

### Expected inputs

For `edf2bids.py`:

- Local raw EDF files in `external/DEMI_EEG_Pipeline/edfs/`.
- EDF filenames matching the project-specific `demi_(\d+)[\s\.]` pattern.
- Either EDF annotations matching numeric trigger strings or a stim channel readable by `mne.find_events()`.

For `eeg_pipeline.py`:

- A local BIDS-like dataset in `external/DEMI_EEG_Pipeline/bids/`.
- Subjects named `sub-XXX`.
- Task name `tracelab`.
- EEG datatype path conventions compatible with `mne_bids.BIDSPath()`.
- Expected event annotations named `stim_on`, `red_on`, `trace_start`, `trace_end`, `accuracy_submit`, and `vividness_submit`.
- Expected channel names for EOG/EMG from the conversion layer.

### Produced outputs

The external README and code indicate these outputs:

- `output/eeg/sub-XXX_eeg_prepped.edf`: final preprocessed EDF+ files for sessions not over the interpolation threshold.
- `output/eeg/bad/too_noisy/sub-XXX_eeg_prepped.edf`: final EDF+ files for sessions with more than 25% interpolated EEG channels.
- `output/eeg/bad/ica_err/sub-XXX_eeg_prepped.fif`: partially processed FIF files when ICA blink detection fails.
- `output/plots/sub_XXX/`: PSD, channel, and ICA diagnostic plots.
- `output/prep_info.csv`: subject-level preprocessing metadata, including annotation counts and bad/interpolated channel summaries.

The current DEMI R EEG import expects copied or relocated versions of these outputs at:

- `_Data/eeg/prep_info.csv`
- `_Data/eeg/edfs/sub-XXX_eeg_prepped.edf`

## Connection To Current R EEG Scripts

The current R EEG path starts after the external preprocessing pipeline has completed. `_Scripts/04_import_eeg.R` hard-fails if `_Data/eeg/prep_info.csv` or `_Data/eeg/edfs/*.edf` is missing.

The connection is:

1. Austin's external pipeline writes preprocessed EDF+ files and `prep_info.csv`.
2. Those files are placed locally under `_Data/eeg/`.
3. `_Scripts/04_import_eeg.R` reads `prep_info.csv` and each EDF with `eeguana::read_edf()`.
4. `_Scripts/04_import_eeg.R` uses `prep_info.csv` to identify remaining bad channels and, optionally, interpolated channels.
5. Because `_Scripts/_settings.R` currently has `using_csd <- TRUE`, `_Scripts/04_import_eeg.R` assumes still-noisy channels have already been dropped during CSD export and does not drop remaining bad channels again.
6. `_Scripts/04_import_eeg.R` separates EOG/EMG from EEG, drops bad EEG channels when relevant, reconstructs trial/event timing, creates `real_trace_start` and `real_trace_end`, epochs baseline/tracing/post-tracing EEG, epochs imagery EMG, downsamples EEG to 100 Hz, and writes `_Scripts/_rds/eeg_imported.Rds`.
7. `_Scripts/05_eeg_prep.R` reads `eeg_imported.Rds`, performs optional EMG-based imagery trial flagging, wavelet decomposition at 1:40 Hz, dB normalization, epoch relabeling to `during`/`after`, aggregation into theta/alpha/beta bands, and joins behavioural covariates.
8. `_Scripts/07_gam_prep.R` onward uses the compact band/epoch/sensor dataset for the prior GAM/maxT analysis.

The external pipeline therefore handled raw-signal preprocessing only. It did not perform DEMI-specific behavioural trial alignment, real movement onset/offset correction, epoch extraction, time-frequency decomposition, trial-level metadata joins, or statistical modelling.

## Decisions To Preserve, Reconsider, Or Replace

### Preserve or port carefully

- Event name vocabulary: `stim_on`, `red_on`, `trace_start`, `trace_end`, `accuracy_submit`, `vividness_submit`, plus derived `real_trace_start` and `real_trace_end`.
- Raw trigger code mapping from `edf2bids.py`, if the raw files still use the same codes.
- Channel type handling for EOG and EMG channels.
- Subject-level preprocessing manifest concept currently represented by `prep_info.csv`.
- Per-subject QC outputs showing PSD/channel traces before and after key steps.
- Bad/interpolated channel accounting.
- Explicit handling of incomplete recordings and ICA failures.
- The downstream event-alignment logic in `_Scripts/_functions/eeg.R`, especially physical-tracing false-start correction and the trial-count mapping to behavioural data.

### Reconsider before porting

- Whether CSD should be primary. The old external pipeline set `perform_csd = True`, and the R settings assume CSD. The new plan treats CSD as an open decision, likely sensitivity-only unless justified.
- The 1 Hz high-pass filter. It may be appropriate for ICA, but the final analysis filter settings should be decided for the planned time-frequency analyses and documented.
- The 50 Hz low-pass filter before time-frequency analysis up to 40 Hz. This likely remains defensible, but the new pipeline should document transition bands and filter design.
- The order of robust reference, filtering, ICA, and CSD. MNE-first implementation should follow current MNE/PyPREP practices rather than automatically matching 2020-era code.
- Fixed ICA settings: 20 components, Picard, decimation by 5, and automatic EOG selection.
- The 25% interpolated-channel threshold. Retain as historical reference, but define a current QC exclusion/interpolation policy before full reprocessing.
- The "complete recording" rule of at least 600 annotations. This is a useful historical check but should become an explicit event-count QC rule tied to expected trial counts.
- The event-repair heuristic that advances mislabeled triggers to the next event label. It may be necessary, but it is risky and should be tested against raw event logs or behavioural timing.
- Exporting preprocessed data to EDF+. For MNE-first work, FIF or BIDS derivatives may preserve richer metadata and avoid custom EDF writing.
- EMG-based imagery trial exclusion. The old R setting currently has `exclude_bad_by_emg <- FALSE`, so the old EMG logic is reference/QC material rather than an active exclusion policy.

### Replace rather than reuse directly

- The old Pipenv/Python 3.7 execution environment.
- The git-pinned PyPREP and Picard dependencies from the old `Pipfile`.
- The custom EDF writer as the primary derivative format.
- The old external preprocessing script as the active raw-to-preprocessed implementation.
- The downstream R wavelet/GAM path as the confirmatory EEG analysis path, while retaining it as historical reference.

## Dependencies, Formats, And Assumptions

### Python dependencies

The old pipeline is built around:

- Python 3.7 via Pipenv.
- MNE `0.21.2`.
- mne-bids `0.6`.
- PyPREP from a git commit (`447634eff3d2d7e49916d8073f46c75d6c9443d3`).
- python-picard from a git commit (`f40b97e600bfe6140eec9e5346934e6032a669ad`).
- matplotlib `3.3.2`.
- PyQt5 for plotting.
- numexpr, statsmodels, psutil, numpy, scipy as transitive/direct dependencies.

These versions are old relative to a new 2026 MNE-first analysis and should not be assumed compatible with modern MNE/PyPREP APIs.

### Data formats

- Raw input: EDF.
- Intermediate organization: BIDS-like EEG directory under `bids/`.
- Final old preprocessing output: EDF+.
- Problem-session output: FIF for ICA failures; EDF+ for too-noisy sessions.
- R downstream cache: RDS files.
- R downstream compact analysis table: `dat_bands_notime.rds`.

### Scientific and technical assumptions

- Raw EDF filenames encode participant IDs using the `demi_###` pattern.
- Raw events are recoverable from EDF annotations or a stim channel.
- Six event types define the basic task timeline.
- Recordings with fewer than 600 annotations are incomplete.
- The standard 1005 montage plus EEGBCI channel-name standardization is appropriate.
- EOG channels are sufficient for automatic blink-component detection.
- EMG/EOG channels can be removed before EEG PREP and re-appended afterward.
- CSD computation after ICA is acceptable and should drop remaining bad channels.
- EDF+ export can safely carry cleaned EEG/CSD, EOG, EMG, and annotations into R.
- `prep_info.csv` is the bridge for downstream channel exclusion decisions.

## Risky Or Unclear Steps

- The BIDS sidecar metadata contains placeholders such as hardware/software filter fields and cap details marked for confirmation.
- The raw-to-BIDS script intentionally ignores split EDF files. If split recordings are real data rather than duplicates, this is a critical limitation.
- The event-repair logic can silently remove duplicate annotations and relabel suspected wrong labels by advancing them to the next expected event. This may be correct for known acquisition issues, but it needs auditable event QC in the new pipeline.
- The `file start` crop is hard-coded around a known duplicate-data issue but generalized to any recording with that annotation.
- The completeness threshold is annotation-count based, not explicitly trial-count based.
- The pipeline modifies `raw_copy._data` directly during line-noise cleaning, which is fragile against MNE internals.
- The custom EDF writer depends on private-ish MNE fields such as `_raw_extras`.
- CSD output uses units that EDF readers and downstream R code may not fully understand semantically.
- The code saves too-noisy subjects as EDF+ under a bad folder but still produces a final-like file, so inclusion/exclusion depends on how outputs are copied into `_Data/eeg/edfs`.
- ICA failure is defined as finding zero EOG components, but there is no manual review gate in the scripted path.
- Plotting uses a Qt backend and may be brittle in headless or modern environments.
- The old pipeline records a random seed, but the main loop uses one fixed seed for all subjects unless `random_seed` is not provided.
- `edf2bids.py` uses an internal mne-bids helper `_write_json`, which is likely to be version-fragile.

## Side-By-Side Pipeline Map

| Old external preprocessing step | Old downstream R step | Proposed MNE-first replacement |
|---|---|---|
| Convert loose EDF files to BIDS-like layout with event/channel metadata | None; R consumes already-prepped EDF files | `analysis/eeg_mne/00_make_manifest.py` and `01_import_raw.py` should discover raw files, validate metadata, and write an auditable manifest/events table |
| Map trigger codes to event names | `_Scripts/_functions/eeg.R::update_events()` assumes named events | Build a tested event-construction module preserving code/name mapping and producing event-count QC |
| Mark HEO/VEO as EOG and EMG-L/EMG-A as EMG | `_Scripts/04_import_eeg.R` separates EOG/EMG from EEG | Use MNE channel typing during import; keep EOG/EMG available for ICA/QC and EMG checks |
| Apply standard montage and channel-name standardization | `_Scripts/05_eeg_prep.R` later joins channel coordinates from `BESA-81.csv` | Validate channel names and montage explicitly; decide authoritative coordinate source before preprocessing |
| Repair duplicate/mislabeled triggers | R event update drops practice/bad-start trials and creates real-trace events | Port repair rules only with event QC tables and tests against behavioural timing |
| Remove line noise with CleanLine-like notch procedure | No downstream line-noise handling | Use current MNE/PyPREP-compatible line-noise strategy; document notch frequencies and diagnostics |
| Robust reference with PyPREP RANSAC | R optionally drops channels listed in `prep_info.csv` | Use current PyPREP/MNE bad-channel detection and interpolation policy; export bad/interpolated channel summaries |
| Interpolate bad channels, flag sessions with >25% interpolated | R can optionally drop interpolated channels but currently does not | Define current bad-channel/interpolation exclusion thresholds before full run |
| Filter 1-50 Hz | R wavelet analysis assumes pre-filtered EDF input | Implement documented filters in MNE; consider separate ICA high-pass copy if needed |
| ICA blink correction using EOG | No downstream blink correction | Use MNE ICA workflow with explicit component-review/QC outputs and reproducible criteria |
| Compute CSD by default | R `using_csd <- TRUE`; drops no remaining bad channels because CSD path already dropped them | Decide whether CSD is primary or sensitivity branch; produce separate derivative if used |
| Write cleaned EDF+ and `prep_info.csv` | `_Scripts/04_import_eeg.R` reads EDF+ and `prep_info.csv` | Prefer MNE FIF/BIDS derivatives plus CSV/JSON QC manifests; only export EDF if needed for interoperability |
| None | `_Scripts/04_import_eeg.R` computes physical real trace start/end from behavioural/tracing data | Port this as `02_build_events.py`; preserve false-start and movement-duration logic with unit tests |
| None | `_Scripts/04_import_eeg.R` epochs baseline, tracing, post-trace EEG and imagery EMG | MNE epoching script with explicit baseline/during/after event definitions and retention QC |
| None | `_Scripts/05_eeg_prep.R` wavelet decomposition, dB normalization, band aggregation | MNE/Python time-frequency extraction and feature tables; R only if a chosen model requires it |
| None | `_Scripts/07_gam_prep.R` onward GAM/maxT inference | New planned primary models and cluster/sensitivity analyses; prior GAM remains historical/reference |

## Data-Locality And Reproducibility Implications

- The external submodule is useful provenance because it records the exact historical preprocessing code and pinned commit.
- Raw EDFs, BIDS conversion outputs, preprocessed EDFs, FIF files, plots, and `prep_info.csv` are local data/derivatives and should remain untracked.
- The current root `.gitignore` already protects `_Data/eeg/prep_info.csv`, `_Data/eeg/edfs/`, and common EEG binary formats.
- The new MNE workflow should keep raw recordings and generated derivatives separate.
- A reproducible manifest should record raw file paths, hashes if feasible, recording metadata, channel sets, event counts, and preprocessing decisions without committing raw data.
- The old pipeline's reproducibility depends on Python 3.7, old MNE/mne-bids versions, git-pinned dependencies, and a custom EDF writer. That makes exact reruns possible in principle but brittle in practice.
- For the new pipeline, source-controlled code plus environment files, local-only data manifests, and generated QC reports will be more reproducible than copying cleaned EDFs between repositories.

## Open Questions For Tony Before Implementation

1. Where are the authoritative raw EEG recordings stored locally, and are they loose EDFs, split EDFs, or already BIDS-like?
2. Should the new pipeline preserve the old trigger-repair heuristics exactly for the first pass, or should it first produce diagnostics comparing raw triggers with behavioural timing?
3. Are split EDF files expected for any valid participant sessions?
4. Is the old `demi_###` filename convention still the best way to identify subjects?
5. Should CSD be a sensitivity derivative rather than the primary signal?
6. Which channel coordinate source should be authoritative: recording metadata, standard 1005 montage, `_Data/eeg/BESA-81.csv`, or `_Scripts/_misc/sensor_latlong_chan_map.csv`?
7. Should a 1 Hz high-pass be used only for ICA, or also for the final time-frequency signal?
8. What is the preferred bad-channel policy: interpolate and retain, drop interpolated channels in sensitivity checks, or exclude sessions above a threshold?
9. Should the historical 25% interpolation threshold remain the exclusion threshold?
10. Should ICA component rejection be fully automated, manually reviewed, or semi-automated with saved reports?
11. Should EMG activity in imagery trials be used as exclusion, QC-only, or sensitivity analysis?
12. What derivative format should be primary: MNE FIF, BIDS derivatives, exported EDF, or feature tables only?
13. Should the new event table preserve both original trigger times and corrected `real_trace_start`/`real_trace_end` times for auditability?
14. Which outputs should be tracked as lightweight reproducibility artifacts, and which should remain local-only derivatives?
15. Should the external submodule remain indefinitely as provenance once the new MNE pipeline is implemented, or eventually move under a labeled archive path?
