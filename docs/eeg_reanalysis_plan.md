# EEG Reanalysis Plan

This repository is being updated for a new analysis of the EEG component of the DEMI study.

The behavioural analysis associated with the published behavioural paper is treated as complete. The current work focuses on reprocessing and reanalyzing the EEG data using a modern Python/MNE-centered workflow, starting from the raw EEG recordings where possible.

The revised EEG work is expected to include:

- local inventory of raw EEG and behavioural source files;
- reproducible EEG preprocessing from raw recordings;
- reconstruction and validation of task events;
- trial-level linkage between EEG data and behavioural metadata;
- time-frequency feature extraction for planned EEG bands and epochs;
- transparent quality-control outputs;
- planned statistical analyses focused on the study hypotheses.

Older EEG preprocessing and analysis code is retained for provenance and comparison while the new workflow is developed. Once the new workflow is sufficiently complete, historical EEG analysis code may be archived or more clearly separated from the active analysis path.

The active code and reproducible preparation order are documented in
[`analysis/eeg_mne/README.md`](../analysis/eeg_mne/README.md). The supporting
scripts in [`analysis/behavior/`](../analysis/behavior/README.md) reconstruct
the frozen task/TraceLab linkage required for EEG event alignment; they do not
constitute a new behavioural reanalysis. The earlier R/GAM workflow under
[`_Scripts/`](../_Scripts/README.md) and the pinned external preprocessing
submodule are provenance sources rather than active processing entry points.

The tracked workflow now separates scripts 00--12, which provide read-only
inventory and parameter evidence, from script 13, which implements production
continuous preprocessing, script 14, which implements the accepted
event/epoch eligibility ledger, and script 15, which constructs the accepted
epoch surface. Script 16 implements the accepted retain-and-diagnose
transient-artifact contract and constructs trial-level Morlet power plus
trial-matched `red_on` dB normalization. Continuous preprocessing is closed: 94
files are complete, with one separately retained accepted historical ICA stop.
The ledger records 8,905 primary, 8,896 strict-clean-only, and 8,798 ordinary
future-ready rows. The completed epoch stage contains 8,798 response-onset,
8,798 response-end, and 8,798 `red_on` support epochs, with 8,789 strict-clean
rows identifiable in each family. The completed time-frequency stage preserves
all 8,798 onset and end trials on a 30-channel, 4--40-Hz, -0.5-to-+1.5-s
surface; artifact flags are diagnostic metadata and do not change eligibility.
Script 17 now provides conventional trial-level fixed-band/window features for
all 30 physical channels and equal-weight predeclared frontal, posterior, and
task-hand-normalized motor ROIs. It retains all 8,798 trials, with dB primary
and unnormalized log power as a sensitivity. No statistical model or
inferential result has yet been produced. Model-ready predictor representation,
statistical modelling, inferential contrasts, CSD sensitivity, and scientific
interpretation remain separate future stages; a whole-scalp GAMM is reserved
only as a possible later spatial/historical sensitivity.

This document is intentionally high-level. Detailed working notes, local data
inventories, accepted event and preprocessing policy records, and temporary
planning materials remain private. Tracked identity, montage, and parameter
contracts plus public tests document the reproducible implementation
boundaries without exposing those private scientific decision records.

The public technical contract for the saved continuous derivatives is
[`eeg_continuous_preprocessing_contract.md`](eeg_continuous_preprocessing_contract.md).
