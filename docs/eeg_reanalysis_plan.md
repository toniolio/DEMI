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

The current tracked workflow stops at read-only inventory, linkage, event
evidence, raw-channel/montage QC, and preprocessing-parameter audits. It has
not yet written production-preprocessed EEG or constructed epochs.

This document is intentionally high-level. Detailed working notes, local data
inventories, accepted event and preprocessing policy records, and temporary
planning materials remain private. Tracked identity, montage, and parameter
contracts plus public tests document the reproducible implementation
boundaries without exposing those private scientific decision records.
