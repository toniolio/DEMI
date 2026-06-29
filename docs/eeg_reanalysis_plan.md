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

This document is intentionally high-level. Detailed working notes, local data inventories, and temporary planning materials are maintained outside the public repository.
