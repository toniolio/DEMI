# EEG Reanalysis Plan

## Purpose

This document defines the planned reanalysis of the EEG component of the DEMI project. The goal is to reprocess the EEG data from raw recordings using a transparent, reproducible, MNE-Python-centered workflow, while preserving the already-published behavioural analysis as the behavioural foundation for the EEG study.

The revised EEG analysis should answer the scientific question directly: whether electrophysiological signatures associated with overt movement error and post-movement evaluation are also present during kinesthetic motor imagery.

## Scope

In scope:

- Reprocess raw EEG recordings using a reproducible Python/MNE pipeline.
- Preserve trial/event definitions from the original experiment wherever possible.
- Join EEG epochs to existing behavioural trial metadata.
- Extract planned time-frequency features for theta, alpha, and beta bands.
- Test planned hypotheses using conventional EEG/statistical approaches.
- Generate auditable QC outputs, including participant-level reports.
- Archive the previous EEG/GAM analysis so it remains available for reference but is not confused with the active analysis.

Out of scope:

- Reworking the published behavioural analysis.
- Reinterpreting behavioural results except as needed to construct EEG trial metadata.
- Treating the prior GAM analysis as confirmatory evidence for the new manuscript.
- Changing hypotheses in response to observed EEG results.

## Analytical Principles

1. The behavioural analysis is treated as frozen unless a technical issue prevents construction of EEG trial metadata.
2. EEG preprocessing should be reproducible, auditable, and based on established EEG/MNE workflows.
3. The primary EEG analysis should use planned hypotheses, bands, epochs, and sensor/ROI definitions.
4. Whole-scalp or time-frequency analyses may be used as supportive/exploratory evidence, but the primary claims should not depend on bespoke modelling machinery.
5. Statistical models should preserve meaningful trial- and participant-level variability where feasible.
6. Simpler inferential machinery is preferred only when it remains scientifically appropriate; “simple” should not mean collapsing away the design.
7. Every major methodological decision should be traceable to prior literature, package documentation, or an explicit documented rationale.

## Previous EEG/GAM Analysis

The previous EEG analysis used hierarchical generalized additive models with spatial smoothing over scalp coordinates and a global max-|T| procedure for planned contrasts. This work should be retained as historical/supplemental reference but removed from the active analysis path.

Recommended handling:

- Move previous GAM-specific scripts and outputs into an archive location.
- Preserve enough documentation to explain what the archived analysis was.
- Mark the archived pipeline as superseded by the raw-MNE reanalysis.
- Do not delete prior work unless redundant generated artifacts are too large or already reproducible.

Suggested archive path:

```text
archive/eeg_gam_analysis/
```

## Proposed Active Workflow

```text
analysis/
  behaviour/
    # Existing/frozen behavioural scripts or wrappers

  eeg_mne/
    00_make_manifest.py
    01_import_raw.py
    02_preprocess.py
    03_epoch.py
    04_qc_reports.py
    05_time_frequency.py
    06_extract_features.py
    07_primary_models.py
    08_cluster_support.py
    09_figures.py

docs/
  eeg_reanalysis_plan.md
  eeg_methods_decisions.md

derivatives/
  behaviour/
    trial_metadata_for_eeg.csv

  eeg/
    manifest/
    preprocessed/
    epochs/
    time_frequency/
    features/
    qc/
    stats/
    figures/
```

Exact folder names may be revised after repository audit.

## Preprocessing Plan

The preprocessing pipeline should be implemented in MNE-Python and should include:

- raw file import;
- channel name and montage validation;
- event extraction and validation;
- filtering;
- bad channel detection and documented handling;
- referencing;
- ICA or equivalent artifact correction for ocular artifacts;
- epoching around task events;
- participant-level QC reporting;
- export of cleaned epochs and QC summaries.

A CSD-transformed branch may be generated as a sensitivity analysis. The main analysis should initially use a conventional referenced EEG signal unless the methods review supports making CSD primary.

## Behaviour-to-EEG Metadata Contract

The EEG pipeline should consume a stable behavioural trial metadata table with one row per trial. Required fields likely include:

- participant ID;
- group;
- block;
- trial number;
- trial type;
- stimulus duration or movement time target;
- complexity/sinuosity;
- self-rated accuracy;
- objective error where available;
- expected performance where available;
- event timing keys needed to align EEG epochs;
- inclusion/exclusion flags.

The expected output is:

```text
derivatives/behaviour/trial_metadata_for_eeg.csv
```

## Planned EEG Features

Primary features should be baseline-corrected time-frequency power in planned bands and epochs:

- theta: approximately 4–8 Hz;
- alpha: approximately 8/9–12 Hz;
- beta: approximately 13–30 Hz.

Planned epochs:

- during movement / imagery, time-locked to response or imagery onset;
- post movement / imagery, time-locked to response or imagery offset.

The exact baseline window, epoch limits, frequency method, and decibel normalization procedure should be finalized before running the full analysis.

## Primary Hypotheses

### H1: Frontal theta and low accuracy

Test whether low-accuracy trials differ from high-accuracy trials in frontal theta power during movement/imagery, and whether that effect differs between overt movement and motor imagery.

### H2: Posterior alpha and low accuracy

Test whether low-accuracy trials differ from high-accuracy trials in posterior alpha power after movement/imagery, and whether that effect differs between overt movement and motor imagery.

### H3: Post-movement beta rebound

Test whether post-movement beta power exceeds during-movement beta power over sensorimotor electrodes, and whether this rebound differs between overt movement and motor imagery.

## Statistical Strategy

The preferred primary analysis should preserve participant- and trial-level structure. Candidate approaches include:

1. Bayesian hierarchical models on extracted trial-level ROI features.
2. Frequentist mixed-effects models if Bayesian models are impractical.
3. Subject-level planned contrasts as a transparent robustness/sensitivity analysis.
4. Cluster-based permutation analyses as supportive whole-scalp/time-frequency evidence.

The analysis should avoid making the paper depend on the prior GAM spatial-smoothing approach.

## Quality Control Outputs

The pipeline should produce:

- participant-level preprocessing reports;
- event count summaries;
- trial retention summaries;
- bad/interpolated channel summaries;
- ICA/artifact correction summaries;
- epoch rejection summaries;
- feature extraction completeness checks;
- final included participant/trial counts by group and condition.

## Open Decisions

Before implementation, decide:

- exact raw EEG file format and local data path conventions;
- whether CSD is primary or sensitivity-only;
- bad-channel detection/interpolation policy;
- ICA component rejection policy;
- exact event definitions for onset/offset epochs;
- exact baseline window;
- exact ROI channel lists;
- whether accuracy is modeled continuously, binned, or both;
- primary statistical model family;
- how archived GAM scripts should be relocated without breaking historical reproducibility.