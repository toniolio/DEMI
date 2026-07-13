# Continuous EEG preprocessing contract

The production namespace is `analysis/eeg_mne/continuous_preprocessing/`, with
the versioned technical configuration in
`analysis/eeg_mne/continuous_preprocessing_config_v1.yaml`. Script 13 runs a
small deterministic validation cohort. It is not authorized to process all 95
recordings.

## Signal contract

The pipeline retains 32 EEG channels, including M1 and M2, plus HEO, VEO,
EMG-L, EMG-A, and Trigger. The fixed 30-channel scalp surface owns PyPREP
global-bad detection and average-reference estimation. M1/M2 are reference
targets but not detector/reference sources or primary interpolation candidates.
The active montage is MNE `standard_1005`.

The analysis branch uses a 0.5--45 Hz zero-phase `firwin` FIR. The line-noise
toggle is fixed for a run and defaults off; its enabled branch uses 60 Hz MNE
`spectrum_fit` removal with a 10 s filter length, 2 Hz multitaper bandwidth,
and p = 0.01. The detector receives a separate native-rate, full-session,
unfiltered, unnotched, unreferenced 30-channel copy.

PyPREP `NoisyChannels` runs twice with seed 20260712, detrending, correlation,
and RANSAC enabled. The authoritative union contains nonfinite, flat,
deviation, high-frequency, correlation, SNR, dropout, and RANSAC criteria.
PSD-only findings remain report-only. Accepted bad scalp channels are excluded
from the explicit average-reference sources, the reference is applied to all
32 EEG targets, and accepted bad scalp channels are then interpolated with MNE
spherical splines. A count of 8/30 or greater is an objective stop.

ICA is fit on a separate 1--45 Hz copy with the same line-noise, reference, and
interpolation surface. The tracked implementation uses rank-aware extended
Infomax, seed 20260712, and fixed 4:1 fitting decimation (effective 250 Hz for
the 1000 Hz DEMI recordings). HEO and VEO are scored independently with MNE's
fixed z-score rule. Zero proposals is valid; more than two proposals stops the
recording. ID 86 always stops at the component-review boundary without an
automatic post-ICA derivative.

## Saved derivatives

Outputs are new files below the ignored directory:

```text
_Data/eeg/mne_preprocessing/continuous_validation_v1/
```

Each published recording directory is atomic and self-contained. Ordinary
terminal results retain one canonical post-ICA continuous FIF and the ICA
object. Review stops retain a pre-ICA continuous FIF and technically valid ICA
evidence where available, but no automatic post-ICA derivative. This avoids a
second full continuous copy for every recording; the alternate branch remains
reconstructible from the immutable EDF, tracked code/configuration, hashes,
and ICA exclusion ledger.

FIF precision is explicitly `single`, not an implicit default. The validation
run records observed precision error and projected storage before this becomes
the long-term full-run contract. Detector tables, component rows, stage
ledgers, source/config/code/environment hashes, runtime, memory, byte counts,
and output checksums accompany each signal derivative.

## Run and resume

From the repository root:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py
```

Use `--max-files 2` for a bounded interruption/resume check. On a later run,
complete and stopped results are skipped only when source, configuration, code,
environment, and every recorded output hash still match. Stale, corrupt, or
incomplete directories are preserved in local history and recomputed. Use an
exact cohort filename with `--recording` and add `--force` for an explicit
focused recomputation; the prior result is preserved rather than silently
overwritten.

Processing is sequential. A per-recording stop or failure is recorded and the
driver continues to the next EDF. Raw EDF hashes, sizes, and modification times
are checked after processing. The driver has no all-recording mode.

## Current boundary

The pipeline does not repair EDF events, construct epochs, run AutoReject,
compute CSD, alter accepted event-policy eligibility, or make participant
inclusion/exclusion decisions. Successful continuous preprocessing does not
imply event availability or epoch eligibility. Epoch construction will be a
separate future stage.
