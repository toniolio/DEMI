# Continuous preprocessing parameter audit

The final parameter evidence pass is separate from production preprocessing. It reads raw EDF recordings into memory, writes only local CSV/JSON/Markdown evidence, and never writes filtered, referenced, interpolated, ICA-cleaned, CSD, or epoched signals.

The machine-readable contract is `analysis/eeg_mne/preprocessing_parameter_audit_contract_v1.yaml`. Participant/file selection and participant-level interpretation remain private.

The accepted source-pool policy uses the 30 scalp channels for global detection and average-reference estimation. M1/M2 remain recorded provenance channels and reference targets, are excluded from the primary theta/alpha/beta feature surface, and remain available for a later mastoid sensitivity analysis. Their acquisition-reference role is not recoverable from the available primary sources.

## Common detector input

MNE LOF, PyPREP `NoisyChannels`, and full PREP receive independent copies of the same full-session input:

- explicit HEO/VEO EOG, EMG-L/EMG-A EMG, and Trigger stim typing;
- EEG channels only;
- MNE `standard_1005` montage;
- audit-only resampling to 250 Hz, retaining 60 Hz;
- no reference, notch, or analysis filter before detector comparison.

LOF records channel scores at its documented primary setting and fixed sensitivity settings. `NoisyChannels` records flat, nonfinite, deviation, high-frequency-noise, correlation, dropout, SNR, PSD, and RANSAC findings separately. Full PREP records the integrated line-removal, iterative robust-reference, detection, interpolation, rereference, and post-interpolation states separately even though PyPREP couples them internally.

Repeated runs use a fixed seed. Thresholds are fixed across recordings; the audit does not tune participants individually. Historical channel records, when present, are comparison evidence rather than ground truth.

AutoReject is not used as a continuous global detector. `autoreject.Ransac` requires an Epochs input, so running it here would require constructing artificial epochs and duplicate PREP RANSAC without a concrete validation purpose. Full AutoReject remains a possible later epoch-level transient-channel/epoch repair step after the accepted epoch surface exists.

## Line-noise/filter comparison

The compact comparison evaluates no explicit 60-Hz removal versus MNE `spectrum_fit` removal at 60 Hz, followed in both cases by the same proposed 0.5--45 Hz zero-phase FIR. It reports the filter response, edge exclusion, residual 60-Hz PSD, theta/alpha/beta PSD, and final waveform differences. The cutoffs remain proposals until the private policy is accepted and production preprocessing is separately implemented.

Run the evidence scripts in order:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/09_audit_historical_prep_and_reference.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/10_compare_global_bad_channel_methods.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/11_compare_line_noise_filters.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/12_audit_m1_m2_detector_sensitivity.py
```

Generated evidence belongs under `_Data/eeg/mne_preprocessing/preprocessing_parameter_audit_v1/` and is not committed.
