# Raw-channel QC and montage contract

The active MNE reanalysis keeps raw-channel QC descriptive and read-only. It does not automatically mark bad channels or perform filtering, referencing, interpolation, ICA, CSD, or epoch construction.

## Channel units

The EDF signal header is the authority for original physical dimensions and digital/physical calibration ranges. In all inspected DEMI EDFs, the 32 scalp channels declare microvolts. MNE returns those values in volts. HEO, VEO, EMG-L, EMG-A, and Trigger have blank physical-dimension fields.

Changing HEO/VEO to EOG and EMG-L/EMG-A to EMG in MNE is a channel-type correction, not evidence of a voltage calibration. Their absolute amplitudes and PSD magnitudes therefore remain unknown acquisition units and must not be compared across files as calibrated physical measurements. Within-channel/file morphology, frequency shape and ratios, flatness, dropout, step, and digital-rail metrics remain useful QC evidence.

Trigger is a digital event channel. It is summarized through states and transitions, separately from continuous amplitude and PSD interpretation.

Script 07 reads every EDF one file at a time, scans the full duration in deterministic chunks, and writes local per-channel evidence under `_Data/eeg/mne_preprocessing/channel_qc_v1/`:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/07_run_unit_aware_channel_qc.py
```

Candidate rows and ranks are descriptive review aids, not bad-channel assignments.

## Historical coordinate source

The tracked `_Data/eeg/BESA-81.csv` is an exact ordered label and six-decimal coordinate match to Robert Oostenveld's [`besa_81.txt`](https://robertoostenveld.nl/download/electrode/besa_81.txt). The [source page](https://robertoostenveld.nl/electrode/) describes it as Cartesian positions according to BESA on a unit sphere.

The table contains 112 rows, including electrodes, accessory locations, legacy aliases, and template fiducials. Its orientation is x toward the right, y anterior toward the nasion, and z superior toward Cz. NAS, LPA, and RPA are present; Nz, T9, and T10 duplicate those template locations.

This provenance does not establish participant electrode digitization, participant head shape, a physical head radius, or an acquisition-specific transform. The recovered historical R path used the table for spherical plotting and spatial model coordinates. The historical external Python preprocessing used MNE `standard_1005` instead.

Script 08 validates the tracked source contract and creates a documented comparison with `standard_1005`:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/08_audit_montage_provenance.py
```

The comparison is audit-only. A fitted scale is not a recovered DEMI head radius.

## Use boundary

- Channel-name validation and historical unit-sphere visualization are supported.
- The historical table is not approved as-is for interpolation, CSD, acquisition-geometry claims, or quantitative spatial inference.
- The active preprocessing montage is MNE `standard_1005`.
- The historical BESA unit-sphere surface is retained only for historical GAM/visualization provenance until its references are audited during later documentation cleanup.

The machine-readable provenance and safety guard are in `analysis/eeg_mne/montage_coordinate_contract_v1.yaml` and `analysis/eeg_mne/montage_contract.py`.
