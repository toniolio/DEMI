# Active MNE EEG preparation path

The active EEG reanalysis starts from the raw EDF recordings and keeps event
evidence separate from preprocessing and epoch construction.

## Identity and event-source evidence

The canonical tracked identity surface is
`eeg_behavior_identity_contract.csv`. It contains one explicit row per raw EEG
recording ID. Code must not infer behavioural identity from numeric equality
when a row is missing. The reviewed exception maps EEG recording 13 to
behavioural participant 11; EEG recording 11 is explicitly unmapped.

`event_source_contract.py` contains the pure, tested identity, source
comparison, and conservative selection rules.

Run the source and corrected-alignment evidence in this order:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/05_inventory_event_sources.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/06_build_corrected_event_evidence.py
```

Script 05 inventories EDF annotations and physical Trigger transitions for
every raw file. It writes local versioned evidence under
`_Data/eeg/event_source_inventory_v1/`.

Script 06 applies the explicit identity mapping and selected source, then
reuses the historical in-memory cleanup and timing-linkage helpers. It writes
local versioned evidence under `_Data/eeg/event_evidence_v1/`.

Scripts 01–03 and their existing output directories preserve the earlier
annotation-only audit history. They are not silently overwritten and should
not be treated as the corrected event-evidence surface.

## Other preparation scripts

- Script 00 creates the raw EDF header/annotation manifest.
- Script 04 is the bounded continuous-readability/raw-QC driver. It does not
  perform the identity/source remediation.
- Script 07 performs full-duration, unit-aware, per-channel descriptive QC. It
  writes local versioned evidence under
  `_Data/eeg/mne_preprocessing/channel_qc_v1/` and never assigns bad channels.
- Script 08 validates the tracked BESA unit-sphere provenance contract and
  compares the 32 scalp labels with MNE `standard_1005`.
- Scripts 09-11 perform the read-only historical/reference, automated global
  bad-channel, and line-noise/filter parameter audits. They write only local
  tabular/Markdown evidence and discard all in-memory signal changes.

Run the two preprocessing-foundation audits with:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/07_run_unit_aware_channel_qc.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/08_audit_montage_provenance.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/09_audit_historical_prep_and_reference.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/10_compare_global_bad_channel_methods.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/11_compare_line_noise_filters.py
```

The public unit/montage and parameter-audit boundaries are documented in
`docs/eeg_raw_qc_and_montage_contract.md` and
`docs/eeg_preprocessing_parameter_audit.md`.

None of these scripts filters EEG, applies a reference, interpolates channels,
runs ICA, computes CSD, writes cleaned EEG, makes inclusion decisions, or
constructs final epochs.

Focused public tests are in `tests/test_event_source_contract.py`,
`tests/test_channel_qc.py`, `tests/test_montage_contract.py`, and
`tests/test_preprocessing_parameter_audit.py`.
