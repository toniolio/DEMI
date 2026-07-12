# Active MNE EEG preparation path

This is the active EEG reanalysis path. It starts from raw EDF recordings and
currently provides inventory, frozen-behaviour linkage, event evidence,
descriptive raw-channel QC, montage provenance, and preprocessing-parameter
audits. The earlier EEG/GAM code in [`../../_Scripts/`](../../_Scripts/README.md)
and the pinned external preprocessing submodule are historical evidence, not
active processing entry points.

No script in this directory currently writes production-preprocessed EEG or
constructs epochs. Scripts 00–12 are inspection/evidence scripts. Script 04's
historical filename notwithstanding, it is a config-driven raw-QC driver and
does not preprocess signals.

## Environment and local inputs

From the repository root:

```sh
python3 -m venv .venv
./.venv/bin/python -m pip install -r requirements-eeg.txt
```

Raw and behavioural source data are not distributed with the repository.
Expected local placement is:

```text
_Data/eeg/raw/     raw EDF recordings
_Data/task/        behavioural task text files
_Data/figure/      TraceLab figure/tracing archives
```

All generated evidence stays below ignored `_Data/behavior/` or `_Data/eeg/`
paths. Private configurations, inventories, reviews, and accepted scientific
policy records stay below ignored `_Private/`. A collaborator who needs to run
config-dependent audits must obtain the reviewed local configuration from the
project owner; public code must not guess or replace it.

## Active preparation order

### 1. Reconstruct the frozen behavioural/tracing linkage

Run scripts 00–04 in [`../behavior/`](../behavior/README.md). They inventory
task/TraceLab files, reconstruct tracing tables and historical filters, and
write the old-compatible event-offset input under
`_Data/behavior/event_offsets/`. This is EEG linkage infrastructure, not a new
behavioural analysis.

### 2. Inventory raw EEG

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/00_make_raw_eeg_manifest.py
```

This writes header and annotation manifests under `_Data/eeg/manifest/`.

### 3. Preserve the annotation-only audit lineage when required

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/01_compare_raw_annotations_to_offsets.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/02_audit_raw_event_sequences.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/03_audit_event_alignment_special_cases.py
```

Scripts 01–03 created the first annotation-only alignment evidence. Their
outputs are retained for provenance and for optional regression comparison,
but they are superseded as the current event-evidence surface by scripts 05–06.
They must not be treated as epoch inputs.

### 4. Run bounded raw-QC inspection

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/04_preprocess_continuous_raw.py
```

Script 04 requires the reviewed private raw-QC configuration. It opens EDFs
read-only, validates headers/channel types, and writes sampled descriptive
evidence under `_Data/eeg/mne_preprocessing/raw_qc/`. It does not filter,
reference, interpolate, run ICA, repair events, or write cleaned EEG.

### 5. Build the corrected event-evidence surface

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/05_inventory_event_sources.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/06_build_corrected_event_evidence.py
```

Script 05 inventories EDF annotations and physical Trigger transitions under
`_Data/eeg/event_source_inventory_v1/`. Script 06 applies the tracked explicit
identity contract and conservative source-selection rules, reuses the
historically grounded in-memory cleanup/linkage helpers, and writes the current
corrected evidence under `_Data/eeg/event_evidence_v1/`.

### 6. Run channel-QC and montage audits

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/07_run_unit_aware_channel_qc.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/08_audit_montage_provenance.py
```

Script 07 scans full recordings descriptively without assigning bad channels.
Script 08 validates the historical Oostenveld/BESA unit-sphere table and
compares its labels with the approved active `standard_1005` template. Outputs
belong under `_Data/eeg/mne_preprocessing/channel_qc_v1/`. See
[`../../docs/eeg_raw_qc_and_montage_contract.md`](../../docs/eeg_raw_qc_and_montage_contract.md).

### 7. Run preprocessing-parameter evidence audits

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/09_audit_historical_prep_and_reference.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/10_compare_global_bad_channel_methods.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/11_compare_line_noise_filters.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/12_audit_m1_m2_detector_sensitivity.py
```

These scripts audit the historical reference/PREP path, compare fixed
bad-channel methods, compare line-noise/filter branches, and isolate the M1/M2
detector-pool effect. They write local evidence under
`_Data/eeg/mne_preprocessing/preprocessing_parameter_audit_v1/` and discard
all in-memory signal changes. See
[`../../docs/eeg_preprocessing_parameter_audit.md`](../../docs/eeg_preprocessing_parameter_audit.md).

## Script index

| Script | Current role | Status |
| --- | --- | --- |
| 00 | Raw EDF header/annotation manifest | Active inventory |
| 01 | Raw annotation versus event-offset comparison | Superseded evidence lineage |
| 02 | Read-only annotation sequence cleanup/alignment audit | Superseded evidence lineage; helpers reused by 06 |
| 03 | Split, concatenated, and raw-only special-case audit | Superseded evidence lineage |
| 04 | Config-driven continuous raw-QC inspection | Active inspection; not preprocessing |
| 05 | Dual-source annotation/physical-Trigger inventory | Active event evidence |
| 06 | Identity-corrected, source-selected alignment evidence | Current event-evidence surface |
| 07 | Full-duration unit-aware channel QC | Active descriptive evidence |
| 08 | Historical coordinate provenance and montage audit | Active contract audit |
| 09 | Historical PREP/reference recovery | Active parameter evidence |
| 10 | Global bad-channel method comparison | Active parameter evidence |
| 11 | Line-noise/filter branch comparison | Active parameter evidence |
| 12 | M1/M2 detector-pool sensitivity | Active parameter evidence |

## Scientific-policy boundary and stopping point

Accepted event policy v1.0 and preprocessing parameter policy v1.0 are private
scientific decision records. Superseded drafts are retained privately for
provenance. Public, reproducible implementation boundaries remain tracked in:

- `eeg_behavior_identity_contract.csv` and `event_source_contract.py`;
- `montage_coordinate_contract_v1.yaml` and `montage_contract.py`;
- `preprocessing_parameter_audit_contract_v1.yaml` and
  `preprocessing_parameter_audit.py`;
- the public contract documents under [`../../docs/`](../../docs/).

The next separately authorized implementation step is a no-write dry run of
the accepted production continuous-preprocessing stages with complete
provenance. The present workflow stops before that implementation. It does not
claim cleaned continuous derivatives, an epoch ledger, epochs, time-frequency
features, or new EEG results.

## Validation

The public test suite does not require private raw EEG:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 -m pytest
python3 tools/check_repo_safety.py
git status --short
```

Focused contracts are covered by `tests/test_event_source_contract.py`,
`tests/test_channel_qc.py`, `tests/test_montage_contract.py`,
`tests/test_preprocessing_parameter_audit.py`, and
`tests/test_raw_qc_guardrails.py`.
