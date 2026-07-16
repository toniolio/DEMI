# Active MNE EEG preparation path

This is the active EEG reanalysis path. It starts from raw EDF recordings and
provides inventory, frozen-behaviour linkage, event evidence, descriptive
raw-channel QC, preprocessing-parameter audits, completed versioned continuous
preprocessing, an accepted-policy event/epoch eligibility ledger, and the
completed accepted epoch surface. The
earlier EEG/GAM code in
[`../../_Scripts/`](../../_Scripts/README.md)
and the pinned external preprocessing submodule are historical evidence, not
active processing entry points.

Scripts 00--12 are inspection/evidence scripts. Script 04's historical filename
notwithstanding, it is a config-driven raw-QC driver and does not preprocess
signals. Script 13 is the completed production continuous-preprocessing driver.
Script 14 joins the accepted event-policy surface to continuous provenance and
records eligibility/readiness statuses without opening signal data. Script 15
constructs and validates the accepted response-onset, response-end, and
`red_on` epoch families without applying artifact rejection or spectral
processing.

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

### 8. Run saved continuous preprocessing

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py
```

The tracked contract is `continuous_preprocessing_config_v1.yaml`; focused
production modules are under `continuous_preprocessing/`. The default command
writes the separate validation cohort atomically below
`_Data/eeg/mne_preprocessing/continuous_validation_v1/`.

Use `--max-files 2` for a bounded interruption/resume check. Valid terminal
recordings are skipped only after source/config/code/environment and artifact
hashes match. Use `--recording 'demi_01 Data.edf' --force` for an explicit
focused recomputation; the old result is preserved in local history.

The authoritative v2 all-readable-EDF production run uses a self-locating
launcher that prevents fragile pasted multiline commands, keeps the Mac awake,
streams progress, and saves a timestamped terminal log:

```sh
tools/run_continuous_v2.sh
```

It selects all 95 script-00 inventory rows with `read_status=ok` in stable
recording-ID/filename order, processes split parts independently, keeps the
60-Hz branch off, and writes only below
`_Data/eeg/mne_preprocessing/continuous_v2/`. The validation and v1 production
roots remain untouched. It requires at least 40 GiB free, checks that all roots
are ignored and separate, hashes the complete raw surface before and after each
invocation, and flushes and prints progress after every recording.

Use `--all-recordings --recording 'demi_01 Data.edf'` for a one-file production
smoke. Add `--force` only with one exact `--recording`. To independently reopen
and rescan every current production FIF/ICA artifact after a run:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py --all-recordings --verify-current
```

Rerun the full command unchanged to prove cache skipping; the final lines print
invocation and current terminal counts. A successful surface is expected to
contain 94 complete results, including files 49 and 54_1 marked
`complete_with_qc_warning`, plus the accepted ID-86 historical ICA stop.

The legacy-v1 historical-ICA routing repair reuses the saved pre-ICA FIF, ICA
object, rank, scores, and provenance from results stopped by the superseded
over-two guardrail. It does not rerun earlier preprocessing or ICA fitting:

```bash
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py --all-recordings --repair-historical-ica-routing
```

The active EOG rule reproduces the historical executable behavior: retain the
score-ranked, deduplicated result returned across HEO and VEO and apply every
detected component without a numeric cap. ID 86 remains on its accepted
historical-stop route.

The historical interpolation boundary is also preserved. More than 25% of the
30 scalp channels produces a prominent QC warning and continues through
reference, interpolation, filtering, ICA, and derivative writing. This is not
an event, epoch, participant, or analytic eligibility decision.

Ordinary results retain one post-ICA continuous FIF plus the ICA object. A
component-review stop retains the pre-ICA continuous FIF and available ICA
evidence but no automatic post-ICA file. See
[`../../docs/eeg_continuous_preprocessing_contract.md`](../../docs/eeg_continuous_preprocessing_contract.md).

### 9. Build the accepted event/epoch eligibility ledger

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/14_build_event_epoch_eligibility_ledger.py
```

Script 14 consumes the accepted local event policy, corrected selected-event
evidence, explicit identity/source contracts, and continuous-v2 manifests. It
writes a versioned local row-level ledger and review summaries below
`_Data/eeg/event_epoch_eligibility_v1/`. The ledger reproduces the 8,905-row
primary and 8,896-row strict-clean-only event-policy surfaces and keeps event
candidacy, continuous availability, post-ICA availability, and future epoch
readiness distinct. It does not create epochs or decide analytic inclusion.

### 10. Construct the accepted epoch families

```sh
tools/run_epochs_v1.sh
```

Script 15 consumes exactly the 8,798 ordinary future-ready ledger rows and only
their linked `continuous_v2` post-ICA FIFs. It writes bounded, resumable
per-recording Epochs shards below `_Data/eeg/epochs_v1/` for response onset,
response end, and separate `red_on` spectral-reference support. Task epochs use
-1.5 to +2.5 s; `red_on` support uses -1.5 to +0.8 s. All retain native 1000 Hz
sampling, include both endpoints, and use `baseline=None`.

The same canonical ledger ordering is preserved in every family. The 8,789
strict-clean rows are provided as metadata/indices over the primary signal
surface, while the nine accepted duration-warning rows remain included and
flagged. File 49 remains included with its continuous QC warning; file 54_1 and
ID 86 produce no epochs under their accepted event/derivative routes. The
stage reopens every shard and verifies source-FIF immutability. It does not run
AutoReject, amplitude rejection, CSD, resampling, TFRs, band-power extraction,
spectral normalization, ROI selection, participant exclusion, or statistics.

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
| 13 | Saved production continuous preprocessing | Complete: 94 ordinary completions plus the accepted ID-86 stop |
| 14 | Accepted event/epoch eligibility ledger | Complete; no epochs constructed |
| 15 | Accepted response-onset, response-end, and `red_on` Epochs | Complete: 8,798 epochs per family |

Scripts 00--12 remain evidence/audit programs. Script 13 owns production
continuous preprocessing, script 14 owns the policy ledger, and script 15 owns
the versioned accepted epoch namespace.

## Scientific-policy boundary and stopping point

Accepted event policy v1.0 and preprocessing parameter policy v1.0 are private
scientific decision records. Superseded drafts are retained privately for
provenance. Public, reproducible implementation boundaries remain tracked in:

- `eeg_behavior_identity_contract.csv` and `event_source_contract.py`;
- `montage_coordinate_contract_v1.yaml` and `montage_contract.py`;
- `preprocessing_parameter_audit_contract_v1.yaml` and
  `preprocessing_parameter_audit.py`;
- the public contract documents under [`../../docs/`](../../docs/).

The accepted continuous-preprocessing, eligibility-ledger, and epoch stages
are complete. Continuous success remains independent of event-source
availability, the accepted 8,905-row candidate surface, and later analytic
inclusion. The current workflow now contains constructed voltage epochs but no
time-frequency features or new EEG results. Transient artifact policy, TFR
construction and normalization, ROI/CSD choices, participant inclusion, and
statistical analysis remain separate future stages.

## Validation

The public test suite does not require private raw EEG:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 -m pytest
python3 tools/check_repo_safety.py
git status --short
```

Focused contracts are covered by `tests/test_event_source_contract.py`,
`tests/test_event_epoch_eligibility.py`,
`tests/test_epoch_construction.py`,
`tests/test_channel_qc.py`, `tests/test_montage_contract.py`,
`tests/test_preprocessing_parameter_audit.py`, and
`tests/test_raw_qc_guardrails.py`. Production continuous contracts, ICA routing,
atomic derivative reopening, and resumability are covered by
`tests/test_continuous_preprocessing_contract.py` and
`tests/test_continuous_preprocessing_storage.py`.
