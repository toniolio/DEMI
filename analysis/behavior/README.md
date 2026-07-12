# Behavioural and tracing linkage for the EEG reanalysis

This directory reconstructs the frozen behavioural/task and TraceLab timing
surface needed to link raw EEG events to trials. It does not revise the
published behavioural analysis, refit its models, or generate new behavioural
results. The historical R workflow under [`../../_Scripts/`](../../_Scripts/README.md)
remains the authority source that these scripts port deliberately.

The directory name is retained because the scripts and their local output
contracts are already used throughout the event-evidence workflow. Moving it
would create broad path churn without clarifying the scientific boundary more
than this guide does.

## Run order

Run from the repository root after placing local task files under
`_Data/task/` and TraceLab figure archives under `_Data/figure/`:

```sh
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/00_make_behavioral_file_manifest.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/01_inspect_task_and_figure_linkage.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/02_reconstruct_tracing_tables.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/03_apply_tracing_filters.py
PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/04_build_event_offsets.py
```

| Script | Role | Local output |
| --- | --- | --- |
| 00 | Inventory task and TraceLab source files. | `_Data/behavior/manifest/` |
| 01 | Inspect task-table and TraceLab filename structure. | `_Data/behavior/manifest/` |
| 02 | Reconstruct TraceLab points, segments, frames, and tracings. | `_Data/behavior/tracing_tables/` |
| 03 | Port the historical tracing filters and timing summaries. | `_Data/behavior/tracing_filters/` |
| 04 | Build the old-compatible task-to-tracing event-offset surface used by EEG alignment. | `_Data/behavior/event_offsets/` |

All outputs are generated, local-only evidence. Some exact historical parity
checks use optional private/local RDS artifacts; the scripts document those
inputs explicitly. Participant-level decisions and the accepted event policy
remain outside this directory.

Continue with the active EEG guide in
[`../eeg_mne/README.md`](../eeg_mne/README.md).
