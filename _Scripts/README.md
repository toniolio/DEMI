# Published behaviour and historical EEG analysis

This directory contains two related provenance surfaces:

- scripts 00–03 and their helpers support the behavioural workflow associated
  with the published behavioural paper;
- scripts 04–13 contain the earlier EEG import, preparation, plotting, and
  hierarchical-GAM workflow used during prior analysis development.

The behavioural analysis is published and frozen. The EEG/GAM scripts are
historical reference material, not the active confirmatory EEG pipeline. They
remain useful for reconstructing task/tracing linkage, event handling, earlier
preprocessing assumptions, and result provenance, and should not be deleted or
silently modernized.

New EEG work starts in [`../analysis/eeg_mne/`](../analysis/eeg_mne/README.md).
The Python linkage reconstruction in
[`../analysis/behavior/`](../analysis/behavior/README.md) ports only the
frozen timing/linkage logic required by that reanalysis; it does not reopen the
published behavioural analysis.

Generated RDS files, model fits, and plot dumps are local-only. The tracked
`../_Data/eeg/BESA-81.csv` file is the historical Oostenveld/BESA unit-sphere
coordinate table expected by scripts 05 and 06; it is not an acquisition
channel map or participant digitization.
