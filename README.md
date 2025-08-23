# Detecting Error in Motor Imagery (DEMI)

Analysis code for the DEMI project. This repo currently contains the **behavior-only** analysis (scripts **00–03**) as used in the published study below. The **EEG + hierarchical GAM** work (scripts 04+) will be added in a later pass.

## Layout
- `R/00–03` — behavior analysis (unchanged)
- `legacy/dissertation/` — archival dissertation-era repo (read-only import)
- `data/` — local only; not tracked
- `outputs/` — small tables/figures (if any)

## References
- Ingram, T. G. J., Hurst, A. J., Solomon, J. P., Stratas, A., & Boe, S. G. (2022).
  *Imagined Movement Accuracy Is Strongly Associated With Drivers of Overt Movement Error and Weakly Associated With Imagery Vividness.*
  **JEP: Human Perception & Performance, 48(12), 1362–1372.** https://doi.org/10.1037/xhp0001064
- Task code: TraceLab — https://github.com/LBRF/TraceLab 
