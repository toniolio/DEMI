# Error Signatures in Motor Imagery — Behavior-Only Pass (00–04)

This repository hosts the **analysis code** for a motor imagery vs. overt movement task. The first pass here focuses on the **published behavior-only paper** and preserves the original scripts **unchanged**:

- `R/00_import_behav_dat.R` (renamed from your legacy `0_import.R`)
- `R/01_preprocessing.R`
- `R/02_imagery_clean.R`
- `R/03_behav_analysis.R`

> **Note:** For full EEG+GAM work, see the project snapshot in the issue tracker or the repo wiki — that adds scripts `04+`, per-participant Parquet, and the GAM + max-|T| inference pipeline. That will be layered on after this behavior-only pass.

## What’s in this pass

- A clean directory layout (`R/`, `data/`, `outputs/`, etc.).
- No large data in the repo. Place your local data under `data/raw/` and derived artifacts in `data/derived/`.
- Small shareable outputs (CSV tables/PNGs) will live in `outputs/`.
- A lightweight installer to get required R packages for the behavior scripts.

## Behavior scripts (kept as-is)

- **00_import_behav_dat.R** — imports the task/behavior data and writes compact RDS/CSV to `data/derived/behav/` (or legacy `_Scripts/_rds/`).  
- **01_preprocessing.R** — preprocessing for the behavior-only paper (unchanged).  
- **02_imagery_clean.R** — imagery cleanup (TraceLab issues), uses a small Stan model to flag outliers and back-fill values.  
- **03_behav_analysis.R** — published analyses (unchanged).  

> These scripts may `source()` helper files from a legacy `_Scripts/` tree (e.g., `_Scripts/_functions/*.R`, `_Scripts/_stan/outlier_mixture.stan`). Either include that tree here (preferred) or update the `source()` paths after we migrate helpers into `R/helpers/`.

## Installing dependencies (R)

From an R session at the repo root:

```r
# Option A: use renv (recommended)
install.packages("renv")
renv::init()         # creates renv.lock from your current R library
source("scripts/install_deps.R")  # installs CRAN/GitHub packages needed by 00–03
renv::snapshot()     # pins exact versions

# Option B: install without renv
source("scripts/install_deps.R")
```

### System toolchains / Stan

- `cmdstanr` (for `02_imagery_clean.R`) requires a working C++ toolchain and a CmdStan install:  
  ```r
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  cmdstanr::install_cmdstan()  # one-time
  ```
- `brms` (used in `03_behav_analysis.R`) defaults to **rstan**; make sure `rstan` is installed and can compile models. You can switch to CmdStan by passing `backend = "cmdstanr"` to `brms::brm()` if you prefer a single backend.

## Data policy

- Raw/large data are **not** distributed.
- Place your local data under `data/raw/`. Derived artifacts go to `data/derived/`.
- Small tables and figures (OK to track) live under `outputs/`.

## Citation

If you use this code, please cite the behavior paper and this repository. The behavior paper is already published in *JEP: Human Perception & Performance* (2019–2022 era). Add the full citation and DOI here once you populate this section.

## License

Add your preferred license (e.g., MIT) here once confirmed with co-authors.
