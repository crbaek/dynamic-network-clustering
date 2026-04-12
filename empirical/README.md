# Empirical folder

This folder contains the two empirical applications from the manuscript in a GitHub-ready form.

## Files

- `scbm_empirical_library.R`
  - shared empirical-analysis library
  - supports PVAR and VHAR fitting, clustering, and PisCES smoothing

- `run_scbm_pvar_empirical_application.R`
  - U.S. payroll PVAR application
  - expects the bundled `data/pvar_payroll_1990_2020.RData`
  - uses the processed quarterly log-difference series from the `.RData` object
  - keeps fit and clustering objects in memory
  - saves one workspace `.RData` before plotting
  - figure saving is optional through `SAVE_PLOTS`

- `run_scbm_vhar_empirical_application.R`
  - global realized-volatility VHAR application
  - expects the bundled `data/rk_mat2000-2022.Rdata`
  - applies the same date filtering, log transform, interpolation, and index removal used in the working code
  - keeps fit and clustering objects in memory
  - saves one workspace `.RData` before plotting
  - figure saving is optional through `SAVE_PLOTS`

- `data/`
  - empirical inputs

- `output/`
  - default location for the saved workspace `.RData`
  - figure files are written only when `SAVE_PLOTS <- TRUE`

## Default empirical settings

### PVAR application

- estimator: lasso
- seasonal structure: `s = 4`, `p = 1`
- lambda selection: `cv_c` with `c_lambda = 0.20`
- clustering path counts: `c(2, 2, 2, 3, 3, 3, 3, 2)`
- PisCES alpha: selected by CV
- no fit/cluster csv export

### VHAR application

- estimator: lasso
- horizon lengths: `bw = 5`, `bm = 22`
- lambda selection: `cv_c` with `c_lambda = 0.25`
- aggregated community counts: `c(3, 3, 3)` in the order short / middle / long
- PisCES alpha: selected by CV
- no fit/cluster csv export

Edit the user configuration block at the top of each script if you want to change these defaults.
