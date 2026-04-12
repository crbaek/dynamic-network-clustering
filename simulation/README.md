# ScBM simulation folder

This folder contains the GitHub-facing entry points for the simulation part of the paper.

## Files

- `scbm_pvar_simulation_library.R`
  - standalone PVAR simulation library
  - current working sparse fixed-design generator
  - OLS / lasso first-stage estimation
  - optional PisCES smoothing
  - public entry point: `scbm_pvar_run_simulation()`

- `run_scbm_pvar_simulation.R`
  - batch runner for the PVAR designs used in the manuscript
  - writes `.RData` and `.csv` summaries in the working directory

- `scbm_vhar_simulation_library.R`
  - standalone generalized VHAR simulation library
  - current working sparse fixed-design generator
  - OLS / lasso first-stage estimation
  - optional PisCES smoothing
  - public entry point: `scbm_vhar_run_simulation()`

- `run_scbm_vhar_simulation.R`
  - batch runner for the generalized VHAR designs used in the manuscript
  - writes `.RData` and `.csv` summaries in the working directory

## Design principles

1. The run scripts use relative paths only; there is no hard-coded `setwd()`.
2. The public file names are cleaned for GitHub, but the underlying working logic from the current project files is preserved.
3. Output objects are saved as `.RData`, matching the preferred workflow used elsewhere in the project.
4. The folder is meant to sit next to an empirical folder that will use the same `.RData` convention.

## Typical usage

From R, set the working directory to this folder and run one of the batch scripts:

```r
source("run_scbm_pvar_simulation.R")
source("run_scbm_vhar_simulation.R")
```

For custom runs, source the relevant library and call the public wrapper directly.
