# Empirical folder: console-only runners

This note describes the console-only versions of the empirical runners.

## Files

* `run_scbm_pvar_empirical_application.R`

  * U.S. payroll PVAR application
  * no workspace save
  * no plot creation
  * no file export
  * prints the flow tables for both effective paths `2 -> 3 -> 3 -> 2` and `2 -> 2 -> 3 -> 2`
* `run_scbm_vhar_empirical_application.R`

  * global realized-volatility VHAR application
  * no workspace save
  * no plot creation
  * no file export
  * prints the aggregate flow tables for both `3-3-3` and `3-4-3` configurations
* `scbm_empirical_library.R`

  * shared empirical-analysis library used by both scripts

## Default settings

### PVAR console runner

* estimator: lasso
* seasonal structure: `s = 4`, `p = 1`
* lambda selection: `cv_c` with `c_lambda = 0.20`
* centering: seasonal mean centering within each quarter
* PisCES alpha: selected by CV
* printed outputs: flow tables for `2332` and `2232`

### VHAR console runner

* estimator: lasso
* horizon lengths: `bw = 5`, `bm = 22`
* lambda selection: `cv_c` with `c_lambda = 0.25`
* centering: global series-wise demeaning after log transform and interpolation
* PisCES alpha: selected by CV
* printed outputs: aggregate flow tables for `333` and `343`

Edit the user configuration block at the top of each script if you want to change these defaults.

