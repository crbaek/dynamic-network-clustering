# Latent community paths in VAR-type models

Minimal reproduction code for the main ScBM--PVAR and generalized ScBM--VHAR results.
The repository contains no legacy patch chain and no `R/internal` directory.

## Files

```text
R/scbm_core.R     common spectral clustering, PisCES, CV, metrics and plotting
R/scbm_pvar.R     PVAR simulation and payroll application
R/scbm_vhar.R     VHAR simulation and volatility application
run_simulation.R  main-text lasso simulations
run_empirical.R   two baseline empirical applications
data/              empirical input files (retained from the existing repository)
packages/          bundled sparseVAR source package (retained from the existing repository)
output/            generated results
```

## Methodological convention

At every stage, the directed matrix analyzed by SVD is `M = t(Phi)`.
The left singular space represents sending communities and the right singular
space represents receiving communities. Both sides retain
`r = min(Ky, Kz)` singular directions; `k`-means still uses `Ky` and `Kz`
clusters. VHAR paths are ordered Long -> Medium -> Short.

The PisCES smoothing parameter is selected by repeated uniform off-diagonal
projector reconstruction CV with masking fraction 0.20, 10 masks, signed
held-out error, and an alpha grid ending at `0.99 * alpha_max`.

## Simulation

Quick verification:

```r
source("run_simulation.R")
```

Full main-text grids (200 replications):

```r
Sys.setenv(SCBM_MODE = "paper", SCBM_NCORE = "50")
source("run_simulation.R")
```

The PVAR grid uses `q = 18, 36, 60`, three paths, two interaction types, and
`T = 200, 500, 1000, 2000`. The VHAR grid uses `q = 18, 24, 36`, the same
three paths and two types, and `T = 500, 1000, 2000, 3000` under the manuscript
signal setting. Only lasso/PisCES main-text results are produced.

## Empirical applications

Place the existing data files in `data/`, then run:

```r
source("run_empirical.R")
```

To run one application:

```r
Sys.setenv(SCBM_APPLICATION = "pvar")
source("run_empirical.R")
```

or replace `pvar` by `vhar`.

The PVAR baseline is `(2,2; 2,3; 3,3; 3,2)` with common ranks `(2,2,3,2)`.
The VHAR baseline is Long--Medium--Short `3--3--3` with common ranks `(3,3,3)`.
The PVAR and VHAR scree figures use 80% and 75% reference levels, respectively;
blue solid curves are cumulative proportions and red dotted curves are
individual proportions. The Sankey figures reproduce the manuscript layout:
pastel community boxes list the constituent series alphabetically and
source-colored ribbons connect the successive states.

Both applications use the manuscript seed `12345`. The empirical alpha-CV and
final clustering seeds follow the production seed-mixing rules, so the reported
community labels and the displayed figures are reproducible rather than merely
partition-equivalent.

## Required packages

`MASS`, `sparseVAR`, and `imputeTS` are required. `ggplot2` and `ggsankey` are
optional and are used only for Sankey figures.
