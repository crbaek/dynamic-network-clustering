# Main-text ScBM--PVAR and generalized ScBM--VHAR simulations.
# Robustness and supplementary OLS grids are intentionally excluded.

get_this_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  ff <- grep("^--file=", args, value = TRUE)
  if (length(ff)) return(dirname(normalizePath(sub("^--file=", "", ff[1L]), winslash = "/", mustWork = FALSE)))
  frames <- sys.frames()
  if (length(frames)) for (ii in rev(seq_along(frames))) {
    of <- frames[[ii]]$ofile
    if (!is.null(of) && nzchar(of)) return(dirname(normalizePath(of, winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

ROOT <- get_this_dir()
source(file.path(ROOT, "R", "scbm_core.R"))
source(file.path(ROOT, "R", "scbm_pvar.R"))
source(file.path(ROOT, "R", "scbm_vhar.R"))

MODE <- tolower(Sys.getenv("SCBM_MODE", "quick"))
if (!MODE %in% c("quick", "paper")) stop("SCBM_MODE must be quick or paper.")
RUN_PVAR <- tolower(Sys.getenv("SCBM_RUN_PVAR", "true")) %in% c("1", "true", "yes")
RUN_VHAR <- tolower(Sys.getenv("SCBM_RUN_VHAR", "true")) %in% c("1", "true", "yes")
NCORE <- as.integer(Sys.getenv("SCBM_NCORE", "2"))
NREP_TXT <- Sys.getenv("SCBM_NREP", "")
NREP <- if (nzchar(NREP_TXT)) as.integer(NREP_TXT) else NULL
SEED <- as.integer(Sys.getenv("SCBM_SEED", "12345"))

if (!RUN_PVAR && !RUN_VHAR) stop("At least one model must be selected.")
cat("Mode:", MODE, "\n")
cat("Convention: M=t(Phi), L=sender, R=receiver, r=min(Ky,Kz).\n")
cat("Alpha CV: uniform off-diagonal masks, rho=0.20, R=10.\n\n")

if (RUN_PVAR) scbm_run_pvar_simulation(MODE, NREP, NCORE, SEED + 1000L, ROOT)
if (RUN_VHAR) scbm_run_vhar_simulation(MODE, NREP, NCORE, SEED + 2000L, ROOT)

cat("\nResults:", file.path(ROOT, "output", "simulation"), "\n")
