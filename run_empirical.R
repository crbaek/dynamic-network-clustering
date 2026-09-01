# Baseline PVAR payroll and generalized VHAR volatility applications.
# Alternative empirical specifications are intentionally excluded.

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

APPLICATION <- tolower(Sys.getenv("SCBM_APPLICATION", "both"))
if (!APPLICATION %in% c("pvar", "vhar", "both")) {
  stop("SCBM_APPLICATION must be pvar, vhar, or both.")
}
SEED <- as.integer(Sys.getenv("SCBM_SEED", "12345"))

if (APPLICATION %in% c("pvar", "both")) scbm_run_pvar_empirical(ROOT, SEED)
if (APPLICATION %in% c("vhar", "both")) scbm_run_vhar_empirical(ROOT, SEED)

cat("\nResults:", file.path(ROOT, "output", "empirical"), "\n")
