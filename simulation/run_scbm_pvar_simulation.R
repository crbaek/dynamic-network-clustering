######################################################
## Lead Code Developer and Maintainer: Changryong Baek
## Set up your working directory properly.
######################################################

rm(list = ls())

# ============================================================
# run_scbm_pvar_simulation.R
# GitHub-ready batch runner for the ScBM-PVAR simulation.
# Settings are aligned with the final manuscript:
#   q in {18, 36, 60}
#   T in {200, 500, 1000, 2000}
#   path in {path1, path2, path3}
#   type in {type1, type2}
#   estimator in {ols, lasso}
# ============================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)))
  }
  if (!is.null(sys.frames()[[1L]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1L]]$ofile, winslash = "/", mustWork = FALSE)))
  }
  getwd()
}
THIS_DIR <- get_script_dir()
source(file.path(THIS_DIR, "scbm_pvar_simulation_library.R"))

# ------------------------------
# User-editable configuration
# ------------------------------
q_grid <- c(18L, 36L, 60L)
gridT <- c(200L, 500L, 1000L, 2000L)
paths <- c("path1", "path2", "path3")
link_types <- c("type1", "type2")
methods <- c("ols", "lasso")

nrep <- 200L
ncore <- 50L
seed_base <- 20260328L

# design settings
sigma_diag <- 0.50
target_cycle_sv <- 0.90
burn <- 500L

# lasso settings
lasso_lambda_mode <- "cv_c"     # one of: "theory", "cv", "cv_c"
lasso_c_lambda <- 0.20
lasso_c_grid <- scbm_pvar_default_c_grid()
lasso_fold <- 10L
lasso_nlambda <- 50L
lasso_max_iter <- 1000L
lasso_tol <- 1e-6

# PisCES settings
alpha_fold <- 5L
alpha_grid <- scbm_pvar_default_alpha_grid()
alpha_criterion <- "holdout"

save_outputs <- TRUE
save_stub <- "scbm_pvar_simulation"

# ------------------------------
# Helper
# ------------------------------
run_one_case <- function(path, q, link_type, estimator) {
  scbm_pvar_run_simulation(
    path = path,
    method = estimator,
    q = q,
    T = gridT,
    nrep = nrep,
    ncore = ncore,
    link_type = link_type,
    target_cycle_sv = target_cycle_sv,
    sigma_diag = sigma_diag,
    burn = burn,
    lasso_lambda_mode = lasso_lambda_mode,
    lasso_c_lambda = lasso_c_lambda,
    lasso_c_grid = lasso_c_grid,
    lasso_fold = lasso_fold,
    lasso_nlambda = lasso_nlambda,
    lasso_max_iter = lasso_max_iter,
    lasso_tol = lasso_tol,
    alpha_fold = alpha_fold,
    alpha_grid = alpha_grid,
    alpha_criterion = alpha_criterion,
    base_seed = seed_base,
    libfile = file.path(THIS_DIR, "scbm_pvar_simulation_library.R")
  )
}

# ------------------------------
# Run all configurations
# ------------------------------
manifest <- expand.grid(
  estimator = methods,
  q = q_grid,
  path = paths,
  type = link_types,
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$estimator, manifest$q, manifest$path, manifest$type), ]
manifest$object <- paste0("out", seq_len(nrow(manifest)))
rownames(manifest) <- NULL

out_all <- vector("list", nrow(manifest))
names(out_all) <- manifest$object

for (ii in seq_len(nrow(manifest))) {
  out_obj <- run_one_case(
    path = manifest$path[ii],
    q = manifest$q[ii],
    link_type = manifest$type[ii],
    estimator = manifest$estimator[ii]
  )

  assign(manifest$object[ii], out_obj, envir = .GlobalEnv)
  out_all[[ii]] <- out_obj

  cat(sprintf(
    "[%d/%d] done: estimator=%s, q=%d, path=%s, type=%s\n",
    ii, nrow(manifest),
    manifest$estimator[ii], manifest$q[ii], manifest$path[ii], manifest$type[ii]
  ))
}

assign("out_manifest", manifest, envir = .GlobalEnv)
assign("out_all", out_all, envir = .GlobalEnv)

# ------------------------------
# Summaries
# ------------------------------
summary_all <- do.call(scbm_pvar_simulation_summary, out_all)
summary_ols <- do.call(scbm_pvar_simulation_summary, out_all[manifest$estimator == "ols"])
summary_lasso <- do.call(scbm_pvar_simulation_summary, out_all[manifest$estimator == "lasso"])
raw_all <- do.call(rbind, lapply(out_all, function(z) z$raw))
summary_none <- subset(summary_all, smoothing == "none")
summary_pisces <- subset(summary_all, smoothing == "pisces")

assign("summary_all", summary_all, envir = .GlobalEnv)
assign("summary_ols", summary_ols, envir = .GlobalEnv)
assign("summary_lasso", summary_lasso, envir = .GlobalEnv)
assign("summary_none", summary_none, envir = .GlobalEnv)
assign("summary_pisces", summary_pisces, envir = .GlobalEnv)
assign("raw_all", raw_all, envir = .GlobalEnv)

# ------------------------------
# Save
# ------------------------------
if (isTRUE(save_outputs)) {
  save_list <- c(
    manifest$object,
    "out_manifest", "out_all", "raw_all",
    "summary_all", "summary_ols", "summary_lasso",
    "summary_none", "summary_pisces"
  )
  save(list = save_list, file = paste0(save_stub, ".RData"))
  utils::write.csv(out_manifest, paste0(save_stub, "-manifest.csv"), row.names = FALSE)
  utils::write.csv(raw_all, paste0(save_stub, "-raw-all.csv"), row.names = FALSE)
  utils::write.csv(summary_all, paste0(save_stub, "-summary-all.csv"), row.names = FALSE)
  utils::write.csv(summary_ols, paste0(save_stub, "-summary-ols.csv"), row.names = FALSE)
  utils::write.csv(summary_lasso, paste0(save_stub, "-summary-lasso.csv"), row.names = FALSE)
  utils::write.csv(summary_none, paste0(save_stub, "-summary-none.csv"), row.names = FALSE)
  utils::write.csv(summary_pisces, paste0(save_stub, "-summary-pisces.csv"), row.names = FALSE)
}

cat("Done. Objects out1, out2, ..., out", nrow(manifest), " and summaries are in the workspace.\n", sep = "")
