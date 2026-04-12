######################################################
## Lead Code Developer and Maintainer: Changryong Baek
## Set up your working directory properly.
######################################################

rm(list = ls())

# ============================================================
# run_scbm_vhar_simulation.R
# GitHub-ready batch runner for the generalized ScBM-VHAR
# simulation. Settings are aligned with the final manuscript:
#   q in {18, 24, 36}
#   T in {500, 1000, 2000, 3000}
#   setting in {sg1, sg2, sg3}
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
source(file.path(THIS_DIR, "scbm_vhar_simulation_library.R"))

# ------------------------------
# User-editable configuration
# ------------------------------
q_grid <- c(18L, 24L, 36L)
gridT <- c(500L, 1000L, 2000L, 3000L)
settings <- c("sg1", "sg2", "sg3")
paths <- c("path1", "path2", "path3")
link_types <- c("type1", "type2")
methods <- c("ols", "lasso")

nrep <- 100L
ncore <- 50L
seed_base <- 20260328L

# generalized VHAR design
bM <- 3L
bL <- 10L
target_var <- 0.5
burn <- 300L

# lasso settings
lasso_lambda_mode <- "cv_c"   # one of: "theory", "cv", "cv_c"
lasso_c_lambda <- 0.25
lasso_c_grid <- scbm_vhar_default_c_grid()
lasso_fold <- 10L
lasso_nlambda <- 100L
lasso_max_iter <- 1000L
lasso_tol <- 1e-6

# PisCES settings
alpha_fold <- 5L
alpha_grid <- scbm_vhar_default_alpha_grid()
alpha_criterion <- "holdout"

save_outputs <- TRUE
save_stub <- "scbm_vhar_simulation"

# ------------------------------
# Helper
# ------------------------------
run_one_case <- function(setting, path, q, link_type, estimator) {
  scbm_vhar_run_simulation(
    case_name = setting,
    path = path,
    method = estimator,
    q = q,
    T = gridT,
    nrep = nrep,
    ncore = ncore,
    link_type = link_type,
    bd = bM,
    bm = bL,
    burn = burn,
    target_var = target_var,
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
    libfile = file.path(THIS_DIR, "scbm_vhar_simulation_library.R")
  )
}

# ------------------------------
# Run all configurations
# ------------------------------
manifest <- expand.grid(
  estimator = methods,
  setting = settings,
  path = paths,
  q = q_grid,
  type = link_types,
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$estimator, manifest$q, manifest$setting, manifest$path, manifest$type), ]
manifest$object <- paste0("out", seq_len(nrow(manifest)))
rownames(manifest) <- NULL

out_all <- vector("list", nrow(manifest))
names(out_all) <- manifest$object
error_log <- data.frame(
  ii = integer(0),
  estimator = character(0),
  q = integer(0),
  setting = character(0),
  path = character(0),
  type = character(0),
  error_message = character(0),
  stringsAsFactors = FALSE
)

for (ii in seq_len(nrow(manifest))) {
  out_obj <- tryCatch(
    run_one_case(
      setting = manifest$setting[ii],
      path = manifest$path[ii],
      q = manifest$q[ii],
      link_type = manifest$type[ii],
      estimator = manifest$estimator[ii]
    ),
    error = function(e) e
  )

  if (inherits(out_obj, "error") || inherits(out_obj, "simpleError")) {
    cat(sprintf(
      "[%d/%d] failed: estimator=%s, q=%d, setting=%s, path=%s, type=%s :: %s\n",
      ii, nrow(manifest),
      manifest$estimator[ii], manifest$q[ii], manifest$setting[ii],
      manifest$path[ii], manifest$type[ii], conditionMessage(out_obj)
    ))
    error_log <- rbind(
      error_log,
      data.frame(
        ii = ii,
        estimator = manifest$estimator[ii],
        q = manifest$q[ii],
        setting = manifest$setting[ii],
        path = manifest$path[ii],
        type = manifest$type[ii],
        error_message = conditionMessage(out_obj),
        stringsAsFactors = FALSE
      )
    )
    next
  }

  assign(manifest$object[ii], out_obj, envir = .GlobalEnv)
  out_all[[ii]] <- out_obj

  cat(sprintf(
    "[%d/%d] done: estimator=%s, q=%d, setting=%s, path=%s, type=%s\n",
    ii, nrow(manifest),
    manifest$estimator[ii], manifest$q[ii], manifest$setting[ii],
    manifest$path[ii], manifest$type[ii]
  ))
}

assign("out_manifest", manifest, envir = .GlobalEnv)
assign("out_all", out_all, envir = .GlobalEnv)
assign("error_log", error_log, envir = .GlobalEnv)

valid_out <- Filter(Negate(is.null), out_all)
if (length(valid_out)) {
  raw_all <- do.call(rbind, lapply(valid_out, function(z) z$raw))
  summary_all <- do.call(scbm_vhar_simulation_summary, valid_out)
} else {
  raw_all <- data.frame()
  summary_all <- data.frame()
}

summary_none <- if (nrow(summary_all)) summary_all[summary_all$smoothing == "none", , drop = FALSE] else data.frame()
summary_pisces <- if (nrow(summary_all)) summary_all[summary_all$smoothing == "pisces", , drop = FALSE] else data.frame()
summary_ols <- if (nrow(summary_all)) summary_all[summary_all$estimator == "ols", , drop = FALSE] else data.frame()
summary_lasso <- if (nrow(summary_all)) summary_all[summary_all$estimator == "lasso", , drop = FALSE] else data.frame()

assign("raw_all", raw_all, envir = .GlobalEnv)
assign("summary_all", summary_all, envir = .GlobalEnv)
assign("summary_none", summary_none, envir = .GlobalEnv)
assign("summary_pisces", summary_pisces, envir = .GlobalEnv)
assign("summary_ols", summary_ols, envir = .GlobalEnv)
assign("summary_lasso", summary_lasso, envir = .GlobalEnv)

# ------------------------------
# Save
# ------------------------------
if (isTRUE(save_outputs)) {
  save_list <- c(
    manifest$object[manifest$object %in% ls(envir = .GlobalEnv)],
    "out_manifest", "out_all", "error_log", "raw_all",
    "summary_all", "summary_none", "summary_pisces",
    "summary_ols", "summary_lasso"
  )
  save(list = save_list, file = paste0(save_stub, ".RData"))
  utils::write.csv(out_manifest, paste0(save_stub, "-manifest.csv"), row.names = FALSE)
  utils::write.csv(error_log, paste0(save_stub, "-error-log.csv"), row.names = FALSE)
  utils::write.csv(raw_all, paste0(save_stub, "-raw-all.csv"), row.names = FALSE)
  utils::write.csv(summary_all, paste0(save_stub, "-summary-all.csv"), row.names = FALSE)
  utils::write.csv(summary_none, paste0(save_stub, "-summary-none.csv"), row.names = FALSE)
  utils::write.csv(summary_pisces, paste0(save_stub, "-summary-pisces.csv"), row.names = FALSE)
  utils::write.csv(summary_ols, paste0(save_stub, "-summary-ols.csv"), row.names = FALSE)
  utils::write.csv(summary_lasso, paste0(save_stub, "-summary-lasso.csv"), row.names = FALSE)
}

cat("Done. Objects out1, out2, ..., out", nrow(manifest), " and summaries are in the workspace.\n", sep = "")
