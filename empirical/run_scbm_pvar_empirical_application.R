######################################################
## Lead Code Developer and Maintainer: Changryong Baek
## PVAR empirical application
## - prints aggregate flow tables only
## - setup your workspace/path properly
######################################################

rm(list = ls())

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

load_rdata_object <- function(data_file, object_name = NULL) {
  if (!file.exists(data_file)) stop("Missing data file: ", data_file, call. = FALSE)
  ee <- new.env(parent = emptyenv())
  load(data_file, envir = ee)
  objs <- ls(envir = ee)
  if (!length(objs)) stop("No objects found in: ", data_file, call. = FALSE)
  if (is.null(object_name)) object_name <- objs[1L]
  if (!object_name %in% objs) {
    stop("Object '", object_name, "' not found in: ", data_file, call. = FALSE)
  }
  get(object_name, envir = ee)
}

first_numeric_matrix <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) return(as.matrix(x))
  if (is.list(x)) {
    for (nm in c("Y_diff", "Y", "Yt", "data", "logrv", "rv", "matrix")) {
      if (!is.null(x[[nm]]) && (is.matrix(x[[nm]]) || is.data.frame(x[[nm]]))) {
        return(as.matrix(x[[nm]]))
      }
    }
    for (elt in x) {
      if (is.matrix(elt) || is.data.frame(elt)) return(as.matrix(elt))
    }
  }
  stop("Could not find a matrix/data.frame inside the loaded object.", call. = FALSE)
}

THIS_DIR <- get_script_dir()
source(file.path(THIS_DIR, "scbm_empirical_library.R"))

# ------------------------------------------------------------
# User-editable configuration
# ------------------------------------------------------------
DATA_FILE <- file.path(THIS_DIR, "data", "pvar_payroll_1990_2020.RData")
DATA_OBJECT <- NULL

# model specification
s <- 4L
p <- 1L
estimator <- "lasso"

# lasso settings
lambda_mode <- "cv_c"
c_lambda <- 0.20
c_grid <- scbm_default_c_grid()
lasso_fold <- 10L
lasso_max_iter <- 1000L
lasso_tol <- 1e-6

diagTF <- TRUE
centerTF <- TRUE
updateSigma <- TRUE
sigma_diag_only <- TRUE

# PisCES settings
alpha_mode <- "cv"
alpha <- NULL
alpha_fold <- 5L
alpha_grid <- scbm_default_alpha_grid()
alpha_criterion <- "paper"

# clustering settings
nstart <- 50L
seed <- 12345L

# ------------------------------------------------------------
# Load and parse data
# ------------------------------------------------------------
obj <- load_rdata_object(DATA_FILE, object_name = DATA_OBJECT)

if (is.list(obj) && !is.null(obj$Y_diff)) {
  payroll_mat <- as.matrix(obj$Y_diff)
  series_names <- if (!is.null(obj$series_codes)) as.character(obj$series_codes) else colnames(payroll_mat)
  time_labels <- if (!is.null(obj$quarterly_labels)) as.character(obj$quarterly_labels) else rownames(payroll_mat)
} else {
  payroll_mat <- first_numeric_matrix(obj)
  series_names <- colnames(payroll_mat)
  time_labels <- rownames(payroll_mat)
}

if (is.null(series_names)) series_names <- paste0("V", seq_len(ncol(payroll_mat)))
colnames(payroll_mat) <- series_names
if (!is.null(time_labels) && length(time_labels) == nrow(payroll_mat)) rownames(payroll_mat) <- time_labels

SERIES_IN_ROWS <- FALSE

# ------------------------------------------------------------
# Step 1: first-stage fit
# ------------------------------------------------------------
fit_pvar <- scbm_pvar_empirical_fit(
  data = payroll_mat,
  s = s,
  p = p,
  series_in_rows = SERIES_IN_ROWS,
  series_names = series_names,
  estimator = estimator,
  lambda_mode = lambda_mode,
  c_lambda = c_lambda,
  c_grid = c_grid,
  lasso_fold = lasso_fold,
  lasso_max_iter = lasso_max_iter,
  lasso_tol = lasso_tol,
  diagTF = diagTF,
  centerTF = centerTF,
  updateSigma = updateSigma,
  sigma_diag_only = sigma_diag_only
)

# ------------------------------------------------------------
# Step 2: clustering and console output only
# ------------------------------------------------------------

## Main configuration: 2 -> 3 -> 3 -> 2
n_comm_2332 <- c(2L, 3L, 3L, 3L, 3L, 2L, 2L, 2L)

cluster_pvar_2332 <- scbm_pvar_empirical_cluster(
  fit_pvar,
  s = s,
  p = p,
  series_names = series_names,
  n_comm = n_comm_2332,
  alpha_mode = alpha_mode,
  alpha = alpha,
  alpha_fold = alpha_fold,
  alpha_grid = alpha_grid,
  alpha_criterion = alpha_criterion,
  nstart = nstart,
  seed = seed
)

cat("\n========================================\n")
cat("PVAR flow table: effective path 2 -> 3 -> 3 -> 2\n")
cat("========================================\n")
print(cluster_pvar_2332$flow_table)

## Alternative configuration: 2 -> 2 -> 3 -> 2
n_comm_2232 <- c(2L, 2L, 2L, 3L, 3L, 3L, 3L, 2L)

cluster_pvar_2232 <- scbm_pvar_empirical_cluster(
  fit_pvar,
  s = s,
  p = p,
  series_names = series_names,
  n_comm = n_comm_2232,
  alpha_mode = alpha_mode,
  alpha = alpha,
  alpha_fold = alpha_fold,
  alpha_grid = alpha_grid,
  alpha_criterion = alpha_criterion,
  nstart = nstart,
  seed = seed
)

cat("\n========================================\n")
cat("PVAR flow table: effective path 2 -> 2 -> 3 -> 2\n")
cat("========================================\n")
print(cluster_pvar_2232$flow_table)
