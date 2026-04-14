######################################################
## Lead Code Developer and Maintainer: Changryong Baek
## HAR empirical application
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
LIB_FILE <- file.path(THIS_DIR, "scbm_empirical_library.R")
if (!file.exists(LIB_FILE)) {
  stop("Missing library file: ", LIB_FILE, call. = FALSE)
}
source(LIB_FILE)
suppressPackageStartupMessages(library(imputeTS))

# ------------------------------------------------------------
# User-editable configuration
# ------------------------------------------------------------
DATA_FILE <- file.path(THIS_DIR, "rk_mat2000-2022.Rdata")
if (!file.exists(DATA_FILE)) DATA_FILE <- file.path(THIS_DIR, "data", "rk_mat2000-2022.Rdata")
DATA_OBJECT <- NULL

START_DATE <- as.Date("2010-01-01")
END_DATE <- as.Date("2019-12-31")
DROP_COLUMNS <- c(5L, 30L)

# model specification
bw <- 5L
bm <- 22L
horizon_names <- c("daily", "weekly", "monthly")
estimator <- "lasso"
CENTER_VHAR <- TRUE

# lasso settings
lambda_mode <- "cv_c"
c_lambda <- 0.25
c_grid <- scbm_default_c_grid()
lasso_fold <- 10L
lasso_nlambda <- 100L
lasso_max_iter <- 1000L
lasso_tol <- 1e-6

diagTF <- TRUE
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
# Load and preprocess data
# ------------------------------------------------------------
obj <- load_rdata_object(DATA_FILE, object_name = DATA_OBJECT)
raw_df <- first_numeric_matrix(obj)
raw_df <- as.data.frame(raw_df, check.names = FALSE)

if (is.null(rownames(raw_df))) {
  stop("The realized-volatility input must have row names that can be parsed as dates.", call. = FALSE)
}
dates_all <- as.Date(rownames(raw_df))
if (any(is.na(dates_all))) {
  stop("Could not parse row names of the realized-volatility input as dates.", call. = FALSE)
}

idx_time <- dates_all >= START_DATE & dates_all <= END_DATE
rv_df <- raw_df[idx_time, , drop = FALSE]
dates <- dates_all[idx_time]

imputedata <- matrix(NA_real_, nrow = nrow(rv_df), ncol = ncol(rv_df))
for (j in seq_len(ncol(rv_df))) {
  x <- as.numeric(rv_df[[j]])
  x[x <= 0] <- NA_real_
  imputedata[, j] <- imputeTS::na_interpolation(log(x))
}

logrv <- imputedata[, -DROP_COLUMNS, drop = FALSE]
series_names <- colnames(rv_df[, -DROP_COLUMNS, drop = FALSE])
series_names <- sub("^\\.", "", series_names)
colnames(logrv) <- series_names
rownames(logrv) <- as.character(dates)

logrv_raw <- logrv
logrv_col_means <- colMeans(logrv_raw, na.rm = TRUE)
if (isTRUE(CENTER_VHAR)) {
  logrv <- sweep(logrv_raw, 2L, logrv_col_means, FUN = "-")
} else {
  logrv <- logrv_raw
}

SERIES_IN_ROWS <- FALSE

# ------------------------------------------------------------
# Step 1: first-stage fit
# ------------------------------------------------------------
fit_vhar <- scbm_vhar_empirical_fit(
  data = logrv,
  bw = bw,
  bm = bm,
  horizon_names = horizon_names,
  estimator = estimator,
  centerTF = FALSE,
  series_in_rows = SERIES_IN_ROWS,
  series_names = series_names,
  lambda_mode = lambda_mode,
  c_lambda = c_lambda,
  c_grid = c_grid,
  lasso_fold = lasso_fold,
  lasso_nlambda = lasso_nlambda,
  lasso_max_iter = lasso_max_iter,
  lasso_tol = lasso_tol,
  diagTF = diagTF,
  updateSigma = updateSigma,
  sigma_diag_only = sigma_diag_only
)

# ------------------------------------------------------------
# Step 2: clustering and console output only
# ------------------------------------------------------------

## 3-3-3 configuration
n_comm_333 <- c(3L, 3L, 3L)

cluster_vhar_333 <- scbm_vhar_empirical_cluster(
  fit_vhar,
  bw = bw,
  bm = bm,
  horizon_names = horizon_names,
  series_names = series_names,
  n_comm = n_comm_333,
  alpha_mode = alpha_mode,
  alpha = alpha,
  alpha_fold = alpha_fold,
  alpha_grid = alpha_grid,
  alpha_criterion = alpha_criterion,
  nstart = nstart,
  seed = seed
)

cat("\n========================================\n")
cat("VHAR aggregate flow table: 3-3-3 configuration\n")
cat("========================================\n")
print(cluster_vhar_333$flow_table)

## 3-4-3 configuration
n_comm_343 <- c(3L, 4L, 3L)

cluster_vhar_343 <- scbm_vhar_empirical_cluster(
  fit_vhar,
  bw = bw,
  bm = bm,
  horizon_names = horizon_names,
  series_names = series_names,
  n_comm = n_comm_343,
  alpha_mode = alpha_mode,
  alpha = alpha,
  alpha_fold = alpha_fold,
  alpha_grid = alpha_grid,
  alpha_criterion = alpha_criterion,
  nstart = nstart,
  seed = seed
)

cat("\n========================================\n")
cat("VHAR aggregate flow table: 3-4-3 configuration\n")
cat("========================================\n")
print(cluster_vhar_343$flow_table)
