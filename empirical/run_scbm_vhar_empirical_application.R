######################################################
## Lead Code Developer and Maintainer: Changryong Baek
## Set up your working directory properly.
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

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

run_plot_block <- function(save_plot, file, width, height, code) {
  if (isTRUE(save_plot)) {
    grDevices::pdf(file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(code)
  } else if (interactive()) {
    force(code)
  }
  invisible(NULL)
}

THIS_DIR <- get_script_dir()
source(file.path(THIS_DIR, "scbm_empirical_library.R"))
suppressPackageStartupMessages(library(imputeTS))

# ------------------------------------------------------------
# User-editable configuration
# ------------------------------------------------------------
DATA_FILE <- file.path(THIS_DIR, "data", "rk_mat2000-2022.Rdata")
DATA_OBJECT <- NULL
OUTPUT_DIR <- file.path(THIS_DIR, "output")
OUTPUT_STUB <- "scbm_vhar_empirical"

START_DATE <- as.Date("2010-01-01")
END_DATE <- as.Date("2019-12-31")
DROP_COLUMNS <- c(5L, 30L)

# model specification
bw <- 5L
bm <- 22L
horizon_names <- c("daily", "weekly", "monthly")
estimator <- "lasso"

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

# aggregated community counts in the order short / middle / long
n_comm <- c(3L, 3L, 3L)

# PisCES settings
alpha_mode <- "cv"
alpha <- NULL
alpha_fold <- 5L
alpha_grid <- scbm_default_alpha_grid()
alpha_criterion <- "holdout"

# clustering settings
nstart <- 50L
seed <- 12345L

# plot / workspace options
SAVE_WORKSPACE <- TRUE
#SAVE_PLOTS <- FALSE
MAKE_SANKEY <- TRUE
FIGURE_DIR <- file.path(OUTPUT_DIR, OUTPUT_STUB, "figures")
#SAVE_PLOTS <- TRUE

# ------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------
ensure_dir(OUTPUT_DIR)
APP_OUTPUT_DIR <- file.path(OUTPUT_DIR, OUTPUT_STUB)
ensure_dir(APP_OUTPUT_DIR)
ensure_dir(FIGURE_DIR)

obj <- load_rdata_object(DATA_FILE, object_name = DATA_OBJECT)
raw_df <- first_numeric_matrix(obj)

rv_df <- as.data.frame(raw_df, check.names = FALSE)
if (is.null(rownames(rv_df))) {
  stop("The realized-volatility input must have row names that can be parsed as dates.", call. = FALSE)
}
dates <- as.Date(rownames(rv_df))
if (any(is.na(dates))) {
  stop("Could not parse row names of the realized-volatility input as dates.", call. = FALSE)
}

idx_time <- dates >= START_DATE & dates <= END_DATE
rv_df <- rv_df[idx_time, , drop = FALSE]
dates <- dates[idx_time]

if (length(DROP_COLUMNS)) {
  keep_idx <- setdiff(seq_len(ncol(rv_df)), DROP_COLUMNS)
  rv_df <- rv_df[, keep_idx, drop = FALSE]
}

series_names <- colnames(rv_df)
series_names <- sub("^\\.", "", series_names)

logrv <- matrix(NA_real_, nrow = nrow(rv_df), ncol = ncol(rv_df))
for (j in seq_len(ncol(rv_df))) {
  x <- as.numeric(rv_df[[j]])
  x[x <= 0] <- NA_real_
  logrv[, j] <- imputeTS::na_interpolation(log(x))
}
colnames(logrv) <- series_names
rownames(logrv) <- as.character(dates)

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
# Step 2: clustering and PisCES
# ------------------------------------------------------------
cluster_vhar <- scbm_vhar_empirical_cluster(
  fit_vhar,
  bw = bw,
  bm = bm,
  horizon_names = horizon_names,
  series_names = series_names,
  n_comm = n_comm,
  alpha_mode = alpha_mode,
  alpha = alpha,
  alpha_fold = alpha_fold,
  alpha_grid = alpha_grid,
  alpha_criterion = alpha_criterion,
  nstart = nstart,
  seed = seed
)

# ------------------------------------------------------------
# Save workspace before plots
# ------------------------------------------------------------
if (isTRUE(SAVE_WORKSPACE)) {
  save.image(file = file.path(APP_OUTPUT_DIR, "workspace_before_plots.RData"))
}

# ------------------------------------------------------------
# Optional Sankey plot
# ------------------------------------------------------------
if (isTRUE(MAKE_SANKEY)) {
  need_pkgs <- c("ggplot2", "ggsankey", "dplyr", "tidyr")
  has_pkgs <- vapply(need_pkgs, requireNamespace, logical(1L), quietly = TRUE)

  if (!all(has_pkgs)) {
    message("Skipping VHAR Sankey plot because required packages are missing: ",
            paste(need_pkgs[!has_pkgs], collapse = ", "))
  } else {
    suppressPackageStartupMessages({
      library(ggplot2)
      library(ggsankey)
      library(dplyr)
      library(tidyr)
    })

    flow <- as.data.frame(cluster_vhar$flow_table[, c("long.cluster", "middle.cluster", "short.cluster"), drop = FALSE])
    mat <- as.data.frame(flow)
    rownames(mat) <- rownames(cluster_vhar$flow_table)

    long_box_cols <- c("1" = "#FDE2D5", "2" = "#DDF1E4", "3" = "#D9EAF7")
    middle_box_cols <- c("1" = "#F6B8A6", "2" = "#A9D4B5", "3" = "#9EC5E5")
    short_box_cols <- c("1" = "#E78B74", "2" = "#6D9CCF", "3" = "#74B68C")

    flow_cols <- c(
      "L1" = unname(long_box_cols["1"]),
      "L2" = unname(long_box_cols["2"]),
      "L3" = unname(long_box_cols["3"]),
      "M1" = unname(middle_box_cols["1"]),
      "M2" = unname(middle_box_cols["2"]),
      "M3" = unname(middle_box_cols["3"])
    )

    stages <- c("long.cluster", "middle.cluster", "short.cluster")
    flow_wide <- data.frame(series = rownames(mat), mat, check.names = FALSE)

    data_flow <- flow_wide %>%
      make_long(long.cluster, middle.cluster, short.cluster) %>%
      mutate(
        x_str = as.character(x),
        node_str = as.character(node),
        flow_group = case_when(
          x_str == "long.cluster"   ~ paste0("L", node_str),
          x_str == "middle.cluster" ~ paste0("M", node_str),
          x_str == "short.cluster"  ~ paste0("S", node_str),
          TRUE ~ NA_character_
        )
      ) %>%
      mutate(
        x = factor(x, levels = stages),
        next_x = factor(next_x, levels = stages),
        node = factor(node, levels = as.character(1:4)),
        next_node = factor(next_node, levels = as.character(1:4))
      ) %>%
      filter(!is.na(node))

    make_stage_layout <- function(v, x_pos, box_col_map, group_levels = 1:4, gap = 2, inner_pad = 0.55, box_halfwidth = 0.16) {
      cnt <- table(factor(v, levels = group_levels))
      present <- as.integer(names(cnt)[cnt > 0])
      total_height <- sum(cnt[cnt > 0]) + gap * (length(present) - 1)
      cur <- - total_height / 2
      y_ranges <- setNames(vector("list", length(group_levels)), as.character(group_levels))
      box_df <- data.frame()
      for (g in present) {
        n_g <- as.integer(cnt[as.character(g)])
        ymin <- cur
        ymax <- cur + n_g
        y_ranges[[as.character(g)]] <- c(ymin + inner_pad, ymax - inner_pad)
        box_df <- bind_rows(box_df, data.frame(
          x = x_pos, xmin = x_pos - box_halfwidth, xmax = x_pos + box_halfwidth,
          ymin = ymin, ymax = ymax, fill = unname(box_col_map[as.character(g)]),
          stringsAsFactors = FALSE
        ))
        cur <- ymax + gap
      }
      list(y_ranges = y_ranges, box_df = box_df)
    }

    make_labels <- function(x_pos, col_idx, y_ranges, mat, txt_size = 3.7) {
      out <- list()
      for (g in names(y_ranges)) {
        rng <- y_ranges[[g]]
        if (is.null(rng)) next
        idx <- which(as.integer(mat[, col_idx]) == as.integer(g))
        if (!length(idx)) next
        nms <- sort(as.character(rownames(mat)[idx]))
        n <- length(nms)
        ys <- if (n == 1) mean(rng) else seq(from = rng[2], to = rng[1], length.out = n)
        out <- c(out, list(
          annotate("text", x = x_pos, y = ys, label = nms,
                   colour = "black", size = txt_size, fontface = 2, hjust = 0.5)
        ))
      }
      out
    }

    lay_long   <- make_stage_layout(mat[, 1], x_pos = 1, box_col_map = long_box_cols, gap = 2.2)
    lay_middle <- make_stage_layout(mat[, 2], x_pos = 2, box_col_map = middle_box_cols, gap = 2.5)
    lay_short  <- make_stage_layout(mat[, 3], x_pos = 3, box_col_map = short_box_cols, gap = 2.5)
    box_df <- bind_rows(lay_long$box_df, lay_middle$box_df, lay_short$box_df)

    vhar_sankey_plot <- ggplot(data_flow, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = flow_group)) +
      geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = box_df$fill, color = NA) +
      geom_sankey(flow.alpha = 0.55, node.color = NA, width = 0.001, show.legend = FALSE) +
      geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = NA, color = "grey40", linewidth = 0.4) +
      scale_fill_manual(values = flow_cols, guide = "none", na.value = "transparent") +
      scale_x_discrete(
        labels = c(
          "long.cluster" = "Long (Monthly)",
          "middle.cluster" = "Medium (Weekly)",
          "short.cluster" = "Short (Daily)"
        )
      ) +
      theme_sankey(base_size = 14) +
      labs(x = NULL, y = NULL, fill = NULL) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold"),
        plot.margin = margin(10, 40, 10, 40)
      ) +
      make_labels(x_pos = 1, col_idx = 1, y_ranges = lay_long$y_ranges, mat = mat, txt_size = 3.7) +
      make_labels(x_pos = 2, col_idx = 2, y_ranges = lay_middle$y_ranges, mat = mat, txt_size = 3.7) +
      make_labels(x_pos = 3, col_idx = 3, y_ranges = lay_short$y_ranges, mat = mat, txt_size = 3.7)

    run_plot_block(
      save_plot = SAVE_PLOTS,
      file = file.path(FIGURE_DIR, "vhar_sankey.pdf"),
      width = 12,
      height = 8,
      code = {
        print(vhar_sankey_plot)
      }
    )
  }
}

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------
cat("\n==============================\n")
cat("ScBM-VHAR empirical analysis done\n")
cat("Estimator      : ", fit_vhar$estimator, "\n", sep = "")
cat("q              : ", fit_vhar$q, "\n", sep = "")
cat("T              : ", fit_vhar$TT, "\n", sep = "")
cat("bw             : ", fit_vhar$bw, "\n", sep = "")
cat("bm             : ", fit_vhar$bm, "\n", sep = "")
cat("alpha selected : ", round(cluster_vhar$alpha_sel, 6), "\n", sep = "")
if (!is.null(fit_vhar$first_stage$c_lambda_sel)) {
  cat("c_lambda       : ", round(fit_vhar$first_stage$c_lambda_sel, 6), "\n", sep = "")
}
cat("==============================\n\n")

cat("Aggregate community flow table:\n")
print(cluster_vhar$flow_table)
