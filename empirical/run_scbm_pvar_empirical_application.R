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

# ------------------------------------------------------------
# User-editable configuration
# ------------------------------------------------------------
DATA_FILE <- file.path(THIS_DIR, "data", "pvar_payroll_1990_2020.RData")
DATA_OBJECT <- NULL
OUTPUT_DIR <- file.path(THIS_DIR, "output")
OUTPUT_STUB <- "scbm_pvar_empirical"

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
centerTF <- FALSE
updateSigma <- TRUE
sigma_diag_only <- TRUE

# community counts: (Ky1, Kz1, Ky2, Kz2, Ky3, Kz3, Ky4, Kz4)
n_comm <- c(2L, 2L, 2L, 3L, 3L, 3L, 3L, 2L)

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
SAVE_PLOTS <- FALSE
MAKE_SANKEY <- TRUE
FIGURE_DIR <- file.path(OUTPUT_DIR, OUTPUT_STUB, "figures")
# SAVE_PLOTS <- TRUE

# ------------------------------------------------------------
# Load and parse data
# ------------------------------------------------------------
ensure_dir(OUTPUT_DIR)
APP_OUTPUT_DIR <- file.path(OUTPUT_DIR, OUTPUT_STUB)
ensure_dir(APP_OUTPUT_DIR)
ensure_dir(FIGURE_DIR)

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
# Step 2: clustering and PisCES
# ------------------------------------------------------------
cluster_pvar <- scbm_pvar_empirical_cluster(
  fit_pvar,
  s = s,
  p = p,
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
    message("Skipping PVAR Sankey plot because required packages are missing: ",
            paste(need_pkgs[!has_pkgs], collapse = ", "))
  } else {
    suppressPackageStartupMessages({
      library(ggplot2)
      library(ggsankey)
      library(dplyr)
      library(tidyr)
    })

    flow <- as.data.frame(cluster_pvar$flow_table[, c("S1.sending", "S2.sending", "S3.sending", "S4.sending"), drop = FALSE])
    mat <- as.data.frame(flow)
    rownames(mat) <- rownames(cluster_pvar$flow_table)
    for (j in seq_len(ncol(mat))) mat[[j]] <- as.integer(mat[[j]])

    stages <- c("S1.sending", "S2.sending", "S3.sending", "S4.sending")
    flow_wide <- data.frame(series = rownames(mat), mat, check.names = FALSE)

    s1_box_cols <- c("1" = "#E3EDF6", "2" = "#FDE8E8")
    s2_box_cols <- c("2" = "#C6DBEF", "1" = "#FBCFCE")
    s3_box_cols <- c("1" = "#9ECAE1", "2" = "#F8AFA8", "3" = "#CCEBC5")
    s4_box_cols <- c("1" = "#F08080", "2" = "#D9CE9B", "3" = "#6BAED6") 
    

    flow_cols <- c(
      "Stg1_1" = unname(s1_box_cols["1"]), "Stg1_2" = unname(s1_box_cols["2"]),
      "Stg2_1" = unname(s2_box_cols["1"]), "Stg2_2" = unname(s2_box_cols["2"]),
      "Stg3_1" = unname(s3_box_cols["1"]), "Stg3_2" = unname(s3_box_cols["2"]), "Stg3_3" = unname(s3_box_cols["3"]),
      "Stg4_1" = unname(s4_box_cols["1"]), "Stg4_2" = unname(s4_box_cols["2"]), "Stg4_3" = unname(s4_box_cols["3"])
    )

    data_flow <- flow_wide %>%
      make_long(S1.sending, S2.sending, S3.sending, S4.sending) %>%
      mutate(
        x_str = as.character(x),
        node_str = as.character(node),
        flow_group = case_when(
          x_str == stages[1L] ~ paste0("Stg1_", node_str),
          x_str == stages[2L] ~ paste0("Stg2_", node_str),
          x_str == stages[3L] ~ paste0("Stg3_", node_str),
          x_str == stages[4L] ~ paste0("Stg4_", node_str),
          TRUE ~ NA_character_
        )
      ) %>%
      mutate(
        x = factor(x, levels = stages),
        next_x = factor(next_x, levels = stages),
        node = factor(node, levels = as.character(1:3)),
        next_node = factor(next_node, levels = as.character(1:3))
      ) %>%
      filter(!is.na(node))

    make_stage_layout <- function(v, x_pos, box_col_map, group_levels = 1:3,
                                  gap = 2.7, inner_pad = 0.55, box_halfwidth = 0.25) {
      cnt <- table(factor(v, levels = group_levels))
      present <- as.integer(names(cnt)[cnt > 0L])
      total_height <- sum(cnt[cnt > 0L]) + gap * (length(present) - 1L)
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

    make_labels <- function(x_pos, col_idx, y_ranges, mat, txt_size = 3.2) {
      out <- list()
      for (g in names(y_ranges)) {
        rng <- y_ranges[[g]]
        if (is.null(rng)) next
        idx <- which(as.integer(mat[, col_idx]) == as.integer(g))
        if (!length(idx)) next
        nms <- sort(as.character(rownames(mat)[idx]))
        ys <- if (length(nms) == 1L) mean(rng) else seq(from = rng[2L], to = rng[1L], length.out = length(nms))
        out <- c(out, list(
          annotate("text", x = x_pos, y = ys, label = nms,
                   colour = "black", size = txt_size, fontface = 2, hjust = 0.5)
        ))
      }
      out
    }

    lay_s1 <- make_stage_layout(mat[, 1L], x_pos = 1, box_col_map = s1_box_cols)
    lay_s2 <- make_stage_layout(mat[, 2L], x_pos = 2, box_col_map = s2_box_cols)
    lay_s3 <- make_stage_layout(mat[, 3L], x_pos = 3, box_col_map = s3_box_cols)
    lay_s4 <- make_stage_layout(mat[, 4L], x_pos = 4, box_col_map = s4_box_cols)
    box_df <- bind_rows(lay_s1$box_df, lay_s2$box_df, lay_s3$box_df, lay_s4$box_df)

    pvar_sankey_plot <- ggplot(data_flow, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = flow_group)) +
      geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = box_df$fill, color = NA) +
      geom_sankey(flow.alpha = 0.5, node.color = NA, width = 0.001, show.legend = FALSE) +
      geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = NA, color = "grey40", linewidth = 0.4) +
      scale_fill_manual(values = flow_cols, guide = "none", na.value = "transparent") +
      scale_x_discrete(labels = c("S1.sending" = "Q1", "S2.sending" = "Q2", "S3.sending" = "Q3", "S4.sending" = "Q4")) +
      theme_sankey(base_size = 14) +
      labs(x = NULL, y = NULL) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(10, 40, 10, 40)
      ) +
      make_labels(x_pos = 1, col_idx = 1, y_ranges = lay_s1$y_ranges, mat = mat) +
      make_labels(x_pos = 2, col_idx = 2, y_ranges = lay_s2$y_ranges, mat = mat) +
      make_labels(x_pos = 3, col_idx = 3, y_ranges = lay_s3$y_ranges, mat = mat) +
      make_labels(x_pos = 4, col_idx = 4, y_ranges = lay_s4$y_ranges, mat = mat)

    run_plot_block(
      save_plot = SAVE_PLOTS,
      file = file.path(FIGURE_DIR, "pvar_sankey.pdf"),
      width = 12,
      height = 7,
      code = {
        print(pvar_sankey_plot)
      }
    )
  }
}

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------
cat("\n==============================\n")
cat("ScBM-PVAR empirical analysis done\n")
cat("Estimator      : ", fit_pvar$estimator, "\n", sep = "")
cat("q              : ", fit_pvar$q, "\n", sep = "")
cat("T              : ", fit_pvar$TT, "\n", sep = "")
cat("s              : ", fit_pvar$s, "\n", sep = "")
cat("p              : ", fit_pvar$p, "\n", sep = "")
cat("alpha selected : ", round(cluster_pvar$alpha_sel, 6), "\n", sep = "")
if (!is.null(fit_pvar$first_stage$c_lambda_sel)) {
  cat("c_lambda       : ", round(fit_pvar$first_stage$c_lambda_sel, 6), "\n", sep = "")
}
cat("==============================\n\n")

cat("Community flow table:\n")
print(cluster_pvar$flow_table)

