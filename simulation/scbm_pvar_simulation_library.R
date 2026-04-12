# ============================================================
# scbm_pvar_simulation_library.R
# GitHub-ready PVAR simulation library for
# 'Dynamic spectral co-clustering of directed networks to unveil
# latent community paths in VAR-type models'.
#
# This file keeps the current working PVAR simulation logic,
# including sparse fixed-design generation, OLS / lasso first-stage
# estimation, and optional PisCES smoothing for community recovery.
#
# Public wrappers added at the end of the file:
#   - scbm_pvar_run_simulation()
#   - scbm_pvar_simulation_summary()
#   - scbm_pvar_default_c_grid()
#   - scbm_pvar_default_alpha_grid()
# ============================================================


suppressPackageStartupMessages({
  library(MASS)
  library(gtools)
  library(foreach)
  library(doParallel)
  library(sparseVAR)
})

# ------------------------------
# Basic utilities
# ------------------------------
pvar2_safe_matrix <- function(X) {
  if (is.null(X)) return(matrix(NA_real_, 1L, 1L))
  X <- tryCatch(as.matrix(X), error = function(e) NULL)
  if (is.null(X)) return(matrix(NA_real_, 1L, 1L))
  if (length(dim(X)) != 2L) X <- matrix(X, ncol = 1L)
  X
}

pvar2_mean <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

pvar2_fro_norm <- function(A) {
  A <- pvar2_safe_matrix(A)
  if (any(!is.finite(A))) return(NA_real_)
  sqrt(sum(A^2))
}

pvar2_operator_norm <- function(A) {
  A <- pvar2_safe_matrix(A)
  if (any(!is.finite(A))) return(NA_real_)
  d <- tryCatch(svd(A, nu = 0L, nv = 0L)$d, error = function(e) numeric(0))
  if (!length(d)) return(NA_real_)
  as.numeric(d[1L])
}

pvar2_safe_row_normalize <- function(X) {
  X <- pvar2_safe_matrix(X)
  rn <- sqrt(rowSums(X^2))
  rn[!is.finite(rn) | rn <= .Machine$double.eps] <- 1
  X / rn
}

pvar2_path_id <- function(path) {
  if (is.numeric(path)) {
    pid <- as.integer(path)[1L]
    if (!pid %in% 1:4) stop("path must be 1, 2, 3, 4 or 'path1', 'path2', 'path3', 'path4'.", call. = FALSE)
    return(pid)
  }
  x <- tolower(as.character(path)[1L])
  if (!x %in% c("path1", "path2", "path3", "path4")) {
    stop("path must be one of 'path1', 'path2', 'path3', 'path4'.", call. = FALSE)
  }
  switch(x, path1 = 1L, path2 = 2L, path3 = 3L, path4 = 4L)
}

pvar2_path_name <- function(path) paste0("path", pvar2_path_id(path))

pvar2_link_type <- function(link_type) {
  x <- tolower(as.character(link_type)[1L])
  if (x %in% c("type1", "1", "t1")) return("type1")
  if (x %in% c("type2", "2", "t2")) return("type2")
  stop("link_type must be one of 'type1' or 'type2'.", call. = FALSE)
}

pvar2_method <- function(method) {
  x <- tolower(as.character(method)[1L])
  if (!x %in% c("ols", "lasso")) stop("method must be 'ols' or 'lasso'.", call. = FALSE)
  x
}

pvar2_lambda_mode <- function(lambda_mode = c("theory", "cv", "cv_c")) {
  match.arg(lambda_mode)
}

pvar2_default_c_grid <- function() seq(from = .1, to = 1, by = .05)

pvar2_round_df <- function(df, digits = 3L) {
  out <- df
  is_num <- vapply(out, is.numeric, logical(1L))
  out[is_num] <- lapply(out[is_num], function(z) round(z, digits))
  out
}

pvar2_normalize_seed <- function(seed, default = 1L) {
  mod <- 2147483647
  x <- suppressWarnings(as.numeric(seed)[1L])
  if (!is.finite(x) || is.na(x)) return(as.integer(default))
  x <- floor(abs(x)) %% mod
  if (x <= 0) x <- as.numeric(default)
  as.integer(x)
}

pvar2_make_seed <- function(base_seed, ..., default = 1L) {
  mod <- 2147483647
  x <- as.numeric(pvar2_normalize_seed(base_seed, default = default))
  vals <- unlist(list(...), use.names = FALSE)
  if (!length(vals)) return(as.integer(x))
  for (ii in seq_along(vals)) {
    v <- suppressWarnings(as.numeric(vals[ii]))
    if (!is.finite(v) || is.na(v)) v <- 0
    v <- floor(abs(v)) %% mod
    x <- ((x * 1000003) + v + 97 * ii) %% mod
    if (x <= 0) x <- as.numeric(default)
  }
  as.integer(x)
}


# ------------------------------
# Community path specification
# ------------------------------
pvar2_balanced_labels <- function(q, K) {
  q <- as.integer(q)
  K <- as.integer(K)
  sz <- rep(q %/% K, K)
  rem <- q %% K
  if (rem > 0L) sz[seq_len(rem)] <- sz[seq_len(rem)] + 1L
  rep(seq_len(K), times = sz)
}

pvar2_split_label_contiguous <- function(labels, label_to_split, new_labels) {
  labels <- as.integer(labels)
  new_labels <- as.integer(new_labels)
  idx <- which(labels == as.integer(label_to_split)[1L])
  if (!length(idx)) return(labels)
  repl <- pvar2_balanced_labels(length(idx), length(new_labels))
  out <- labels
  out[idx] <- new_labels[repl]
  out
}

pvar2_half_swap <- function(labels) {
  labels <- as.integer(labels)
  if (max(labels) != 2L) stop("pvar2_half_swap currently implemented for K=2 only.", call. = FALSE)
  idx1 <- which(labels == 1L)
  idx2 <- which(labels == 2L)
  nswap <- min(length(idx1), length(idx2)) %/% 2L
  if (nswap <= 0L) return(labels)
  take1 <- tail(idx1, nswap)
  take2 <- head(idx2, nswap)
  out <- labels
  out[take1] <- 2L
  out[take2] <- 1L
  out
}

pvar2_path_spec <- function(path) {
  pid <- pvar2_path_id(path)
  if (pid == 1L) return(c(4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L))
  if (pid == 4L) return(c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L))
  if (pid == 3L) return(c(2L, 2L, 2L, 3L, 3L, 4L, 4L, 2L))
  if (pid == 2L) return(c(2L, 3L, 3L, 3L, 3L, 2L, 2L, 2L))
  stop("Unsupported path.", call. = FALSE)
}

pvar2_make_persistent_paths <- function(q, path) {
  path_name <- pvar2_path_name(path)

  if (path_name == "path1") {
    y1 <- pvar2_balanced_labels(q, 4L)
    y2 <- y1
    y3 <- y1
    y4 <- y1
  } else if (path_name == "path4") {
    y1 <- pvar2_balanced_labels(q, 2L)
    y2 <- pvar2_half_swap(y1)
    y3 <- y2
    y4 <- y2
  } else if (path_name == "path3") {
    y1 <- pvar2_balanced_labels(q, 2L)
    y2 <- y1
    y3 <- pvar2_split_label_contiguous(y2, label_to_split = 2L, new_labels = c(2L, 3L))
    y4 <- y3
    idx1 <- which(y3 == 1L)
    if (length(idx1)) {
      repl <- pvar2_balanced_labels(length(idx1), 2L)
      y4[idx1] <- c(1L, 2L)[repl]
    }
    y4[y3 == 2L] <- 3L
    y4[y3 == 3L] <- 4L
  } else if (path_name == "path2") {
    y1 <- pvar2_balanced_labels(q, 2L)
    idx1 <- which(y1 == 1L)
    idx2 <- which(y1 == 2L)

    n11 <- floor((2 * length(idx1)) / 3)
    n21 <- floor(length(idx2) / 3)
    if (n11 <= 0L || n21 <= 0L) stop("q is too small for path2 construction.", call. = FALSE)

    y2 <- integer(q)
    # S1:C1 -> S2:C1 (upper 2/3), S2:C2 (lower 1/3)
    y2[idx1[seq_len(n11)]] <- 1L
    y2[idx1[(n11 + 1L):length(idx1)]] <- 2L
    # S1:C2 -> S2:C2 (upper 1/3), S2:C3 (lower 2/3)
    y2[idx2[seq_len(n21)]] <- 2L
    y2[idx2[(n21 + 1L):length(idx2)]] <- 3L

    y3 <- y2

    y4 <- integer(q)
    idx_c1 <- which(y3 == 1L)
    idx_c2 <- which(y3 == 2L)
    idx_c3 <- which(y3 == 3L)
    nmid <- floor(length(idx_c2) / 2)
    if (nmid <= 0L) stop("Middle block too small for path4 merge step.", call. = FALSE)

    # S3:C1 -> S4:C1; S3:C2 split half-half; S3:C3 -> S4:C2
    y4[idx_c1] <- 1L
    y4[idx_c2[seq_len(nmid)]] <- 1L
    y4[idx_c2[(nmid + 1L):length(idx_c2)]] <- 2L
    y4[idx_c3] <- 2L
  } else {
    stop("Unsupported path.", call. = FALSE)
  }

  y_list <- list(as.integer(y1), as.integer(y2), as.integer(y3), as.integer(y4))
  z_list <- list(y_list[[2L]], y_list[[3L]], y_list[[4L]], y_list[[1L]])
  list(y = y_list, z = z_list)
}

pvar2_generate_path_labels <- function(q, n_comm, path = "path1") {
  path_obj <- pvar2_make_persistent_paths(q = q, path = path)
  s <- length(path_obj$y)
  if (length(n_comm) != 2L * s) stop("n_comm length mismatch.", call. = FALSE)
  for (mm in seq_len(s)) {
    Ky <- n_comm[2L * mm - 1L]
    Kz <- n_comm[2L * mm]
    if (max(path_obj$y[[mm]]) != Ky) stop("Sender community count mismatch in season ", mm, ".", call. = FALSE)
    if (max(path_obj$z[[mm]]) != Kz) stop("Receiver community count mismatch in season ", mm, ".", call. = FALSE)
  }
  path_obj
}

# ------------------------------
# Sparse benchmark DGP generation
# ------------------------------
pvar2_make_B <- function(Ky, Kz, link_type = "type1") {
  link_type <- pvar2_link_type(link_type)
  bu <- 0.05
  bl <- if (link_type == "type1") 0.10 else 0.15
  B <- matrix(bl, nrow = Ky, ncol = Kz)
  for (ii in seq_len(Ky)) {
    for (jj in seq_len(Kz)) {
      if (ii < jj) B[ii, jj] <- bu
    }
  }
  dd <- min(Ky, Kz)
  for (ii in seq_len(dd)) B[ii, ii] <- 0.50
  B
}

pvar2_block_pair_type <- function(a, b, Ky, Kz) {
  dd <- min(Ky, Kz)
  if (a <= dd && b <= dd && a == b) return("diag")
  if (a < b) return("upper")
  "lower"
}

pvar2_sparse_support_targets <- function(link_type = "type1") {
  link_type <- pvar2_link_type(link_type)
  if (link_type == "type1") {
    return(c(diag = 2.6, upper = 0.6, lower = 0.9))
  }
  c(diag = 2.2, upper = 0.7, lower = 1.3)
}

pvar2_sparse_magnitudes <- function(link_type = "type1") {
  link_type <- pvar2_link_type(link_type)
  if (link_type == "type1") {
    return(c(self = 0.30, diag = 0.14, upper = 0.04, lower = 0.06))
  }
  c(self = 0.28, diag = 0.12, upper = 0.05, lower = 0.07)
}

pvar2_make_sparse_prob <- function(y_labels, z_labels, link_type = "type1", cap = 0.95) {
  y_labels <- as.integer(y_labels)
  z_labels <- as.integer(z_labels)
  q <- length(y_labels)
  Ky <- max(y_labels)
  Kz <- max(z_labels)
  z_sizes <- as.numeric(tabulate(z_labels, nbins = Kz))
  tg <- pvar2_sparse_support_targets(link_type)

  P <- matrix(0, nrow = q, ncol = q)
  for (ii in seq_len(q)) {
    a <- y_labels[ii]
    for (jj in seq_len(q)) {
      if (ii == jj) next
      b <- z_labels[jj]
      pair_type <- pvar2_block_pair_type(a, b, Ky, Kz)
      target <- tg[pair_type]
      den <- max(1, z_sizes[b] - as.integer(a == b))
      P[ii, jj] <- min(cap, as.numeric(target) / den)
    }
  }
  P
}

pvar2_sample_sparse_support <- function(y_labels, z_labels, link_type = "type1") {
  Prob <- pvar2_make_sparse_prob(y_labels = y_labels, z_labels = z_labels, link_type = link_type)
  q <- nrow(Prob)
  A <- matrix(rbinom(q * q, 1L, as.vector(Prob)), nrow = q, ncol = q)
  diag(A) <- 0L
  list(A = A, Prob = Prob)
}

pvar2_build_sparse_phi <- function(A, y_labels, z_labels, link_type = "type1", self_value = NULL) {
  A <- pvar2_safe_matrix(A)
  q <- nrow(A)
  Ky <- max(y_labels)
  Kz <- max(z_labels)
  mg <- pvar2_sparse_magnitudes(link_type)
  if (is.null(self_value)) self_value <- mg[["self"]]

  Phi <- matrix(0, nrow = q, ncol = q)
  for (ii in seq_len(q)) {
    Phi[ii, ii] <- as.numeric(self_value) * runif(1L, 0.95, 1.05)
    a <- y_labels[ii]
    for (jj in seq_len(q)) {
      if (ii == jj || A[ii, jj] == 0) next
      b <- z_labels[jj]
      pair_type <- pvar2_block_pair_type(a, b, Ky, Kz)
      base <- as.numeric(mg[pair_type])
      Phi[ii, jj] <- base * runif(1L, 0.90, 1.10)
    }
  }
  Phi
}

pvar2_largest_singular <- function(A) {
  A <- pvar2_safe_matrix(A)
  d <- tryCatch(svd(A, nu = 0L, nv = 0L)$d, error = function(e) numeric(0))
  if (!length(d)) return(0)
  max(d)
}

pvar2_spectral_radius <- function(A) {
  ev <- tryCatch(eigen(pvar2_safe_matrix(A), only.values = TRUE)$values, error = function(e) complex(0))
  if (!length(ev)) return(NA_real_)
  max(Mod(ev))
}

pvar2_cycle_radius <- function(Phi_list) {
  prod_mat <- Reduce("%*%", rev(Phi_list))
  pvar2_spectral_radius(prod_mat)
}

pvar2_generate_benchmark_pvar <- function(q, path, link_type,
                                          target_cycle_sv = 0.90,
                                          sparse_prob_cap = 0.95) {
  q <- as.integer(q)
  n_comm <- pvar2_path_spec(path)
  path_obj <- pvar2_generate_path_labels(q = q, n_comm = n_comm, path = path)
  s <- length(path_obj$y)

  B_list <- vector("list", s)
  A_list <- vector("list", s)
  Prob_list <- vector("list", s)
  Phi_raw <- vector("list", s)
  nnz_row_mean <- numeric(s)
  density <- numeric(s)

  for (mm in seq_len(s)) {
    Ky <- n_comm[2L * mm - 1L]
    Kz <- n_comm[2L * mm]
    B_list[[mm]] <- pvar2_make_B(Ky = Ky, Kz = Kz, link_type = link_type)
    samp <- pvar2_sample_sparse_support(
      y_labels = path_obj$y[[mm]],
      z_labels = path_obj$z[[mm]],
      link_type = link_type
    )
    A_list[[mm]] <- samp$A
    Prob_list[[mm]] <- samp$Prob
    Phi_raw[[mm]] <- pvar2_build_sparse_phi(
      A = samp$A,
      y_labels = path_obj$y[[mm]],
      z_labels = path_obj$z[[mm]],
      link_type = link_type
    )
    nnz_row_mean[mm] <- mean(rowSums(samp$A != 0) + 1)
    density[mm] <- mean(Phi_raw[[mm]] != 0)
  }

  raw_prod <- Reduce("%*%", rev(Phi_raw))
  sv_raw <- pvar2_largest_singular(raw_prod)
  if (!is.finite(sv_raw) || sv_raw <= 0) stop("Invalid cycle singular value.", call. = FALSE)
  scale <- (as.numeric(target_cycle_sv) / sv_raw)^(1 / s)
  Phi_list <- lapply(Phi_raw, function(P) scale * P)

  list(
    y = path_obj$y,
    z = path_obj$z,
    B = B_list,
    Prob = Prob_list,
    A = A_list,
    Phi_raw = Phi_raw,
    Phi = Phi_list,
    scale = scale,
    rho = pvar2_cycle_radius(Phi_list),
    target_cycle_sv = as.numeric(target_cycle_sv),
    nnz_row_mean = nnz_row_mean,
    density = density,
    same_phi_across_rep = TRUE,
    same_phi_across_T = TRUE
  )
}

pvar2_simulate <- function(TT, Phi_list, sigma_diag = 0.5, burn = 500L, seed_data = NULL) {
  if (!is.null(seed_data)) set.seed(pvar2_normalize_seed(seed_data))
  s <- length(Phi_list)
  q <- nrow(Phi_list[[1L]])
  Sigma <- diag(as.numeric(sigma_diag), q)
  innovations <- MASS::mvrnorm(n = TT + burn, mu = rep(0, q), Sigma = Sigma)

  Y <- matrix(0, nrow = TT + burn, ncol = q)
  Y[1L, ] <- innovations[1L, ]
  if (TT + burn >= 2L) {
    for (tt in 2:(TT + burn)) {
      mm <- ((tt - 1L) %% s) + 1L
      Y[tt, ] <- drop(Phi_list[[mm]] %*% Y[tt - 1L, ]) + innovations[tt, ]
    }
  }
  t(Y[-seq_len(burn), , drop = FALSE])
}


# ------------------------------
# OLS and lasso estimation
# ------------------------------
# Note:
# sparseVAR::sPVAR_adalasso_fista currently performs its own internal
# lambda selection. For cv_c, we therefore keep the same sparseVAR FISTA
# back-end but refit season by season at fixed lambdas determined by a
# common c_lambda grid, which mirrors the library32.R strategy used for
# the sparse VHAR case.
pvar2_split_blocks <- function(Phi_hat, s = 4L) {
  Phi_hat <- pvar2_safe_matrix(Phi_hat)
  q <- nrow(Phi_hat)
  lapply(seq_len(s), function(mm) Phi_hat[, ((mm - 1L) * q + 1L):(mm * q), drop = FALSE])
}

pvar2_ols_internal <- function(Yt, s = 4L) {
  TT <- ncol(Yt)
  q <- nrow(Yt)
  Phi_hat <- matrix(NA_real_, nrow = q, ncol = q * s)
  for (mm in seq_len(s)) {
    idx <- which(((seq_len(TT) - 1L) %% s) + 1L == mm)
    idx <- idx[idx > 1L]
    Z <- Yt[, idx, drop = FALSE]
    X <- Yt[, idx - 1L, drop = FALSE]
    XXt <- X %*% t(X) + 1e-8 * diag(1, q)
    Phi_hat[, ((mm - 1L) * q + 1L):(mm * q)] <- solve(XXt, X %*% t(Z))
  }
  list(Phi_hat = Phi_hat)
}

pvar2_fit_ols <- function(Yt, s = 4L) {
  fit <- tryCatch(sparseVAR::PVAR_ols(Yt = Yt, s = s, p = 1), error = function(e) NULL)
  if (is.null(fit)) fit <- pvar2_ols_internal(Yt, s = s)
  list(
    Phi_hat = pvar2_safe_matrix(fit$Phi_hat),
    Phi_blocks = pvar2_split_blocks(fit$Phi_hat, s = s),
    lambda_pkg = rep(NA_real_, s),
    lambda_text = rep(NA_real_, s),
    c_lambda_sel = NA_real_,
    lambda_mode = "ols"
  )
}

pvar2_center_series <- function(Yt, s = 4L, centerTF = TRUE) {
  if (!isTRUE(centerTF)) return(list(Yc = pvar2_safe_matrix(Yt), mu = NULL))
  q <- nrow(Yt)
  TT <- ncol(Yt)
  mu <- matrix(0, q, s)
  for (j in seq_len(s)) mu[, j] <- rowMeans(Yt[, seq(j, TT, by = s), drop = FALSE])
  season_id <- ((seq_len(TT) - 1L) %% s) + 1L
  list(Yc = Yt - mu[, season_id, drop = FALSE], mu = mu)
}

pvar2_get_season_design <- function(Yc, s = 4L, season = 1L) {
  TT <- ncol(Yc)
  ids <- seq(from = as.integer(season), to = TT, by = s)
  ids <- ids[ids > 1L]
  n <- length(ids)
  if (n <= 1L) stop("Not enough observations for season ", season, ".", call. = FALSE)
  X <- t(Yc[, ids - 1L, drop = FALSE])
  Y <- t(Yc[, ids, drop = FALSE])
  list(X = X, Y = Y, id = ids)
}

pvar2_season_sample_sizes <- function(TT, s = 4L) {
  vapply(seq_len(s), function(mm) {
    ids <- which(((seq_len(TT) - 1L) %% s) + 1L == mm)
    ids <- ids[ids > 1L]
    length(ids)
  }, integer(1L))
}

pvar2_lambda_textbook <- function(q, s, n_eff, theory_const = 1) {
  as.numeric(theory_const) * sqrt(log(as.numeric(s) * q^2) / as.numeric(n_eff))
}

pvar2_lambda_pkg_from_text <- function(lambda_text, n_eff) {
  0.5 * as.numeric(n_eff) * as.numeric(lambda_text)
}

pvar2_lambda_text_from_pkg <- function(lambda_pkg, n_eff) {
  if (any(!is.finite(lambda_pkg)) || any(!is.finite(n_eff)) || any(n_eff <= 0)) {
    out <- rep(NA_real_, max(length(lambda_pkg), length(n_eff)))
    return(out)
  }
  2 * as.numeric(lambda_pkg) / as.numeric(n_eff)
}

pvar2_weight_matrix <- function(q, diagTF = TRUE) {
  W <- matrix(1, q, q)
  if (isTRUE(diagTF)) diag(W) <- 0
  W
}

pvar2_sigma_inv_from_ols <- function(X, Y, sigma_diag_only = TRUE) {
  eps <- 1e-8
  B0 <- solve(crossprod(X) + diag(eps, ncol(X)), crossprod(X, Y))
  R0 <- Y - X %*% B0
  Sig <- crossprod(R0) / nrow(R0)
  if (isTRUE(sigma_diag_only)) {
    diag(1 / pmax(diag(Sig), 1e-12), ncol(Y), ncol(Y))
  } else {
    sv <- svd(Sig)
    dd <- pmax(sv$d, 1e-12)
    sv$u %*% diag(1 / dd, ncol(Y)) %*% t(sv$u)
  }
}

pvar2_fit_lasso_fixed <- function(Yt, s = 4L,
                                  lambda_pkg,
                                  diagTF = TRUE,
                                  centerTF = TRUE,
                                  updateSigma = FALSE,
                                  sigma_diag_only = TRUE,
                                  max_iter = 1000L,
                                  tol = 1e-6) {
  q <- nrow(Yt)
  TT <- ncol(Yt)
  lambda_pkg <- as.numeric(lambda_pkg)
  if (length(lambda_pkg) == 1L) lambda_pkg <- rep(lambda_pkg, s)
  if (length(lambda_pkg) != s) stop("lambda_pkg must have length 1 or s.", call. = FALSE)

  cen <- pvar2_center_series(Yt, s = s, centerTF = centerTF)
  Yc <- cen$Yc
  Phi_hat <- matrix(0, q, q * s)
  W0 <- pvar2_weight_matrix(q = q, diagTF = diagTF)

  for (mm in seq_len(s)) {
    des <- pvar2_get_season_design(Yc, s = s, season = mm)
    X <- as.matrix(des$X)
    Y <- as.matrix(des$Y)
    SigmaInv <- NULL
    if (isTRUE(updateSigma)) {
      SigmaInv <- pvar2_sigma_inv_from_ols(X, Y, sigma_diag_only = sigma_diag_only)
    }
    Bhat <- sparseVAR:::fista_lasso_multi_cpp(
      X = X,
      Y = Y,
      weights = W0,
      lambda = lambda_pkg[mm],
      SigmaInv = if (is.null(SigmaInv)) matrix(0, 0, 0) else SigmaInv,
      max_iter = as.integer(max_iter),
      tol = tol
    )
    cols <- ((mm - 1L) * q + 1L):(mm * q)
    Phi_hat[, cols] <- t(Bhat)
  }

  n_eff <- pvar2_season_sample_sizes(TT, s = s)
  lambda_text <- pvar2_lambda_text_from_pkg(lambda_pkg, n_eff)
  list(
    Phi_hat = Phi_hat,
    Phi_blocks = pvar2_split_blocks(Phi_hat, s = s),
    lambda_pkg = lambda_pkg,
    lambda_text = lambda_text
  )
}

pvar2_cv_score_one_c <- function(Yc, s = 4L, c_value,
                                 q,
                                 fold = 10L,
                                 diagTF = TRUE,
                                 updateSigma = FALSE,
                                 sigma_diag_only = TRUE,
                                 max_iter = 1000L,
                                 tol = 1e-6) {
  total <- 0
  W0 <- pvar2_weight_matrix(q = q, diagTF = diagTF)

  for (mm in seq_len(s)) {
    des <- pvar2_get_season_design(Yc, s = s, season = mm)
    X <- as.matrix(des$X)
    Y <- as.matrix(des$Y)
    n_m <- nrow(X)
    lambda_text_base <- pvar2_lambda_textbook(q = q, s = s, n_eff = n_m, theory_const = 1)
    lambda_pkg_base <- pvar2_lambda_pkg_from_text(lambda_text_base, n_eff = n_m)
    lambda_use <- as.numeric(c_value) * as.numeric(lambda_pkg_base)

    SigmaInv <- NULL
    if (isTRUE(updateSigma)) {
      SigmaInv <- pvar2_sigma_inv_from_ols(X, Y, sigma_diag_only = sigma_diag_only)
    }

    fold_use <- max(2L, min(as.integer(fold), n_m))
    folds <- rep(seq_len(fold_use), each = floor(n_m / fold_use) + 1L)[1L:n_m]
    for (ff in seq_len(fold_use)) {
      tr <- which(folds != ff)
      va <- which(folds == ff)
      Bhat <- sparseVAR:::fista_lasso_multi_cpp(
        X = X[tr, , drop = FALSE],
        Y = Y[tr, , drop = FALSE],
        weights = W0,
        lambda = lambda_use,
        SigmaInv = if (is.null(SigmaInv)) matrix(0, 0, 0) else SigmaInv,
        max_iter = as.integer(max_iter),
        tol = tol
      )
      Rv <- Y[va, , drop = FALSE] - X[va, , drop = FALSE] %*% Bhat
      if (is.null(SigmaInv)) {
        total <- total + sum(Rv^2)
      } else {
        total <- total + sum((Rv %*% SigmaInv) * Rv)
      }
    }
  }
  total
}

pvar2_select_c_lambda_blockcv <- function(Yt, s = 4L,
                                          c_grid = pvar2_default_c_grid(),
                                          fold = 10L,
                                          diagTF = TRUE,
                                          centerTF = TRUE,
                                          updateSigma = FALSE,
                                          sigma_diag_only = TRUE,
                                          max_iter = 1000L,
                                          tol = 1e-6) {
  q <- nrow(Yt)
  c_grid <- unique(as.numeric(c_grid))
  c_grid <- c_grid[is.finite(c_grid) & c_grid > 0]
  if (!length(c_grid)) stop("c_grid must contain positive finite values.", call. = FALSE)
  c_grid <- sort(c_grid)

  cen <- pvar2_center_series(Yt, s = s, centerTF = centerTF)
  Yc <- cen$Yc
  scores <- vapply(c_grid, function(cc) {
    pvar2_cv_score_one_c(
      Yc = Yc, s = s, c_value = cc, q = q,
      fold = fold, diagTF = diagTF,
      updateSigma = updateSigma, sigma_diag_only = sigma_diag_only,
      max_iter = max_iter, tol = tol
    )
  }, numeric(1L))

  id_best <- which.min(scores)
  c_sel <- c_grid[id_best]
  n_eff <- pvar2_season_sample_sizes(ncol(Yt), s = s)
  lambda_text_base <- pvar2_lambda_textbook(q = q, s = s, n_eff = n_eff, theory_const = 1)
  lambda_pkg_base <- pvar2_lambda_pkg_from_text(lambda_text_base, n_eff = n_eff)

  list(
    c_lambda_sel = as.numeric(c_sel),
    score = scores,
    c_grid = c_grid,
    lambda_pkg_sel = as.numeric(c_sel) * lambda_pkg_base,
    lambda_text_sel = as.numeric(c_sel) * lambda_text_base,
    lambda_pkg_base = lambda_pkg_base,
    lambda_text_base = lambda_text_base
  )
}

pvar2_fit_lasso <- function(Yt, s = 4L,
                            lambda_mode = c("theory", "cv", "cv_c"),
                            c_lambda = 0.20,
                            c_grid = pvar2_default_c_grid(),
                            fold = 10L,
                            nlambda = 50L,
                            diagTF = TRUE,
                            centerTF = TRUE,
                            updateSigma = FALSE,
                            sigma_diag_only = TRUE,
                            max_outer = 5L,
                            tol_outer = 1e-2,
                            max_iter = 1000L,
                            tol = 1e-6) {
  lambda_mode <- pvar2_lambda_mode(lambda_mode)
  q <- nrow(Yt)
  TT <- ncol(Yt)
  n_eff <- pvar2_season_sample_sizes(TT, s = s)
  lambda_text_base <- pvar2_lambda_textbook(q = q, s = s, n_eff = n_eff, theory_const = 1)
  lambda_pkg_base <- pvar2_lambda_pkg_from_text(lambda_text_base, n_eff = n_eff)

  if (lambda_mode == "cv") {
    fit <- sparseVAR::sPVAR_adalasso_fista(
      Yt = Yt,
      s = s,
      p = 1,
      type = "lasso",
      fold = as.integer(fold),
      diagTF = isTRUE(diagTF),
      centerTF = isTRUE(centerTF),
      updateSigma = isTRUE(updateSigma),
      sigma_diag_only = isTRUE(sigma_diag_only),
      nlambda = as.integer(nlambda),
      max_outer = as.integer(max_outer),
      tol_outer = tol_outer,
      max_iter = as.integer(max_iter),
      tol = tol
    )
    lambda_pkg <- as.numeric(fit$lambda_lasso)
    if (length(lambda_pkg) == 1L) lambda_pkg <- rep(lambda_pkg, s)
    lambda_text <- pvar2_lambda_text_from_pkg(lambda_pkg, n_eff)
    c_vec <- lambda_pkg / lambda_pkg_base
    c_sel <- pvar2_mean(c_vec)
    Phi_hat <- pvar2_safe_matrix(fit$Phi_hat_lasso)
  } else {
    c_sel <- as.numeric(c_lambda)[1L]
    if (lambda_mode == "cv_c") {
      sel <- pvar2_select_c_lambda_blockcv(
        Yt = Yt, s = s, c_grid = c_grid,
        fold = fold, diagTF = diagTF, centerTF = centerTF,
        updateSigma = updateSigma, sigma_diag_only = sigma_diag_only,
        max_iter = max_iter, tol = tol
      )
      c_sel <- sel$c_lambda_sel
    }
    lambda_pkg <- c_sel * lambda_pkg_base
    fixed_fit <- pvar2_fit_lasso_fixed(
      Yt = Yt, s = s, lambda_pkg = lambda_pkg,
      diagTF = diagTF, centerTF = centerTF,
      updateSigma = updateSigma, sigma_diag_only = sigma_diag_only,
      max_iter = max_iter, tol = tol
    )
    Phi_hat <- fixed_fit$Phi_hat
    lambda_text <- fixed_fit$lambda_text
  }

  list(
    Phi_hat = Phi_hat,
    Phi_blocks = pvar2_split_blocks(Phi_hat, s = s),
    lambda_pkg = as.numeric(lambda_pkg),
    lambda_text = as.numeric(lambda_text),
    c_lambda_sel = as.numeric(c_sel),
    lambda_mode = lambda_mode
  )
}

# ------------------------------
# Spectral co-clustering and scoring
# ------------------------------
pvar2_embeddings_from_blocks <- function(blocks, n_comm) {
  s <- length(blocks)
  emb_L <- vector("list", s)
  emb_R <- vector("list", s)

  for (mm in seq_len(s)) {
    Ky <- n_comm[2L * mm - 1L]
    Kz <- n_comm[2L * mm]
    Ahat <- t(blocks[[mm]])
    sv <- svd(Ahat, nu = Ky, nv = Kz)
    XL <- sv$u[, seq_len(Ky), drop = FALSE]
    XR <- sv$v[, seq_len(Kz), drop = FALSE]
    emb_L[[mm]] <- pvar2_safe_row_normalize(XL)
    emb_R[[mm]] <- pvar2_safe_row_normalize(XR)
  }
  list(emb_L = emb_L, emb_R = emb_R)
}

pvar2_kmeans_labels <- function(X, K, nstart = 20L) {
  if (K <= 1L) return(rep(1L, nrow(X)))
  X <- pvar2_safe_matrix(X)
  ux_n <- nrow(unique(round(X, 12)))
  if (ux_n < K) {
    set.seed(1L)
    X <- X + matrix(rnorm(length(X), sd = 1e-8), nrow = nrow(X), ncol = ncol(X))
  }
  fit <- stats::kmeans(X, centers = K, nstart = as.integer(nstart), iter.max = 100L)
  as.integer(fit$cluster)
}

pvar2_cluster <- function(blocks, n_comm) {
  s <- length(blocks)
  tmp <- pvar2_embeddings_from_blocks(blocks, n_comm)
  emb_L <- tmp$emb_L
  emb_R <- tmp$emb_R

  group_L <- vector("list", s)
  group_R <- vector("list", s)

  Ky1 <- n_comm[1L]
  Kzs <- n_comm[2L * s]
  if (Ky1 == Kzs) {
    lab <- pvar2_kmeans_labels(cbind(emb_R[[s]], emb_L[[1L]]), Ky1)
    group_R[[s]] <- lab
    group_L[[1L]] <- lab
  } else {
    group_R[[s]] <- pvar2_kmeans_labels(emb_R[[s]], Kzs)
    group_L[[1L]] <- pvar2_kmeans_labels(emb_L[[1L]], Ky1)
  }

  if (s >= 2L) {
    for (mm in 2:s) {
      Ky <- n_comm[2L * mm - 1L]
      Kzprev <- n_comm[2L * (mm - 1L)]
      if (Ky == Kzprev) {
        lab <- pvar2_kmeans_labels(cbind(emb_R[[mm - 1L]], emb_L[[mm]]), Ky)
        group_R[[mm - 1L]] <- lab
        group_L[[mm]] <- lab
      } else {
        group_R[[mm - 1L]] <- pvar2_kmeans_labels(emb_R[[mm - 1L]], Kzprev)
        group_L[[mm]] <- pvar2_kmeans_labels(emb_L[[mm]], Ky)
      }
    }
  }

  list(group_L = group_L, group_R = group_R)
}

pvar2_all_perms <- function(x) {
  x <- as.integer(x)
  if (length(x) <= 1L) return(matrix(x, nrow = 1L))
  gtools::permutations(n = length(x), r = length(x), v = x)
}

pvar2_adjusted_rand_index <- function(labels_true, labels_est) {
  x <- as.integer(labels_true)
  y <- as.integer(labels_est)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n <= 1L) return(1)
  tab <- table(x, y)
  nij2 <- sum(tab * (tab - 1) / 2)
  ai <- rowSums(tab)
  bj <- colSums(tab)
  ai2 <- sum(ai * (ai - 1) / 2)
  bj2 <- sum(bj * (bj - 1) / 2)
  n2 <- n * (n - 1) / 2
  exp_idx <- ai2 * bj2 / n2
  max_idx <- 0.5 * (ai2 + bj2)
  den <- max_idx - exp_idx
  if (abs(den) <= .Machine$double.eps) return(0)
  as.numeric((nij2 - exp_idx) / den)
}

pvar2_match_accuracy <- function(true_labels, est_labels) {
  tru <- as.integer(true_labels)
  est <- as.integer(est_labels)
  ok <- is.finite(tru) & is.finite(est)
  tru <- tru[ok]
  est <- est[ok]
  if (!length(tru)) return(list(acc = NA_real_, ari = NA_real_, matched = integer(0)))
  K <- max(max(tru), max(est))
  perms <- pvar2_all_perms(seq_len(K))
  best_acc <- -Inf
  best_lab <- est
  for (ii in seq_len(nrow(perms))) {
    mp <- perms[ii, ]
    est2 <- mp[est]
    acc <- mean(tru == est2)
    if (acc > best_acc) {
      best_acc <- acc
      best_lab <- est2
    }
  }
  list(acc = as.numeric(best_acc), ari = pvar2_adjusted_rand_index(tru, best_lab), matched = best_lab)
}

pvar2_score <- function(true_y, fit_obj) {
  s <- length(true_y)
  acc <- ari <- numeric(s)
  for (mm in seq_len(s)) {
    met <- pvar2_match_accuracy(true_y[[mm]], fit_obj$group_L[[mm]])
    acc[mm] <- met$acc
    ari[mm] <- met$ari
  }
  names(acc) <- paste0("season", seq_len(s))
  names(ari) <- paste0("season", seq_len(s))
  list(acc = acc, ari = ari, acc_overall = mean(acc), ari_overall = mean(ari))
}

pvar2_fro_errors <- function(Phi_true, Phi_hat) {
  s <- length(Phi_true)
  out <- numeric(s)
  names(out) <- paste0("fro_s", seq_len(s))
  for (mm in seq_len(s)) out[mm] <- pvar2_fro_norm(Phi_hat[[mm]] - Phi_true[[mm]])
  out
}

pvar2_op_errors <- function(Phi_true, Phi_hat) {
  s <- length(Phi_true)
  out <- numeric(s)
  names(out) <- paste0("op_s", seq_len(s))
  for (mm in seq_len(s)) out[mm] <- pvar2_operator_norm(Phi_hat[[mm]] - Phi_true[[mm]])
  out
}

# ------------------------------
# Result rows and runner
# ------------------------------
pvar2_make_result_row <- function(path, q, link_type, estimator, T, rep_id,
                                  seed, design_id,
                                  score = NULL, frob = NULL, opb = NULL,
                                  c_lambda_sel = NA_real_,
                                  lambda_pkg = rep(NA_real_, 4L),
                                  lambda_text = rep(NA_real_, 4L),
                                  lambda_mode = NA_character_,
                                  error_message = "") {
  if (is.null(score)) score <- list(acc = rep(NA_real_, 4L), ari = rep(NA_real_, 4L), acc_overall = NA_real_, ari_overall = NA_real_)
  if (is.null(frob)) {
    frob <- rep(NA_real_, 4L)
    names(frob) <- paste0("fro_s", 1:4)
  }
  if (is.null(opb)) {
    opb <- rep(NA_real_, 4L)
    names(opb) <- paste0("op_s", 1:4)
  }
  lambda_pkg <- as.numeric(lambda_pkg)
  lambda_text <- as.numeric(lambda_text)
  if (length(lambda_pkg) == 1L) lambda_pkg <- rep(lambda_pkg, 4L)
  if (length(lambda_text) == 1L) lambda_text <- rep(lambda_text, 4L)

  data.frame(
    path = as.character(path),
    q = as.integer(q),
    type = as.character(link_type),
    estimator = as.character(estimator),
    T = as.integer(T),
    rep_id = as.integer(rep_id),
    seed = as.integer(seed),
    design_id = as.integer(design_id),
    fro_s1 = as.numeric(frob[1L]),
    fro_s2 = as.numeric(frob[2L]),
    fro_s3 = as.numeric(frob[3L]),
    fro_s4 = as.numeric(frob[4L]),
    op_s1 = as.numeric(opb[1L]),
    op_s2 = as.numeric(opb[2L]),
    op_s3 = as.numeric(opb[3L]),
    op_s4 = as.numeric(opb[4L]),
    Fro = pvar2_mean(frob),
    Spec = pvar2_mean(opb),
    acc_s1 = as.numeric(score$acc[1L]),
    acc_s2 = as.numeric(score$acc[2L]),
    acc_s3 = as.numeric(score$acc[3L]),
    acc_s4 = as.numeric(score$acc[4L]),
    ari_s1 = as.numeric(score$ari[1L]),
    ari_s2 = as.numeric(score$ari[2L]),
    ari_s3 = as.numeric(score$ari[3L]),
    ari_s4 = as.numeric(score$ari[4L]),
    acc_overall = as.numeric(score$acc_overall),
    ari_overall = as.numeric(score$ari_overall),
    lambda_pkg_s1 = as.numeric(lambda_pkg[1L]),
    lambda_pkg_s2 = as.numeric(lambda_pkg[2L]),
    lambda_pkg_s3 = as.numeric(lambda_pkg[3L]),
    lambda_pkg_s4 = as.numeric(lambda_pkg[4L]),
    lambda_text_s1 = as.numeric(lambda_text[1L]),
    lambda_text_s2 = as.numeric(lambda_text[2L]),
    lambda_text_s3 = as.numeric(lambda_text[3L]),
    lambda_text_s4 = as.numeric(lambda_text[4L]),
    lambda_pkg_mean = pvar2_mean(lambda_pkg),
    lambda_text_mean = pvar2_mean(lambda_text),
    c_lambda_sel = as.numeric(c_lambda_sel),
    lambda_mode = as.character(lambda_mode),
    error_message = as.character(error_message),
    stringsAsFactors = FALSE
  )
}

pvar2_pick_c_lambda <- function(c_map, TT, fallback = 0.20) {
  if (is.null(c_map)) return(as.numeric(fallback)[1L])
  if (is.list(c_map) && !is.null(c_map[[as.character(TT)]])) {
    return(as.numeric(c_map[[as.character(TT)]])[1L])
  }
  if (!is.list(c_map)) {
    if (!is.null(names(c_map)) && as.character(TT) %in% names(c_map)) {
      return(as.numeric(c_map[as.character(TT)])[1L])
    }
    if (length(c_map) == 1L) return(as.numeric(c_map)[1L])
  }
  as.numeric(fallback)[1L]
}

pvar2_design_key <- function(q, path, link_type) {
  q <- as.integer(q)
  pid <- pvar2_path_id(path)
  tid <- if (pvar2_link_type(link_type) == "type1") 1L else 2L
  as.integer(q * 1000L + pid * 10L + tid)
}


pvar2_run_one_rep <- function(rep_id,
                              fixed_gen,
                              q,
                              T,
                              path,
                              link_type,
                              method,
                              seed_data,
                              design_id = NA_integer_,
                              sigma_diag = 0.50,
                              burn = 500L,
                              lasso_lambda_mode = "cv_c",
                              lasso_c_lambda = 0.20,
                              lasso_c_grid = pvar2_default_c_grid(),
                              lasso_fold = 10L,
                              lasso_nlambda = 50L,
                              lasso_max_iter = 1000L,
                              lasso_tol = 1e-6,
                              c_map = NULL) {
  method <- pvar2_method(method)
  path_name <- pvar2_path_name(path)
  type_name <- pvar2_link_type(link_type)

  Ylong <- pvar2_simulate(
    TT = max(T),
    Phi_list = fixed_gen$Phi,
    sigma_diag = sigma_diag,
    burn = burn,
    seed_data = seed_data
  )
  n_comm <- pvar2_path_spec(path_name)

  rows <- lapply(T, function(TT) {
    Ysub <- Ylong[, seq_len(TT), drop = FALSE]
    fit <- tryCatch({
      if (method == "ols") {
        pvar2_fit_ols(Ysub, s = 4L)
      } else {
        lambda_mode_use <- if (lasso_lambda_mode == "cv_c") "theory" else lasso_lambda_mode
        c_use <- pvar2_pick_c_lambda(c_map, TT, fallback = lasso_c_lambda)
        pvar2_fit_lasso(
          Yt = Ysub,
          s = 4L,
          lambda_mode = lambda_mode_use,
          c_lambda = c_use,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          nlambda = lasso_nlambda,
          diagTF = TRUE,
          centerTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE,
          max_iter = lasso_max_iter,
          tol = lasso_tol
        )
      }
    }, error = function(e) e)

    if (inherits(fit, "error") || inherits(fit, "simpleError")) {
      pvar2_make_result_row(
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        c_lambda_sel = if (method == "lasso") pvar2_pick_c_lambda(c_map, TT, fallback = lasso_c_lambda) else NA_real_,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        error_message = conditionMessage(fit)
      )
    } else {
      blocks_hat <- pvar2_split_blocks(fit$Phi_hat, s = 4L)
      clfit <- pvar2_cluster(blocks_hat, n_comm = n_comm)
      score <- pvar2_score(true_y = fixed_gen$y, fit_obj = clfit)
      frob <- pvar2_fro_errors(fixed_gen$Phi, blocks_hat)
      opb <- pvar2_op_errors(fixed_gen$Phi, blocks_hat)
      pvar2_make_result_row(
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        score = score,
        frob = frob,
        opb = opb,
        c_lambda_sel = if (method == "lasso") fit$c_lambda_sel else NA_real_,
        lambda_pkg = fit$lambda_pkg,
        lambda_text = fit$lambda_text,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        error_message = ""
      )
    }
  })

  do.call(rbind, rows)
}

run_dosim_pvar2 <- function(path,
                            method = c("ols", "lasso"),
                            q = 36L,
                            T = c(200L, 500L, 1000L, 2000L),
                            nrep = 100L,
                            ncore = 25L,
                            link_type = "type1",
                            target_cycle_sv = 0.90,
                            sigma_diag = 0.50,
                            diag_add = 0.5,
                            burn = 500L,
                            lasso_lambda_mode = c("theory", "cv", "cv_c"),
                            lasso_c_lambda = 0.20,
                            lasso_c_grid = pvar2_default_c_grid(),
                            lasso_fold = 10L,
                            lasso_nlambda = 50L,
                            lasso_max_iter = 1000L,
                            lasso_tol = 1e-6,
                            base_seed = 20260328L,
                            libfile = "pvar2-sparsefixed.R") {
  method <- pvar2_method(method)
  path_name <- pvar2_path_name(path)
  type_name <- pvar2_link_type(link_type)
  lasso_lambda_mode <- pvar2_lambda_mode(lasso_lambda_mode)
  q <- as.integer(q)
  T <- sort(unique(as.integer(T)))
  nrep <- as.integer(nrep)
  ncore <- as.integer(ncore)
  design_key <- pvar2_design_key(q = q, path = path_name, link_type = type_name)

  phi_seed <- pvar2_make_seed(base_seed, 900000L, design_key)
  set.seed(phi_seed)
  fixed_gen <- pvar2_generate_benchmark_pvar(
    q = q,
    path = path_name,
    link_type = type_name,
    target_cycle_sv = target_cycle_sv
  )

  c_map <- lasso_c_lambda
  cv_c_info <- NULL
  if (method == "lasso" && identical(lasso_lambda_mode, "cv_c")) {
    pilot_seed_data <- pvar2_make_seed(base_seed, 700000L, design_key)
    Ypilot_long <- pvar2_simulate(
      TT = max(T),
      Phi_list = fixed_gen$Phi,
      sigma_diag = sigma_diag,
      burn = burn,
      seed_data = pilot_seed_data
    )
    c_map <- setNames(rep(NA_real_, length(T)), as.character(T))
    cv_c_info <- vector("list", length(T))
    names(cv_c_info) <- as.character(T)
    for (TT in T) {
      sel <- tryCatch(
        pvar2_select_c_lambda_blockcv(
          Yt = Ypilot_long[, seq_len(TT), drop = FALSE],
          s = 4L,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          diagTF = TRUE,
          centerTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE,
          max_iter = lasso_max_iter,
          tol = lasso_tol
        ),
        error = function(e) NULL
      )
      if (is.null(sel) || !is.finite(sel$c_lambda_sel)) {
        c_map[as.character(TT)] <- as.numeric(lasso_c_lambda)[1L]
        cv_c_info[[as.character(TT)]] <- list(error = TRUE, c_lambda_sel = as.numeric(lasso_c_lambda)[1L])
      } else {
        c_map[as.character(TT)] <- sel$c_lambda_sel
        cv_c_info[[as.character(TT)]] <- sel
      }
    }
  }

  cl <- parallel::makeCluster(ncore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  libfile <- as.character(libfile)[1L]
  raw_list <- foreach::foreach(
    rep_id = seq_len(nrep),
    .packages = c("MASS", "gtools", "sparseVAR")
  ) %dopar% {
    source(libfile, local = TRUE)
    seed_i <- pvar2_make_seed(base_seed, 100000L * design_key, rep_id)
    pvar2_run_one_rep(
      rep_id = rep_id,
      fixed_gen = fixed_gen,
      q = q,
      T = T,
      path = path_name,
      link_type = type_name,
      method = method,
      seed_data = seed_i,
      design_id = design_key,
      sigma_diag = sigma_diag,
      burn = burn,
      lasso_lambda_mode = lasso_lambda_mode,
      lasso_c_lambda = lasso_c_lambda,
      lasso_c_grid = lasso_c_grid,
      lasso_fold = lasso_fold,
      lasso_nlambda = lasso_nlambda,
      lasso_max_iter = lasso_max_iter,
      lasso_tol = lasso_tol,
      c_map = c_map
    )
  }

  raw <- do.call(rbind, raw_list)
  raw <- raw[order(raw$rep_id, raw$T), , drop = FALSE]
  rownames(raw) <- NULL

  out <- list(
    method = method,
    estimator = method,
    path = path_name,
    q = q,
    type = type_name,
    T = T,
    nrep = nrep,
    ncore = ncore,
    target_cycle_sv = target_cycle_sv,
    sigma_diag = sigma_diag,
    burn = burn,
    same_phi_across_rep = TRUE,
    same_phi_across_T = TRUE,
    fixed_gen = fixed_gen,
    lasso_lambda_mode = if (method == "lasso") lasso_lambda_mode else NA_character_,
    lasso_c_lambda = if (method == "lasso") as.numeric(lasso_c_lambda) else NA_real_,
    lasso_c_grid = if (method == "lasso") as.numeric(lasso_c_grid) else NA_real_,
    c_lambda_map = if (method == "lasso") c_map else NULL,
    cv_c_info = cv_c_info,
    raw = raw
  )
  class(out) <- c("pvar_dosim1", class(out))
  out
}


# ------------------------------
# Summary tables
# ------------------------------
table_summary_pvar2 <- function(...) {
  xs <- list(...)
  if (length(xs) == 1L && is.list(xs[[1L]]) && !inherits(xs[[1L]], "pvar_dosim1")) {
    xs <- xs[[1L]]
  }
  if (!length(xs)) stop("Provide at least one pvar_dosim1 object.", call. = FALSE)
  ok <- vapply(xs, function(z) inherits(z, "pvar_dosim1"), logical(1L))
  if (!all(ok)) stop("All inputs must be pvar_dosim1 objects.", call. = FALSE)

  raw <- do.call(rbind, lapply(xs, function(z) z$raw))
  raw_ok <- raw[is.na(raw$error_message) | raw$error_message == "", , drop = FALSE]
  if (!nrow(raw_ok)) return(data.frame())

  agg <- aggregate(
    cbind(Fro, Spec, acc_overall, ari_overall) ~ path + q + type + estimator + T,
    data = raw_ok,
    FUN = function(z) mean(z, na.rm = TRUE)
  )

  keys <- unique(agg[, c("path", "q", "type", "estimator"), drop = FALSE])
  keys <- keys[order(keys$path, keys$q, keys$type, keys$estimator), , drop = FALSE]
  out <- keys
  key_out <- paste(out$path, out$q, out$type, out$estimator, sep = "__")

  for (TT in sort(unique(agg$T))) {
    sub <- agg[agg$T == TT, , drop = FALSE]
    key_sub <- paste(sub$path, sub$q, sub$type, sub$estimator, sep = "__")
    out[[paste0("Fro_T", TT)]] <- as.numeric(setNames(sub$Fro, key_sub)[key_out])
    out[[paste0("Spec_T", TT)]] <- as.numeric(setNames(sub$Spec, key_sub)[key_out])
    out[[paste0("acc_T", TT)]] <- as.numeric(setNames(sub$acc_overall, key_sub)[key_out])
    out[[paste0("ari_T", TT)]] <- as.numeric(setNames(sub$ari_overall, key_sub)[key_out])
  }

  rownames(out) <- NULL
  pvar2_round_df(out, digits = 3L)
}



pvar2_alpha_max <- function() 1 / (4 * sqrt(2) + 2)

pvar2_default_alpha_grid <- function() {
  amax <- pvar2_alpha_max()
  c(0, exp(seq(log(0.01 * amax), log(amax), length.out = 20L)))
}

# ------------------------------
# PisCES helpers
# ------------------------------
pvar2_labels_to_membership <- function(labels, K) {
  labels <- as.integer(labels)
  q <- length(labels)
  out <- matrix(0, nrow = q, ncol = as.integer(K))
  for (kk in seq_len(as.integer(K))) out[labels == kk, kk] <- 1
  out
}

pvar2_embeddings_with_projectors <- function(blocks, n_comm) {
  s <- length(blocks)
  emb_L <- vector("list", s)
  emb_R <- vector("list", s)
  proj_L <- vector("list", s)
  proj_R <- vector("list", s)

  for (mm in seq_len(s)) {
    Ky <- n_comm[2L * mm - 1L]
    Kz <- n_comm[2L * mm]
    Ahat <- t(blocks[[mm]])
    sv <- svd(Ahat, nu = Ky, nv = Kz)
    XL <- sv$u[, seq_len(Ky), drop = FALSE]
    XR <- sv$v[, seq_len(Kz), drop = FALSE]
    emb_L[[mm]] <- pvar2_safe_row_normalize(XL)
    emb_R[[mm]] <- pvar2_safe_row_normalize(XR)
    proj_L[[mm]] <- XL %*% t(XL)
    proj_R[[mm]] <- XR %*% t(XR)
  }

  list(emb_L = emb_L, emb_R = emb_R, proj_L = proj_L, proj_R = proj_R)
}

pvar2_projector_topK <- function(M, K) {
  ee <- eigen(M, symmetric = TRUE)
  U <- ee$vectors[, seq_len(as.integer(K)), drop = FALSE]
  U %*% t(U)
}

pvar2_pisces_smooth <- function(proj_list, ranks, alpha,
                                tol = 1e-5, max_iter = 1000L) {
  s <- length(proj_list)
  cur <- proj_list
  nxt <- vector("list", s)
  diff <- Inf
  iter <- 0L

  while (iter < as.integer(max_iter) && diff > tol) {
    iter <- iter + 1L
    diff <- 0
    for (mm in seq_len(s)) {
      M <- proj_list[[mm]]
      if (mm == 1L) {
        S <- M + alpha * cur[[2L]]
      } else if (mm == s) {
        S <- alpha * cur[[s - 1L]] + M
      } else {
        S <- alpha * cur[[mm - 1L]] + M + alpha * cur[[mm + 1L]]
      }
      nxt[[mm]] <- pvar2_projector_topK(S, ranks[mm])
      diff <- diff + sqrt(sum((nxt[[mm]] - cur[[mm]])^2))
    }
    cur <- nxt
  }
  list(projectors = cur, iter = iter)
}

pvar2_embeddings_from_projectors <- function(proj_list, ranks) {
  s <- length(proj_list)
  emb <- vector("list", s)
  for (mm in seq_len(s)) {
    ee <- eigen(proj_list[[mm]], symmetric = TRUE)
    U <- ee$vectors[, seq_len(as.integer(ranks[mm])), drop = FALSE]
    emb[[mm]] <- pvar2_safe_row_normalize(U)
  }
  emb
}

pvar2_cluster_pisces <- function(blocks, n_comm, alpha = 0) {
  s <- length(blocks)
  tmp <- pvar2_embeddings_with_projectors(blocks, n_comm)

  if (alpha > 0) {
    ranks_L <- n_comm[seq(1L, 2L * s, by = 2L)]
    ranks_R <- n_comm[seq(2L, 2L * s, by = 2L)]
    smL <- pvar2_pisces_smooth(tmp$proj_L, ranks = ranks_L, alpha = alpha)
    smR <- pvar2_pisces_smooth(tmp$proj_R, ranks = ranks_R, alpha = alpha)
    emb_L <- pvar2_embeddings_from_projectors(smL$projectors, ranks_L)
    emb_R <- pvar2_embeddings_from_projectors(smR$projectors, ranks_R)
  } else {
    emb_L <- tmp$emb_L
    emb_R <- tmp$emb_R
  }

  group_L <- vector("list", s)
  group_R <- vector("list", s)

  Ky1 <- n_comm[1L]
  Kzs <- n_comm[2L * s]
  if (Ky1 == Kzs) {
    lab <- pvar2_kmeans_labels(cbind(emb_R[[s]], emb_L[[1L]]), Ky1)
    group_R[[s]] <- lab
    group_L[[1L]] <- lab
  } else {
    group_R[[s]] <- pvar2_kmeans_labels(emb_R[[s]], Kzs)
    group_L[[1L]] <- pvar2_kmeans_labels(emb_L[[1L]], Ky1)
  }

  if (s >= 2L) {
    for (mm in 2:s) {
      Ky <- n_comm[2L * mm - 1L]
      Kzprev <- n_comm[2L * (mm - 1L)]
      if (Ky == Kzprev) {
        lab <- pvar2_kmeans_labels(cbind(emb_R[[mm - 1L]], emb_L[[mm]]), Ky)
        group_R[[mm - 1L]] <- lab
        group_L[[mm]] <- lab
      } else {
        group_R[[mm - 1L]] <- pvar2_kmeans_labels(emb_R[[mm - 1L]], Kzprev)
        group_L[[mm]] <- pvar2_kmeans_labels(emb_L[[mm]], Ky)
      }
    }
  }

  list(group_L = group_L, group_R = group_R, emb_L = emb_L, emb_R = emb_R, alpha = alpha)
}

pvar2_completed_from_mask <- function(M, keep_mask, rankK) {
  M0 <- M * keep_mask
  sv <- svd(M0, nu = as.integer(rankK), nv = as.integer(rankK))
  sv$u[, seq_len(as.integer(rankK)), drop = FALSE] %*%
    diag(sv$d[seq_len(as.integer(rankK))], as.integer(rankK)) %*%
    t(sv$v[, seq_len(as.integer(rankK)), drop = FALSE])
}

pvar2_predict_from_labels <- function(M, yhat, zhat) {
  Mpos <- pmax(M, 0)
  q <- nrow(Mpos)
  Ky <- max(yhat)
  Kz <- max(zhat)
  dy <- pmax(rowSums(Mpos), 1e-8)
  dz <- pmax(colSums(Mpos), 1e-8)
  Y <- pvar2_labels_to_membership(yhat, Ky)
  Z <- pvar2_labels_to_membership(zhat, Kz)
  Bhat <- matrix(0, nrow = Ky, ncol = Kz)
  for (ky in seq_len(Ky)) {
    idx_y <- which(yhat == ky)
    for (kz in seq_len(Kz)) {
      idx_z <- which(zhat == kz)
      num <- sum(Mpos[idx_y, idx_z, drop = FALSE])
      den <- sum(outer(dy[idx_y], dz[idx_z]))
      Bhat[ky, kz] <- if (den > 0) num / den else 0
    }
  }
  diag(dy, q) %*% Y %*% Bhat %*% t(Z) %*% diag(dz, q)
}

pvar2_select_alpha <- function(blocks, n_comm,
                               folds = 5L,
                               alpha_grid = pvar2_default_alpha_grid(),
                               criterion = c("holdout", "paper"),
                               seed = NULL) {
  criterion <- match.arg(criterion)
  if (!is.null(seed)) set.seed(as.integer(seed))
  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  alpha_grid <- alpha_grid[is.finite(alpha_grid) & alpha_grid >= 0]
  if (!length(alpha_grid)) stop("alpha_grid must contain nonnegative finite values.", call. = FALSE)

  s <- length(blocks)
  q <- nrow(blocks[[1L]])
  offdiag <- which(row(matrix(0, q, q)) != col(matrix(0, q, q)))
  folds <- as.integer(folds)
  folds <- max(2L, folds)

  heldout_masks <- vector("list", folds)
  comp_blocks <- vector("list", folds)

  for (ff in seq_len(folds)) {
    heldout_masks[[ff]] <- vector("list", s)
    comp_blocks[[ff]] <- vector("list", s)
    for (mm in seq_len(s)) {
      rankK <- min(n_comm[2L * mm - 1L], n_comm[2L * mm])
      hold_idx <- sample(offdiag, size = floor(length(offdiag) / 2), replace = FALSE)
      keep_mask <- matrix(1, q, q)
      keep_mask[hold_idx] <- 0
      heldout_masks[[ff]][[mm]] <- keep_mask
      comp_blocks[[ff]][[mm]] <- pvar2_completed_from_mask(blocks[[mm]], keep_mask, rankK)
    }
  }

  score_alpha <- numeric(length(alpha_grid))
  for (aa in seq_along(alpha_grid)) {
    alpha <- alpha_grid[aa]
    score_fold <- numeric(folds)
    for (ff in seq_len(folds)) {
      fit <- pvar2_cluster_pisces(comp_blocks[[ff]], n_comm = n_comm, alpha = alpha)
      sc <- 0
      for (mm in seq_len(s)) {
        Mcomp <- comp_blocks[[ff]][[mm]]
        Mhat <- pvar2_predict_from_labels(Mcomp, fit$group_L[[mm]], fit$group_R[[mm]])
        hold_mask <- 1 - heldout_masks[[ff]][[mm]]
        if (criterion == "holdout") {
          diff <- (blocks[[mm]] - Mhat) * hold_mask
          sc <- sc + mean(diff[hold_mask == 1]^2)
        } else {
          trM <- sum(diag(pmax(Mcomp, 0))) / q
          trH <- sum(diag(Mhat)) / q
          sc <- sc + trM * (1 - trH)
        }
      }
      score_fold[ff] <- sc
    }
    score_alpha[aa] <- mean(score_fold)
  }

  best <- which.min(score_alpha)
  list(alpha = alpha_grid[best], alpha_grid = alpha_grid, score = score_alpha)
}

# ------------------------------
# Row wrappers with smoothing info
# ------------------------------
pvar2_make_result_row_pisces <- function(..., smoothing = c("none", "pisces"), alpha_sel = NA_real_) {
  smoothing <- match.arg(smoothing)
  out <- pvar2_make_result_row(...)
  out$smoothing <- smoothing
  out$alpha_sel <- as.numeric(alpha_sel)
  out
}

# ------------------------------
# One-rep runner with none vs PisCES comparison
# ------------------------------
pvar2_run_one_rep_pisces <- function(rep_id,
                                     fixed_gen,
                                     q,
                                     T,
                                     path,
                                     link_type,
                                     method,
                                     seed_data,
                                     design_id = NA_integer_,
                                     sigma_diag = 0.50,
                                     burn = 500L,
                                     lasso_lambda_mode = "cv_c",
                                     lasso_c_lambda = 0.20,
                                     lasso_c_grid = pvar2_default_c_grid(),
                                     lasso_fold = 10L,
                                     lasso_nlambda = 50L,
                                     lasso_max_iter = 1000L,
                                     lasso_tol = 1e-6,
                                     c_map = NULL,
                                     alpha_map = NULL) {
  method <- pvar2_method(method)
  path_name <- pvar2_path_name(path)
  type_name <- pvar2_link_type(link_type)

  Ylong <- pvar2_simulate(
    TT = max(T),
    Phi_list = fixed_gen$Phi,
    sigma_diag = sigma_diag,
    burn = burn,
    seed_data = seed_data
  )
  n_comm <- pvar2_path_spec(path_name)

  rows <- lapply(T, function(TT) {
    Ysub <- Ylong[, seq_len(TT), drop = FALSE]
    fit <- tryCatch({
      if (method == "ols") {
        pvar2_fit_ols(Ysub, s = 4L)
      } else {
        lambda_mode_use <- if (lasso_lambda_mode == "cv_c") "theory" else lasso_lambda_mode
        c_use <- pvar2_pick_c_lambda(c_map, TT, fallback = lasso_c_lambda)
        pvar2_fit_lasso(
          Yt = Ysub,
          s = 4L,
          lambda_mode = lambda_mode_use,
          c_lambda = c_use,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          nlambda = lasso_nlambda,
          diagTF = TRUE,
          centerTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE,
          max_iter = lasso_max_iter,
          tol = lasso_tol
        )
      }
    }, error = function(e) e)

    alpha_use <- pvar2_pick_c_lambda(alpha_map, TT, fallback = 0)

    if (inherits(fit, "error") || inherits(fit, "simpleError")) {
      r1 <- pvar2_make_result_row_pisces(
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        c_lambda_sel = if (method == "lasso") pvar2_pick_c_lambda(c_map, TT, fallback = lasso_c_lambda) else NA_real_,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        error_message = conditionMessage(fit),
        smoothing = "none",
        alpha_sel = 0
      )
      r2 <- pvar2_make_result_row_pisces(
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        c_lambda_sel = if (method == "lasso") pvar2_pick_c_lambda(c_map, TT, fallback = lasso_c_lambda) else NA_real_,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        error_message = conditionMessage(fit),
        smoothing = "pisces",
        alpha_sel = alpha_use
      )
      return(rbind(r1, r2))
    }

    blocks_hat <- pvar2_split_blocks(fit$Phi_hat, s = 4L)

    clfit_none <- pvar2_cluster_pisces(blocks_hat, n_comm = n_comm, alpha = 0)
    score_none <- pvar2_score(true_y = fixed_gen$y, fit_obj = clfit_none)
    frob <- pvar2_fro_errors(fixed_gen$Phi, blocks_hat)
    opb <- pvar2_op_errors(fixed_gen$Phi, blocks_hat)

    row_none <- pvar2_make_result_row_pisces(
      path = path_name,
      q = q,
      link_type = type_name,
      estimator = method,
      T = TT,
      rep_id = rep_id,
      seed = seed_data,
      design_id = design_id,
      score = score_none,
      frob = frob,
      opb = opb,
      c_lambda_sel = if (method == "lasso") fit$c_lambda_sel else NA_real_,
      lambda_pkg = fit$lambda_pkg,
      lambda_text = fit$lambda_text,
      lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
      error_message = "",
      smoothing = "none",
      alpha_sel = 0
    )

    clfit_pisces <- pvar2_cluster_pisces(blocks_hat, n_comm = n_comm, alpha = alpha_use)
    score_pisces <- pvar2_score(true_y = fixed_gen$y, fit_obj = clfit_pisces)
    row_pisces <- pvar2_make_result_row_pisces(
      path = path_name,
      q = q,
      link_type = type_name,
      estimator = method,
      T = TT,
      rep_id = rep_id,
      seed = seed_data,
      design_id = design_id,
      score = score_pisces,
      frob = frob,
      opb = opb,
      c_lambda_sel = if (method == "lasso") fit$c_lambda_sel else NA_real_,
      lambda_pkg = fit$lambda_pkg,
      lambda_text = fit$lambda_text,
      lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
      error_message = "",
      smoothing = "pisces",
      alpha_sel = alpha_use
    )

    rbind(row_none, row_pisces)
  })

  do.call(rbind, rows)
}

# ------------------------------
# Main runner with alpha CV on pilot data
# ------------------------------
run_dosim_pvar2_pisces <- function(path,
                                   method = c("ols", "lasso"),
                                   q = 36L,
                                   T = c(200L, 500L, 1000L, 2000L),
                                   nrep = 100L,
                                   ncore = 25L,
                                   link_type = "type1",
                                   target_cycle_sv = 0.90,
                                   sigma_diag = 0.50,
                                   diag_add = 0.5,
                                   burn = 500L,
                                   lasso_lambda_mode = c("theory", "cv", "cv_c"),
                                   lasso_c_lambda = 0.20,
                                   lasso_c_grid = pvar2_default_c_grid(),
                                   lasso_fold = 10L,
                                   lasso_nlambda = 50L,
                                   lasso_max_iter = 1000L,
                                   lasso_tol = 1e-6,
                                   alpha_fold = 5L,
                                   alpha_grid = pvar2_default_alpha_grid(),
                                   alpha_criterion = c("holdout", "paper"),
                                   base_seed = 20260328L,
                                   libfile = "pvar2-sparsefixed.R") {
  method <- pvar2_method(method)
  path_name <- pvar2_path_name(path)
  type_name <- pvar2_link_type(link_type)
  lasso_lambda_mode <- pvar2_lambda_mode(lasso_lambda_mode)
  alpha_criterion <- match.arg(alpha_criterion)
  q <- as.integer(q)
  T <- sort(unique(as.integer(T)))
  nrep <- as.integer(nrep)
  ncore <- as.integer(ncore)
  design_key <- pvar2_design_key(q = q, path = path_name, link_type = type_name)

  phi_seed <- pvar2_make_seed(base_seed, 900000L, design_key)
  set.seed(phi_seed)
  fixed_gen <- pvar2_generate_benchmark_pvar(
    q = q,
    path = path_name,
    link_type = type_name,
    target_cycle_sv = target_cycle_sv
  )

  c_map <- lasso_c_lambda
  cv_c_info <- NULL
  if (method == "lasso" && identical(lasso_lambda_mode, "cv_c")) {
    pilot_seed_data <- pvar2_make_seed(base_seed, 700000L, design_key)
    Ypilot_long <- pvar2_simulate(
      TT = max(T),
      Phi_list = fixed_gen$Phi,
      sigma_diag = sigma_diag,
      burn = burn,
      seed_data = pilot_seed_data
    )
    c_map <- setNames(rep(NA_real_, length(T)), as.character(T))
    cv_c_info <- vector("list", length(T))
    names(cv_c_info) <- as.character(T)
    for (TT in T) {
      sel <- tryCatch(
        pvar2_select_c_lambda_blockcv(
          Yt = Ypilot_long[, seq_len(TT), drop = FALSE],
          s = 4L,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          diagTF = TRUE,
          centerTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE,
          max_iter = lasso_max_iter,
          tol = lasso_tol
        ),
        error = function(e) NULL
      )
      if (is.null(sel) || !is.finite(sel$c_lambda_sel)) {
        c_map[as.character(TT)] <- as.numeric(lasso_c_lambda)[1L]
        cv_c_info[[as.character(TT)]] <- list(error = TRUE, c_lambda_sel = as.numeric(lasso_c_lambda)[1L])
      } else {
        c_map[as.character(TT)] <- sel$c_lambda_sel
        cv_c_info[[as.character(TT)]] <- sel
      }
    }
  }

  # alpha selection on pilot data, using the same method and c_lambda map
  pilot_seed_data2 <- pvar2_make_seed(base_seed, 800000L, design_key)
  Ypilot_long2 <- pvar2_simulate(
    TT = max(T),
    Phi_list = fixed_gen$Phi,
    sigma_diag = sigma_diag,
    burn = burn,
    seed_data = pilot_seed_data2
  )
  alpha_map <- setNames(rep(0, length(T)), as.character(T))
  alpha_cv_info <- vector("list", length(T))
  names(alpha_cv_info) <- as.character(T)
  n_comm <- pvar2_path_spec(path_name)

  for (TT in T) {
    Ypilot <- Ypilot_long2[, seq_len(TT), drop = FALSE]
    fit_pilot <- tryCatch({
      if (method == "ols") {
        pvar2_fit_ols(Ypilot, s = 4L)
      } else {
        lambda_mode_use <- if (lasso_lambda_mode == "cv_c") "theory" else lasso_lambda_mode
        c_use <- pvar2_pick_c_lambda(c_map, TT, fallback = lasso_c_lambda)
        pvar2_fit_lasso(
          Yt = Ypilot,
          s = 4L,
          lambda_mode = lambda_mode_use,
          c_lambda = c_use,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          nlambda = lasso_nlambda,
          diagTF = TRUE,
          centerTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE,
          max_iter = lasso_max_iter,
          tol = lasso_tol
        )
      }
    }, error = function(e) NULL)

    if (is.null(fit_pilot)) {
      alpha_map[as.character(TT)] <- 0
      alpha_cv_info[[as.character(TT)]] <- list(error = TRUE, alpha = 0)
    } else {
      blocks_hat <- pvar2_split_blocks(fit_pilot$Phi_hat, s = 4L)
      sel_alpha <- tryCatch(
        pvar2_select_alpha(
          blocks = blocks_hat,
          n_comm = n_comm,
          folds = alpha_fold,
          alpha_grid = alpha_grid,
          criterion = alpha_criterion,
          seed = pvar2_make_seed(base_seed, 810000L, design_key, TT)
        ),
        error = function(e) NULL
      )
      if (is.null(sel_alpha) || !is.finite(sel_alpha$alpha)) {
        alpha_map[as.character(TT)] <- 0
        alpha_cv_info[[as.character(TT)]] <- list(error = TRUE, alpha = 0)
      } else {
        alpha_map[as.character(TT)] <- sel_alpha$alpha
        alpha_cv_info[[as.character(TT)]] <- sel_alpha
      }
    }
  }

  cl <- parallel::makeCluster(ncore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  libfile <- as.character(libfile)[1L]
  raw_list <- foreach::foreach(
    rep_id = seq_len(nrep),
    .packages = c("MASS", "gtools", "sparseVAR")
  ) %dopar% {
    source(libfile, local = TRUE)
    seed_i <- pvar2_make_seed(base_seed, 100000L * design_key, rep_id)
    pvar2_run_one_rep_pisces(
      rep_id = rep_id,
      fixed_gen = fixed_gen,
      q = q,
      T = T,
      path = path_name,
      link_type = type_name,
      method = method,
      seed_data = seed_i,
      design_id = design_key,
      sigma_diag = sigma_diag,
      burn = burn,
      lasso_lambda_mode = lasso_lambda_mode,
      lasso_c_lambda = lasso_c_lambda,
      lasso_c_grid = lasso_c_grid,
      lasso_fold = lasso_fold,
      lasso_nlambda = lasso_nlambda,
      lasso_max_iter = lasso_max_iter,
      lasso_tol = lasso_tol,
      c_map = c_map,
      alpha_map = alpha_map
    )
  }

  raw <- do.call(rbind, raw_list)
  raw <- raw[order(raw$rep_id, raw$T, raw$smoothing), , drop = FALSE]
  rownames(raw) <- NULL

  out <- list(
    method = method,
    estimator = method,
    path = path_name,
    q = q,
    type = type_name,
    T = T,
    nrep = nrep,
    ncore = ncore,
    target_cycle_sv = target_cycle_sv,
    sigma_diag = sigma_diag,
    burn = burn,
    same_phi_across_rep = TRUE,
    same_phi_across_T = TRUE,
    fixed_gen = fixed_gen,
    lasso_lambda_mode = if (method == "lasso") lasso_lambda_mode else NA_character_,
    lasso_c_lambda = if (method == "lasso") as.numeric(lasso_c_lambda) else NA_real_,
    lasso_c_grid = if (method == "lasso") as.numeric(lasso_c_grid) else NA_real_,
    c_lambda_map = if (method == "lasso") c_map else NULL,
    cv_c_info = cv_c_info,
    alpha_map = alpha_map,
    alpha_cv_info = alpha_cv_info,
    raw = raw
  )
  class(out) <- c("pvar_dosim2_pisces", class(out))
  out
}

# ------------------------------
# Summary tables with smoothing
# ------------------------------
table_summary_pvar2_pisces <- function(...) {
  xs <- list(...)
  if (length(xs) == 1L && is.list(xs[[1L]]) && !inherits(xs[[1L]], "pvar_dosim2_pisces")) {
    xs <- xs[[1L]]
  }
  if (!length(xs)) stop("Provide at least one pvar_dosim2_pisces object.", call. = FALSE)
  ok <- vapply(xs, function(z) inherits(z, "pvar_dosim2_pisces"), logical(1L))
  if (!all(ok)) stop("All inputs must be pvar_dosim2_pisces objects.", call. = FALSE)

  raw <- do.call(rbind, lapply(xs, function(z) z$raw))
  raw_ok <- raw[is.na(raw$error_message) | raw$error_message == "", , drop = FALSE]
  if (!nrow(raw_ok)) return(data.frame())

  agg <- aggregate(
    cbind(Fro, Spec, acc_overall, ari_overall, alpha_sel) ~ path + q + type + estimator + smoothing + T,
    data = raw_ok,
    FUN = function(z) mean(z, na.rm = TRUE)
  )

  keys <- unique(agg[, c("path", "q", "type", "estimator", "smoothing"), drop = FALSE])
  keys <- keys[order(keys$path, keys$q, keys$type, keys$estimator, keys$smoothing), , drop = FALSE]
  out <- keys
  key_out <- paste(out$path, out$q, out$type, out$estimator, out$smoothing, sep = "__")

  for (TT in sort(unique(agg$T))) {
    sub <- agg[agg$T == TT, , drop = FALSE]
    key_sub <- paste(sub$path, sub$q, sub$type, sub$estimator, sub$smoothing, sep = "__")
    out[[paste0("Fro_T", TT)]] <- as.numeric(setNames(sub$Fro, key_sub)[key_out])
    out[[paste0("Spec_T", TT)]] <- as.numeric(setNames(sub$Spec, key_sub)[key_out])
    out[[paste0("acc_T", TT)]] <- as.numeric(setNames(sub$acc_overall, key_sub)[key_out])
    out[[paste0("ari_T", TT)]] <- as.numeric(setNames(sub$ari_overall, key_sub)[key_out])
    out[[paste0("alpha_T", TT)]] <- as.numeric(setNames(sub$alpha_sel, key_sub)[key_out])
  }

  rownames(out) <- NULL
  pvar2_round_df(out, digits = 3L)
}


# backward-compatible convenient aliases
run_dosim_pvar2 <- run_dosim_pvar2_pisces
table_summary_pvar2 <- table_summary_pvar2_pisces

# ------------------------------
# GitHub-facing aliases
# ------------------------------
scbm_pvar_default_c_grid <- pvar2_default_c_grid
scbm_pvar_default_alpha_grid <- pvar2_default_alpha_grid
scbm_pvar_run_simulation <- run_dosim_pvar2_pisces
scbm_pvar_simulation_summary <- table_summary_pvar2_pisces
