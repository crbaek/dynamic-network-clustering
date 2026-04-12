# ============================================================
# scbm_vhar_simulation_library.R
# GitHub-ready VHAR simulation library for
# 'Dynamic spectral co-clustering of directed networks to unveil
# latent community paths in VAR-type models'.
#
# This file keeps the current working generalized VHAR simulation
# logic, including fixed-design generation, OLS / lasso first-stage
# estimation, and optional PisCES smoothing for community recovery.
#
# Public wrappers added at the end of the file:
#   - scbm_vhar_run_simulation()
#   - scbm_vhar_simulation_summary()
#   - scbm_vhar_default_c_grid()
#   - scbm_vhar_default_alpha_grid()
# ============================================================


suppressPackageStartupMessages({
  library(MASS)
  library(gtools)
  library(foreach)
  library(doParallel)
  library(sparseVAR)
})

# ------------------------------
# Basic helpers
# ------------------------------
vhar3_safe_matrix <- function(X) {
  if (is.null(X)) return(matrix(NA_real_, 1L, 1L))
  X <- tryCatch(as.matrix(X), error = function(e) NULL)
  if (is.null(X)) return(matrix(NA_real_, 1L, 1L))
  if (length(dim(X)) != 2L) X <- matrix(X, ncol = 1L)
  X
}

vhar3_mean <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

vhar3_fro_norm <- function(A) {
  A <- vhar3_safe_matrix(A)
  if (any(!is.finite(A))) return(NA_real_)
  sqrt(sum(A^2))
}

vhar3_operator_norm <- function(A) {
  A <- vhar3_safe_matrix(A)
  if (any(!is.finite(A))) return(NA_real_)
  d <- tryCatch(svd(A, nu = 0L, nv = 0L)$d, error = function(e) numeric(0))
  if (!length(d)) return(NA_real_)
  as.numeric(d[1L])
}

vhar3_row_normalize <- function(X) {
  X <- vhar3_safe_matrix(X)
  rn <- sqrt(rowSums(X^2))
  rn[!is.finite(rn) | rn <= .Machine$double.eps] <- 1
  X / rn
}

vhar3_round_df <- function(df, digits = 3L) {
  out <- df
  is_num <- vapply(out, is.numeric, logical(1L))
  out[is_num] <- lapply(out[is_num], function(z) round(z, digits))
  out
}

vhar3_normalize_seed <- function(seed, default = 1L) {
  mod <- 2147483647
  x <- suppressWarnings(as.numeric(seed)[1L])
  if (!is.finite(x) || is.na(x)) return(as.integer(default))
  x <- floor(abs(x)) %% mod
  if (x <= 0) x <- as.numeric(default)
  as.integer(x)
}

vhar3_make_seed <- function(base_seed, ..., default = 1L) {
  mod <- 2147483647
  x <- as.numeric(vhar3_normalize_seed(base_seed, default = default))
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

vhar3_default_c_grid <- function() seq(from = 0.1, to = 1, by = 0.05)

vhar3_alpha_max <- function() 1 / (4 * sqrt(2) + 2)

vhar3_default_alpha_grid <- function() {
  amax <- vhar3_alpha_max()
  c(0, exp(seq(log(0.01 * amax), log(amax), length.out = 20L)))
}

# ------------------------------
# Matchers
# ------------------------------
vhar3_setting_name <- function(setting_id) {
  if (is.numeric(setting_id)) {
    sid <- as.integer(setting_id)[1L]
    if (!sid %in% 1:3) stop("setting_id must be 1, 2, 3 or 'sg1', 'sg2', 'sg3'.", call. = FALSE)
    return(paste0("sg", sid))
  }
  x <- tolower(as.character(setting_id)[1L])
  if (!x %in% c("sg1", "sg2", "sg3")) {
    stop("setting_id must be one of 'sg1', 'sg2', 'sg3'.", call. = FALSE)
  }
  x
}

vhar3_path_name <- function(path) {
  if (is.numeric(path)) {
    pid <- as.integer(path)[1L]
    if (!pid %in% 1:3) stop("path must be 1, 2, 3 or 'path1', 'path2', 'path3'.", call. = FALSE)
    return(paste0("path", pid))
  }
  x <- tolower(as.character(path)[1L])
  if (!x %in% c("path1", "path2", "path3")) {
    stop("path must be one of 'path1', 'path2', 'path3'.", call. = FALSE)
  }
  x
}

vhar3_link_type <- function(link_type) {
  x <- tolower(as.character(link_type)[1L])
  if (x %in% c("type1", "1", "t1")) return("type1")
  if (x %in% c("type2", "2", "t2")) return("type2")
  stop("link_type must be one of 'type1' or 'type2'.", call. = FALSE)
}

vhar3_method <- function(method) {
  x <- tolower(as.character(method)[1L])
  if (!x %in% c("ols", "lasso")) stop("method must be 'ols' or 'lasso'.", call. = FALSE)
  x
}

vhar3_lambda_mode <- function(lambda_mode = c("theory", "cv", "cv_c")) {
  match.arg(lambda_mode)
}

# ------------------------------
# Design settings
# ------------------------------
vhar3_setting_parms <- function(setting_id = "sg1") {
  x <- vhar3_setting_name(setting_id)
  switch(
    x,
    sg1 = list(
      coeff = c(d = 0.34, w = 0.28, m = 0.24),
      sv_target = 0.90,
      pdiag = c(d = 0.95, w = 0.93, m = 0.91),
      poff = c(d = 0.02, w = 0.03, m = 0.04)
    ),
    sg2 = list(
      coeff = c(d = 0.30, w = 0.24, m = 0.20),
      sv_target = 0.84,
      pdiag = c(d = 0.84, w = 0.81, m = 0.78),
      poff = c(d = 0.06, w = 0.08, m = 0.10)
    ),
    sg3 = list(
      coeff = c(d = 0.26, w = 0.21, m = 0.17),
      sv_target = 0.78,
      pdiag = c(d = 0.72, w = 0.69, m = 0.66),
      poff = c(d = 0.11, w = 0.13, m = 0.15)
    )
  )
}

vhar3_make_B <- function(y_labels, z_labels,
                         link_type = "type1",
                         horizon = c("daily", "weekly", "monthly"),
                         setting_id = "sg1") {
  link_type <- vhar3_link_type(link_type)
  horizon <- match.arg(horizon)
  pars <- vhar3_setting_parms(setting_id)
  key <- substr(horizon, 1L, 1L)
  diag_p <- unname(pars$pdiag[key])
  off_p <- unname(pars$poff[key])
  Ky <- max(y_labels)
  Kz <- max(z_labels)

  if (link_type == "type1") {
    B <- matrix(off_p, nrow = Ky, ncol = Kz)
    dd <- min(Ky, Kz)
    B[cbind(seq_len(dd), seq_len(dd))] <- diag_p
  } else {
    upper_p <- min(0.995, off_p + 0.03)
    lower_p <- min(0.995, off_p + 0.06)
    B <- matrix(lower_p, nrow = Ky, ncol = Kz)
    for (ii in seq_len(Ky)) {
      for (jj in seq_len(Kz)) {
        if (ii < jj) B[ii, jj] <- upper_p
      }
    }
    dd <- min(Ky, Kz)
    B[cbind(seq_len(dd), seq_len(dd))] <- pmax(0.001, diag_p - 0.03)
  }

  B[B < 0.001] <- 0.001
  B[B > 0.995] <- 0.995
  B
}

vhar3_sample_graph <- function(y_labels, z_labels, B, allow_self = TRUE) {
  q <- length(y_labels)
  P <- matrix(0, nrow = q, ncol = q)
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      P[i, j] <- B[y_labels[i], z_labels[j]]
    }
  }
  A <- matrix(rbinom(q * q, 1L, as.vector(P)), nrow = q, ncol = q)
  if (!allow_self) diag(A) <- 0
  A
}

vhar3_scale_by_sender_group <- function(A, y_labels) {
  A <- vhar3_safe_matrix(A)
  grp_sizes <- as.numeric(table(factor(y_labels, levels = seq_len(max(y_labels)))))
  div <- grp_sizes[y_labels]
  div[!is.finite(div) | div <= .Machine$double.eps] <- 1
  A / div
}

# ------------------------------
# Path construction
# ------------------------------
vhar3_balanced_labels <- function(q, K) {
  q <- as.integer(q)
  K <- as.integer(K)
  sz <- rep(q %/% K, K)
  rem <- q %% K
  if (rem > 0L) sz[seq_len(rem)] <- sz[seq_len(rem)] + 1L
  rep(seq_len(K), times = sz)
}

vhar3_perm_labels <- function(labels, perm) perm[as.integer(labels)]

vhar3_path2_weekly_labels <- function(q) {
  q <- as.integer(q)
  if (q %% 6L != 0L) stop("path2 requires q divisible by 6.", call. = FALSE)
  a <- q %/% 6L
  out <- integer(q)
  out[seq_len(2L * a)] <- 1L
  out[(2L * a + 1L):(3L * a)] <- 2L
  out[(3L * a + 1L):(4L * a)] <- 2L
  out[(4L * a + 1L):(6L * a)] <- 3L
  out
}

vhar3_make_partitions <- function(q) {
  C2 <- vhar3_balanced_labels(q, 2L)
  C2_swap <- vhar3_perm_labels(C2, c(2L, 1L))
  C3 <- vhar3_balanced_labels(q, 3L)
  C3_231 <- vhar3_perm_labels(C3, c(2L, 3L, 1L))
  list(
    C2 = C2,
    C2_swap = C2_swap,
    C3 = C3,
    C3_231 = C3_231,
    P2_weekly = vhar3_path2_weekly_labels(q)
  )
}

vhar3_generate_path_labels <- function(q, path = "path1") {
  P <- vhar3_make_partitions(q)
  path_name <- vhar3_path_name(path)

  if (path_name == "path1") {
    y1 <- P$C3; z1 <- P$C3
    y2 <- z1;   z2 <- P$C3
    y3 <- z2;   z3 <- P$C3
  } else if (path_name == "path2") {
    # daily 2 -> weekly 3 -> monthly 2
    # daily c1 = 1/3 -> weekly c1, 1/6 -> weekly c2
    # daily c2 = 1/6 -> weekly c2, remaining 1/3 -> weekly c3
    # weekly -> monthly is the exact reverse, so monthly = daily partition
    y1 <- P$C2
    z1 <- P$P2_weekly
    y2 <- z1
    z2 <- P$C2
    y3 <- z2
    z3 <- P$C2
  } else {
    # keep original 2 => 2 => 3 => 3 refinement path
    y1 <- P$C2
    z1 <- P$C2_swap
    y2 <- z1
    z2 <- P$C3_231
    y3 <- z2
    z3 <- P$C3_231
  }

  list(
    y = list(as.integer(y1), as.integer(y2), as.integer(y3)),
    z = list(as.integer(z1), as.integer(z2), as.integer(z3))
  )
}

vhar3_path_spec <- function(path = "path1") {
  path_name <- vhar3_path_name(path)
  if (path_name == "path1") return(c(3L, 3L, 3L, 3L, 3L, 3L))
  if (path_name == "path2") return(c(2L, 3L, 3L, 2L, 2L, 2L))
  c(2L, 2L, 2L, 3L, 3L, 3L)
}

# ------------------------------
# Generation and scaling
# ------------------------------
vhar3_effective_snr_multiplier <- function(bw = 3L, bm = 10L) {
  c(1, sqrt(as.numeric(bw)), sqrt(as.numeric(bm)))
}

vhar3_scale_phi_to_target_sv <- function(Phi_list, sv_target = 0.90) {
  sv0 <- vhar3_operator_norm(Reduce(`+`, Phi_list))
  fac <- if (!is.finite(sv0) || sv0 <= 0) 1 else as.numeric(sv_target) / as.numeric(sv0)
  Phi2 <- lapply(Phi_list, function(P) fac * P)
  list(
    Phi = Phi2,
    scale_fac = fac,
    sv0 = sv0,
    sv = vhar3_operator_norm(Reduce(`+`, Phi2))
  )
}

vhar3_vhar_to_varcat <- function(Phi_h, bw = 3L, bm = 10L) {
  q <- nrow(Phi_h[[1L]])
  A_cat <- matrix(0, nrow = q, ncol = q * bm)
  A_cat[, 1:q] <- Phi_h[[1L]] + Phi_h[[2L]] / bw + Phi_h[[3L]] / bm
  if (bw >= 2L) {
    for (hh in 2:bw) {
      idx <- ((hh - 1L) * q + 1L):(hh * q)
      A_cat[, idx] <- Phi_h[[2L]] / bw + Phi_h[[3L]] / bm
    }
  }
  if (bm > bw) {
    for (hh in (bw + 1L):bm) {
      idx <- ((hh - 1L) * q + 1L):(hh * q)
      A_cat[, idx] <- Phi_h[[3L]] / bm
    }
  }
  A_cat
}

vhar3_companion_matrix <- function(Phi_cat, q, p) {
  F <- matrix(0, nrow = q * p, ncol = q * p)
  F[1:q, 1:(q * p)] <- Phi_cat
  if (p >= 2L) F[(q + 1L):(q * p), 1:(q * (p - 1L))] <- diag(q * (p - 1L))
  F
}

vhar3_companion_radius <- function(Phi_h, bw = 3L, bm = 10L) {
  A_cat <- vhar3_vhar_to_varcat(Phi_h, bw = bw, bm = bm)
  q <- nrow(Phi_h[[1L]])
  F <- vhar3_companion_matrix(A_cat, q, bm)
  val <- tryCatch(max(Mod(eigen(F, only.values = TRUE)$values)), error = function(e) NA_real_)
  as.numeric(val)
}

vhar3_stationary_variance_unit <- function(Phi_h, bw = 3L, bm = 10L,
                                           tol = 1e-10, max_iter = 5000L) {
  q <- nrow(Phi_h[[1L]])
  A_cat <- vhar3_vhar_to_varcat(Phi_h, bw = bw, bm = bm)
  F <- vhar3_companion_matrix(A_cat, q, bm)
  m <- nrow(F)
  G <- matrix(0, nrow = m, ncol = q)
  G[1:q, 1:q] <- diag(q)

  M <- G
  Sigma_y <- diag(q)
  contrib_norm <- vhar3_operator_norm(diag(q))
  iter <- 0L
  repeat {
    iter <- iter + 1L
    M <- F %*% M
    B <- M[1:q, , drop = FALSE]
    Sigma_add <- B %*% t(B)
    Sigma_y <- Sigma_y + Sigma_add
    contrib_norm <- vhar3_operator_norm(Sigma_add)
    if (!is.finite(contrib_norm) || contrib_norm < tol || iter >= max_iter) break
  }
  list(
    Sigma_y = Sigma_y,
    mean_var = mean(diag(Sigma_y)),
    min_var = min(diag(Sigma_y)),
    max_var = max(diag(Sigma_y)),
    iter = iter,
    last_contrib_op = contrib_norm
  )
}

vhar3_noise_scale_for_target_var <- function(Phi_h, bw = 3L, bm = 10L,
                                             target_var = 0.5,
                                             tol = 1e-10, max_iter = 5000L) {
  unit <- vhar3_stationary_variance_unit(Phi_h, bw = bw, bm = bm, tol = tol, max_iter = max_iter)
  mu <- as.numeric(unit$mean_var)
  if (!is.finite(mu) || mu <= 0) {
    sd <- 1.0
    achieved <- NA_real_
  } else {
    sd <- sqrt(as.numeric(target_var) / mu)
    achieved <- (sd ^ 2) * mu
  }
  c(unit, list(
    noise_sd = sd,
    target_var = as.numeric(target_var),
    achieved_mean_var = achieved,
    achieved_min_var = (sd ^ 2) * as.numeric(unit$min_var),
    achieved_max_var = (sd ^ 2) * as.numeric(unit$max_var)
  ))
}

vhar3_generate_design <- function(q = 18L,
                                  setting_id = "sg1",
                                  path = "path1",
                                  link_type = "type1",
                                  bw = 3L,
                                  bm = 10L,
                                  target_var = 0.5,
                                  noise_sd = NULL) {
  q <- as.integer(q)
  setting_name <- vhar3_setting_name(setting_id)
  path_name <- vhar3_path_name(path)
  type_name <- vhar3_link_type(link_type)

  path_obj <- vhar3_generate_path_labels(q = q, path = path_name)
  pars <- vhar3_setting_parms(setting_name)
  eff <- vhar3_effective_snr_multiplier(bw = bw, bm = bm)
  horizons <- c("daily", "weekly", "monthly")

  A_list <- vector("list", 3L)
  B_list <- vector("list", 3L)
  Phi_raw <- vector("list", 3L)

  for (hh in 1:3) {
    B_list[[hh]] <- vhar3_make_B(
      y_labels = path_obj$y[[hh]],
      z_labels = path_obj$z[[hh]],
      link_type = type_name,
      horizon = horizons[hh],
      setting_id = setting_name
    )
    A_list[[hh]] <- vhar3_sample_graph(
      y_labels = path_obj$y[[hh]],
      z_labels = path_obj$z[[hh]],
      B = B_list[[hh]],
      allow_self = TRUE
    )
    Phi_raw[[hh]] <- vhar3_scale_by_sender_group(
      A = A_list[[hh]],
      y_labels = path_obj$y[[hh]]
    )
  }

  Phi_raw[[1L]] <- pars$coeff[["d"]] * eff[1L] * Phi_raw[[1L]]
  Phi_raw[[2L]] <- pars$coeff[["w"]] * eff[2L] * Phi_raw[[2L]]
  Phi_raw[[3L]] <- pars$coeff[["m"]] * eff[3L] * Phi_raw[[3L]]

  scaled <- vhar3_scale_phi_to_target_sv(Phi_raw, sv_target = pars$sv_target)

  out <- list(
    Phi = lapply(scaled$Phi, vhar3_safe_matrix),
    A = A_list,
    B = B_list,
    y = path_obj$y,
    z = path_obj$z,
    scale_fac = scaled$scale_fac,
    sv = scaled$sv,
    rho = vhar3_companion_radius(scaled$Phi, bw = bw, bm = bm),
    setting_id = setting_name,
    path = path_name,
    type = type_name,
    bw = as.integer(bw),
    bm = as.integer(bm)
  )

  if (is.null(noise_sd) && length(target_var) == 1L && is.finite(target_var)) {
    ni <- vhar3_noise_scale_for_target_var(out$Phi, bw = bw, bm = bm, target_var = target_var)
    out$noise_sd <- as.numeric(ni$noise_sd)
    out$target_var <- as.numeric(target_var)
    out$noise_info <- ni
  } else {
    out$noise_sd <- if (is.null(noise_sd)) sqrt(0.5) else as.numeric(noise_sd)
    out$target_var <- if (length(target_var) == 1L && is.finite(target_var)) as.numeric(target_var) else NA_real_
    out$noise_info <- NULL
  }
  out
}

# ------------------------------
# Simulation
# ------------------------------
vhar3_simulate <- function(Phi, T_max, bw = 3L, bm = 10L,
                           burn = 300L, noise_sd = 0.5,
                           seed_data = NULL) {
  if (!is.null(seed_data)) set.seed(as.integer(seed_data))
  q <- nrow(Phi[[1L]])
  total <- as.integer(T_max + burn + bm)
  Y <- matrix(0, nrow = q, ncol = total)
  E <- as.numeric(noise_sd) * matrix(rnorm(q * total), nrow = q, ncol = total)

  for (tt in (bm + 1L):total) {
    yd <- Y[, tt - 1L]
    yw <- rowMeans(Y[, (tt - bw):(tt - 1L), drop = FALSE])
    ym <- rowMeans(Y[, (tt - bm):(tt - 1L), drop = FALSE])
    Y[, tt] <- Phi[[1L]] %*% yd + Phi[[2L]] %*% yw + Phi[[3L]] %*% ym + E[, tt]
  }
  Y[, (burn + bm + 1L):total, drop = FALSE]
}

# ------------------------------
# Fitting
# ------------------------------
vhar3_split_phi_blocks <- function(Bhat, q) {
  Bhat <- vhar3_safe_matrix(Bhat)
  list(
    vhar3_safe_matrix(Bhat[, 1:q, drop = FALSE]),
    vhar3_safe_matrix(Bhat[, (q + 1L):(2L * q), drop = FALSE]),
    vhar3_safe_matrix(Bhat[, (2L * q + 1L):(3L * q), drop = FALSE])
  )
}

vhar3_fit_ols <- function(Yt, bw = 3L, bm = 10L) {
  fit <- sparseVAR::VHAR_ols(Yt = Yt, bd = as.integer(bw), bm = as.integer(bm))
  Bhat <- vhar3_safe_matrix(fit$Phi_hat)
  q <- nrow(Yt)
  list(
    Phi_hat = vhar3_split_phi_blocks(Bhat, q = q),
    Bhat = Bhat,
    lambda_pkg = rep(NA_real_, 3L),
    lambda_text = rep(NA_real_, 3L),
    c_lambda_sel = NA_real_,
    lambda_mode = "ols"
  )
}

vhar3_lambda_textbook <- function(q, N_eff, theory_const = 1) {
  as.numeric(theory_const) * sqrt(log(3 * q^2) / as.numeric(N_eff))
}

vhar3_lambda_pkg_from_text <- function(lambda_text, N_eff) {
  0.5 * as.numeric(N_eff) * as.numeric(lambda_text)
}

vhar3_lambda_text_from_pkg <- function(lambda_pkg, N_eff) {
  if (any(!is.finite(lambda_pkg)) || !is.finite(N_eff) || N_eff <= 0) {
    out <- rep(NA_real_, length(lambda_pkg))
    return(out)
  }
  2 * as.numeric(lambda_pkg) / as.numeric(N_eff)
}

vhar3_select_c_lambda_blockcv <- function(Yt, bw = 3L, bm = 10L,
                                          c_grid = vhar3_default_c_grid(),
                                          fold = 10L,
                                          max_iter = 1000L,
                                          tol = 1e-6,
                                          diagTF = TRUE,
                                          updateSigma = FALSE,
                                          sigma_diag_only = TRUE) {
  q <- nrow(Yt)
  N_eff <- ncol(Yt) - as.integer(bm)
  if (N_eff <= 0L) stop("Need ncol(Yt) > bm for c-lambda CV.", call. = FALSE)

  c_grid <- unique(as.numeric(c_grid))
  c_grid <- c_grid[is.finite(c_grid) & c_grid > 0]
  if (!length(c_grid)) stop("c_grid must contain positive finite values.", call. = FALSE)
  c_grid <- sort(c_grid)

  lambda_text_base <- vhar3_lambda_textbook(q = q, N_eff = N_eff, theory_const = 1)
  lambda_pkg_base <- vhar3_lambda_pkg_from_text(lambda_text_base, N_eff = N_eff)
  lambda_seq <- sort(unique(as.numeric(c_grid) * as.numeric(lambda_pkg_base)), decreasing = TRUE)

  fit <- sparseVAR::VHAR_adalasso_fista(
    Yt = Yt,
    type = "lasso",
    fold = as.integer(fold),
    lambda_seq = as.numeric(lambda_seq),
    nlambda = as.integer(length(lambda_seq)),
    lambda_min_ratio = NULL,
    diagTF = isTRUE(diagTF),
    updateSigma = isTRUE(updateSigma),
    sigma_diag_only = isTRUE(sigma_diag_only),
    max_iter = as.integer(max_iter),
    tol = tol,
    bd = as.integer(bw),
    bm = as.integer(bm)
  )

  lambda_sel <- if (!is.null(fit$lambda_lasso)) as.numeric(fit$lambda_lasso) else NA_real_
  if (!is.finite(lambda_sel)) {
    lambda_sel <- lambda_seq[which.min(abs(lambda_seq - stats::median(lambda_seq)))]
  }

  c_sel <- lambda_sel / lambda_pkg_base
  c_sel <- c_grid[which.min(abs(c_grid - c_sel))]
  lambda_sel <- lambda_pkg_base * c_sel
  lambda_text_sel <- vhar3_lambda_text_from_pkg(lambda_sel, N_eff = N_eff)

  list(
    c_lambda_sel = as.numeric(c_sel),
    lambda_pkg_sel = as.numeric(lambda_sel),
    lambda_text_sel = as.numeric(lambda_text_sel),
    lambda_pkg_base = as.numeric(lambda_pkg_base),
    lambda_text_base = as.numeric(lambda_text_base),
    c_grid = c_grid,
    lambda_seq = lambda_seq,
    fit = fit
  )
}

vhar3_fit_lasso <- function(Yt, bw = 3L, bm = 10L,
                            lambda_mode = c("theory", "cv", "cv_c"),
                            c_lambda = 0.25,
                            c_grid = vhar3_default_c_grid(),
                            fold = 10L,
                            nlambda = 100L,
                            max_iter = 1000L,
                            tol = 1e-6,
                            diagTF = TRUE,
                            updateSigma = FALSE,
                            sigma_diag_only = TRUE) {
  lambda_mode <- vhar3_lambda_mode(lambda_mode)
  q <- nrow(Yt)
  N_eff <- ncol(Yt) - as.integer(bm)
  if (N_eff <= 0L) stop("Need ncol(Yt) > bm.", call. = FALSE)

  if (lambda_mode == "theory") {
    lambda_text <- vhar3_lambda_textbook(q = q, N_eff = N_eff, theory_const = c_lambda)
    lambda_pkg <- vhar3_lambda_pkg_from_text(lambda_text, N_eff = N_eff)
    fit <- sparseVAR::VHAR_adalasso_fista(
      Yt = Yt,
      type = "lasso",
      lambda = as.numeric(lambda_pkg),
      diagTF = isTRUE(diagTF),
      updateSigma = isTRUE(updateSigma),
      sigma_diag_only = isTRUE(sigma_diag_only),
      max_iter = as.integer(max_iter),
      tol = tol,
      bd = as.integer(bw),
      bm = as.integer(bm)
    )
    c_sel <- as.numeric(c_lambda)[1L]
  } else if (lambda_mode == "cv") {
    fit <- sparseVAR::VHAR_adalasso_fista(
      Yt = Yt,
      type = "lasso",
      fold = as.integer(fold),
      nlambda = as.integer(nlambda),
      lambda_min_ratio = NULL,
      diagTF = isTRUE(diagTF),
      updateSigma = isTRUE(updateSigma),
      sigma_diag_only = isTRUE(sigma_diag_only),
      max_iter = as.integer(max_iter),
      tol = tol,
      bd = as.integer(bw),
      bm = as.integer(bm)
    )
    lambda_pkg <- if (!is.null(fit$lambda_lasso)) as.numeric(fit$lambda_lasso) else NA_real_
    if (length(lambda_pkg) == 1L) lambda_pkg <- rep(lambda_pkg, 3L)
    lambda_text <- vhar3_lambda_text_from_pkg(lambda_pkg, N_eff = N_eff)
    lambda_text_base <- vhar3_lambda_textbook(q = q, N_eff = N_eff, theory_const = 1)
    c_sel <- vhar3_mean(lambda_text / lambda_text_base)
  } else {
    sel <- vhar3_select_c_lambda_blockcv(
      Yt = Yt,
      bw = bw,
      bm = bm,
      c_grid = c_grid,
      fold = fold,
      max_iter = max_iter,
      tol = tol,
      diagTF = diagTF,
      updateSigma = updateSigma,
      sigma_diag_only = sigma_diag_only
    )
    fit <- sparseVAR::VHAR_adalasso_fista(
      Yt = Yt,
      type = "lasso",
      lambda = as.numeric(sel$lambda_pkg_sel),
      diagTF = isTRUE(diagTF),
      updateSigma = isTRUE(updateSigma),
      sigma_diag_only = isTRUE(sigma_diag_only),
      max_iter = as.integer(max_iter),
      tol = tol,
      bd = as.integer(bw),
      bm = as.integer(bm)
    )
    lambda_pkg <- as.numeric(sel$lambda_pkg_sel)
    lambda_text <- as.numeric(sel$lambda_text_sel)
    c_sel <- as.numeric(sel$c_lambda_sel)
  }

  Braw <- fit$Phi_hat_lasso
  if (is.null(Braw)) Braw <- fit$Phi_lasso_hat
  Braw <- vhar3_safe_matrix(Braw)

  list(
    Phi_hat = vhar3_split_phi_blocks(Braw, q = q),
    Bhat = Braw,
    lambda_pkg = rep(as.numeric(lambda_pkg), 3L),
    lambda_text = rep(as.numeric(lambda_text), 3L),
    c_lambda_sel = as.numeric(c_sel),
    lambda_mode = lambda_mode
  )
}

# ------------------------------
# PisCES and clustering
# ------------------------------
vhar3_embeddings_from_blocks <- function(blocks, n_comm) {
  s <- length(blocks)
  emb_L <- vector("list", s)
  emb_R <- vector("list", s)
  proj_L <- vector("list", s)
  proj_R <- vector("list", s)

  for (mm in seq_len(s)) {
    Ky <- n_comm[2L * mm - 1L]
    Kz <- n_comm[2L * mm]
    Ahat <- t(blocks[[mm]])
    sv <- svd(Ahat, nu = as.integer(Ky), nv = as.integer(Kz))
    XL <- sv$u[, seq_len(as.integer(Ky)), drop = FALSE]
    XR <- sv$v[, seq_len(as.integer(Kz)), drop = FALSE]
    emb_L[[mm]] <- vhar3_row_normalize(XL)
    emb_R[[mm]] <- vhar3_row_normalize(XR)
    proj_L[[mm]] <- XL %*% t(XL)
    proj_R[[mm]] <- XR %*% t(XR)
  }
  list(emb_L = emb_L, emb_R = emb_R, proj_L = proj_L, proj_R = proj_R)
}

vhar3_projector_topK <- function(M, K) {
  ee <- eigen(M, symmetric = TRUE)
  U <- ee$vectors[, seq_len(as.integer(K)), drop = FALSE]
  U %*% t(U)
}

vhar3_pisces_smooth <- function(proj_list, ranks, alpha,
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
      nxt[[mm]] <- vhar3_projector_topK(S, ranks[mm])
      diff <- diff + sqrt(sum((nxt[[mm]] - cur[[mm]])^2))
    }
    cur <- nxt
  }
  list(projectors = cur, iter = iter)
}

vhar3_embeddings_from_projectors <- function(proj_list, ranks) {
  s <- length(proj_list)
  emb <- vector("list", s)
  for (mm in seq_len(s)) {
    ee <- eigen(proj_list[[mm]], symmetric = TRUE)
    U <- ee$vectors[, seq_len(as.integer(ranks[mm])), drop = FALSE]
    emb[[mm]] <- vhar3_row_normalize(U)
  }
  emb
}

vhar3_kmeans_labels <- function(X, K, nstart = 50L) {
  K <- as.integer(K)
  if (K <= 1L) return(rep(1L, nrow(X)))
  X <- vhar3_safe_matrix(X)
  if (nrow(X) < K) stop("Number of rows must be at least K.", call. = FALSE)

  run_once <- function(Xin, algorithm = "Hartigan-Wong", nstart_use = nstart) {
    stats::kmeans(Xin, centers = K, nstart = as.integer(nstart_use), iter.max = 200L, algorithm = algorithm)$cluster
  }

  ux_n <- nrow(unique(round(X, 12)))
  if (ux_n < K) {
    set.seed(1L)
    X <- X + matrix(rnorm(length(X), sd = 1e-8), nrow = nrow(X), ncol = ncol(X))
  }

  out <- tryCatch(run_once(X, algorithm = "Hartigan-Wong", nstart_use = nstart), error = function(e) NULL)
  if (!is.null(out)) return(as.integer(out))

  set.seed(11L)
  Xj <- X + matrix(rnorm(length(X), sd = 1e-7), nrow = nrow(X), ncol = ncol(X))
  out <- tryCatch(run_once(Xj, algorithm = "Lloyd", nstart_use = nstart), error = function(e) NULL)
  if (!is.null(out)) return(as.integer(out))

  out <- tryCatch(run_once(Xj, algorithm = "MacQueen", nstart_use = 1L), error = function(e) NULL)
  if (!is.null(out)) return(as.integer(out))

  hc <- hclust(dist(Xj), method = "ward.D2")
  as.integer(cutree(hc, k = K))
}

vhar3_cluster <- function(blocks, n_comm, alpha = 0) {
  s <- 3L
  tmp <- vhar3_embeddings_from_blocks(blocks, n_comm)

  if (alpha > 0) {
    ranks_L <- n_comm[c(1L, 3L, 5L)]
    ranks_R <- n_comm[c(2L, 4L, 6L)]
    smL <- vhar3_pisces_smooth(tmp$proj_L, ranks = ranks_L, alpha = alpha)
    smR <- vhar3_pisces_smooth(tmp$proj_R, ranks = ranks_R, alpha = alpha)
    emb_L <- vhar3_embeddings_from_projectors(smL$projectors, ranks_L)
    emb_R <- vhar3_embeddings_from_projectors(smR$projectors, ranks_R)
  } else {
    emb_L <- tmp$emb_L
    emb_R <- tmp$emb_R
  }

  group_L <- vector("list", s)
  group_R <- vector("list", s)

  group_L[[1L]] <- vhar3_kmeans_labels(emb_L[[1L]], n_comm[1L])
  group_R[[3L]] <- vhar3_kmeans_labels(emb_R[[3L]], n_comm[6L])

  if (n_comm[2L] == n_comm[3L]) {
    lab12 <- vhar3_kmeans_labels(cbind(emb_R[[1L]], emb_L[[2L]]), n_comm[2L])
    group_R[[1L]] <- lab12
    group_L[[2L]] <- lab12
  } else {
    group_R[[1L]] <- vhar3_kmeans_labels(emb_R[[1L]], n_comm[2L])
    group_L[[2L]] <- vhar3_kmeans_labels(emb_L[[2L]], n_comm[3L])
  }

  if (n_comm[4L] == n_comm[5L]) {
    lab23 <- vhar3_kmeans_labels(cbind(emb_R[[2L]], emb_L[[3L]]), n_comm[4L])
    group_R[[2L]] <- lab23
    group_L[[3L]] <- lab23
  } else {
    group_R[[2L]] <- vhar3_kmeans_labels(emb_R[[2L]], n_comm[4L])
    group_L[[3L]] <- vhar3_kmeans_labels(emb_L[[3L]], n_comm[5L])
  }

  list(group_L = group_L, group_R = group_R, emb_L = emb_L, emb_R = emb_R, alpha = alpha)
}

vhar3_labels_to_membership <- function(labels, K) {
  labels <- as.integer(labels)
  q <- length(labels)
  out <- matrix(0, nrow = q, ncol = as.integer(K))
  for (kk in seq_len(as.integer(K))) out[labels == kk, kk] <- 1
  out
}

vhar3_completed_from_mask <- function(M, keep_mask, rankK) {
  M0 <- M * keep_mask
  sv <- svd(M0, nu = as.integer(rankK), nv = as.integer(rankK))
  sv$u[, seq_len(as.integer(rankK)), drop = FALSE] %*%
    diag(sv$d[seq_len(as.integer(rankK))], as.integer(rankK)) %*%
    t(sv$v[, seq_len(as.integer(rankK)), drop = FALSE])
}

vhar3_predict_from_labels <- function(M, yhat, zhat) {
  Mpos <- pmax(M, 0)
  q <- nrow(Mpos)
  Ky <- max(yhat)
  Kz <- max(zhat)
  dy <- pmax(rowSums(Mpos), 1e-8)
  dz <- pmax(colSums(Mpos), 1e-8)
  Y <- vhar3_labels_to_membership(yhat, Ky)
  Z <- vhar3_labels_to_membership(zhat, Kz)
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

vhar3_select_alpha <- function(blocks, n_comm,
                               folds = 5L,
                               alpha_grid = vhar3_default_alpha_grid(),
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
  folds <- max(2L, as.integer(folds))

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
      comp_blocks[[ff]][[mm]] <- vhar3_completed_from_mask(blocks[[mm]], keep_mask, rankK)
    }
  }

  score_alpha <- numeric(length(alpha_grid))
  for (aa in seq_along(alpha_grid)) {
    alpha <- alpha_grid[aa]
    score_fold <- numeric(folds)
    for (ff in seq_len(folds)) {
      fit <- vhar3_cluster(comp_blocks[[ff]], n_comm = n_comm, alpha = alpha)
      sc <- 0
      for (mm in seq_len(s)) {
        Mcomp <- comp_blocks[[ff]][[mm]]
        Mhat <- vhar3_predict_from_labels(Mcomp, fit$group_L[[mm]], fit$group_R[[mm]])
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
# Accuracy and errors
# ------------------------------
vhar3_all_perms <- function(x) {
  x <- as.integer(x)
  if (length(x) <= 1L) return(matrix(x, nrow = 1L))
  gtools::permutations(n = length(x), r = length(x), v = x)
}

vhar3_adjusted_rand_index <- function(labels_true, labels_est) {
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

vhar3_match_accuracy <- function(true_labels, est_labels) {
  tru <- as.integer(true_labels)
  est <- as.integer(est_labels)
  ok <- is.finite(tru) & is.finite(est)
  tru <- tru[ok]
  est <- est[ok]
  if (!length(tru)) return(list(acc = NA_real_, ari = NA_real_, matched = integer(0)))
  K <- max(max(tru), max(est))
  perms <- vhar3_all_perms(seq_len(K))
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
  list(acc = as.numeric(best_acc), ari = vhar3_adjusted_rand_index(tru, best_lab), matched = best_lab)
}

vhar3_score <- function(true_y, true_z, fit_obj) {
  met1 <- vhar3_match_accuracy(true_y[[1L]], fit_obj$group_L[[1L]])
  met2 <- vhar3_match_accuracy(true_z[[1L]], fit_obj$group_R[[1L]])
  met3 <- vhar3_match_accuracy(true_z[[2L]], fit_obj$group_R[[2L]])
  met4 <- vhar3_match_accuracy(true_z[[3L]], fit_obj$group_R[[3L]])

  acc <- c(met1$acc, met2$acc, met3$acc, met4$acc)
  ari <- c(met1$ari, met2$ari, met3$ari, met4$ari)
  names(acc) <- c("daily_send", "daily_recv_weekly_send", "weekly_recv_monthly_send", "monthly_recv")
  names(ari) <- names(acc)
  list(acc = acc, ari = ari, acc_overall = mean(acc), ari_overall = mean(ari))
}

vhar3_fro_errors <- function(Phi_true, Phi_hat) {
  c(
    fro_d = vhar3_fro_norm(Phi_hat[[1L]] - Phi_true[[1L]]),
    fro_w = vhar3_fro_norm(Phi_hat[[2L]] - Phi_true[[2L]]),
    fro_m = vhar3_fro_norm(Phi_hat[[3L]] - Phi_true[[3L]])
  )
}

vhar3_op_errors <- function(Phi_true, Phi_hat) {
  c(
    op_d = vhar3_operator_norm(Phi_hat[[1L]] - Phi_true[[1L]]),
    op_w = vhar3_operator_norm(Phi_hat[[2L]] - Phi_true[[2L]]),
    op_m = vhar3_operator_norm(Phi_hat[[3L]] - Phi_true[[3L]])
  )
}

# ------------------------------
# Result rows and helpers
# ------------------------------
vhar3_make_result_row <- function(setting, path, q, link_type, estimator, T, rep_id,
                                  seed, design_id,
                                  score = NULL, frob = NULL, opb = NULL,
                                  c_lambda_sel = NA_real_,
                                  lambda_pkg = rep(NA_real_, 3L),
                                  lambda_text = rep(NA_real_, 3L),
                                  lambda_mode = NA_character_,
                                  smoothing = c("none", "pisces"),
                                  alpha_sel = NA_real_,
                                  error_message = "") {
  smoothing <- match.arg(smoothing)
  if (is.null(score)) {
    score <- list(acc = rep(NA_real_, 4L), ari = rep(NA_real_, 4L), acc_overall = NA_real_, ari_overall = NA_real_)
  }
  if (is.null(frob)) {
    frob <- c(fro_d = NA_real_, fro_w = NA_real_, fro_m = NA_real_)
  }
  if (is.null(opb)) {
    opb <- c(op_d = NA_real_, op_w = NA_real_, op_m = NA_real_)
  }
  lambda_pkg <- as.numeric(lambda_pkg)
  lambda_text <- as.numeric(lambda_text)
  if (length(lambda_pkg) == 1L) lambda_pkg <- rep(lambda_pkg, 3L)
  if (length(lambda_text) == 1L) lambda_text <- rep(lambda_text, 3L)

  data.frame(
    setting = as.character(setting),
    path = as.character(path),
    q = as.integer(q),
    type = as.character(link_type),
    estimator = as.character(estimator),
    smoothing = as.character(smoothing),
    T = as.integer(T),
    rep_id = as.integer(rep_id),
    seed = as.integer(seed),
    design_id = as.integer(design_id),
    fro_d = as.numeric(frob["fro_d"]),
    fro_w = as.numeric(frob["fro_w"]),
    fro_m = as.numeric(frob["fro_m"]),
    op_d = as.numeric(opb["op_d"]),
    op_w = as.numeric(opb["op_w"]),
    op_m = as.numeric(opb["op_m"]),
    Fro = vhar3_mean(frob),
    Spec = vhar3_mean(opb),
    acc_daily_send = as.numeric(score$acc[1L]),
    acc_daily_recv_weekly_send = as.numeric(score$acc[2L]),
    acc_weekly_recv_monthly_send = as.numeric(score$acc[3L]),
    acc_monthly_recv = as.numeric(score$acc[4L]),
    ari_daily_send = as.numeric(score$ari[1L]),
    ari_daily_recv_weekly_send = as.numeric(score$ari[2L]),
    ari_weekly_recv_monthly_send = as.numeric(score$ari[3L]),
    ari_monthly_recv = as.numeric(score$ari[4L]),
    acc_overall = as.numeric(score$acc_overall),
    ari_overall = as.numeric(score$ari_overall),
    lambda_pkg_d = as.numeric(lambda_pkg[1L]),
    lambda_pkg_w = as.numeric(lambda_pkg[2L]),
    lambda_pkg_m = as.numeric(lambda_pkg[3L]),
    lambda_text_d = as.numeric(lambda_text[1L]),
    lambda_text_w = as.numeric(lambda_text[2L]),
    lambda_text_m = as.numeric(lambda_text[3L]),
    lambda_pkg_mean = vhar3_mean(lambda_pkg),
    lambda_text_mean = vhar3_mean(lambda_text),
    c_lambda_sel = as.numeric(c_lambda_sel),
    lambda_mode = as.character(lambda_mode),
    alpha_sel = as.numeric(alpha_sel),
    error_message = as.character(error_message),
    stringsAsFactors = FALSE
  )
}

vhar3_pick_map_value <- function(map_obj, key, fallback = 0) {
  if (is.null(map_obj)) return(as.numeric(fallback)[1L])
  if (!is.null(names(map_obj)) && as.character(key) %in% names(map_obj)) {
    return(as.numeric(map_obj[as.character(key)])[1L])
  }
  if (length(map_obj) == 1L) return(as.numeric(map_obj)[1L])
  as.numeric(fallback)[1L]
}

vhar3_design_key <- function(q, setting_id, path, link_type) {
  q <- as.integer(q)
  sid <- as.integer(sub("sg", "", vhar3_setting_name(setting_id)))
  pid <- as.integer(sub("path", "", vhar3_path_name(path)))
  tid <- if (vhar3_link_type(link_type) == "type1") 1L else 2L
  as.integer(q * 10000L + sid * 100L + pid * 10L + tid)
}

# ------------------------------
# One-rep simulation
# ------------------------------
vhar3_run_one_rep <- function(rep_id,
                              fixed_gen,
                              setting_id,
                              path,
                              method,
                              q,
                              T,
                              link_type,
                              bw,
                              bm,
                              burn,
                              seed_data,
                              design_id = NA_integer_,
                              lasso_lambda_mode = "cv_c",
                              lasso_c_lambda = 0.25,
                              lasso_c_grid = vhar3_default_c_grid(),
                              lasso_fold = 10L,
                              lasso_nlambda = 100L,
                              lasso_max_iter = 1000L,
                              lasso_tol = 1e-6,
                              alpha_map = NULL,
                              c_map = NULL) {
  method <- vhar3_method(method)
  setting_name <- vhar3_setting_name(setting_id)
  path_name <- vhar3_path_name(path)
  type_name <- vhar3_link_type(link_type)

  Ylong <- vhar3_simulate(
    Phi = fixed_gen$Phi,
    T_max = max(T),
    bw = bw,
    bm = bm,
    burn = burn,
    noise_sd = fixed_gen$noise_sd,
    seed_data = seed_data
  )
  n_comm <- vhar3_path_spec(path_name)

  rows <- lapply(T, function(TT) {
    alpha_use <- vhar3_pick_map_value(alpha_map, TT, fallback = 0)

    fit <- tryCatch({
      if (method == "ols") {
        vhar3_fit_ols(Ylong[, seq_len(TT), drop = FALSE], bw = bw, bm = bm)
      } else {
        lambda_mode_use <- if (lasso_lambda_mode == "cv_c") "theory" else lasso_lambda_mode
        c_use <- vhar3_pick_map_value(c_map, TT, fallback = lasso_c_lambda)
        vhar3_fit_lasso(
          Yt = Ylong[, seq_len(TT), drop = FALSE],
          bw = bw,
          bm = bm,
          lambda_mode = lambda_mode_use,
          c_lambda = c_use,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          nlambda = lasso_nlambda,
          max_iter = lasso_max_iter,
          tol = lasso_tol,
          diagTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE
        )
      }
    }, error = function(e) e)

    if (inherits(fit, "error") || inherits(fit, "simpleError")) {
      row_none <- vhar3_make_result_row(
        setting = setting_name,
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        c_lambda_sel = if (method == "lasso") vhar3_pick_map_value(c_map, TT, fallback = lasso_c_lambda) else NA_real_,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        smoothing = "none",
        alpha_sel = 0,
        error_message = conditionMessage(fit)
      )
      row_pisces <- vhar3_make_result_row(
        setting = setting_name,
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        c_lambda_sel = if (method == "lasso") vhar3_pick_map_value(c_map, TT, fallback = lasso_c_lambda) else NA_real_,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        smoothing = "pisces",
        alpha_sel = alpha_use,
        error_message = conditionMessage(fit)
      )
      return(rbind(row_none, row_pisces))
    }

    blocks_hat <- fit$Phi_hat
    frob <- vhar3_fro_errors(fixed_gen$Phi, blocks_hat)
    opb <- vhar3_op_errors(fixed_gen$Phi, blocks_hat)

    clfit_none <- tryCatch(vhar3_cluster(blocks_hat, n_comm = n_comm, alpha = 0), error = function(e) e)
    if (inherits(clfit_none, "error") || inherits(clfit_none, "simpleError")) {
      row_none <- vhar3_make_result_row(
        setting = setting_name,
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        frob = frob,
        opb = opb,
        c_lambda_sel = if (method == "lasso") fit$c_lambda_sel else NA_real_,
        lambda_pkg = fit$lambda_pkg,
        lambda_text = fit$lambda_text,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        smoothing = "none",
        alpha_sel = 0,
        error_message = conditionMessage(clfit_none)
      )
    } else {
      score_none <- vhar3_score(true_y = fixed_gen$y, true_z = fixed_gen$z, fit_obj = clfit_none)
      row_none <- vhar3_make_result_row(
        setting = setting_name,
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
        smoothing = "none",
        alpha_sel = 0,
        error_message = ""
      )
    }

    clfit_pisces <- tryCatch(vhar3_cluster(blocks_hat, n_comm = n_comm, alpha = alpha_use), error = function(e) e)
    if (inherits(clfit_pisces, "error") || inherits(clfit_pisces, "simpleError")) {
      row_pisces <- vhar3_make_result_row(
        setting = setting_name,
        path = path_name,
        q = q,
        link_type = type_name,
        estimator = method,
        T = TT,
        rep_id = rep_id,
        seed = seed_data,
        design_id = design_id,
        frob = frob,
        opb = opb,
        c_lambda_sel = if (method == "lasso") fit$c_lambda_sel else NA_real_,
        lambda_pkg = fit$lambda_pkg,
        lambda_text = fit$lambda_text,
        lambda_mode = if (method == "lasso") lasso_lambda_mode else "ols",
        smoothing = "pisces",
        alpha_sel = alpha_use,
        error_message = conditionMessage(clfit_pisces)
      )
    } else {
      score_pisces <- vhar3_score(true_y = fixed_gen$y, true_z = fixed_gen$z, fit_obj = clfit_pisces)
      row_pisces <- vhar3_make_result_row(
        setting = setting_name,
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
        smoothing = "pisces",
        alpha_sel = alpha_use,
        error_message = ""
      )
    }

    rbind(row_none, row_pisces)
  })

  do.call(rbind, rows)
}

# ------------------------------
# Main simulation runner
# ------------------------------
vhar3_run_dosim <- function(case_name,
                            path,
                            method = c("ols", "lasso"),
                            q = 18L,
                            T = c(500L, 1000L, 2000L, 3000L, 4000L),
                            nrep = 100L,
                            ncore = 50L,
                            link_type = "type1",
                            bd = 3L,
                            bm = 10L,
                            burn = 300L,
                            target_var = 0.5,
                            noise_sd = NULL,
                            lasso_lambda_mode = c("theory", "cv", "cv_c"),
                            lasso_c_lambda = 0.25,
                            lasso_c_grid = vhar3_default_c_grid(),
                            lasso_fold = 10L,
                            lasso_nlambda = 100L,
                            lasso_max_iter = 1000L,
                            lasso_tol = 1e-6,
                            alpha_fold = 5L,
                            alpha_grid = vhar3_default_alpha_grid(),
                            alpha_criterion = c("holdout", "paper"),
                            base_seed = 20260401L,
                            libfile = "vhar3sparsefixed.R") {
  method <- vhar3_method(method)
  setting_name <- vhar3_setting_name(case_name)
  path_name <- vhar3_path_name(path)
  type_name <- vhar3_link_type(link_type)
  lasso_lambda_mode <- vhar3_lambda_mode(lasso_lambda_mode)
  alpha_criterion <- match.arg(alpha_criterion)

  q <- as.integer(q)
  T <- sort(unique(as.integer(T)))
  nrep <- as.integer(nrep)
  ncore <- as.integer(ncore)
  design_key <- vhar3_design_key(q = q, setting_id = setting_name, path = path_name, link_type = type_name)

  phi_seed <- vhar3_make_seed(base_seed, 900000L, design_key)
  set.seed(phi_seed)
  fixed_gen <- vhar3_generate_design(
    q = q,
    setting_id = setting_name,
    path = path_name,
    link_type = type_name,
    bw = bd,
    bm = bm,
    target_var = target_var,
    noise_sd = noise_sd
  )

  # select c_lambda for each T if requested
  c_map <- lasso_c_lambda
  cv_c_info <- NULL
  if (method == "lasso" && identical(lasso_lambda_mode, "cv_c")) {
    pilot_seed <- vhar3_make_seed(base_seed, 700000L, design_key)
    Ypilot_long <- vhar3_simulate(
      Phi = fixed_gen$Phi,
      T_max = max(T),
      bw = bd,
      bm = bm,
      burn = burn,
      noise_sd = fixed_gen$noise_sd,
      seed_data = pilot_seed
    )
    c_map <- setNames(rep(NA_real_, length(T)), as.character(T))
    cv_c_info <- vector("list", length(T))
    names(cv_c_info) <- as.character(T)

    for (TT in T) {
      sel <- tryCatch(
        vhar3_select_c_lambda_blockcv(
          Yt = Ypilot_long[, seq_len(TT), drop = FALSE],
          bw = bd,
          bm = bm,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          max_iter = lasso_max_iter,
          tol = lasso_tol,
          diagTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE
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

  # select alpha for each T from pilot fits
  alpha_map <- setNames(rep(0, length(T)), as.character(T))
  alpha_cv_info <- vector("list", length(T))
  names(alpha_cv_info) <- as.character(T)
  n_comm <- vhar3_path_spec(path_name)
  pilot_seed2 <- vhar3_make_seed(base_seed, 800000L, design_key)
  Ypilot_long2 <- vhar3_simulate(
    Phi = fixed_gen$Phi,
    T_max = max(T),
    bw = bd,
    bm = bm,
    burn = burn,
    noise_sd = fixed_gen$noise_sd,
    seed_data = pilot_seed2
  )

  for (TT in T) {
    Ypilot <- Ypilot_long2[, seq_len(TT), drop = FALSE]
    fit_pilot <- tryCatch({
      if (method == "ols") {
        vhar3_fit_ols(Ypilot, bw = bd, bm = bm)
      } else {
        lambda_mode_use <- if (lasso_lambda_mode == "cv_c") "theory" else lasso_lambda_mode
        c_use <- vhar3_pick_map_value(c_map, TT, fallback = lasso_c_lambda)
        vhar3_fit_lasso(
          Yt = Ypilot,
          bw = bd,
          bm = bm,
          lambda_mode = lambda_mode_use,
          c_lambda = c_use,
          c_grid = lasso_c_grid,
          fold = lasso_fold,
          nlambda = lasso_nlambda,
          max_iter = lasso_max_iter,
          tol = lasso_tol,
          diagTF = TRUE,
          updateSigma = FALSE,
          sigma_diag_only = TRUE
        )
      }
    }, error = function(e) NULL)

    if (is.null(fit_pilot)) {
      alpha_map[as.character(TT)] <- 0
      alpha_cv_info[[as.character(TT)]] <- list(error = TRUE, alpha = 0)
    } else {
      sel_alpha <- tryCatch(
        vhar3_select_alpha(
          blocks = fit_pilot$Phi_hat,
          n_comm = n_comm,
          folds = alpha_fold,
          alpha_grid = alpha_grid,
          criterion = alpha_criterion,
          seed = vhar3_make_seed(base_seed, 810000L, design_key, TT)
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
    seed_i <- vhar3_make_seed(base_seed, 100000L * design_key, rep_id)
    vhar3_run_one_rep(
      rep_id = rep_id,
      fixed_gen = fixed_gen,
      setting_id = setting_name,
      path = path_name,
      method = method,
      q = q,
      T = T,
      link_type = type_name,
      bw = bd,
      bm = bm,
      burn = burn,
      seed_data = seed_i,
      design_id = design_key,
      lasso_lambda_mode = lasso_lambda_mode,
      lasso_c_lambda = lasso_c_lambda,
      lasso_c_grid = lasso_c_grid,
      lasso_fold = lasso_fold,
      lasso_nlambda = lasso_nlambda,
      lasso_max_iter = lasso_max_iter,
      lasso_tol = lasso_tol,
      alpha_map = alpha_map,
      c_map = c_map
    )
  }

  raw <- do.call(rbind, raw_list)
  raw <- raw[order(raw$rep_id, raw$T, raw$smoothing), , drop = FALSE]
  rownames(raw) <- NULL

  out <- list(
    setting = setting_name,
    path = path_name,
    method = method,
    estimator = method,
    q = q,
    type = type_name,
    T = T,
    nrep = nrep,
    ncore = ncore,
    bw = as.integer(bd),
    bm = as.integer(bm),
    burn = as.integer(burn),
    target_var = if (length(target_var) == 1L && is.finite(target_var)) as.numeric(target_var) else NA_real_,
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
  class(out) <- c("vhar3_dosim", class(out))
  out
}

# ------------------------------
# Summary table
# ------------------------------
vhar3_table_summary <- function(...) {
  xs <- list(...)
  if (length(xs) == 1L && is.list(xs[[1L]]) && !inherits(xs[[1L]], "vhar3_dosim")) {
    xs <- xs[[1L]]
  }
  if (!length(xs)) stop("Provide at least one vhar3_dosim object.", call. = FALSE)
  ok <- vapply(xs, function(z) inherits(z, "vhar3_dosim"), logical(1L))
  if (!all(ok)) stop("All inputs must be vhar3_dosim objects.", call. = FALSE)

  raw <- do.call(rbind, lapply(xs, function(z) z$raw))
  raw_ok <- raw[is.na(raw$error_message) | raw$error_message == "", , drop = FALSE]
  if (!nrow(raw_ok)) return(data.frame())

  agg <- aggregate(
    cbind(Fro, Spec, acc_overall, ari_overall, alpha_sel) ~ setting + path + q + type + estimator + smoothing + T,
    data = raw_ok,
    FUN = function(z) mean(z, na.rm = TRUE)
  )

  keys <- unique(agg[, c("setting", "path", "q", "type", "estimator", "smoothing"), drop = FALSE])
  keys <- keys[order(keys$setting, keys$path, keys$q, keys$type, keys$estimator, keys$smoothing), , drop = FALSE]
  out <- keys
  key_out <- paste(out$setting, out$path, out$q, out$type, out$estimator, out$smoothing, sep = "__")

  for (TT in sort(unique(agg$T))) {
    sub <- agg[agg$T == TT, , drop = FALSE]
    key_sub <- paste(sub$setting, sub$path, sub$q, sub$type, sub$estimator, sub$smoothing, sep = "__")
    out[[paste0("Fro_T", TT)]] <- as.numeric(setNames(sub$Fro, key_sub)[key_out])
    out[[paste0("Spec_T", TT)]] <- as.numeric(setNames(sub$Spec, key_sub)[key_out])
    out[[paste0("acc_T", TT)]] <- as.numeric(setNames(sub$acc_overall, key_sub)[key_out])
    out[[paste0("ari_T", TT)]] <- as.numeric(setNames(sub$ari_overall, key_sub)[key_out])
    out[[paste0("alpha_T", TT)]] <- as.numeric(setNames(sub$alpha_sel, key_sub)[key_out])
  }

  rownames(out) <- NULL
  vhar3_round_df(out, digits = 3L)
}

# ------------------------------
# GitHub-facing aliases
# ------------------------------
scbm_vhar_default_c_grid <- vhar3_default_c_grid
scbm_vhar_default_alpha_grid <- vhar3_default_alpha_grid
scbm_vhar_run_simulation <- vhar3_run_dosim
scbm_vhar_simulation_summary <- vhar3_table_summary
