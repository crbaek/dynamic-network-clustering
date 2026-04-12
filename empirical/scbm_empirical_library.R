# ============================================================
# library3.R
# Empirical analysis library for ScBM-PVAR and ScBM-VHAR
# - assumes the data matrix is given
# - user supplies lag order / horizon lengths and community counts
# - supports OLS / lasso first-stage estimation
# - supports optional PisCES smoothing by CV or fixed alpha
# - returns coefficient estimates, community assignments,
#   flow tables, and exportable summaries
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("sparseVAR", quietly = TRUE)) {
    stop("Package 'sparseVAR' is required but not installed.", call. = FALSE)
  }
})

# ------------------------------
# Basic helpers
# ------------------------------
scbm_safe_matrix <- function(X) {
  if (is.null(X)) return(matrix(NA_real_, 1L, 1L))
  X <- tryCatch(as.matrix(X), error = function(e) NULL)
  if (is.null(X)) return(matrix(NA_real_, 1L, 1L))
  if (length(dim(X)) != 2L) X <- matrix(X, ncol = 1L)
  storage.mode(X) <- "double"
  X
}

scbm_mean <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

scbm_row_normalize <- function(X) {
  X <- scbm_safe_matrix(X)
  rn <- sqrt(rowSums(X^2))
  rn[!is.finite(rn) | rn <= .Machine$double.eps] <- 1
  X / rn
}

scbm_operator_norm <- function(A) {
  A <- scbm_safe_matrix(A)
  if (any(!is.finite(A))) return(NA_real_)
  d <- tryCatch(svd(A, nu = 0L, nv = 0L)$d, error = function(e) numeric(0))
  if (!length(d)) return(NA_real_)
  as.numeric(d[1L])
}

scbm_fro_norm <- function(A) {
  A <- scbm_safe_matrix(A)
  if (any(!is.finite(A))) return(NA_real_)
  sqrt(sum(A^2))
}

scbm_round_df <- function(df, digits = 4L) {
  out <- df
  is_num <- vapply(out, is.numeric, logical(1L))
  out[is_num] <- lapply(out[is_num], function(z) round(z, digits))
  out
}

scbm_default_c_grid <- function() seq(0.10, 1.00, by = 0.05)
scbm_alpha_max <- function() 1 / (4 * sqrt(2) + 2)
scbm_default_alpha_grid <- function() {
  amax <- scbm_alpha_max()
  c(0, exp(seq(log(0.01 * amax), log(amax), length.out = 20L)))
}

scbm_as_Yt <- function(data, series_in_rows = FALSE, drop_non_numeric = TRUE) {
  if (is.data.frame(data) && isTRUE(drop_non_numeric)) {
    keep <- vapply(data, is.numeric, logical(1L))
    data <- data[, keep, drop = FALSE]
  }
  Y <- scbm_safe_matrix(data)
  if (!all(is.finite(Y))) stop("Input data must be numeric and finite.", call. = FALSE)
  if (!isTRUE(series_in_rows)) Y <- t(Y)
  Y
}

scbm_block_scree <- function(blocks, block_names = NULL) {
  s <- length(blocks)
  if (is.null(block_names)) block_names <- paste0("block", seq_len(s))
  out <- vector("list", s)
  for (mm in seq_len(s)) {
    B <- scbm_safe_matrix(blocks[[mm]])
    sv <- svd(B, nu = 0L, nv = 0L)$d
    if (!length(sv)) sv <- rep(0, nrow(B))
    tot <- sum(sv^2)
    prop <- if (tot <= 0) rep(0, length(sv)) else (sv^2) / tot
    out[[mm]] <- data.frame(
      block = block_names[mm],
      index = seq_along(sv),
      singular_value = sv,
      prop = prop,
      cumprop = cumsum(prop),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

scbm_validate_n_comm <- function(n_comm, s) {
  n_comm <- as.integer(n_comm)
  if (length(n_comm) != 2L * as.integer(s)) {
    stop("n_comm must have length 2 * s.", call. = FALSE)
  }
  if (any(!is.finite(n_comm)) || any(n_comm < 1L)) {
    stop("n_comm must contain positive integers.", call. = FALSE)
  }
  n_comm
}

scbm_make_stage_names <- function(send_names, receive_names) {
  if (length(send_names) != length(receive_names)) {
    stop("send_names and receive_names must have the same length.", call. = FALSE)
  }
  as.vector(rbind(paste0(send_names, ".sending"), paste0(receive_names, ".receiving")))
}

# ------------------------------
# Generic spectral co-clustering and PisCES helpers
# ------------------------------
scbm_labels_to_membership <- function(labels, K) {
  labels <- as.integer(labels)
  q <- length(labels)
  out <- matrix(0, nrow = q, ncol = as.integer(K))
  for (kk in seq_len(as.integer(K))) out[labels == kk, kk] <- 1
  out
}

scbm_predict_from_labels <- function(M, yhat, zhat) {
  Mpos <- pmax(scbm_safe_matrix(M), 0)
  q <- nrow(Mpos)
  Ky <- max(as.integer(yhat))
  Kz <- max(as.integer(zhat))
  dy <- pmax(rowSums(Mpos), 1e-8)
  dz <- pmax(colSums(Mpos), 1e-8)
  Y <- scbm_labels_to_membership(yhat, Ky)
  Z <- scbm_labels_to_membership(zhat, Kz)
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

scbm_completed_from_mask <- function(M, keep_mask, rankK) {
  M0 <- scbm_safe_matrix(M) * keep_mask
  sv <- svd(M0, nu = as.integer(rankK), nv = as.integer(rankK))
  sv$u[, seq_len(as.integer(rankK)), drop = FALSE] %*%
    diag(sv$d[seq_len(as.integer(rankK))], as.integer(rankK)) %*%
    t(sv$v[, seq_len(as.integer(rankK)), drop = FALSE])
}

scbm_projector_topK <- function(M, K) {
  ee <- eigen(scbm_safe_matrix(M), symmetric = TRUE)
  U <- ee$vectors[, seq_len(as.integer(K)), drop = FALSE]
  U %*% t(U)
}

scbm_embeddings_with_projectors <- function(blocks, n_comm) {
  s <- length(blocks)
  emb_L <- vector("list", s)
  emb_R <- vector("list", s)
  proj_L <- vector("list", s)
  proj_R <- vector("list", s)

  for (mm in seq_len(s)) {
    Ky <- as.integer(n_comm[2L * mm - 1L])
    Kz <- as.integer(n_comm[2L * mm])
    Ahat <- t(scbm_safe_matrix(blocks[[mm]]))
    sv <- svd(Ahat, nu = Ky, nv = Kz)
    XL <- sv$u[, seq_len(Ky), drop = FALSE]
    XR <- sv$v[, seq_len(Kz), drop = FALSE]
    emb_L[[mm]] <- scbm_row_normalize(XL)
    emb_R[[mm]] <- scbm_row_normalize(XR)
    proj_L[[mm]] <- XL %*% t(XL)
    proj_R[[mm]] <- XR %*% t(XR)
  }

  list(emb_L = emb_L, emb_R = emb_R, proj_L = proj_L, proj_R = proj_R)
}

scbm_embeddings_from_projectors <- function(proj_list, ranks) {
  s <- length(proj_list)
  emb <- vector("list", s)
  for (mm in seq_len(s)) {
    ee <- eigen(scbm_safe_matrix(proj_list[[mm]]), symmetric = TRUE)
    U <- ee$vectors[, seq_len(as.integer(ranks[mm])), drop = FALSE]
    emb[[mm]] <- scbm_row_normalize(U)
  }
  emb
}

scbm_pisces_smooth <- function(proj_list, ranks, alpha, tol = 1e-5, max_iter = 1000L) {
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
      nxt[[mm]] <- scbm_projector_topK(S, ranks[mm])
      diff <- diff + sqrt(sum((nxt[[mm]] - cur[[mm]])^2))
    }
    cur <- nxt
  }

  list(projectors = cur, iter = iter)
}

scbm_kmeans_labels <- function(X, K, nstart = 50L) {
  K <- as.integer(K)
  X <- scbm_safe_matrix(X)
  if (K <= 1L) return(rep(1L, nrow(X)))
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

scbm_cluster_blocks <- function(blocks,
                                n_comm,
                                alpha = 0,
                                circular = FALSE,
                                nstart = 50L,
                                tol = 1e-5,
                                max_iter = 1000L) {
  s <- length(blocks)
  n_comm <- scbm_validate_n_comm(n_comm, s)
  tmp <- scbm_embeddings_with_projectors(blocks, n_comm)

  if (alpha > 0) {
    ranks_L <- n_comm[seq(1L, 2L * s, by = 2L)]
    ranks_R <- n_comm[seq(2L, 2L * s, by = 2L)]
    smL <- scbm_pisces_smooth(tmp$proj_L, ranks = ranks_L, alpha = alpha, tol = tol, max_iter = max_iter)
    smR <- scbm_pisces_smooth(tmp$proj_R, ranks = ranks_R, alpha = alpha, tol = tol, max_iter = max_iter)
    emb_L <- scbm_embeddings_from_projectors(smL$projectors, ranks_L)
    emb_R <- scbm_embeddings_from_projectors(smR$projectors, ranks_R)
  } else {
    emb_L <- tmp$emb_L
    emb_R <- tmp$emb_R
  }

  group_L <- vector("list", s)
  group_R <- vector("list", s)

  if (isTRUE(circular)) {
    Ky1 <- n_comm[1L]
    Kzs <- n_comm[2L * s]
    if (Ky1 == Kzs) {
      lab <- scbm_kmeans_labels(cbind(emb_R[[s]], emb_L[[1L]]), Ky1, nstart = nstart)
      group_R[[s]] <- lab
      group_L[[1L]] <- lab
    } else {
      group_R[[s]] <- scbm_kmeans_labels(emb_R[[s]], Kzs, nstart = nstart)
      group_L[[1L]] <- scbm_kmeans_labels(emb_L[[1L]], Ky1, nstart = nstart)
    }
  } else {
    group_L[[1L]] <- scbm_kmeans_labels(emb_L[[1L]], n_comm[1L], nstart = nstart)
    group_R[[s]] <- scbm_kmeans_labels(emb_R[[s]], n_comm[2L * s], nstart = nstart)
  }

  if (s >= 2L) {
    for (mm in 2:s) {
      Ky <- n_comm[2L * mm - 1L]
      Kzprev <- n_comm[2L * (mm - 1L)]
      if (Ky == Kzprev) {
        lab <- scbm_kmeans_labels(cbind(emb_R[[mm - 1L]], emb_L[[mm]]), Ky, nstart = nstart)
        group_R[[mm - 1L]] <- lab
        group_L[[mm]] <- lab
      } else {
        group_R[[mm - 1L]] <- scbm_kmeans_labels(emb_R[[mm - 1L]], Kzprev, nstart = nstart)
        group_L[[mm]] <- scbm_kmeans_labels(emb_L[[mm]], Ky, nstart = nstart)
      }
    }
  }

  list(
    group_L = group_L,
    group_R = group_R,
    emb_L = emb_L,
    emb_R = emb_R,
    alpha = as.numeric(alpha),
    circular = isTRUE(circular)
  )
}

scbm_select_alpha <- function(blocks,
                              n_comm,
                              circular = FALSE,
                              folds = 5L,
                              alpha_grid = scbm_default_alpha_grid(),
                              criterion = c("holdout", "paper"),
                              nstart = 50L,
                              seed = NULL) {
  criterion <- match.arg(criterion)
  if (!is.null(seed)) set.seed(as.integer(seed))
  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  alpha_grid <- alpha_grid[is.finite(alpha_grid) & alpha_grid >= 0]
  if (!length(alpha_grid)) stop("alpha_grid must contain nonnegative finite values.", call. = FALSE)

  s <- length(blocks)
  n_comm <- scbm_validate_n_comm(n_comm, s)
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
      comp_blocks[[ff]][[mm]] <- scbm_completed_from_mask(blocks[[mm]], keep_mask, rankK)
    }
  }

  score_alpha <- numeric(length(alpha_grid))
  for (aa in seq_along(alpha_grid)) {
    alpha_use <- alpha_grid[aa]
    score_fold <- numeric(folds)
    for (ff in seq_len(folds)) {
      fit <- scbm_cluster_blocks(
        comp_blocks[[ff]],
        n_comm = n_comm,
        alpha = alpha_use,
        circular = circular,
        nstart = nstart
      )
      sc <- 0
      for (mm in seq_len(s)) {
        Mcomp <- comp_blocks[[ff]][[mm]]
        Mhat <- scbm_predict_from_labels(Mcomp, fit$group_L[[mm]], fit$group_R[[mm]])
        hold_mask <- 1 - heldout_masks[[ff]][[mm]]
        if (criterion == "holdout") {
          diff <- (blocks[[mm]] - Mhat) * hold_mask
          sc <- sc + mean(diff[hold_mask == 1]^2, na.rm = TRUE)
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
# Flow / export helpers
# ------------------------------
scbm_make_flow_table <- function(group_L, group_R, send_names, receive_names, row_names = NULL) {
  s <- length(group_L)
  if (length(group_R) != s) stop("group_L and group_R must have the same length.", call. = FALSE)
  if (length(send_names) != s || length(receive_names) != s) {
    stop("send_names and receive_names must both have length s.", call. = FALSE)
  }
  q <- length(group_L[[1L]])
  out <- matrix(NA_integer_, nrow = q, ncol = 2L * s)
  cols <- character(2L * s)
  for (mm in seq_len(s)) {
    out[, 2L * mm - 1L] <- as.integer(group_L[[mm]])
    out[, 2L * mm] <- as.integer(group_R[[mm]])
    cols[2L * mm - 1L] <- send_names[mm]
    cols[2L * mm] <- receive_names[mm]
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- cols
  rn <- if (is.null(row_names)) paste0("V", seq_len(q)) else as.character(row_names)
  rownames(out) <- rn
  out
}

scbm_make_path_summary <- function(flow_table) {
  x <- as.data.frame(flow_table, stringsAsFactors = FALSE)
  nm <- rownames(x)
  labs <- as.matrix(x)
  n_comms <- apply(labs, 1L, function(z) length(unique(z[!is.na(z)])))
  n_switch <- apply(labs, 1L, function(z) {
    z <- z[!is.na(z)]
    if (length(z) <= 1L) return(0L)
    sum(z[-1L] != z[-length(z)])
  })
  data.frame(series = nm, n_communities = n_comms, n_switches = n_switch, stringsAsFactors = FALSE)
}

scbm_make_sankey_df <- function(flow_table, cycle = FALSE) {
  x <- as.data.frame(flow_table, stringsAsFactors = FALSE)
  stage_names <- colnames(x)
  q <- nrow(x)
  ids <- rownames(x)
  pieces <- vector("list", if (isTRUE(cycle)) length(stage_names) else (length(stage_names) - 1L))
  for (jj in seq_along(pieces)) {
    j2 <- if (jj < length(stage_names)) (jj + 1L) else 1L
    pieces[[jj]] <- data.frame(
      series = ids,
      x = stage_names[jj],
      node = x[[jj]],
      next_x = stage_names[j2],
      next_node = x[[j2]],
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

scbm_reorder_block <- function(block, send_labels, receive_labels, row_names = NULL, col_names = NULL) {
  B <- scbm_safe_matrix(block)
  r_ord <- order(as.integer(send_labels), seq_along(send_labels))
  c_ord <- order(as.integer(receive_labels), seq_along(receive_labels))
  out <- B[r_ord, c_ord, drop = FALSE]
  if (!is.null(row_names)) rownames(out) <- as.character(row_names)[r_ord]
  if (!is.null(col_names)) colnames(out) <- as.character(col_names)[c_ord]
  list(matrix = out, row_order = r_ord, col_order = c_ord)
}

# ------------------------------
# PVAR: fitting helpers
# ------------------------------
scbm_pvar_center_series <- function(Yt, s = 4L, centerTF = TRUE) {
  Yt <- scbm_safe_matrix(Yt)
  if (!isTRUE(centerTF)) return(list(Yc = Yt, mu = NULL))
  q <- nrow(Yt)
  TT <- ncol(Yt)
  mu <- matrix(0, q, s)
  for (j in seq_len(s)) {
    idx <- seq(j, TT, by = s)
    mu[, j] <- rowMeans(Yt[, idx, drop = FALSE])
  }
  season_id <- ((seq_len(TT) - 1L) %% s) + 1L
  list(Yc = Yt - mu[, season_id, drop = FALSE], mu = mu)
}

scbm_pvar_get_season_design <- function(Yc, s = 4L, p = 1L, season = 1L) {
  Yc <- scbm_safe_matrix(Yc)
  q <- nrow(Yc)
  TT <- ncol(Yc)
  ids <- which(((seq_len(TT) - 1L) %% s) + 1L == as.integer(season))
  ids <- ids[ids > as.integer(p)]
  if (length(ids) <= 1L) {
    stop("Not enough observations for season ", season, ".", call. = FALSE)
  }
  X <- t(vapply(
    ids,
    FUN.VALUE = numeric(q * p),
    FUN = function(tt) as.vector(Yc[, (tt - 1L):(tt - p), drop = FALSE])
  ))
  Y <- t(Yc[, ids, drop = FALSE])
  list(X = X, Y = Y, id = ids, n = length(ids))
}

scbm_pvar_season_sample_sizes <- function(TT, s = 4L, p = 1L) {
  vapply(seq_len(s), function(mm) {
    ids <- which(((seq_len(TT) - 1L) %% s) + 1L == mm)
    ids <- ids[ids > as.integer(p)]
    length(ids)
  }, integer(1L))
}

scbm_pvar_weight_matrix <- function(q, p, diagTF = TRUE) {
  W <- matrix(1, nrow = q * p, ncol = q)
  if (isTRUE(diagTF)) {
    for (hh in seq_len(p)) {
      idx <- ((hh - 1L) * q + 1L):(hh * q)
      W[idx, ] <- 1
      W[cbind(idx, seq_len(q))] <- 0
    }
  }
  W
}

scbm_pvar_sigma_inv_from_ols <- function(X, Y, sigma_diag_only = TRUE) {
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

scbm_pvar_lambda_textbook <- function(q, s, p, n_eff, theory_const = 1) {
  as.numeric(theory_const) * sqrt(log(as.numeric(s) * as.numeric(p) * q^2) / as.numeric(n_eff))
}

scbm_pvar_lambda_pkg_from_text <- function(lambda_text, n_eff) {
  0.5 * as.numeric(n_eff) * as.numeric(lambda_text)
}

scbm_pvar_lambda_text_from_pkg <- function(lambda_pkg, n_eff) {
  if (any(!is.finite(lambda_pkg)) || any(!is.finite(n_eff)) || any(n_eff <= 0)) {
    out <- rep(NA_real_, max(length(lambda_pkg), length(n_eff)))
    return(out)
  }
  2 * as.numeric(lambda_pkg) / as.numeric(n_eff)
}

scbm_pvar_fit_ols <- function(Yt, s = 4L, p = 1L, centerTF = TRUE, ridge = 1e-8) {
  Yt <- scbm_safe_matrix(Yt)
  q <- nrow(Yt)
  TT <- ncol(Yt)
  cen <- scbm_pvar_center_series(Yt, s = s, centerTF = centerTF)
  Yc <- cen$Yc

  lag_blocks <- vector("list", s)
  blocks <- vector("list", s)
  Phi_hat_full <- matrix(NA_real_, nrow = q, ncol = q * p * s)

  for (mm in seq_len(s)) {
    des <- scbm_pvar_get_season_design(Yc, s = s, p = p, season = mm)
    B <- solve(crossprod(des$X) + diag(ridge, q * p), crossprod(des$X, des$Y))
    Phi_season <- t(B)
    lag_blocks[[mm]] <- lapply(seq_len(p), function(hh) {
      Phi_season[, ((hh - 1L) * q + 1L):(hh * q), drop = FALSE]
    })
    blocks[[mm]] <- Reduce(`+`, lag_blocks[[mm]])
    Phi_hat_full[, ((mm - 1L) * q * p + 1L):(mm * q * p)] <- Phi_season
  }

  list(
    Phi_hat = Phi_hat_full,
    Bhat = Phi_hat_full,
    lag_blocks = lag_blocks,
    blocks = blocks,
    mu = cen$mu,
    lambda_pkg = rep(NA_real_, s),
    lambda_text = rep(NA_real_, s),
    c_lambda_sel = NA_real_,
    lambda_mode = "ols"
  )
}

scbm_pvar_fit_lasso_fixed <- function(Yt,
                                      s = 4L,
                                      p = 1L,
                                      lambda_pkg,
                                      diagTF = TRUE,
                                      centerTF = TRUE,
                                      updateSigma = FALSE,
                                      sigma_diag_only = TRUE,
                                      max_iter = 1000L,
                                      tol = 1e-6) {
  Yt <- scbm_safe_matrix(Yt)
  q <- nrow(Yt)
  TT <- ncol(Yt)
  lambda_pkg <- as.numeric(lambda_pkg)
  if (length(lambda_pkg) == 1L) lambda_pkg <- rep(lambda_pkg, s)
  if (length(lambda_pkg) != s) stop("lambda_pkg must have length 1 or s.", call. = FALSE)

  cen <- scbm_pvar_center_series(Yt, s = s, centerTF = centerTF)
  Yc <- cen$Yc
  W0 <- scbm_pvar_weight_matrix(q = q, p = p, diagTF = diagTF)
  Phi_hat_full <- matrix(0, nrow = q, ncol = q * p * s)
  lag_blocks <- vector("list", s)
  blocks <- vector("list", s)

  for (mm in seq_len(s)) {
    des <- scbm_pvar_get_season_design(Yc, s = s, p = p, season = mm)
    SigmaInv <- NULL
    if (isTRUE(updateSigma)) {
      SigmaInv <- scbm_pvar_sigma_inv_from_ols(des$X, des$Y, sigma_diag_only = sigma_diag_only)
    }
    Bhat <- sparseVAR:::fista_lasso_multi_cpp(
      X = des$X,
      Y = des$Y,
      weights = W0,
      lambda = lambda_pkg[mm],
      SigmaInv = if (is.null(SigmaInv)) matrix(0, 0, 0) else SigmaInv,
      max_iter = as.integer(max_iter),
      tol = tol
    )
    Phi_season <- t(Bhat)
    lag_blocks[[mm]] <- lapply(seq_len(p), function(hh) {
      Phi_season[, ((hh - 1L) * q + 1L):(hh * q), drop = FALSE]
    })
    blocks[[mm]] <- Reduce(`+`, lag_blocks[[mm]])
    Phi_hat_full[, ((mm - 1L) * q * p + 1L):(mm * q * p)] <- Phi_season
  }

  n_eff <- scbm_pvar_season_sample_sizes(TT, s = s, p = p)
  lambda_text <- scbm_pvar_lambda_text_from_pkg(lambda_pkg, n_eff)
  list(
    Phi_hat = Phi_hat_full,
    Bhat = Phi_hat_full,
    lag_blocks = lag_blocks,
    blocks = blocks,
    mu = cen$mu,
    lambda_pkg = lambda_pkg,
    lambda_text = lambda_text
  )
}

scbm_pvar_cv_score_one_c <- function(Yc,
                                     s = 4L,
                                     p = 1L,
                                     c_value,
                                     q,
                                     fold = 10L,
                                     diagTF = TRUE,
                                     updateSigma = FALSE,
                                     sigma_diag_only = TRUE,
                                     max_iter = 1000L,
                                     tol = 1e-6) {
  total <- 0
  W0 <- scbm_pvar_weight_matrix(q = q, p = p, diagTF = diagTF)

  for (mm in seq_len(s)) {
    des <- scbm_pvar_get_season_design(Yc, s = s, p = p, season = mm)
    n_m <- nrow(des$X)
    lambda_text_base <- scbm_pvar_lambda_textbook(q = q, s = s, p = p, n_eff = n_m, theory_const = 1)
    lambda_pkg_base <- scbm_pvar_lambda_pkg_from_text(lambda_text_base, n_eff = n_m)
    lambda_use <- as.numeric(c_value) * as.numeric(lambda_pkg_base)

    SigmaInv <- NULL
    if (isTRUE(updateSigma)) {
      SigmaInv <- scbm_pvar_sigma_inv_from_ols(des$X, des$Y, sigma_diag_only = sigma_diag_only)
    }

    fold_use <- max(2L, min(as.integer(fold), n_m))
    folds <- rep(seq_len(fold_use), each = floor(n_m / fold_use) + 1L)[1L:n_m]
    for (ff in seq_len(fold_use)) {
      tr <- which(folds != ff)
      va <- which(folds == ff)
      Bhat <- sparseVAR:::fista_lasso_multi_cpp(
        X = des$X[tr, , drop = FALSE],
        Y = des$Y[tr, , drop = FALSE],
        weights = W0,
        lambda = lambda_use,
        SigmaInv = if (is.null(SigmaInv)) matrix(0, 0, 0) else SigmaInv,
        max_iter = as.integer(max_iter),
        tol = tol
      )
      Rv <- des$Y[va, , drop = FALSE] - des$X[va, , drop = FALSE] %*% Bhat
      if (is.null(SigmaInv)) {
        total <- total + sum(Rv^2)
      } else {
        total <- total + sum((Rv %*% SigmaInv) * Rv)
      }
    }
  }

  total
}

scbm_pvar_select_c_lambda_blockcv <- function(Yt,
                                              s = 4L,
                                              p = 1L,
                                              c_grid = scbm_default_c_grid(),
                                              fold = 10L,
                                              diagTF = TRUE,
                                              centerTF = TRUE,
                                              updateSigma = FALSE,
                                              sigma_diag_only = TRUE,
                                              max_iter = 1000L,
                                              tol = 1e-6) {
  Yt <- scbm_safe_matrix(Yt)
  q <- nrow(Yt)
  c_grid <- unique(as.numeric(c_grid))
  c_grid <- c_grid[is.finite(c_grid) & c_grid > 0]
  if (!length(c_grid)) stop("c_grid must contain positive finite values.", call. = FALSE)
  c_grid <- sort(c_grid)

  cen <- scbm_pvar_center_series(Yt, s = s, centerTF = centerTF)
  scores <- vapply(c_grid, function(cc) {
    scbm_pvar_cv_score_one_c(
      Yc = cen$Yc,
      s = s,
      p = p,
      c_value = cc,
      q = q,
      fold = fold,
      diagTF = diagTF,
      updateSigma = updateSigma,
      sigma_diag_only = sigma_diag_only,
      max_iter = max_iter,
      tol = tol
    )
  }, numeric(1L))

  id_best <- which.min(scores)
  c_sel <- c_grid[id_best]
  n_eff <- scbm_pvar_season_sample_sizes(ncol(Yt), s = s, p = p)
  lambda_text_base <- scbm_pvar_lambda_textbook(q = q, s = s, p = p, n_eff = n_eff, theory_const = 1)
  lambda_pkg_base <- scbm_pvar_lambda_pkg_from_text(lambda_text_base, n_eff = n_eff)

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

scbm_pvar_fit_lasso <- function(Yt,
                                s = 4L,
                                p = 1L,
                                lambda_mode = c("theory", "cv_c"),
                                c_lambda = 0.20,
                                c_grid = scbm_default_c_grid(),
                                fold = 10L,
                                diagTF = TRUE,
                                centerTF = TRUE,
                                updateSigma = FALSE,
                                sigma_diag_only = TRUE,
                                max_iter = 1000L,
                                tol = 1e-6) {
  lambda_mode <- match.arg(lambda_mode)
  Yt <- scbm_safe_matrix(Yt)
  q <- nrow(Yt)
  TT <- ncol(Yt)
  n_eff <- scbm_pvar_season_sample_sizes(TT, s = s, p = p)

  if (any(n_eff <= 0)) stop("Need enough observations in every season after lag p.", call. = FALSE)

  if (lambda_mode == "cv_c") {
    sel <- scbm_pvar_select_c_lambda_blockcv(
      Yt = Yt,
      s = s,
      p = p,
      c_grid = c_grid,
      fold = fold,
      diagTF = diagTF,
      centerTF = centerTF,
      updateSigma = updateSigma,
      sigma_diag_only = sigma_diag_only,
      max_iter = max_iter,
      tol = tol
    )
    c_sel <- sel$c_lambda_sel
    lambda_pkg <- sel$lambda_pkg_sel
  } else {
    c_sel <- as.numeric(c_lambda)[1L]
    lambda_text_base <- scbm_pvar_lambda_textbook(q = q, s = s, p = p, n_eff = n_eff, theory_const = 1)
    lambda_pkg_base <- scbm_pvar_lambda_pkg_from_text(lambda_text_base, n_eff = n_eff)
    lambda_pkg <- as.numeric(c_sel) * lambda_pkg_base
  }

  fixed_fit <- scbm_pvar_fit_lasso_fixed(
    Yt = Yt,
    s = s,
    p = p,
    lambda_pkg = lambda_pkg,
    diagTF = diagTF,
    centerTF = centerTF,
    updateSigma = updateSigma,
    sigma_diag_only = sigma_diag_only,
    max_iter = max_iter,
    tol = tol
  )

  list(
    Phi_hat = fixed_fit$Phi_hat,
    Bhat = fixed_fit$Bhat,
    lag_blocks = fixed_fit$lag_blocks,
    blocks = fixed_fit$blocks,
    mu = fixed_fit$mu,
    lambda_pkg = as.numeric(fixed_fit$lambda_pkg),
    lambda_text = as.numeric(fixed_fit$lambda_text),
    c_lambda_sel = as.numeric(c_sel),
    lambda_mode = lambda_mode
  )
}

# ------------------------------
# VHAR: fitting helpers
# ------------------------------
scbm_vhar_split_phi_blocks <- function(Bhat, q) {
  Bhat <- scbm_safe_matrix(Bhat)
  list(
    scbm_safe_matrix(Bhat[, 1:q, drop = FALSE]),
    scbm_safe_matrix(Bhat[, (q + 1L):(2L * q), drop = FALSE]),
    scbm_safe_matrix(Bhat[, (2L * q + 1L):(3L * q), drop = FALSE])
  )
}

scbm_vhar_lambda_textbook <- function(q, N_eff, theory_const = 1) {
  as.numeric(theory_const) * sqrt(log(3 * q^2) / as.numeric(N_eff))
}

scbm_vhar_lambda_pkg_from_text <- function(lambda_text, N_eff) {
  0.5 * as.numeric(N_eff) * as.numeric(lambda_text)
}

scbm_vhar_lambda_text_from_pkg <- function(lambda_pkg, N_eff) {
  if (any(!is.finite(lambda_pkg)) || !is.finite(N_eff) || N_eff <= 0) {
    out <- rep(NA_real_, length(lambda_pkg))
    return(out)
  }
  2 * as.numeric(lambda_pkg) / as.numeric(N_eff)
}

scbm_vhar_fit_ols <- function(Yt, bw = 5L, bm = 22L) {
  Yt <- scbm_safe_matrix(Yt)
  fit <- sparseVAR::VHAR_ols(Yt = Yt, bd = as.integer(bw), bm = as.integer(bm))
  Bhat <- scbm_safe_matrix(fit$Phi_hat)
  q <- nrow(Yt)
  list(
    Phi_hat = Bhat,
    Bhat = Bhat,
    blocks = scbm_vhar_split_phi_blocks(Bhat, q = q),
    lambda_pkg = rep(NA_real_, 3L),
    lambda_text = rep(NA_real_, 3L),
    c_lambda_sel = NA_real_,
    lambda_mode = "ols"
  )
}

scbm_vhar_select_c_lambda_blockcv <- function(Yt,
                                              bw = 5L,
                                              bm = 22L,
                                              c_grid = scbm_default_c_grid(),
                                              fold = 10L,
                                              max_iter = 1000L,
                                              tol = 1e-6,
                                              diagTF = TRUE,
                                              updateSigma = FALSE,
                                              sigma_diag_only = TRUE) {
  Yt <- scbm_safe_matrix(Yt)
  q <- nrow(Yt)
  N_eff <- ncol(Yt) - as.integer(bm)
  if (N_eff <= 0L) stop("Need ncol(Yt) > bm for c-lambda CV.", call. = FALSE)

  c_grid <- unique(as.numeric(c_grid))
  c_grid <- c_grid[is.finite(c_grid) & c_grid > 0]
  if (!length(c_grid)) stop("c_grid must contain positive finite values.", call. = FALSE)
  c_grid <- sort(c_grid)

  lambda_text_base <- scbm_vhar_lambda_textbook(q = q, N_eff = N_eff, theory_const = 1)
  lambda_pkg_base <- scbm_vhar_lambda_pkg_from_text(lambda_text_base, N_eff = N_eff)
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
  lambda_text_sel <- scbm_vhar_lambda_text_from_pkg(lambda_sel, N_eff = N_eff)

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

scbm_vhar_fit_lasso <- function(Yt,
                                bw = 5L,
                                bm = 22L,
                                lambda_mode = c("theory", "cv", "cv_c"),
                                c_lambda = 0.25,
                                c_grid = scbm_default_c_grid(),
                                fold = 10L,
                                nlambda = 100L,
                                max_iter = 1000L,
                                tol = 1e-6,
                                diagTF = TRUE,
                                updateSigma = FALSE,
                                sigma_diag_only = TRUE) {
  lambda_mode <- match.arg(lambda_mode)
  Yt <- scbm_safe_matrix(Yt)
  q <- nrow(Yt)
  N_eff <- ncol(Yt) - as.integer(bm)
  if (N_eff <= 0L) stop("Need ncol(Yt) > bm.", call. = FALSE)

  if (lambda_mode == "theory") {
    lambda_text <- scbm_vhar_lambda_textbook(q = q, N_eff = N_eff, theory_const = c_lambda)
    lambda_pkg <- scbm_vhar_lambda_pkg_from_text(lambda_text, N_eff = N_eff)
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
    lambda_text <- scbm_vhar_lambda_text_from_pkg(lambda_pkg, N_eff = N_eff)
    lambda_text_base <- scbm_vhar_lambda_textbook(q = q, N_eff = N_eff, theory_const = 1)
    c_sel <- scbm_mean(lambda_text / lambda_text_base)
  } else {
    sel <- scbm_vhar_select_c_lambda_blockcv(
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
  Braw <- scbm_safe_matrix(Braw)

  list(
    Phi_hat = Braw,
    Bhat = Braw,
    blocks = scbm_vhar_split_phi_blocks(Braw, q = q),
    lambda_pkg = rep(as.numeric(lambda_pkg), 3L),
    lambda_text = rep(as.numeric(lambda_text), 3L),
    c_lambda_sel = as.numeric(c_sel),
    lambda_mode = lambda_mode
  )
}

# ------------------------------
# VHAR aggregate clustering helpers
# ------------------------------
scbm_vhar_parse_aggregate_n_comm <- function(n_comm) {
  n_comm <- as.integer(n_comm)
  if (length(n_comm) == 3L) {
    if (any(!is.finite(n_comm)) || any(n_comm < 1L)) {
      stop("For aggregated VHAR clustering, n_comm must contain positive integers.", call. = FALSE)
    }
    k_short <- n_comm[1L]
    k_middle <- n_comm[2L]
    k_long <- n_comm[3L]
  } else if (length(n_comm) == 6L) {
    if (any(!is.finite(n_comm)) || any(n_comm < 1L)) {
      stop("For aggregated VHAR clustering, n_comm must contain positive integers.", call. = FALSE)
    }
    if (n_comm[2L] != n_comm[3L]) {
      stop("For aggregated VHAR clustering, n_comm[2] and n_comm[3] must be equal.", call. = FALSE)
    }
    if (length(unique(n_comm[4L:6L])) != 1L) {
      stop("For aggregated VHAR clustering, n_comm[4:6] must be equal.", call. = FALSE)
    }
    k_short <- n_comm[1L]
    k_middle <- n_comm[2L]
    k_long <- n_comm[4L]
  } else {
    stop(
      "For VHAR, n_comm must have length 3 in the order c(DS, DR_WS, WR_MS_MR), or a consistent length-6 vector.",
      call. = FALSE
    )
  }

  expanded <- as.integer(c(k_short, k_middle, k_middle, k_long, k_long, k_long))
  list(
    short = as.integer(k_short),
    middle = as.integer(k_middle),
    long = as.integer(k_long),
    input = as.integer(c(k_short, k_middle, k_long)),
    expanded = expanded
  )
}

scbm_cluster_vhar_aggregated <- function(blocks,
                                         n_comm,
                                         alpha = 0,
                                         nstart = 50L,
                                         tol = 1e-5,
                                         max_iter = 1000L) {
  if (length(blocks) != 3L) stop("VHAR aggregation expects exactly 3 horizon blocks.", call. = FALSE)
  nc <- scbm_vhar_parse_aggregate_n_comm(n_comm)
  tmp <- scbm_embeddings_with_projectors(blocks, nc$expanded)

  if (alpha > 0) {
    ranks_L <- nc$expanded[c(1L, 3L, 5L)]
    ranks_R <- nc$expanded[c(2L, 4L, 6L)]
    smL <- scbm_pisces_smooth(tmp$proj_L, ranks = ranks_L, alpha = alpha, tol = tol, max_iter = max_iter)
    smR <- scbm_pisces_smooth(tmp$proj_R, ranks = ranks_R, alpha = alpha, tol = tol, max_iter = max_iter)
    emb_L <- scbm_embeddings_from_projectors(smL$projectors, ranks_L)
    emb_R <- scbm_embeddings_from_projectors(smR$projectors, ranks_R)
  } else {
    emb_L <- tmp$emb_L
    emb_R <- tmp$emb_R
  }

  lab_short <- scbm_kmeans_labels(emb_L[[1L]], nc$short, nstart = nstart)
  lab_middle <- scbm_kmeans_labels(cbind(emb_R[[1L]], emb_L[[2L]]), nc$middle, nstart = nstart)
  lab_long <- scbm_kmeans_labels(cbind(emb_R[[2L]], emb_L[[3L]], emb_R[[3L]]), nc$long, nstart = nstart)

  group_L <- list(lab_short, lab_middle, lab_long)
  group_R <- list(lab_middle, lab_long, lab_long)

  list(
    group_L = group_L,
    group_R = group_R,
    emb_L = emb_L,
    emb_R = emb_R,
    alpha = as.numeric(alpha),
    circular = FALSE,
    aggregate_labels = list(short = lab_short, middle = lab_middle, long = lab_long),
    n_comm_agg = nc$input,
    n_comm_expanded = nc$expanded
  )
}

scbm_make_vhar_flow_outputs <- function(community_fit,
                                        horizon_names = c("daily", "weekly", "monthly"),
                                        row_names = NULL) {
  if (length(horizon_names) != 3L) stop("horizon_names must have length 3.", call. = FALSE)
  q <- length(community_fit$group_L[[1L]])
  rn <- if (is.null(row_names)) paste0("V", seq_len(q)) else as.character(row_names)

  flow_table <- data.frame(
    long.cluster = as.integer(community_fit$aggregate_labels$long),
    middle.cluster = as.integer(community_fit$aggregate_labels$middle),
    short.cluster = as.integer(community_fit$aggregate_labels$short),
    stringsAsFactors = FALSE,
    row.names = rn
  )

  flow_table_full <- data.frame(
    stringsAsFactors = FALSE,
    row.names = rn,
    check.names = FALSE
  )
  flow_table_full[[paste0(horizon_names[3L], ".receiving")]] <- as.integer(community_fit$group_R[[3L]])
  flow_table_full[[paste0(horizon_names[3L], ".sending")]] <- as.integer(community_fit$group_L[[3L]])
  flow_table_full[[paste0(horizon_names[2L], ".receiving")]] <- as.integer(community_fit$group_R[[2L]])
  flow_table_full[[paste0(horizon_names[2L], ".sending")]] <- as.integer(community_fit$group_L[[2L]])
  flow_table_full[[paste0(horizon_names[1L], ".receiving")]] <- as.integer(community_fit$group_R[[1L]])
  flow_table_full[[paste0(horizon_names[1L], ".sending")]] <- as.integer(community_fit$group_L[[1L]])

  aggregate_definition <- data.frame(
    aggregate_stage = c("long.cluster", "middle.cluster", "short.cluster"),
    original_stages = c(
      paste0(horizon_names[3L], ".receiving + ", horizon_names[3L], ".sending + ", horizon_names[2L], ".receiving"),
      paste0(horizon_names[2L], ".sending + ", horizon_names[1L], ".receiving"),
      paste0(horizon_names[1L], ".sending")
    ),
    stringsAsFactors = FALSE
  )

  list(
    flow_table = flow_table,
    flow_table_full = flow_table_full,
    aggregate_definition = aggregate_definition
  )
}

# ------------------------------
# Top-level empirical wrappers
# ------------------------------
scbm_pvar_empirical_fit <- function(data,
                                s = 4L,
                                p = 1L,
                                series_in_rows = FALSE,
                                series_names = NULL,
                                estimator = c("ols", "lasso"),
                                lambda_mode = c("cv_c", "theory"),
                                c_lambda = 0.20,
                                c_grid = scbm_default_c_grid(),
                                lasso_fold = 10L,
                                lasso_max_iter = 1000L,
                                lasso_tol = 1e-6,
                                diagTF = TRUE,
                                centerTF = TRUE,
                                updateSigma = FALSE,
                                sigma_diag_only = TRUE) {
  estimator <- match.arg(estimator)
  lambda_mode <- match.arg(lambda_mode)

  Yt <- scbm_as_Yt(data, series_in_rows = series_in_rows)
  q <- nrow(Yt)
  TT <- ncol(Yt)
  s <- as.integer(s)
  p <- as.integer(p)
  if (TT <= p) stop("Need ncol(Yt) > p.", call. = FALSE)
  if (is.null(series_names)) series_names <- colnames(t(Yt))
  if (is.null(series_names)) series_names <- paste0("V", seq_len(q))

  fit_first <- if (estimator == "ols") {
    scbm_pvar_fit_ols(Yt, s = s, p = p, centerTF = centerTF)
  } else {
    scbm_pvar_fit_lasso(
      Yt,
      s = s,
      p = p,
      lambda_mode = lambda_mode,
      c_lambda = c_lambda,
      c_grid = c_grid,
      fold = lasso_fold,
      diagTF = diagTF,
      centerTF = centerTF,
      updateSigma = updateSigma,
      sigma_diag_only = sigma_diag_only,
      max_iter = lasso_max_iter,
      tol = lasso_tol
    )
  }
  names(fit_first$blocks) <- paste0("season", seq_len(s))
  if (!is.null(fit_first$lag_blocks)) names(fit_first$lag_blocks) <- paste0("season", seq_len(s))

  out <- list(
    model = "PVAR",
    estimator = estimator,
    Yt = Yt,
    q = q,
    TT = TT,
    s = s,
    p = p,
    first_stage = fit_first,
    blocks = fit_first$blocks,
    lag_blocks = fit_first$lag_blocks,
    series_names = as.character(series_names),
    scree = scbm_block_scree(fit_first$blocks, block_names = paste0("season", seq_len(s)))
  )
  class(out) <- c("scbm_pvar_empirical", class(out))
  out
}

scbm_pvar_empirical_cluster <- function(fit,
                                        s = 4L,
                                        p = 1L,
                                        series_names = NULL,
                                        n_comm,
                                        alpha_mode = c("cv", "fixed", "none"),
                                        alpha = NULL,
                                        alpha_fold = 5L,
                                        alpha_grid = scbm_default_alpha_grid(),
                                        alpha_criterion = c("holdout", "paper"),
                                        nstart = 50L,
                                        seed = NULL) {

  alpha_mode <- match.arg(alpha_mode)
  alpha_criterion <- match.arg(alpha_criterion)
  if (!is.null(seed)) set.seed(as.integer(seed))

  q <- fit$q
  s <- as.integer(s)
  p <- as.integer(p)
  n_comm <- scbm_validate_n_comm(n_comm, s)
  if (is.null(series_names)) series_names <- fit$series_names
  if (is.null(series_names)) series_names <- paste0("V", seq_len(q))

  names(fit$first_stage$blocks) <- paste0("season", seq_len(s))
  if (!is.null(fit$first_stage$lag_blocks)) names(fit$first_stage$lag_blocks) <- paste0("season", seq_len(s))

  alpha_sel <- 0
  alpha_cv <- NULL
  if (alpha_mode == "fixed") {
    if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
      stop("For alpha_mode = 'fixed', alpha must be one nonnegative number.", call. = FALSE)
    }
    alpha_sel <- as.numeric(alpha)
  } else if (alpha_mode == "cv") {
    alpha_cv <- scbm_select_alpha(
      blocks = fit$first_stage$blocks,
      n_comm = n_comm,
      circular = TRUE,
      folds = alpha_fold,
      alpha_grid = alpha_grid,
      criterion = alpha_criterion,
      nstart = nstart,
      seed = if (is.null(seed)) NULL else as.integer(seed) + 101L
    )
    alpha_sel <- alpha_cv$alpha
  }

  community_fit <- scbm_cluster_blocks(
    fit$first_stage$blocks,
    n_comm = n_comm,
    alpha = alpha_sel,
    circular = TRUE,
    nstart = nstart
  )

  send_names <- paste0("S", seq_len(s))
  receive_names <- paste0("S", seq_len(s))
  flow_table <- scbm_make_flow_table(
    community_fit$group_L,
    community_fit$group_R,
    send_names = paste0(send_names, ".sending"),
    receive_names = paste0(receive_names, ".receiving"),
    row_names = fit$series_names
  )

  reordered_blocks <- lapply(seq_len(s), function(mm) {
    scbm_reorder_block(
      fit$first_stage$blocks[[mm]],
      community_fit$group_L[[mm]],
      community_fit$group_R[[mm]],
      row_names = series_names,
      col_names = series_names
    )
  })

  out <- list(
    model = "PVAR",
    n_comm = n_comm,
    alpha_mode = alpha_mode,
    alpha_sel = as.numeric(alpha_sel),
    alpha_cv = alpha_cv,
    community_fit = community_fit,
    flow_table = flow_table,
    path_summary = scbm_make_path_summary(flow_table),
    sankey_df = scbm_make_sankey_df(flow_table, cycle = FALSE),
    reordered_blocks = reordered_blocks,
    series_names = as.character(series_names)
  )
  class(out) <- c("scbm_pvar_empirical", class(out))
  out
}


scbm_vhar_empirical_fit <- function(data,
                                bw = 5L,
                                bm = 22L,
                                horizon_names = c("daily", "weekly", "monthly"),
                                estimator = c("ols", "lasso"),
                                series_in_rows = FALSE,
                                series_names = NULL,
                                lambda_mode = c("cv_c", "cv", "theory"),
                                c_lambda = 0.25,
                                c_grid = scbm_default_c_grid(),
                                lasso_fold = 10L,
                                lasso_nlambda = 100L,
                                lasso_max_iter = 1000L,
                                lasso_tol = 1e-6,
                                diagTF = TRUE,
                                updateSigma = FALSE,
                                sigma_diag_only = TRUE) {
  
  estimator <- match.arg(estimator)
  lambda_mode <- match.arg(lambda_mode)

  Yt <- scbm_as_Yt(data, series_in_rows = series_in_rows)
  q <- nrow(Yt)
  TT <- ncol(Yt)
  bw <- as.integer(bw)
  bm <- as.integer(bm)
  if (bw <= 1L || bm <= bw) stop("Need 1 < bw < bm.", call. = FALSE)
  if (TT <= bm) stop("Need ncol(Yt) > bm.", call. = FALSE)
  if (length(horizon_names) != 3L) stop("horizon_names must have length 3.", call. = FALSE)
  if (is.null(series_names)) series_names <- colnames(t(Yt))
  if (is.null(series_names)) series_names <- paste0("V", seq_len(q))

  fit_first <- if (estimator == "ols") {
    scbm_vhar_fit_ols(Yt, bw = bw, bm = bm)
  } else {
    scbm_vhar_fit_lasso(
      Yt,
      bw = bw,
      bm = bm,
      lambda_mode = lambda_mode,
      c_lambda = c_lambda,
      c_grid = c_grid,
      fold = lasso_fold,
      nlambda = lasso_nlambda,
      max_iter = lasso_max_iter,
      tol = lasso_tol,
      diagTF = diagTF,
      updateSigma = updateSigma,
      sigma_diag_only = sigma_diag_only
    )
  }
  names(fit_first$blocks) <- as.character(horizon_names)

  out <- list(
    model = "VHAR",
    estimator = estimator,
    Yt = Yt,
    q = q,
    TT = TT,
    bw = bw,
    bm = bm,
    first_stage = fit_first,
    blocks = fit_first$blocks,
    scree = scbm_block_scree(fit_first$blocks, block_names = as.character(horizon_names)),
    horizon_names = as.character(horizon_names),
    series_names = as.character(series_names)
  )
  class(out) <- c("scbm_vhar_empirical", class(out))
  out
}


scbm_vhar_empirical_cluster <- function(fit,
                                    bw = 5L,
                                    bm = 22L,
                                    horizon_names = c("daily", "weekly", "monthly"),
                                    series_names = NULL,
                                    n_comm,
                                    alpha_mode = c("cv", "fixed", "none"),
                                    alpha = NULL,
                                    alpha_fold = 5L,
                                    alpha_grid = scbm_default_alpha_grid(),
                                    alpha_criterion = c("holdout", "paper"),
                                    nstart = 50L,
                                    seed = 12345) {

  alpha_mode <- match.arg(alpha_mode)
  alpha_criterion <- match.arg(alpha_criterion)
  if (!is.null(seed)) set.seed(as.integer(seed))

  q <- fit$q
  bw <- as.integer(if (missing(bw)) fit$bw else bw)
  bm <- as.integer(if (missing(bm)) fit$bm else bm)
  if (bw <= 1L || bm <= bw) stop("Need 1 < bw < bm.", call. = FALSE)
  if (length(horizon_names) != 3L) stop("horizon_names must have length 3.", call. = FALSE)
  if (is.null(series_names)) series_names <- fit$series_names
  if (is.null(series_names)) series_names <- paste0("V", seq_len(q))

  n_comm_info <- scbm_vhar_parse_aggregate_n_comm(n_comm)
  n_comm_expanded <- n_comm_info$expanded

  alpha_sel <- 0
  alpha_cv <- NULL
  if (alpha_mode == "fixed") {
    if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
      stop("For alpha_mode = 'fixed', alpha must be one nonnegative number.", call. = FALSE)
    }
    alpha_sel <- as.numeric(alpha)
  } else if (alpha_mode == "cv") {
    alpha_cv <- scbm_select_alpha(
      blocks = fit$first_stage$blocks,
      n_comm = n_comm_expanded,
      circular = FALSE,
      folds = alpha_fold,
      alpha_grid = alpha_grid,
      criterion = alpha_criterion,
      nstart = nstart,
      seed = if (is.null(seed)) NULL else as.integer(seed) + 202L
    )
    alpha_sel <- alpha_cv$alpha
  }

  community_fit <- scbm_cluster_vhar_aggregated(
    fit$first_stage$blocks,
    n_comm = n_comm_info$input,
    alpha = alpha_sel,
    nstart = nstart
  )

  flow_outputs <- scbm_make_vhar_flow_outputs(
    community_fit = community_fit,
    horizon_names = horizon_names,
    row_names = series_names
  )

  reordered_blocks <- list(
    monthly = scbm_reorder_block(
      fit$first_stage$blocks[[3L]],
      community_fit$group_L[[3L]],
      community_fit$group_R[[3L]],
      row_names = series_names,
      col_names = series_names
    ),
    weekly = scbm_reorder_block(
      fit$first_stage$blocks[[2L]],
      community_fit$group_L[[2L]],
      community_fit$group_R[[2L]],
      row_names = series_names,
      col_names = series_names
    ),
    daily = scbm_reorder_block(
      fit$first_stage$blocks[[1L]],
      community_fit$group_L[[1L]],
      community_fit$group_R[[1L]],
      row_names = series_names,
      col_names = series_names
    )
  )

  out <- list(
    model = "VHAR",
    bw = bw,
    bm = bm,
    n_comm = as.integer(n_comm_info$input),
    n_comm_expanded = as.integer(n_comm_expanded),
    alpha_mode = alpha_mode,
    alpha_sel = as.numeric(alpha_sel),
    alpha_cv = alpha_cv,
    community_fit = community_fit,
    aggregate_definition = flow_outputs$aggregate_definition,
    flow_table = flow_outputs$flow_table,
    flow_table_full = flow_outputs$flow_table_full,
    path_summary = scbm_make_path_summary(flow_outputs$flow_table),
    sankey_df = scbm_make_sankey_df(flow_outputs$flow_table, cycle = FALSE),
    reordered_blocks = reordered_blocks,
    series_names = as.character(series_names),
    horizon_names = as.character(horizon_names)
  )
  class(out) <- c("scbm_vhar_empirical", class(out))
  out
}


# ------------------------------
# Export helper
# ------------------------------
scbm_fit_empirical_export <- function(result, prefix = "scbm_empirical_output") {
  if (missing(result) || is.null(result)) stop("Provide a fitted scbm_empirical_* object.", call. = FALSE)
  prefix <- as.character(prefix)[1L]

  save(list = c("result"), file = paste0(prefix, ".RData"))
  utils::write.csv(result$scree, paste0(prefix, "_scree.csv"), row.names = FALSE)

  for (mm in seq_along(result$blocks)) {
    nm <- if (!is.null(names(result$blocks)[mm]) && nzchar(names(result$blocks)[mm])) names(result$blocks)[mm] else paste0("block", mm)
    utils::write.csv(result$blocks[[mm]], paste0(prefix, "_", nm, "_block.csv"), row.names = TRUE)
  }

  if (!is.null(result$lag_blocks)) {
    for (mm in seq_along(result$lag_blocks)) {
      for (hh in seq_along(result$lag_blocks[[mm]])) {
        utils::write.csv(
          result$lag_blocks[[mm]][[hh]],
          paste0(prefix, "_season", mm, "_lag", hh, ".csv"),
          row.names = TRUE
        )
      }
    }
  }

  invisible(prefix)
}


scbm_cluster_empirical_export <- function(result, prefix = "scbm_empirical_output") {
  if (missing(result) || is.null(result)) stop("Provide a fitted scbm_empirical_* object.", call. = FALSE)
  prefix <- as.character(prefix)[1L]

  save(list = c("result"), file = paste0(prefix, ".RData"))
  utils::write.csv(result$flow_table, paste0(prefix, "_community_flow.csv"), row.names = TRUE)
  if (!is.null(result$flow_table_full)) {
    utils::write.csv(result$flow_table_full, paste0(prefix, "_community_flow_full.csv"), row.names = TRUE)
  }
  if (!is.null(result$aggregate_definition)) {
    utils::write.csv(result$aggregate_definition, paste0(prefix, "_community_definition.csv"), row.names = FALSE)
  }
  utils::write.csv(result$path_summary, paste0(prefix, "_community_summary.csv"), row.names = FALSE)
  utils::write.csv(result$sankey_df, paste0(prefix, "_sankey.csv"), row.names = FALSE)

  invisible(prefix)
}
