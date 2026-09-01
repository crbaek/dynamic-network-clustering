# Generalized ScBM--VHAR simulation and global-volatility application.
# Coefficient blocks are stored short/medium/long, while clustering and
# interpretation use the canonical long -> medium -> short path.

scbm_vhar_path_name <- function(path) {
  if (is.numeric(path)) path <- paste0("path", as.integer(path)[1L])
  path <- tolower(as.character(path)[1L])
  if (!path %in% paste0("path", 1:3)) stop("path must be path1, path2, or path3.")
  path
}

scbm_vhar_link_type <- function(type) {
  type <- tolower(as.character(type)[1L])
  if (type %in% c("1", "t1")) type <- "type1"
  if (type %in% c("2", "t2")) type <- "type2"
  if (!type %in% c("type1", "type2")) stop("type must be type1 or type2.")
  type
}

scbm_vhar_partitions <- function(q) {
  c2 <- scbm_balanced_labels(q, 2L)
  c2_swap <- c(2L, 1L)[c2]
  c3 <- scbm_balanced_labels(q, 3L)
  c3_231 <- c(2L, 3L, 1L)[c3]
  if (q %% 6L != 0L) stop("VHAR designs require q divisible by 6.")
  a <- q %/% 6L
  weekly <- integer(q)
  weekly[seq_len(2L * a)] <- 1L
  weekly[(2L * a + 1L):(4L * a)] <- 2L
  weekly[(4L * a + 1L):(6L * a)] <- 3L
  list(c2 = c2, c2_swap = c2_swap, c3 = c3,
       c3_231 = c3_231, weekly = weekly)
}

scbm_vhar_states <- function(q, path) {
  path <- scbm_vhar_path_name(path)
  p <- scbm_vhar_partitions(q)
  states <- switch(
    path,
    path1 = list(p$c3, p$c3, p$c3, p$c3),
    path2 = list(p$c2, p$c2, p$weekly, p$c2),
    path3 = list(p$c3_231, p$c3_231, p$c2_swap, p$c2)
  )
  names(states) <- c("long.L", "long.R", "medium.R", "short.R")
  pairs <- list(
    long = list(y = states[[1L]], z = states[[2L]]),
    medium = list(y = states[[2L]], z = states[[3L]]),
    short = list(y = states[[3L]], z = states[[4L]])
  )
  n_comm <- as.integer(unlist(lapply(pairs, function(z) c(max(z$y), max(z$z)))))
  list(states = lapply(states, as.integer), pairs = pairs, n_comm = n_comm)
}

scbm_vhar_setting <- function() {
  list(
    coeff = c(short = 0.34, medium = 0.28, long = 0.24),
    sv_target = 0.90,
    pdiag = c(short = 0.95, medium = 0.93, long = 0.91),
    poff = c(short = 0.02, medium = 0.03, long = 0.04)
  )
}

scbm_vhar_B <- function(y, z, type, horizon) {
  type <- scbm_vhar_link_type(type)
  pars <- scbm_vhar_setting()
  ky <- max(y); kz <- max(z)
  diag_p <- pars$pdiag[horizon]; off_p <- pars$poff[horizon]
  if (type == "type1") {
    b <- matrix(off_p, ky, kz)
    d <- min(ky, kz); b[cbind(seq_len(d), seq_len(d))] <- diag_p
  } else {
    b <- matrix(min(0.995, off_p + 0.06), ky, kz)
    for (i in seq_len(ky)) for (j in seq_len(kz)) {
      if (i < j) b[i, j] <- min(0.995, off_p + 0.03)
    }
    d <- min(ky, kz); b[cbind(seq_len(d), seq_len(d))] <- max(0.001, diag_p - 0.03)
  }
  b
}

scbm_vhar_sender_matrix <- function(y, z, type, horizon) {
  b <- scbm_vhar_B(y, z, type, horizon)
  q <- length(y)
  prob <- outer(y, z, Vectorize(function(a, c) b[a, c]))
  a <- matrix(rbinom(q * q, 1L, as.vector(prob)), q, q)
  group_size <- tabulate(y, nbins = max(y))[y]
  list(matrix = a / group_size, support = a, probability = prob, B = b)
}

scbm_vhar_to_var <- function(phi, bw = 3L, bm = 10L) {
  q <- nrow(phi[[1L]])
  out <- matrix(0, q, q * bm)
  out[, 1:q] <- phi[[1L]] + phi[[2L]] / bw + phi[[3L]] / bm
  if (bw >= 2L) for (h in 2:bw) {
    idx <- ((h - 1L) * q + 1L):(h * q)
    out[, idx] <- phi[[2L]] / bw + phi[[3L]] / bm
  }
  if (bm > bw) for (h in (bw + 1L):bm) {
    idx <- ((h - 1L) * q + 1L):(h * q)
    out[, idx] <- phi[[3L]] / bm
  }
  out
}

scbm_companion <- function(a, q, p) {
  f <- matrix(0, q * p, q * p)
  f[seq_len(q), ] <- a
  if (p > 1L) f[(q + 1L):(q * p), seq_len(q * (p - 1L))] <- diag(q * (p - 1L))
  f
}

scbm_vhar_radius <- function(phi, bw = 3L, bm = 10L) {
  f <- scbm_companion(scbm_vhar_to_var(phi, bw, bm), nrow(phi[[1L]]), bm)
  max(Mod(eigen(f, only.values = TRUE)$values))
}

scbm_vhar_unit_variance <- function(phi, bw = 3L, bm = 10L,
                                     tol = 1e-10, max_iter = 5000L) {
  q <- nrow(phi[[1L]])
  f <- scbm_companion(scbm_vhar_to_var(phi, bw, bm), q, bm)
  g <- matrix(0, q * bm, q); g[seq_len(q), ] <- diag(q)
  m <- g; sig <- diag(q)
  for (ii in seq_len(max_iter)) {
    m <- f %*% m
    b <- m[seq_len(q), , drop = FALSE]
    add <- b %*% t(b); sig <- sig + add
    if (scbm_operator_norm(add) < tol) break
  }
  mean(diag(sig))
}

scbm_vhar_generate_design <- function(q, path, type, bw = 3L, bm = 10L,
                                       target_var = 0.5, seed = 1L) {
  set.seed(scbm_normalize_seed(seed))
  path <- scbm_vhar_path_name(path); type <- scbm_vhar_link_type(type)
  structure <- scbm_vhar_states(q, path)
  pars <- scbm_vhar_setting()
  # DGP storage order: short, medium, long.
  pair <- rev(structure$pairs)
  horizons <- c("short", "medium", "long")
  raw <- Map(function(z, h) scbm_vhar_sender_matrix(z$y, z$z, type, h),
             pair, horizons)
  sender <- lapply(raw, `[[`, "matrix")
  mult <- pars$coeff * c(short = 1, medium = sqrt(bw), long = sqrt(bm))
  sender <- Map(function(m, w) w * m, sender, as.list(mult[horizons]))
  sv0 <- scbm_operator_norm(Reduce("+", sender))
  scale <- pars$sv_target / sv0
  sender <- lapply(sender, function(x) scale * x)
  phi <- lapply(sender, t)
  unit_var <- scbm_vhar_unit_variance(phi, bw, bm)
  noise_sd <- sqrt(target_var / unit_var)
  list(
    Phi = phi, sender_receiver = sender,
    support = lapply(raw, `[[`, "support"),
    states = structure$states,
    n_comm = structure$n_comm,       # canonical long, medium, short order
    path = path, type = type, q = as.integer(q), bw = bw, bm = bm,
    scale = scale, rho = scbm_vhar_radius(phi, bw, bm),
    noise_sd = noise_sd, target_var = target_var
  )
}

scbm_vhar_simulate <- function(phi, T, bw = 3L, bm = 10L,
                                burn = 300L, noise_sd = 0.5, seed = 1L) {
  set.seed(scbm_normalize_seed(seed))
  q <- nrow(phi[[1L]]); total <- as.integer(T + burn + bm)
  y <- matrix(0, q, total)
  e <- noise_sd * matrix(rnorm(q * total), q, total)
  for (tt in (bm + 1L):total) {
    short <- y[, tt - 1L]
    medium <- rowMeans(y[, (tt - bw):(tt - 1L), drop = FALSE])
    long <- rowMeans(y[, (tt - bm):(tt - 1L), drop = FALSE])
    y[, tt] <- phi[[1L]] %*% short + phi[[2L]] %*% medium +
      phi[[3L]] %*% long + e[, tt]
  }
  y[, (burn + bm + 1L):total, drop = FALSE]
}

scbm_vhar_split_blocks <- function(phi, q) {
  phi <- scbm_safe_matrix(phi)
  list(phi[, 1:q, drop = FALSE],
       phi[, (q + 1L):(2L * q), drop = FALSE],
       phi[, (2L * q + 1L):(3L * q), drop = FALSE])
}

scbm_vhar_ols <- function(y, bw = 5L, bm = 22L) {
  scbm_require_packages("sparseVAR")
  fit <- sparseVAR::VHAR_ols(Yt = y, bd = as.integer(bw), bm = as.integer(bm))
  phi <- scbm_safe_matrix(fit$Phi_hat)
  list(Phi_hat = phi, blocks = scbm_vhar_split_blocks(phi, nrow(y)), estimator = "ols")
}

scbm_vhar_lambda <- function(q, n, c_value) {
  text <- as.numeric(c_value) * sqrt(log(3 * q^2) / n)
  list(text = text, package = 0.5 * n * text)
}

scbm_vhar_extract_phi <- function(fit) {
  x <- fit$Phi_hat_lasso %||% fit$Phi_lasso_hat %||% fit$Phi_hat
  if (is.null(x)) stop("Could not extract the VHAR lasso coefficient matrix.")
  scbm_safe_matrix(x)
}

scbm_vhar_fit_lasso <- function(y, c_value, bw = 5L, bm = 22L,
                                 update_sigma = FALSE,
                                 max_iter = 1000L, tol = 1e-6) {
  scbm_require_packages("sparseVAR")
  n <- ncol(y) - bm; q <- nrow(y)
  lam <- scbm_vhar_lambda(q, n, c_value)
  fit <- sparseVAR::VHAR_adalasso_fista(
    Yt = y, type = "lasso", lambda = lam$package,
    diagTF = TRUE, updateSigma = isTRUE(update_sigma), sigma_diag_only = TRUE,
    max_iter = as.integer(max_iter), tol = tol,
    bd = as.integer(bw), bm = as.integer(bm)
  )
  phi <- scbm_vhar_extract_phi(fit)
  list(
    Phi_hat = phi, blocks = scbm_vhar_split_blocks(phi, q),
    estimator = "lasso", c_lambda = as.numeric(c_value),
    lambda_pkg = lam$package, lambda_text = lam$text
  )
}

scbm_vhar_cv_c <- function(y, c_grid = seq(0.1, 1, by = 0.05),
                            bw = 5L, bm = 22L, update_sigma = FALSE,
                            folds = 10L, max_iter = 1000L, tol = 1e-6) {
  scbm_require_packages("sparseVAR")
  q <- nrow(y); n <- ncol(y) - bm
  base <- scbm_vhar_lambda(q, n, 1)$package
  lambda_seq <- sort(unique(base * c_grid), decreasing = TRUE)
  fit <- sparseVAR::VHAR_adalasso_fista(
    Yt = y, type = "lasso", fold = as.integer(folds),
    lambda_seq = lambda_seq, nlambda = length(lambda_seq),
    diagTF = TRUE, updateSigma = isTRUE(update_sigma), sigma_diag_only = TRUE,
    max_iter = as.integer(max_iter), tol = tol,
    bd = as.integer(bw), bm = as.integer(bm)
  )
  selected <- as.numeric(fit$lambda_lasso %||% fit$lambda %||% NA_real_)[1L]
  if (!is.finite(selected)) selected <- stats::median(lambda_seq)
  csel <- c_grid[which.min(abs(c_grid - selected / base))]
  list(c = as.numeric(csel),
       curve = data.frame(c = c_grid, lambda = base * c_grid),
       selected_lambda = base * csel)
}

scbm_vhar_spec_error <- function(truth, estimate) {
  mean(unlist(Map(function(a, b) scbm_operator_norm(a - b), truth, estimate)))
}

scbm_vhar_estimated_states <- function(cluster) {
  # cluster stages are ordered long, medium, short.
  list(cluster$group_L[[1L]], cluster$group_R[[1L]],
       cluster$group_R[[2L]], cluster$group_R[[3L]])
}

scbm_vhar_fit_and_score <- function(y, design, c_value, alpha,
                                     seed, T_value) {
  fit <- scbm_vhar_fit_lasso(y, c_value, design$bw, design$bm)
  canonical_blocks <- rev(fit$blocks)
  cluster <- scbm_cluster_chain(
    canonical_blocks, design$n_comm, alpha = alpha, cyclic = FALSE,
    nstart = 50L, seed = scbm_make_seed(seed, T_value, round(alpha * 1e8))
  )
  score <- scbm_score_states(design$states, scbm_vhar_estimated_states(cluster))
  data.frame(
    T = as.integer(T_value),
    Spec = scbm_vhar_spec_error(design$Phi, fit$blocks),
    Accuracy = score$accuracy_overall, ARI = score$ari_overall,
    alpha = as.numeric(alpha), c_lambda = as.numeric(c_value),
    stringsAsFactors = FALSE
  )
}

scbm_vhar_tuning <- function(design, T_grid, seed) {
  maxT <- max(T_grid)
  yc <- scbm_vhar_simulate(design$Phi, maxT, design$bw, design$bm,
                           noise_sd = design$noise_sd, seed = scbm_make_seed(seed, 7001L))
  ya <- scbm_vhar_simulate(design$Phi, maxT, design$bw, design$bm,
                           noise_sd = design$noise_sd, seed = scbm_make_seed(seed, 8001L))
  c_map <- alpha_map <- setNames(numeric(length(T_grid)), as.character(T_grid))
  c_curves <- alpha_curves <- vector("list", length(T_grid)); names(c_curves) <- names(alpha_curves) <- names(c_map)
  for (TT in T_grid) {
    cfit <- scbm_vhar_cv_c(yc[, seq_len(TT), drop = FALSE],
                           bw = design$bw, bm = design$bm)
    c_map[as.character(TT)] <- cfit$c; c_curves[[as.character(TT)]] <- cfit$curve
    fit <- scbm_vhar_fit_lasso(ya[, seq_len(TT), drop = FALSE], cfit$c,
                               design$bw, design$bm)
    afit <- scbm_select_alpha(
      rev(fit$blocks), design$n_comm, rho = 0.20, repeats = 10L,
      seed = scbm_make_seed(seed, 8002L, TT)
    )
    alpha_map[as.character(TT)] <- afit$alpha
    alpha_curves[[as.character(TT)]] <- afit$curve
  }
  list(c = c_map, alpha = alpha_map, c_curve = c_curves, alpha_curve = alpha_curves)
}

scbm_vhar_rep_worker <- function(rep_id, design, T_grid, tuning, base_seed) {
  y <- scbm_vhar_simulate(
    design$Phi, max(T_grid), design$bw, design$bm,
    noise_sd = design$noise_sd,
    seed = scbm_make_seed(base_seed, design$q, rep_id)
  )
  out <- lapply(T_grid, function(TT) {
    tryCatch({
      z <- scbm_vhar_fit_and_score(
        y[, seq_len(TT), drop = FALSE], design,
        tuning$c[as.character(TT)], tuning$alpha[as.character(TT)],
        seed = scbm_make_seed(base_seed, rep_id), T_value = TT
      )
      z$replication <- as.integer(rep_id); z$error <- ""; z
    }, error = function(e) data.frame(
      T = TT, Spec = NA_real_, Accuracy = NA_real_, ARI = NA_real_,
      alpha = tuning$alpha[as.character(TT)], c_lambda = tuning$c[as.character(TT)],
      replication = as.integer(rep_id), error = conditionMessage(e),
      stringsAsFactors = FALSE
    ))
  })
  do.call(rbind, out)
}

scbm_vhar_run_cell <- function(q, path, type, T_grid, nrep, ncore,
                                seed, root = scbm_repo_root()) {
  design_seed <- scbm_make_seed(seed, q, match(path, paste0("path", 1:3)),
                                match(type, c("type1", "type2")))
  design <- scbm_vhar_generate_design(
    q, path, type, bw = 3L, bm = 10L, seed = design_seed
  )
  tuning <- scbm_vhar_tuning(design, T_grid, scbm_make_seed(design_seed, 99L))
  fun <- function(ii) scbm_vhar_rep_worker(ii, design, T_grid, tuning, design_seed)
  rows <- scbm_parallel_lapply(
    seq_len(as.integer(nrep)), fun, ncore = ncore, root = root,
    modules = c("scbm_core.R", "scbm_pvar.R", "scbm_vhar.R"),
    packages = c("MASS", "sparseVAR")
  )
  raw <- do.call(rbind, rows)
  raw$q <- q; raw$path <- path; raw$type <- type
  list(raw = raw, tuning = tuning, design = design)
}

scbm_run_vhar_simulation <- function(mode = c("quick", "paper"),
                                      nrep = NULL, ncore = 2L,
                                      seed = 12345L, root = scbm_repo_root()) {
  mode <- match.arg(mode)
  scbm_require_packages(c("MASS", "sparseVAR"))
  if (is.null(nrep)) nrep <- if (mode == "paper") 200L else 20L
  if (mode == "paper") {
    grid <- expand.grid(
      q = c(18L, 24L, 36L), path = paste0("path", 1:3),
      type = c("type1", "type2"), stringsAsFactors = FALSE
    )
    T_grid <- c(500L, 1000L, 2000L, 3000L)
  } else {
    grid <- data.frame(
      q = c(18L, 24L), path = c("path1", "path2"),
      type = c("type1", "type2"), stringsAsFactors = FALSE
    )
    T_grid <- c(500L, 3000L)
  }
  output <- scbm_ensure_dir(file.path(root, "output", "simulation", "VHAR"))
  cells <- vector("list", nrow(grid))
  for (ii in seq_len(nrow(grid))) {
    message(sprintf("VHAR %d/%d: q=%d, %s, %s", ii, nrow(grid),
                    grid$q[ii], grid$path[ii], grid$type[ii]))
    cells[[ii]] <- scbm_vhar_run_cell(
      grid$q[ii], grid$path[ii], grid$type[ii], T_grid,
      nrep, ncore, scbm_make_seed(seed, ii), root
    )
  }
  raw <- do.call(rbind, lapply(cells, `[[`, "raw"))
  ok <- raw[raw$error == "", , drop = FALSE]
  summary <- aggregate(
    cbind(Spec, Accuracy, ARI, alpha, c_lambda) ~ q + path + type + T,
    data = ok, FUN = mean
  )
  summary <- summary[order(summary$q, summary$path, summary$type, summary$T), ]
  scbm_write_csv(summary, file.path(output, paste0("vhar_", mode, "_summary.csv")))
  saveRDS(raw, file.path(output, paste0("vhar_", mode, "_raw.rds")))
  saveRDS(lapply(cells, `[[`, "tuning"), file.path(output, paste0("vhar_", mode, "_tuning.rds")))
  invisible(list(summary = summary, raw = raw, output_dir = output))
}

scbm_prepare_vhar_data <- function(file) {
  scbm_require_packages("imputeTS")
  obj <- scbm_load_rdata_object(file)
  df <- if (is.matrix(obj) || is.data.frame(obj)) {
    as.data.frame(obj, check.names = FALSE)
  } else as.data.frame(scbm_first_numeric_matrix(obj), check.names = FALSE)
  dates <- as.Date(rownames(df))
  if (anyNA(dates)) stop("Could not parse dates from row names.")
  keep <- dates >= as.Date("2010-01-01") & dates <= as.Date("2019-12-31")
  df <- df[keep, , drop = FALSE]
  x <- matrix(NA_real_, nrow(df), ncol(df), dimnames = dimnames(df))
  for (j in seq_len(ncol(df))) {
    z <- as.numeric(df[[j]]); z[z <= 0] <- NA_real_
    x[, j] <- imputeTS::na_interpolation(log(z))
  }
  x <- x[, -c(5L, 30L), drop = FALSE]
  colnames(x) <- sub("^\\.", "", colnames(x))
  t(scale(x, center = TRUE, scale = FALSE))
}

scbm_run_vhar_empirical <- function(root = scbm_repo_root(), seed = 12345L) {
  scbm_require_packages(c("sparseVAR", "imputeTS"))
  y <- scbm_prepare_vhar_data(file.path(root, "data", "rk_mat2000-2022.Rdata"))
  output <- scbm_ensure_dir(file.path(root, "output", "empirical", "VHAR"))
  fit_ols <- scbm_vhar_ols(y, bw = 5L, bm = 22L)
  cfit <- scbm_vhar_cv_c(y, bw = 5L, bm = 22L, update_sigma = TRUE)
  fit <- scbm_vhar_fit_lasso(y, cfit$c, bw = 5L, bm = 22L, update_sigma = TRUE)
  n_comm_canonical <- rep(3L, 6L)
  afit <- scbm_select_alpha(
    fit$blocks, n_comm_canonical, rho = 0.20, repeats = 10L,
    seed = scbm_make_seed(seed, 202L)
  )
  cluster_seed <- scbm_make_seed(
    seed, 920001L, nrow(y), 5L, 22L, scbm_alpha_key(afit$alpha)
  )
  cluster <- scbm_cluster_vhar_empirical(
    fit$blocks, afit$alpha, state_counts = c(3L, 3L, 3L),
    nstart = 50L, seed = cluster_seed
  )
  stage <- c("Long (Monthly)", "Medium (Weekly)", "Short (Daily)")
  flow <- scbm_flow_table(cluster$states, rownames(y), stage)
  scree <- scbm_scree(fit_ols$blocks, c("daily", "weekly", "monthly"))
  scbm_plot_scree(scree, file.path(output, "vhar_scree.pdf"), 0.75,
                   width = 12, height = 4.5)
  scbm_save_vhar_sankey(flow, file.path(output, "vhar_sankey.pdf"))
  scbm_write_csv(flow, file.path(output, "vhar_membership_path.csv"), row.names = TRUE)
  scbm_write_csv(cluster$aggregate_definition,
                 file.path(output, "vhar_aggregate_definition.csv"))
  scbm_write_csv(afit$curve, file.path(output, "vhar_alpha_cv.csv"))
  scbm_write_csv(cfit$curve, file.path(output, "vhar_lambda_cv.csv"))
  summary <- data.frame(
    model = "VHAR", estimator = "lasso", q = nrow(y), T = ncol(y),
    specification = "Long--Medium--Short: 3--3--3",
    common_ranks = "(3,3,3)", alpha = afit$alpha, c_lambda = cfit$c,
    alpha_at_grid_upper = isTRUE(afit$selected_at_upper),
    cluster_seed = as.integer(cluster_seed)
  )
  scbm_write_csv(summary, file.path(output, "vhar_summary.csv"))
  saveRDS(list(fit = fit, cluster = cluster, alpha_cv = afit, lambda_cv = cfit),
          file.path(output, "vhar_result.rds"))
  invisible(list(summary = summary, fit = fit, cluster = cluster, output_dir = output))
}
