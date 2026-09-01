# ScBM--PVAR simulation and payroll application.
# This file is standalone apart from R/scbm_core.R and the sparseVAR package.

scbm_pvar_path_name <- function(path) {
  if (is.numeric(path)) path <- paste0("path", as.integer(path)[1L])
  path <- tolower(as.character(path)[1L])
  if (!path %in% paste0("path", 1:3)) stop("path must be path1, path2, or path3.")
  path
}

scbm_pvar_path_spec <- function(path) {
  switch(
    scbm_pvar_path_name(path),
    path1 = c(4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L),
    path2 = c(2L, 3L, 3L, 3L, 3L, 2L, 2L, 2L),
    path3 = c(2L, 2L, 2L, 3L, 3L, 4L, 4L, 2L)
  )
}

scbm_pvar_paths <- function(q, path) {
  path <- scbm_pvar_path_name(path)
  if (path == "path1") {
    y1 <- scbm_balanced_labels(q, 4L)
    y2 <- y3 <- y4 <- y1
  } else if (path == "path2") {
    y1 <- scbm_balanced_labels(q, 2L)
    i1 <- which(y1 == 1L); i2 <- which(y1 == 2L)
    n11 <- floor(2 * length(i1) / 3); n21 <- floor(length(i2) / 3)
    y2 <- integer(q)
    y2[i1[seq_len(n11)]] <- 1L
    y2[i1[(n11 + 1L):length(i1)]] <- 2L
    y2[i2[seq_len(n21)]] <- 2L
    y2[i2[(n21 + 1L):length(i2)]] <- 3L
    y3 <- y2
    c1 <- which(y3 == 1L); c2 <- which(y3 == 2L); c3 <- which(y3 == 3L)
    half <- floor(length(c2) / 2)
    y4 <- integer(q)
    y4[c1] <- 1L
    y4[c2[seq_len(half)]] <- 1L
    y4[c2[(half + 1L):length(c2)]] <- 2L
    y4[c3] <- 2L
  } else {
    y1 <- scbm_balanced_labels(q, 2L)
    y2 <- y1
    y3 <- scbm_split_label(y2, 2L, c(2L, 3L))
    y4 <- y3
    i1 <- which(y3 == 1L)
    repl <- scbm_balanced_labels(length(i1), 2L)
    y4[i1] <- c(1L, 2L)[repl]
    y4[y3 == 2L] <- 3L
    y4[y3 == 3L] <- 4L
  }
  y <- list(as.integer(y1), as.integer(y2), as.integer(y3), as.integer(y4))
  z <- list(y[[2L]], y[[3L]], y[[4L]], y[[1L]])
  list(y = y, z = z)
}

scbm_pvar_link_type <- function(type) {
  type <- tolower(as.character(type)[1L])
  if (type %in% c("1", "t1")) type <- "type1"
  if (type %in% c("2", "t2")) type <- "type2"
  if (!type %in% c("type1", "type2")) stop("type must be type1 or type2.")
  type
}

scbm_pvar_pair_class <- function(a, b, ky, kz) {
  d <- min(ky, kz)
  if (a <= d && b <= d && a == b) return("diag")
  if (a < b) "upper" else "lower"
}

scbm_pvar_support_prob <- function(y, z, type, cap = 0.95) {
  type <- scbm_pvar_link_type(type)
  q <- length(y); ky <- max(y); kz <- max(z)
  target <- if (type == "type1") {
    c(diag = 2.6, upper = 0.6, lower = 0.9)
  } else c(diag = 2.2, upper = 0.7, lower = 1.3)
  zsize <- tabulate(z, nbins = kz)
  p <- matrix(0, q, q)
  for (i in seq_len(q)) for (j in seq_len(q)) {
    if (i == j) next
    cls <- scbm_pvar_pair_class(y[i], z[j], ky, kz)
    den <- max(1, zsize[z[j]] - as.integer(y[i] == z[j]))
    p[i, j] <- min(cap, target[cls] / den)
  }
  p
}

scbm_pvar_sender_matrix <- function(y, z, type, cap = 0.95) {
  type <- scbm_pvar_link_type(type)
  p <- scbm_pvar_support_prob(y, z, type, cap)
  q <- length(y); ky <- max(y); kz <- max(z)
  a <- matrix(rbinom(q * q, 1L, as.vector(p)), q, q)
  diag(a) <- 0
  mag <- if (type == "type1") {
    c(self = 0.30, diag = 0.14, upper = 0.04, lower = 0.06)
  } else c(self = 0.28, diag = 0.12, upper = 0.05, lower = 0.07)
  m <- matrix(0, q, q)
  for (i in seq_len(q)) {
    m[i, i] <- mag["self"] * runif(1L, 0.95, 1.05)
    for (j in seq_len(q)) {
      if (i == j || a[i, j] == 0L) next
      cls <- scbm_pvar_pair_class(y[i], z[j], ky, kz)
      m[i, j] <- mag[cls] * runif(1L, 0.90, 1.10)
    }
  }
  list(matrix = m, support = a, probability = p)
}

scbm_pvar_generate_design <- function(q, path, type,
                                       target_cycle_sv = 0.90,
                                       seed = 1L) {
  set.seed(scbm_normalize_seed(seed))
  path <- scbm_pvar_path_name(path)
  type <- scbm_pvar_link_type(type)
  membership <- scbm_pvar_paths(q, path)
  raw <- Map(function(y, z) scbm_pvar_sender_matrix(y, z, type),
             membership$y, membership$z)
  sender <- lapply(raw, `[[`, "matrix")
  transition_raw <- lapply(sender, t)
  product <- Reduce("%*%", rev(transition_raw))
  sv <- scbm_operator_norm(product)
  if (!is.finite(sv) || sv <= 0) stop("Invalid PVAR cycle singular value.")
  scale <- (target_cycle_sv / sv)^(1 / 4)
  phi <- lapply(transition_raw, function(x) scale * x)
  list(
    Phi = phi,
    sender_receiver = sender,
    support = lapply(raw, `[[`, "support"),
    probability = lapply(raw, `[[`, "probability"),
    y = membership$y, z = membership$z,
    n_comm = scbm_pvar_path_spec(path),
    path = path, type = type, q = as.integer(q),
    scale = scale, target_cycle_sv = target_cycle_sv
  )
}

scbm_pvar_simulate <- function(phi, T, sigma = 0.5, burn = 500L, seed = 1L) {
  set.seed(scbm_normalize_seed(seed))
  q <- nrow(phi[[1L]]); s <- length(phi); total <- as.integer(T + burn)
  e <- MASS::mvrnorm(total, mu = rep(0, q), Sigma = diag(sigma, q))
  y <- matrix(0, total, q); y[1L, ] <- e[1L, ]
  for (tt in 2:total) {
    season <- ((tt - 1L) %% s) + 1L
    y[tt, ] <- drop(phi[[season]] %*% y[tt - 1L, ]) + e[tt, ]
  }
  t(y[-seq_len(burn), , drop = FALSE])
}

scbm_pvar_split_blocks <- function(phi, s = 4L) {
  phi <- scbm_safe_matrix(phi); q <- nrow(phi)
  lapply(seq_len(s), function(mm) {
    phi[, ((mm - 1L) * q + 1L):(mm * q), drop = FALSE]
  })
}

scbm_pvar_ols <- function(y, s = 4L) {
  y <- scbm_safe_matrix(y); q <- nrow(y); T <- ncol(y)
  out <- matrix(0, q, q * s)
  for (mm in seq_len(s)) {
    idx <- which(((seq_len(T) - 1L) %% s) + 1L == mm)
    idx <- idx[idx > 1L]
    x <- t(y[, idx - 1L, drop = FALSE]); z <- t(y[, idx, drop = FALSE])
    b <- solve(crossprod(x) + diag(1e-8, q), crossprod(x, z))
    out[, ((mm - 1L) * q + 1L):(mm * q)] <- t(b)
  }
  list(Phi_hat = out, blocks = scbm_pvar_split_blocks(out, s), estimator = "ols")
}

scbm_pvar_center <- function(y, s = 4L) {
  q <- nrow(y); T <- ncol(y)
  season <- ((seq_len(T) - 1L) %% s) + 1L
  mu <- vapply(seq_len(s), function(mm) rowMeans(y[, season == mm, drop = FALSE]),
               numeric(q))
  y - mu[, season, drop = FALSE]
}

scbm_pvar_design_matrix <- function(y, s, season) {
  T <- ncol(y)
  idx <- which(((seq_len(T) - 1L) %% s) + 1L == season)
  idx <- idx[idx > 1L]
  list(X = t(y[, idx - 1L, drop = FALSE]), Y = t(y[, idx, drop = FALSE]))
}

scbm_pvar_lambda <- function(q, s, n, c_value) {
  text <- as.numeric(c_value) * sqrt(log(s * q^2) / n)
  list(text = text, package = 0.5 * n * text)
}

scbm_fista_solver <- function() {
  getFromNamespace("fista_lasso_multi_cpp", "sparseVAR")
}

scbm_precision_from_design <- function(x, y, diagonal = TRUE) {
  b0 <- solve(crossprod(x) + diag(1e-8, ncol(x)), crossprod(x, y))
  r <- y - x %*% b0
  sigma <- crossprod(r) / nrow(r)
  if (diagonal) return(diag(1 / pmax(diag(sigma), 1e-12), ncol(y)))
  sv <- svd(sigma)
  sv$u %*% diag(1 / pmax(sv$d, 1e-12), nrow(sigma)) %*% t(sv$u)
}

scbm_pvar_fit_lasso <- function(y, c_value, s = 4L,
                                 update_sigma = FALSE, sigma_diagonal = TRUE,
                                 max_iter = 1000L, tol = 1e-6) {
  scbm_require_packages("sparseVAR")
  y <- scbm_pvar_center(scbm_safe_matrix(y), s)
  q <- nrow(y); out <- matrix(0, q, q * s)
  lambda_pkg <- lambda_text <- numeric(s)
  weight <- matrix(1, q, q); diag(weight) <- 0
  solver <- scbm_fista_solver()
  for (mm in seq_len(s)) {
    des <- scbm_pvar_design_matrix(y, s, mm)
    lam <- scbm_pvar_lambda(q, s, nrow(des$X), c_value)
    sigma_inv <- if (isTRUE(update_sigma)) {
      scbm_precision_from_design(des$X, des$Y, sigma_diagonal)
    } else matrix(0, 0, 0)
    bhat <- solver(
      X = des$X, Y = des$Y, weights = weight, lambda = lam$package,
      SigmaInv = sigma_inv, max_iter = as.integer(max_iter), tol = tol
    )
    out[, ((mm - 1L) * q + 1L):(mm * q)] <- t(bhat)
    lambda_pkg[mm] <- lam$package; lambda_text[mm] <- lam$text
  }
  list(
    Phi_hat = out, blocks = scbm_pvar_split_blocks(out, s),
    estimator = "lasso", c_lambda = as.numeric(c_value),
    lambda_pkg = lambda_pkg, lambda_text = lambda_text
  )
}

scbm_pvar_cv_c <- function(y, c_grid = seq(0.1, 1, by = 0.05),
                            s = 4L, folds = 10L,
                            update_sigma = FALSE, sigma_diagonal = TRUE,
                            max_iter = 1000L, tol = 1e-6) {
  scbm_require_packages("sparseVAR")
  y <- scbm_pvar_center(scbm_safe_matrix(y), s)
  q <- nrow(y); weight <- matrix(1, q, q); diag(weight) <- 0
  solver <- scbm_fista_solver()
  scores <- vapply(c_grid, function(cc) {
    total <- 0
    for (mm in seq_len(s)) {
      des <- scbm_pvar_design_matrix(y, s, mm)
      n <- nrow(des$X); k <- max(2L, min(as.integer(folds), n))
      fold_id <- rep(seq_len(k), each = floor(n / k) + 1L)[seq_len(n)]
      lam <- scbm_pvar_lambda(q, s, n, cc)$package
      sigma_inv <- if (isTRUE(update_sigma)) {
        scbm_precision_from_design(des$X, des$Y, sigma_diagonal)
      } else NULL
      for (ff in seq_len(k)) {
        tr <- fold_id != ff; va <- !tr
        xtr <- des$X[tr, , drop = FALSE]
        ytr <- des$Y[tr, , drop = FALSE]
        b <- solver(
          X = xtr, Y = ytr, weights = weight, lambda = lam,
          SigmaInv = if (is.null(sigma_inv)) matrix(0, 0, 0) else sigma_inv,
          max_iter = as.integer(max_iter), tol = tol
        )
        r <- des$Y[va, , drop = FALSE] - des$X[va, , drop = FALSE] %*% b
        total <- total + if (is.null(sigma_inv)) {
          sum(r^2)
        } else {
          sum((r %*% sigma_inv) * r)
        }
      }
    }
    total
  }, numeric(1L))
  list(c = as.numeric(c_grid[which.min(scores)]),
       curve = data.frame(c = c_grid, score = scores))
}

scbm_pvar_spec_error <- function(truth, estimate) {
  mean(unlist(Map(function(a, b) scbm_operator_norm(a - b), truth, estimate)))
}

scbm_pvar_fit_and_score <- function(y, design, c_value, alpha,
                                     seed, T_value) {
  fit <- scbm_pvar_fit_lasso(y, c_value)
  cluster <- scbm_cluster_chain(
    fit$blocks, design$n_comm, alpha = alpha, cyclic = TRUE,
    nstart = 50L, seed = scbm_make_seed(seed, T_value, round(alpha * 1e8))
  )
  score <- scbm_score_states(design$y, cluster$group_L)
  data.frame(
    T = as.integer(T_value),
    Spec = scbm_pvar_spec_error(design$Phi, fit$blocks),
    Accuracy = score$accuracy_overall,
    ARI = score$ari_overall,
    alpha = as.numeric(alpha), c_lambda = as.numeric(c_value),
    stringsAsFactors = FALSE
  )
}

scbm_pvar_tuning <- function(design, T_grid, seed) {
  maxT <- max(T_grid)
  y_c <- scbm_pvar_simulate(design$Phi, maxT, seed = scbm_make_seed(seed, 7001L))
  y_a <- scbm_pvar_simulate(design$Phi, maxT, seed = scbm_make_seed(seed, 8001L))
  c_map <- alpha_map <- setNames(numeric(length(T_grid)), as.character(T_grid))
  c_curves <- alpha_curves <- vector("list", length(T_grid)); names(c_curves) <- names(alpha_curves) <- names(c_map)
  for (TT in T_grid) {
    cfit <- scbm_pvar_cv_c(y_c[, seq_len(TT), drop = FALSE])
    c_map[as.character(TT)] <- cfit$c; c_curves[[as.character(TT)]] <- cfit$curve
    fit <- scbm_pvar_fit_lasso(y_a[, seq_len(TT), drop = FALSE], cfit$c)
    afit <- scbm_select_alpha(
      fit$blocks, design$n_comm, rho = 0.20, repeats = 10L,
      seed = scbm_make_seed(seed, 8002L, TT)
    )
    alpha_map[as.character(TT)] <- afit$alpha
    alpha_curves[[as.character(TT)]] <- afit$curve
  }
  list(c = c_map, alpha = alpha_map, c_curve = c_curves, alpha_curve = alpha_curves)
}

scbm_pvar_rep_worker <- function(rep_id, design, T_grid, tuning, base_seed) {
  y <- scbm_pvar_simulate(
    design$Phi, max(T_grid),
    seed = scbm_make_seed(base_seed, design$q, rep_id)
  )
  out <- lapply(T_grid, function(TT) {
    tryCatch({
      z <- scbm_pvar_fit_and_score(
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

scbm_pvar_run_cell <- function(q, path, type, T_grid, nrep, ncore,
                                seed, root = scbm_repo_root()) {
  design_seed <- scbm_make_seed(seed, q, match(path, paste0("path", 1:3)),
                                match(type, c("type1", "type2")))
  design <- scbm_pvar_generate_design(q, path, type, seed = design_seed)
  tuning <- scbm_pvar_tuning(design, T_grid, scbm_make_seed(design_seed, 99L))
  fun <- function(ii) scbm_pvar_rep_worker(ii, design, T_grid, tuning, design_seed)
  rows <- scbm_parallel_lapply(
    seq_len(as.integer(nrep)), fun, ncore = ncore, root = root,
    modules = c("scbm_core.R", "scbm_pvar.R"),
    packages = c("MASS", "sparseVAR")
  )
  raw <- do.call(rbind, rows)
  raw$q <- q; raw$path <- path; raw$type <- type
  list(raw = raw, tuning = tuning, design = design)
}

scbm_run_pvar_simulation <- function(mode = c("quick", "paper"),
                                      nrep = NULL, ncore = 2L,
                                      seed = 12345L, root = scbm_repo_root()) {
  mode <- match.arg(mode)
  scbm_require_packages(c("MASS", "sparseVAR"))
  if (is.null(nrep)) nrep <- if (mode == "paper") 200L else 20L
  if (mode == "paper") {
    grid <- expand.grid(
      q = c(18L, 36L, 60L), path = paste0("path", 1:3),
      type = c("type1", "type2"), stringsAsFactors = FALSE
    )
    T_grid <- c(200L, 500L, 1000L, 2000L)
  } else {
    grid <- data.frame(
      q = c(18L, 36L), path = c("path1", "path2"),
      type = c("type1", "type2"), stringsAsFactors = FALSE
    )
    T_grid <- c(200L, 2000L)
  }
  output <- scbm_ensure_dir(file.path(root, "output", "simulation", "PVAR"))
  cells <- vector("list", nrow(grid))
  for (ii in seq_len(nrow(grid))) {
    message(sprintf("PVAR %d/%d: q=%d, %s, %s", ii, nrow(grid),
                    grid$q[ii], grid$path[ii], grid$type[ii]))
    cells[[ii]] <- scbm_pvar_run_cell(
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
  scbm_write_csv(summary, file.path(output, paste0("pvar_", mode, "_summary.csv")))
  saveRDS(raw, file.path(output, paste0("pvar_", mode, "_raw.rds")))
  saveRDS(lapply(cells, `[[`, "tuning"), file.path(output, paste0("pvar_", mode, "_tuning.rds")))
  invisible(list(summary = summary, raw = raw, output_dir = output))
}

scbm_pvar_load_empirical <- function(file) {
  obj <- scbm_load_rdata_object(file)
  if (is.list(obj) && !is.null(obj$Y_diff)) {
    y <- as.matrix(obj$Y_diff)
    nm <- obj$series_codes %||% colnames(y)
    tm <- obj$quarterly_labels %||% rownames(y)
  } else {
    y <- scbm_first_numeric_matrix(obj); nm <- colnames(y); tm <- rownames(y)
  }
  if (nrow(y) > ncol(y)) y <- t(y)
  if (is.null(nm) || length(nm) != nrow(y)) nm <- paste0("V", seq_len(nrow(y)))
  rownames(y) <- nm
  list(Y = y, names = nm, time = tm)
}

scbm_run_pvar_empirical <- function(root = scbm_repo_root(), seed = 12345L) {
  scbm_require_packages(c("sparseVAR", "MASS"))
  data <- scbm_pvar_load_empirical(file.path(root, "data", "pvar_payroll_1990_2020.RData"))
  output <- scbm_ensure_dir(file.path(root, "output", "empirical", "PVAR"))
  fit_ols <- scbm_pvar_ols(data$Y)
  cfit <- scbm_pvar_cv_c(data$Y, update_sigma = TRUE)
  fit <- scbm_pvar_fit_lasso(data$Y, cfit$c, update_sigma = TRUE)
  n_comm <- c(2L, 2L, 2L, 3L, 3L, 3L, 3L, 2L)
  afit <- scbm_select_alpha(
    fit$blocks, n_comm, rho = 0.20, repeats = 10L,
    seed = scbm_make_seed(seed, 101L)
  )
  cluster_seed <- scbm_make_seed(
    seed, 910001L, nrow(data$Y), 4L, 1L, scbm_alpha_key(afit$alpha)
  )
  cluster <- scbm_cluster_pvar_empirical(
    fit$blocks, n_comm, alpha = afit$alpha,
    nstart = 50L, seed = cluster_seed
  )
  stage <- c("Q4/Q1 (Spring)", "Q1/Q2 (Summer)",
             "Q2/Q3 (Fall)", "Q3/Q4 (Winter)")
  flow <- scbm_flow_table(cluster$group_L, data$names, stage)
  scree <- scbm_scree(fit_ols$blocks, paste0("season", 1:4))
  scbm_plot_scree(scree, file.path(output, "pvar_scree.pdf"), 0.80)
  scbm_save_pvar_sankey(flow, file.path(output, "pvar_sankey.pdf"))
  scbm_write_csv(flow, file.path(output, "pvar_membership_path.csv"), row.names = TRUE)
  scbm_write_csv(afit$curve, file.path(output, "pvar_alpha_cv.csv"))
  scbm_write_csv(cfit$curve, file.path(output, "pvar_lambda_cv.csv"))
  summary <- data.frame(
    model = "PVAR", estimator = "lasso", q = nrow(data$Y), T = ncol(data$Y),
    specification = "(2,2; 2,3; 3,3; 3,2)",
    common_ranks = "(2,2,3,2)", alpha = afit$alpha, c_lambda = cfit$c,
    alpha_at_grid_upper = isTRUE(afit$selected_at_upper),
    cluster_seed = as.integer(cluster_seed)
  )
  scbm_write_csv(summary, file.path(output, "pvar_summary.csv"))
  saveRDS(list(fit = fit, cluster = cluster, alpha_cv = afit, lambda_cv = cfit),
          file.path(output, "pvar_result.rds"))
  invisible(list(summary = summary, fit = fit, cluster = cluster, output_dir = output))
}
