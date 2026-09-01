# Core utilities for ScBM--PVAR and generalized ScBM--VHAR reproduction.
#
# Canonical directed convention used throughout:
#   M = t(Phi), left singular space = sending communities,
#   right singular space = receiving communities, and
#   r_m = min(K_y,m, K_z,m) on both sides.

`%||%` <- function(x, y) if (is.null(x)) y else x

scbm_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  ff <- grep("^--file=", args, value = TRUE)
  if (length(ff)) {
    return(dirname(normalizePath(
      sub("^--file=", "", ff[1L]), winslash = "/", mustWork = FALSE
    )))
  }
  frames <- sys.frames()
  if (length(frames)) {
    for (ii in rev(seq_along(frames))) {
      of <- frames[[ii]]$ofile
      if (!is.null(of) && nzchar(of)) {
        return(dirname(normalizePath(of, winslash = "/", mustWork = FALSE)))
      }
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

scbm_repo_root <- function() {
  here <- scbm_script_dir()
  if (basename(here) == "R") return(dirname(here))
  if (file.exists(file.path(here, "R", "scbm_core.R"))) return(here)
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

scbm_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

scbm_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

scbm_safe_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L) x <- matrix(x, ncol = 1L)
  x
}

scbm_row_normalize <- function(x) {
  x <- scbm_safe_matrix(x)
  rn <- sqrt(rowSums(x^2))
  rn[!is.finite(rn) | rn <= .Machine$double.eps] <- 1
  x / rn
}

scbm_operator_norm <- function(x) {
  x <- scbm_safe_matrix(x)
  if (any(!is.finite(x))) return(NA_real_)
  d <- tryCatch(svd(x, nu = 0L, nv = 0L)$d, error = function(e) numeric())
  if (length(d)) as.numeric(d[1L]) else NA_real_
}

scbm_mean_finite <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(NA_real_)
  mean(x[is.finite(x)])
}

scbm_normalize_seed <- function(seed, default = 1L) {
  mod <- 2147483647
  x <- suppressWarnings(as.numeric(seed)[1L])
  if (!is.finite(x)) x <- default
  x <- floor(abs(x)) %% mod
  if (x <= 0) x <- default
  as.integer(x)
}

scbm_make_seed <- function(base_seed, ..., default = 1L) {
  mod <- 2147483647
  x <- as.numeric(scbm_normalize_seed(base_seed, default))
  vals <- unlist(list(...), use.names = FALSE)
  for (ii in seq_along(vals)) {
    v <- suppressWarnings(as.numeric(vals[ii]))
    if (!is.finite(v)) v <- 0
    x <- ((x * 1000003) + floor(abs(v)) + 97 * ii) %% mod
    if (x <= 0) x <- default
  }
  as.integer(x)
}

scbm_load_rdata_object <- function(file, object_name = NULL) {
  if (!file.exists(file)) stop("Missing data file: ", file, call. = FALSE)
  ee <- new.env(parent = emptyenv())
  load(file, envir = ee)
  nm <- ls(ee)
  if (!length(nm)) stop("No object found in ", file, call. = FALSE)
  object_name <- object_name %||% nm[1L]
  if (!object_name %in% nm) stop("Object not found: ", object_name, call. = FALSE)
  get(object_name, envir = ee)
}

scbm_first_numeric_matrix <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) return(as.matrix(x))
  if (is.list(x)) {
    preferred <- c("Y_diff", "Y", "Yt", "data", "logrv", "rv", "matrix")
    for (nm in preferred) {
      z <- x[[nm]]
      if (is.matrix(z) || is.data.frame(z)) return(as.matrix(z))
    }
    for (z in x) if (is.matrix(z) || is.data.frame(z)) return(as.matrix(z))
  }
  stop("No numeric matrix/data.frame found in the loaded object.", call. = FALSE)
}

scbm_write_csv <- function(x, file, row.names = FALSE) {
  if (is.null(x)) return(invisible(FALSE))
  utils::write.csv(as.data.frame(x), file, row.names = row.names)
  invisible(TRUE)
}


scbm_balanced_labels <- function(q, k) {
  q <- as.integer(q); k <- as.integer(k)
  size <- rep(q %/% k, k)
  if (q %% k) size[seq_len(q %% k)] <- size[seq_len(q %% k)] + 1L
  rep(seq_len(k), times = size)
}

scbm_split_label <- function(labels, old, new) {
  labels <- as.integer(labels)
  idx <- which(labels == old)
  repl <- scbm_balanced_labels(length(idx), length(new))
  labels[idx] <- new[repl]
  labels
}

scbm_pairs <- function(n_comm, n_stage = NULL) {
  n_comm <- as.integer(n_comm)
  if (!length(n_comm) || length(n_comm) %% 2L != 0L || any(n_comm < 1L)) {
    stop("n_comm must contain positive (Ky,Kz) pairs.", call. = FALSE)
  }
  out <- matrix(n_comm, ncol = 2L, byrow = TRUE,
                dimnames = list(NULL, c("Ky", "Kz")))
  if (!is.null(n_stage) && nrow(out) != as.integer(n_stage)) {
    stop("n_comm stage count mismatch.", call. = FALSE)
  }
  out
}

scbm_common_ranks <- function(n_comm) {
  pp <- scbm_pairs(n_comm)
  as.integer(pmin(pp[, "Ky"], pp[, "Kz"]))
}

scbm_alpha_max <- function() 1 / (4 * sqrt(2) + 2)

scbm_alpha_grid <- function(n_positive = 20L,
                            lower_fraction = 0.01,
                            upper_fraction = 0.99) {
  amax <- scbm_alpha_max()
  c(0, exp(seq(
    log(lower_fraction * amax), log(upper_fraction * amax),
    length.out = as.integer(n_positive)
  )))
}

# Integer key used by the production empirical RNG policy.
scbm_alpha_key <- function(alpha) {
  alpha <- as.numeric(alpha)[1L]
  if (!is.finite(alpha) || alpha < 0) {
    stop("alpha must be nonnegative and finite.", call. = FALSE)
  }
  as.integer(round(alpha / scbm_alpha_max() * 1e6))
}

scbm_permutations <- function(x) {
  x <- as.integer(x)
  if (length(x) <= 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(ii) {
    rest <- scbm_permutations(x[-ii])
    cbind(x[ii], rest)
  }))
}

scbm_adjusted_rand <- function(x, y) {
  x <- as.integer(x); y <- as.integer(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n <= 1L) return(1)
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  index <- sum(choose2(tab))
  ai <- sum(choose2(rowSums(tab)))
  bj <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- ai * bj / total
  maximum <- 0.5 * (ai + bj)
  den <- maximum - expected
  if (abs(den) <= .Machine$double.eps) {
    same_partition <- identical(outer(x, x, `==`), outer(y, y, `==`))
    return(as.numeric(same_partition))
  }
  as.numeric((index - expected) / den)
}

scbm_match_labels <- function(truth, estimate) {
  truth <- as.integer(truth); estimate <- as.integer(estimate)
  ok <- is.finite(truth) & is.finite(estimate)
  truth <- truth[ok]; estimate <- estimate[ok]
  if (!length(truth)) return(list(accuracy = NA_real_, ari = NA_real_, labels = integer()))
  k <- max(c(truth, estimate))
  best <- -Inf; best_labels <- estimate
  permutations <- scbm_permutations(seq_len(k))
  for (ii in seq_len(nrow(permutations))) {
    map <- permutations[ii, ]
    mapped <- map[estimate]
    acc <- mean(mapped == truth)
    if (acc > best) { best <- acc; best_labels <- mapped }
  }
  list(
    accuracy = as.numeric(best),
    ari = scbm_adjusted_rand(truth, best_labels),
    labels = best_labels
  )
}

scbm_kmeans <- function(x, k, nstart = 50L, seed = NULL) {
  x <- scbm_safe_matrix(x)
  k <- as.integer(k)[1L]
  if (k <= 1L) return(rep.int(1L, nrow(x)))
  if (nrow(x) < k) stop("Number of rows is smaller than k.", call. = FALSE)
  if (!is.null(seed)) set.seed(scbm_normalize_seed(seed))

  run_once <- function(xin, algorithm = "Hartigan-Wong", nstart_use = nstart) {
    stats::kmeans(
      xin, centers = k, nstart = as.integer(nstart_use),
      iter.max = 200L, algorithm = algorithm
    )$cluster
  }

  if (nrow(unique(round(x, 12))) < k) {
    set.seed(1L)
    x <- x + matrix(rnorm(length(x), sd = 1e-8), nrow(x), ncol(x))
  }

  out <- tryCatch(
    run_once(x, algorithm = "Hartigan-Wong", nstart_use = nstart),
    error = function(e) NULL
  )
  if (!is.null(out)) return(as.integer(out))

  set.seed(11L)
  xj <- x + matrix(rnorm(length(x), sd = 1e-7), nrow(x), ncol(x))
  out <- tryCatch(
    run_once(xj, algorithm = "Lloyd", nstart_use = nstart),
    error = function(e) NULL
  )
  if (!is.null(out)) return(as.integer(out))

  out <- tryCatch(
    run_once(xj, algorithm = "MacQueen", nstart_use = 1L),
    error = function(e) NULL
  )
  if (!is.null(out)) return(as.integer(out))

  as.integer(stats::cutree(stats::hclust(stats::dist(xj), method = "ward.D2"), k = k))
}

scbm_svd_stage <- function(phi, ky, kz) {
  phi <- scbm_safe_matrix(phi)
  m <- t(phi)
  r <- min(as.integer(ky), as.integer(kz))
  if (r < 1L || r > min(dim(m))) stop("Invalid retained rank.", call. = FALSE)
  sv <- svd(m, nu = r, nv = r)
  u <- sv$u[, seq_len(r), drop = FALSE]
  v <- sv$v[, seq_len(r), drop = FALSE]
  list(
    M = m,
    U = u,
    V = v,
    emb_L = scbm_row_normalize(u),
    emb_R = scbm_row_normalize(v),
    proj_L = u %*% t(u),
    proj_R = v %*% t(v),
    d = as.numeric(sv$d),
    rank = r,
    Ky = as.integer(ky),
    Kz = as.integer(kz)
  )
}

scbm_projector_top <- function(x, rank) {
  x <- (scbm_safe_matrix(x) + t(scbm_safe_matrix(x))) / 2
  ee <- eigen(x, symmetric = TRUE)
  u <- ee$vectors[, seq_len(as.integer(rank)), drop = FALSE]
  u %*% t(u)
}

scbm_pisces <- function(projectors, ranks, alpha,
                        tol = 1e-5, max_iter = 1000L) {
  alpha <- as.numeric(alpha)[1L]
  ranks <- as.integer(ranks)
  s <- length(projectors)
  if (length(ranks) != s) stop("ranks length mismatch.", call. = FALSE)
  if (!is.finite(alpha) || alpha < 0) stop("alpha must be nonnegative.", call. = FALSE)
  if (alpha == 0 || s == 1L) {
    return(list(projectors = projectors, iter = 0L, diff = 0))
  }
  current <- projectors
  diff <- Inf; iter <- 0L
  while (iter < as.integer(max_iter) && diff > tol) {
    iter <- iter + 1L
    next_projectors <- vector("list", s)
    diff <- 0
    for (mm in seq_len(s)) {
      z <- projectors[[mm]]
      if (mm > 1L) z <- z + alpha * current[[mm - 1L]]
      if (mm < s)  z <- z + alpha * current[[mm + 1L]]
      next_projectors[[mm]] <- scbm_projector_top(z, ranks[mm])
      diff <- diff + sqrt(sum((next_projectors[[mm]] - current[[mm]])^2))
    }
    current <- next_projectors
  }
  list(projectors = current, iter = iter, diff = diff)
}

scbm_embeddings_from_projectors <- function(projectors, ranks) {
  Map(function(p, r) {
    ee <- eigen((p + t(p)) / 2, symmetric = TRUE)
    scbm_row_normalize(ee$vectors[, seq_len(as.integer(r)), drop = FALSE])
  }, projectors, as.list(as.integer(ranks)))
}

scbm_embeddings <- function(blocks, n_comm, alpha = 0) {
  pp <- scbm_pairs(n_comm, length(blocks))
  raw <- Map(function(phi, ky, kz) scbm_svd_stage(phi, ky, kz),
             blocks, pp[, "Ky"], pp[, "Kz"])
  ranks <- vapply(raw, function(z) as.integer(z$rank), integer(1L))
  if (alpha > 0) {
    sm_l <- scbm_pisces(lapply(raw, `[[`, "proj_L"), ranks, alpha)
    sm_r <- scbm_pisces(lapply(raw, `[[`, "proj_R"), ranks, alpha)
    emb_l <- scbm_embeddings_from_projectors(sm_l$projectors, ranks)
    emb_r <- scbm_embeddings_from_projectors(sm_r$projectors, ranks)
  } else {
    emb_l <- lapply(raw, `[[`, "emb_L")
    emb_r <- lapply(raw, `[[`, "emb_R")
  }
  list(L = emb_l, R = emb_r, raw = raw, ranks = ranks)
}

scbm_cluster_chain <- function(blocks, n_comm, alpha = 0,
                               cyclic = FALSE, nstart = 50L, seed = 1L) {
  pp <- scbm_pairs(n_comm, length(blocks))
  emb <- scbm_embeddings(blocks, n_comm, alpha)
  s <- nrow(pp)
  group_l <- vector("list", s)
  group_r <- vector("list", s)
  seed <- scbm_normalize_seed(seed)

  if (!cyclic) {
    group_l[[1L]] <- scbm_kmeans(emb$L[[1L]], pp[1L, "Ky"], nstart,
                                 scbm_make_seed(seed, 1L, 1L))
    group_r[[s]] <- scbm_kmeans(emb$R[[s]], pp[s, "Kz"], nstart,
                                scbm_make_seed(seed, s, 2L))
  }

  pairs <- if (cyclic) seq_len(s) else seq_len(s - 1L)
  for (mm in pairs) {
    nn <- if (mm == s) 1L else mm + 1L
    kr <- pp[mm, "Kz"]; kl <- pp[nn, "Ky"]
    if (kr == kl) {
      lab <- scbm_kmeans(
        cbind(emb$R[[mm]], emb$L[[nn]]), kr, nstart,
        scbm_make_seed(seed, mm, nn, kr)
      )
      group_r[[mm]] <- lab
      group_l[[nn]] <- lab
    } else {
      group_r[[mm]] <- scbm_kmeans(
        emb$R[[mm]], kr, nstart, scbm_make_seed(seed, mm, 31L)
      )
      group_l[[nn]] <- scbm_kmeans(
        emb$L[[nn]], kl, nstart, scbm_make_seed(seed, nn, 47L)
      )
    }
  }

  list(
    group_L = group_l, group_R = group_r,
    emb_L = emb$L, emb_R = emb$R,
    ranks = emb$ranks, alpha = as.numeric(alpha)
  )
}

# Empirical PVAR clustering in the exact production call order.
# The RNG is initialized once; consecutive k-means calls consume the same
# stream as the manuscript code.
scbm_cluster_pvar_empirical <- function(blocks, n_comm, alpha = 0,
                                         nstart = 50L, seed = 1L) {
  pp <- scbm_pairs(n_comm, length(blocks))
  emb <- scbm_embeddings(blocks, n_comm, alpha)
  s <- nrow(pp)
  group_l <- vector("list", s)
  group_r <- vector("list", s)
  set.seed(scbm_normalize_seed(seed))

  # Cyclic closure first: receiving state in season s and sending state in season 1.
  ky1 <- pp[1L, "Ky"]
  kzs <- pp[s, "Kz"]
  if (ky1 == kzs) {
    lab <- scbm_kmeans(cbind(emb$R[[s]], emb$L[[1L]]), ky1, nstart)
    group_r[[s]] <- lab
    group_l[[1L]] <- lab
  } else {
    group_r[[s]] <- scbm_kmeans(emb$R[[s]], kzs, nstart)
    group_l[[1L]] <- scbm_kmeans(emb$L[[1L]], ky1, nstart)
  }

  # Remaining linked pairs: R_{m-1}=L_m, m=2,...,s.
  if (s >= 2L) {
    for (mm in 2:s) {
      ky <- pp[mm, "Ky"]
      kz_prev <- pp[mm - 1L, "Kz"]
      if (ky == kz_prev) {
        lab <- scbm_kmeans(cbind(emb$R[[mm - 1L]], emb$L[[mm]]), ky, nstart)
        group_r[[mm - 1L]] <- lab
        group_l[[mm]] <- lab
      } else {
        group_r[[mm - 1L]] <- scbm_kmeans(emb$R[[mm - 1L]], kz_prev, nstart)
        group_l[[mm]] <- scbm_kmeans(emb$L[[mm]], ky, nstart)
      }
    }
  }

  list(
    group_L = group_l, group_R = group_r,
    emb_L = emb$L, emb_R = emb$R,
    ranks = emb$ranks, alpha = as.numeric(alpha)
  )
}

scbm_cluster_vhar_empirical <- function(blocks_short_to_long, alpha,
                                         state_counts = c(3L, 3L, 3L),
                                         nstart = 50L, seed = 1L) {
  state_counts <- as.integer(state_counts)
  if (length(blocks_short_to_long) != 3L || length(state_counts) != 3L) {
    stop("VHAR empirical clustering requires three horizons.", call. = FALSE)
  }

  blocks <- rev(blocks_short_to_long)  # long, medium, short
  n_comm <- rep(state_counts, each = 2L)
  emb <- scbm_embeddings(blocks, n_comm, alpha)
  set.seed(scbm_normalize_seed(seed))

  long <- scbm_kmeans(
    cbind(emb$L[[1L]], emb$R[[1L]], emb$L[[2L]]),
    state_counts[1L], nstart
  )
  medium <- scbm_kmeans(
    cbind(emb$R[[2L]], emb$L[[3L]]),
    state_counts[2L], nstart
  )
  short <- scbm_kmeans(emb$R[[3L]], state_counts[3L], nstart)

  list(
    states = list(long, medium, short),
    group_L = emb$L, group_R = emb$R,
    ranks = emb$ranks, alpha = as.numeric(alpha),
    aggregate_definition = data.frame(
      state = c("Long", "Medium", "Short"),
      embeddings = c(
        "monthly.L + monthly.R + weekly.L",
        "weekly.R + daily.L",
        "daily.R"
      ), stringsAsFactors = FALSE
    )
  )
}

scbm_uniform_split <- function(m, rho = 0.20) {
  m <- scbm_safe_matrix(m)
  if (nrow(m) != ncol(m)) stop("Matrix must be square.", call. = FALSE)
  off <- which(row(m) != col(m))
  nh <- min(length(off) - 1L, max(1L, ceiling(rho * length(off))))
  idx <- sample(off, nh, replace = FALSE)
  train <- m; train[idx] <- 0
  list(original = m, train = train, hold = idx)
}

scbm_select_alpha <- function(blocks, n_comm,
                              rho = 0.20, repeats = 10L,
                              alpha_grid = scbm_alpha_grid(), seed = 1L) {
  pp <- scbm_pairs(n_comm, length(blocks))
  ranks <- as.integer(pmin(pp[, "Ky"], pp[, "Kz"]))
  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  alpha_grid <- alpha_grid[is.finite(alpha_grid) & alpha_grid >= 0]
  if (!length(alpha_grid) || max(alpha_grid) >= scbm_alpha_max()) {
    stop("alpha_grid must be nonempty and strictly below alpha_max.", call. = FALSE)
  }
  set.seed(scbm_normalize_seed(seed))
  graph_blocks <- lapply(blocks, function(phi) t(scbm_safe_matrix(phi)))
  bank <- vector("list", as.integer(repeats))
  denominator <- 0
  for (rr in seq_len(as.integer(repeats))) {
    splits <- lapply(graph_blocks, scbm_uniform_split, rho = rho)
    proj <- Map(function(sp, r) {
      sv <- svd(sp$train, nu = r, nv = r)
      u <- sv$u[, seq_len(r), drop = FALSE]
      v <- sv$v[, seq_len(r), drop = FALSE]
      list(L = u %*% t(u), R = v %*% t(v))
    }, splits, as.list(ranks))
    bank[[rr]] <- list(
      split = splits,
      L = lapply(proj, `[[`, "L"),
      R = lapply(proj, `[[`, "R")
    )
    denominator <- denominator + sum(vapply(splits, function(sp) {
      sum(sp$original[sp$hold]^2)
    }, numeric(1L)))
  }
  if (!is.finite(denominator) || denominator <= 1e-14) {
    stop("Held-out energy is too small for alpha selection.", call. = FALSE)
  }
  curve <- lapply(alpha_grid, function(alpha) {
    numerator <- 0
    for (rr in seq_along(bank)) {
      sm_l <- scbm_pisces(bank[[rr]]$L, ranks, alpha)$projectors
      sm_r <- scbm_pisces(bank[[rr]]$R, ranks, alpha)$projectors
      for (mm in seq_along(blocks)) {
        sp <- bank[[rr]]$split[[mm]]
        pred <- sm_l[[mm]] %*% sp$train %*% sm_r[[mm]]
        numerator <- numerator + sum((sp$original[sp$hold] - pred[sp$hold])^2)
      }
    }
    data.frame(alpha = alpha, cv = numerator / denominator,
               numerator = numerator, denominator = denominator)
  })
  curve <- do.call(rbind, curve)
  best <- which.min(curve$cv)
  list(alpha = curve$alpha[best], curve = curve,
       selected_at_zero = best == 1L,
       selected_at_upper = best == nrow(curve),
       rho = rho, repeats = as.integer(repeats))
}

scbm_scree <- function(blocks, labels = NULL) {
  labels <- labels %||% paste0("stage", seq_along(blocks))
  do.call(rbind, Map(function(phi, lab) {
    d <- svd(t(scbm_safe_matrix(phi)), nu = 0L, nv = 0L)$d
    prop <- d^2 / sum(d^2)
    data.frame(
      block = lab, index = seq_along(d), singular_value = d,
      prop = prop, cumprop = cumsum(prop), stringsAsFactors = FALSE
    )
  }, blocks, labels))
}

scbm_plot_scree <- function(scree, file, threshold, width = 11, height = 8) {
  blocks <- unique(as.character(scree$block))
  nr <- if (length(blocks) <= 3L) 1L else 2L
  nc <- ceiling(length(blocks) / nr)
  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2.2, 1))
  for (bb in blocks) {
    d <- scree[scree$block == bb, , drop = FALSE]
    graphics::plot(
      d$index, d$cumprop, type = "b", pch = 19, lty = 1, col = "blue",
      ylim = c(0, 1), xlab = "Mode index", ylab = "Cumulative proportion",
      main = bb
    )
    graphics::lines(d$index, d$prop, type = "b", pch = 1, lty = 3, col = "red")
    graphics::abline(h = threshold, lty = 3, col = "black")
  }
  invisible(file)
}

scbm_flow_table <- function(states, names = NULL, stage_names = NULL) {
  mat <- do.call(cbind, lapply(states, as.integer))
  if (!is.null(stage_names)) colnames(mat) <- stage_names
  if (!is.null(names)) rownames(mat) <- names
  as.data.frame(mat, check.names = FALSE)
}

scbm_sankey_long <- function(flow, stage_keys) {
  flow <- as.data.frame(flow, check.names = FALSE)
  out <- vector("list", nrow(flow) * length(stage_keys))
  ii_out <- 0L
  for (ii in seq_len(nrow(flow))) {
    for (jj in seq_along(stage_keys)) {
      ii_out <- ii_out + 1L
      out[[ii_out]] <- data.frame(
        series = rownames(flow)[ii],
        x = stage_keys[jj],
        next_x = if (jj < length(stage_keys)) stage_keys[jj + 1L] else NA_character_,
        node = as.character(flow[ii, jj]),
        next_node = if (jj < length(stage_keys)) as.character(flow[ii, jj + 1L]) else NA_character_,
        flow_group = if (jj < length(stage_keys)) paste0("Stg", jj, "_", flow[ii, jj]) else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, out)
  kmax <- max(as.integer(as.matrix(flow)))
  out$x <- factor(out$x, levels = stage_keys)
  out$next_x <- factor(out$next_x, levels = stage_keys)
  out$node <- factor(out$node, levels = as.character(seq_len(kmax)))
  out$next_node <- factor(out$next_node, levels = as.character(seq_len(kmax)))
  out
}

scbm_sankey_stage_layout <- function(v, x_pos, box_col_map,
                                      group_levels, gap,
                                      inner_pad = 0.55,
                                      box_halfwidth = 0.16) {
  cnt <- table(factor(v, levels = group_levels))
  present <- as.integer(names(cnt)[cnt > 0])
  total_height <- sum(cnt[cnt > 0]) + gap * (length(present) - 1L)
  current <- -total_height / 2
  y_ranges <- setNames(vector("list", length(group_levels)), as.character(group_levels))
  box_df <- data.frame()

  for (gg in present) {
    ng <- as.integer(cnt[as.character(gg)])
    ymin <- current
    ymax <- current + ng
    y_ranges[[as.character(gg)]] <- c(ymin + inner_pad, ymax - inner_pad)
    fill_g <- unname(box_col_map[as.character(gg)])
    if (!length(fill_g) || is.na(fill_g)) {
      stop("Missing Sankey box color for community ", gg, ".", call. = FALSE)
    }
    box_df <- rbind(
      box_df,
      data.frame(
        x = x_pos, xmin = x_pos - box_halfwidth, xmax = x_pos + box_halfwidth,
        ymin = ymin, ymax = ymax, fill = fill_g, community = gg,
        stringsAsFactors = FALSE
      )
    )
    current <- ymax + gap
  }
  list(y_ranges = y_ranges, box_df = box_df)
}

scbm_sankey_label_layers <- function(x_pos, col_idx, y_ranges, mat, text_size) {
  layers <- list()
  for (gg in names(y_ranges)) {
    rng <- y_ranges[[gg]]
    if (is.null(rng)) next
    idx <- which(as.integer(mat[, col_idx]) == as.integer(gg))
    if (!length(idx)) next
    labels <- sort(as.character(rownames(mat)[idx]))
    yy <- if (length(labels) == 1L) mean(rng) else {
      seq(from = rng[2L], to = rng[1L], length.out = length(labels))
    }
    layers[[length(layers) + 1L]] <- ggplot2::annotate(
      "text", x = x_pos, y = yy, label = labels,
      colour = "black", size = text_size, fontface = 2, hjust = 0.5
    )
  }
  layers
}

scbm_save_manuscript_sankey <- function(
    flow, file, stage_labels, box_palettes, gaps,
    flow_alpha = 0.50, text_size = 3.2,
    width = 12, height = 7, title = NULL) {

  scbm_require_packages(c("ggplot2", "ggsankey"))
  flow <- as.data.frame(flow, check.names = FALSE)
  if (!nrow(flow) || ncol(flow) < 2L) stop("Sankey flow table is empty.", call. = FALSE)
  if (is.null(rownames(flow))) stop("Sankey flow table must have series row names.", call. = FALSE)
  if (length(stage_labels) != ncol(flow) || length(box_palettes) != ncol(flow)) {
    stop("Sankey stage metadata does not match the flow table.", call. = FALSE)
  }
  if (length(gaps) == 1L) gaps <- rep(gaps, ncol(flow))
  if (length(gaps) != ncol(flow)) stop("gaps length mismatch.", call. = FALSE)

  # Alphabetical order makes ribbon stacking and printed labels identical.
  flow <- flow[order(rownames(flow)), , drop = FALSE]
  stage_keys <- paste0("stage", seq_len(ncol(flow)))
  mat <- as.matrix(flow)
  storage.mode(mat) <- "integer"
  colnames(mat) <- stage_keys
  long <- scbm_sankey_long(as.data.frame(mat, check.names = FALSE), stage_keys)

  layouts <- lapply(seq_len(ncol(mat)), function(jj) {
    scbm_sankey_stage_layout(
      mat[, jj], x_pos = jj, box_col_map = box_palettes[[jj]],
      group_levels = seq_len(max(mat[, jj])), gap = gaps[jj]
    )
  })
  box_df <- do.call(rbind, lapply(layouts, `[[`, "box_df"))

  flow_cols <- unlist(lapply(seq_len(ncol(mat) - 1L), function(jj) {
    kk <- seq_len(max(mat[, jj]))
    stats::setNames(
      unname(box_palettes[[jj]][as.character(kk)]),
      paste0("Stg", jj, "_", kk)
    )
  }))

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = x, next_x = next_x, node = node, next_node = next_node,
      fill = flow_group
    )
  ) +
    ggplot2::geom_rect(
      data = box_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = box_df$fill, color = NA
    ) +
    ggsankey::geom_sankey(
      flow.alpha = flow_alpha, node.color = NA, width = 0.001,
      show.legend = FALSE
    ) +
    ggplot2::geom_rect(
      data = box_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = NA, color = "grey40", linewidth = 0.40
    ) +
    ggplot2::scale_fill_manual(
      values = flow_cols, guide = "none", na.value = "transparent"
    ) +
    ggplot2::scale_x_discrete(labels = stats::setNames(stage_labels, stage_keys)) +
    ggsankey::theme_sankey(base_size = 14) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL, title = title) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 13, face = "bold"),
      plot.margin = ggplot2::margin(10, 40, 10, 40)
    )

  for (jj in seq_len(ncol(mat))) {
    p <- p + scbm_sankey_label_layers(
      x_pos = jj, col_idx = jj, y_ranges = layouts[[jj]]$y_ranges,
      mat = mat, text_size = text_size
    )
  }

  device_fun <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  ggplot2::ggsave(
    filename = file, plot = p, width = width, height = height,
    units = "in", device = device_fun
  )
  invisible(p)
}

scbm_save_pvar_sankey <- function(flow, file) {
  flow_plot <- as.data.frame(flow, check.names = FALSE)
  if (ncol(flow_plot) != 4L) {
    stop("The manuscript PVAR Sankey requires four seasonal states.", call. = FALSE)
  }

  palettes <- list(
    c("1" = "#FDE8E8", "2" = "#E3EDF6", "3" = "#E8E8E8", "4" = "#F3E5F5"),
    c("1" = "#FBCFCE", "2" = "#C6DBEF", "3" = "#E2F0D9", "4" = "#E8D5EF"),
    c("1" = "#9ECAE1", "2" = "#F8AFA8", "3" = "#CCEBC5", "4" = "#D7BDE2"),
    c("1" = "#74C476", "2" = "#D9CE9B", "3" = "#CE93D8", "4" = "#9ECAE1")
  )

  # Visualization-only relabeling used in the manuscript figure:
  # the last seasonal state is displayed from its highest-numbered
  # community to its lowest-numbered community.  The estimated partition
  # and the exported membership table are not modified.
  k_last <- max(as.integer(flow_plot[[4L]]))
  old_to_display <- stats::setNames(
    rev(seq_len(k_last)),
    as.character(seq_len(k_last))
  )
  displayed_last <- unname(old_to_display[as.character(flow_plot[[4L]])])
  if (anyNA(displayed_last)) {
    stop("Failed to relabel the final PVAR Sankey state.", call. = FALSE)
  }
  flow_plot[[4L]] <- as.integer(displayed_last)

  old_palette <- palettes[[4L]]
  palettes[[4L]] <- stats::setNames(
    unname(old_palette[as.character(rev(seq_len(k_last)))]),
    as.character(seq_len(k_last))
  )

  scbm_save_manuscript_sankey(
    flow = flow_plot, file = file,
    stage_labels = c(
      "Q4/Q1 (Spring)", "Q1/Q2 (Summer)",
      "Q2/Q3 (Fall)", "Q3/Q4 (Winter)"
    ),
    box_palettes = palettes, gaps = rep(2.7, 4L),
    flow_alpha = 0.50, text_size = 3.2, width = 12, height = 7
  )
}

scbm_save_vhar_sankey <- function(flow, file) {
  palettes <- list(
    c("1" = "#FDE2D5", "2" = "#DDF1E4", "3" = "#D9EAF7"),
    c("1" = "#9EC5E5", "2" = "#A9D4B5", "3" = "#F6B8A6"),
    c("1" = "#6D9CCF", "2" = "#CDB4DB", "3" = "#74B68C")
  )
  scbm_save_manuscript_sankey(
    flow = flow, file = file,
    stage_labels = c("Long (Monthly)", "Medium (Weekly)", "Short (Daily)"),
    box_palettes = palettes, gaps = c(2.2, 2.5, 2.5),
    flow_alpha = 0.55, text_size = 3.7, width = 12, height = 8
  )
}

# Backward-compatible generic entry point.
scbm_save_sankey <- function(flow, file, title = NULL) {
  if (ncol(flow) == 4L) return(scbm_save_pvar_sankey(flow, file))
  if (ncol(flow) == 3L) return(scbm_save_vhar_sankey(flow, file))
  stop("Only the manuscript PVAR (4-stage) and VHAR (3-stage) Sankey layouts are supported.",
       call. = FALSE)
}

scbm_score_states <- function(true_states, estimated_states) {
  if (length(true_states) != length(estimated_states)) stop("State count mismatch.")
  met <- Map(scbm_match_labels, true_states, estimated_states)
  list(
    accuracy = vapply(met, `[[`, numeric(1L), "accuracy"),
    ari = vapply(met, `[[`, numeric(1L), "ari"),
    accuracy_overall = mean(vapply(met, `[[`, numeric(1L), "accuracy")),
    ari_overall = mean(vapply(met, `[[`, numeric(1L), "ari"))
  )
}

scbm_parallel_lapply <- function(x, fun, ncore = 1L, root = scbm_repo_root(),
                                 modules = c("scbm_core.R"), packages = character()) {
  ncore <- max(1L, as.integer(ncore))
  if (ncore == 1L) return(lapply(x, fun))
  cl <- parallel::makeCluster(ncore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(root, modules, packages) {
    for (pkg in packages) requireNamespace(pkg, quietly = TRUE)
    for (ff in modules) source(file.path(root, "R", ff), local = .GlobalEnv)
    NULL
  }, root = root, modules = modules, packages = packages)
  parallel::parLapply(cl, x, fun)
}
