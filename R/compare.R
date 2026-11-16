#' Compare prespecified group treatment effects across estimators
#'
#' Compare subgroup estimates using means, OLS, causal forests, or BART.
#' @export
compare_hte <- function(data, outcome, treatment, group, covariates = NULL,
                        methods = c("means", "ols", "cf"), contrasts = "reference",
                        reference = NULL, propensity = NULL, standardize = TRUE,
                        J = 5, seed = 1, level = 0.95, cf_args = list(),
                        bart_args = list(), min_train_cell = 2) {
  methods <- check_methods(methods)
  if ("cf" %in% methods && !requireNamespace("grf", quietly = TRUE)) {
    stop("Method 'cf' requires the grf package: install.packages(\"grf\")", call. = FALSE)
  }
  if ("bart" %in% methods && !requireNamespace("dbarts", quietly = TRUE)) {
    stop("Method 'bart' requires the dbarts package: install.packages(\"dbarts\")", call. = FALSE)
  }
  check_number(J, "J", lower = 1, whole = TRUE)
  check_number(level, "level", lower = 0, upper = 1)
  check_number(seed, "seed", whole = TRUE)
  check_number(min_train_cell, "min_train_cell", lower = -1, whole = TRUE)
  check_flag(standardize, "standardize")
  if (!is.list(cf_args) || !is.list(bart_args)) {
    stop("`cf_args` and `bart_args` must be lists.", call. = FALSE)
  }

  # Use the same sample for diagnostics and estimation.
  cd <- clean_data(data, group, treatment, outcome, covariates)
  diag <- suppressMessages(diagnose_imbalance(data, group, treatment, outcome,
                                              covariates = covariates,
                                              J = J, alpha = 1 - level,
                                              standardize = standardize))

  # Standardize using the full analysis sample.
  scale_y <- 1
  if (standardize) scale_y <- control_sd(cd$Y, cd$T)

  # Drop groups without enough treated or control observations.
  n1 <- tapply(cd$T == 1, cd$G, sum)
  n0 <- tapply(cd$T == 0, cd$G, sum)
  excluded <- names(n1)[pmin(n1, n0) < 2]
  if (length(excluded)) {
    warning("Excluded from estimation (fewer than 2 respondents in an arm): ",
            paste(excluded, collapse = ", "), call. = FALSE)
  }
  keep <- !(cd$G %in% excluded)
  if (length(setdiff(levels(cd$G), excluded)) < 2) {
    stop("Fewer than two groups have at least 2 treated and 2 control respondents.", call. = FALSE)
  }
  G <- droplevels(cd$G[keep])
  Tr <- cd$T[keep]
  Y <- cd$Y[keep]
  d <- cd$d[keep, , drop = FALSE]
  groups <- levels(G)

  sample_info <- cd$sample
  sample_info$estimated <- sample_info$group %in% groups
  sample_info$n_treated <- as.integer(n1[sample_info$group]); sample_info$n_treated[is.na(sample_info$n_treated)] <- 0L
  sample_info$n_control <- as.integer(n0[sample_info$group]); sample_info$n_control[is.na(sample_info$n_control)] <- 0L

  Y <- Y / scale_y
  Gdum <- stats::model.matrix(~ G)[, -1, drop = FALSE]
  colnames(Gdum) <- paste0("grp_", groups[-1])
  X <- if (length(covariates)) {
    stats::model.matrix(stats::as.formula(paste("~", paste(sprintf("`%s`", covariates),
                                                           collapse = " + "))), data = d)[, -1, drop = FALSE]
  } else matrix(numeric(0), nrow(d), 0)
  Z <- cbind(Gdum, X)

  if (is.null(propensity)) {
    e <- stats::ave(Tr, G, FUN = mean)
    e_note <- "treated share within each group"
  } else if (is.character(propensity)) {
    check_name(propensity, "propensity")
    if (!propensity %in% names(d)) {
      stop("Column `", propensity, "` (propensity) not found in `data`.", call. = FALSE)
    }
    e <- d[[propensity]]
    if (!is.numeric(e)) stop("Column `", propensity, "` (propensity) must be numeric.", call. = FALSE)
    e_note <- paste0("column `", propensity, "`")
  } else {
    check_number(propensity, "propensity", lower = 0, upper = 1)
    e <- rep(propensity, nrow(d)); e_note <- paste0("fixed at ", propensity)
  }
  if (any(!is.finite(e) | e <= 0 | e >= 1)) {
    stop("Treatment probabilities must lie strictly between 0 and 1.", call. = FALSE)
  }

  # Estimates
  dm <- dim_by_group(Y, Tr, G, level)
  gates <- matrix(NA_real_, length(groups), length(methods),
                  dimnames = list(groups, methods))
  gates[, "means"] <- dm$estimate

  model_methods <- setdiff(methods, "means")
  tau <- NULL
  fold_cells <- NULL
  failures <- data.frame(method = character(0), fold = integer(0), groups = character(0),
                         reason = character(0), stringsAsFactors = FALSE)
  if (length(model_methods)) {
    set.seed(seed)
    folds <- make_folds(G, Tr, J)

    # Check the training sample in each fold.
    fold_cells <- do.call(rbind, lapply(seq_len(J), function(k) {
      tr <- folds != k
      data.frame(fold = k, group = groups,
                 train_treated = as.integer(tapply(Tr[tr] == 1, G[tr], sum)[groups]),
                 train_control = as.integer(tapply(Tr[tr] == 0, G[tr], sum)[groups]),
                 test_n = as.integer(table(G[!tr])[groups]),
                 stringsAsFactors = FALSE)
    }))
    fold_cells$train_treated[is.na(fold_cells$train_treated)] <- 0L
    fold_cells$train_control[is.na(fold_cells$train_control)] <- 0L
    thin <- fold_cells[pmin(fold_cells$train_treated, fold_cells$train_control) < min_train_cell, ]
    if (nrow(thin)) {
      thin_g <- unique(thin$group)
      warning("Some training folds have fewer than ", min_train_cell,
              " treated or control respondents for: ", paste(thin_g, collapse = ", "),
              " (J = ", J, "). Model-based estimates for these groups lean heavily on ",
              "other groups and may be unstable or fail; see `$fold_cells`. ",
              "Consider a smaller J.", call. = FALSE)
    }

    tau <- matrix(NA_real_, length(Y), length(model_methods),
                  dimnames = list(NULL, model_methods))
    for (k in seq_len(J)) {
      test <- which(folds == k)
      train <- which(folds != k)
      if (!length(test)) next
      for (m in model_methods) {
        res <- tryCatch(
          switch(m,
            ols = fit_ols(Y, Tr, Z, train, test),
            cf = do.call(fit_cf, c(list(Y, Tr, Z, e, train, test), cf_args)),
            bart = do.call(fit_bart, c(list(Y, Tr, Z, train, test), bart_args))),
          error = function(err) structure(conditionMessage(err), class = "fit_error"))
        if (!inherits(res, "fit_error") && (length(res) != length(test) || any(!is.finite(res)))) {
          res <- structure("returned missing or non-finite predictions", class = "fit_error")
        }
        if (inherits(res, "fit_error")) {
          hit <- paste(sort(unique(as.character(G[test]))), collapse = ", ")
          failures[nrow(failures) + 1, ] <- list(m, k, hit, unclass(res))
          warning(method_labels[m], " failed in fold ", k, " (no estimate for: ", hit, "): ",
                  unclass(res), call. = FALSE)
          res <- rep(NA_real_, length(test))
        }
        tau[test, m] <- res
      }
    }
    for (m in model_methods) {
      # Only report complete group predictions.
      gates[, m] <- tapply(tau[, m], G, function(v) if (anyNA(v)) NA_real_ else mean(v))[groups]
    }
  }

  # Sensitivity summary
  rng <- function(v) if (sum(!is.na(v)) < 2) NA_real_ else diff(range(v, na.rm = TRUE))
  sens <- data.frame(
    group = groups,
    n = as.integer(table(G)[groups]),
    min = apply(gates, 1, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)),
    max = apply(gates, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)),
    range = apply(gates, 1, rng),
    sign_disagree = apply(gates, 1, function(v) any(v > 0, na.rm = TRUE) && any(v < 0, na.rm = TRUE)),
    row.names = NULL
  )
  sens <- sens[order(sens$n), ]
  rownames(sens) <- NULL

  # Contrasts
  pairs <- make_pairs(contrasts, reference, groups, table(G))
  con <- NULL
  if (length(pairs)) {
    con <- do.call(rbind, lapply(pairs, function(p) {
      est <- stats::setNames(as.numeric(gates[p[1], ] - gates[p[2], ]), colnames(gates))
      va <- dm$se[groups == p[1]]^2
      vb <- dm$se[groups == p[2]]^2
      se <- sqrt(va + vb)
      df <- (va + vb)^2 / (va^2 / dm$df[groups == p[1]] + vb^2 / dm$df[groups == p[2]])
      if (!is.finite(df)) df <- NA_real_
      q <- stats::qt(1 - (1 - level) / 2, df)
      data.frame(contrast = paste(p[1], "-", p[2]), group_a = p[1], group_b = p[2],
                 t(est), means_se = se, means_df = df,
                 means_ci_low = est[["means"]] - q * se, means_ci_high = est[["means"]] + q * se,
                 range = rng(est),
                 sign_disagree = any(est > 0, na.rm = TRUE) && any(est < 0, na.rm = TRUE),
                 check.names = FALSE)
    }))
    rownames(con) <- NULL
  }

  structure(list(gates = gates, means = dm, sensitivity = sens, contrasts = con,
                 diagnostics = diag, tau = tau, methods = methods,
                 excluded = excluded, sample = sample_info,
                 n_original = cd$n_original, n_dropped = cd$n_dropped,
                 fold_cells = fold_cells, failures = failures, settings = list(
                   outcome = outcome, treatment = treatment, group = group,
                   covariates = covariates, standardize = standardize,
                   scale_y = scale_y, J = J, seed = seed, level = level,
                   min_train_cell = min_train_cell,
                   propensity = e_note)),
            class = "hte_compare")
}

make_pairs <- function(contrasts, reference, groups, sizes) {
  if (is.null(contrasts) || identical(contrasts, "none")) return(list())
  if (identical(contrasts, "reference")) {
    if (is.null(reference)) reference <- names(sizes)[which.max(sizes)]
    if (!reference %in% groups) {
      stop("`reference` group \"", reference, "\" is not among the estimated groups.", call. = FALSE)
    }
    return(lapply(setdiff(groups, reference), function(g) c(g, reference)))
  }
  if (identical(contrasts, "all")) {
    cb <- utils::combn(groups, 2)
    return(lapply(seq_len(ncol(cb)), function(j) cb[, j]))
  }
  if (is.list(contrasts)) {
    for (p in contrasts) {
      if (length(p) != 2 || !all(p %in% groups)) {
        stop("Each element of `contrasts` must name two estimated groups, e.g. c(\"Asian\", \"White\").",
             call. = FALSE)
      }
    }
    return(lapply(contrasts, as.character))
  }
  stop("`contrasts` must be \"reference\", \"all\", \"none\", or a list of group pairs.", call. = FALSE)
}

# Fit each estimator

fit_ols <- function(Y, Tr, Z, train, test) {
  Zc <- sweep(Z, 2, colMeans(Z[train, , drop = FALSE]))
  M <- cbind(1, Zc[train, , drop = FALSE], Tr[train], Tr[train] * Zc[train, , drop = FALSE])
  b <- stats::lm.fit(M, Y[train])$coefficients
  p <- ncol(Z)
  bT <- b[p + 2]
  bTZ <- if (p > 0) b[(p + 3):(2 * p + 2)] else numeric(0)
  if (is.na(bT)) stop("treatment coefficient is not identified")
  # Group interactions must be identified.
  ng <- sum(grepl("^grp_", colnames(Z)))
  if (ng > 0 && anyNA(bTZ[seq_len(ng)])) stop("a group-by-treatment term is not identified")
  bTZ[is.na(bTZ)] <- 0
  drop(bT + Zc[test, , drop = FALSE] %*% bTZ)
}

fit_cf <- function(Y, Tr, Z, e, train, test, num.trees = 2000, ...) {
  cf <- grf::causal_forest(X = Z[train, , drop = FALSE], Y = Y[train], W = Tr[train],
                           W.hat = e[train], num.trees = num.trees, ...)
  drop(stats::predict(cf, newdata = Z[test, , drop = FALSE])$predictions)
}

fit_bart <- function(Y, Tr, Z, train, test, ntree = 200, ndpost = 500, nskip = 250, ...) {
  xtr <- cbind(Z[train, , drop = FALSE], treat = Tr[train])
  zt <- Z[test, , drop = FALSE]
  xte <- rbind(cbind(zt, treat = 1), cbind(zt, treat = 0))
  fit <- dbarts::bart(x.train = xtr, y.train = Y[train], x.test = xte, ntree = ntree,
                      ndpost = ndpost, nskip = nskip, verbose = FALSE, ...)
  m <- length(test)
  yhat <- colMeans(fit$yhat.test)
  yhat[seq_len(m)] - yhat[m + seq_len(m)]
}

# Print and plot methods

#' @export
print.hte_compare <- function(x, digits = 2, ...) {
  s <- x$settings
  unit <- if (s$standardize) " (control-group SD units)" else ""
  cat("Estimator sensitivity: group treatment effects by estimator", unit, "\n\n", sep = "")
  g <- x$gates[x$sensitivity$group, , drop = FALSE]
  colnames(g) <- method_labels[colnames(g)]
  out <- data.frame(Group = rownames(g), n = x$sensitivity$n, round(g, digits),
                    Range = round(x$sensitivity$range, digits),
                    `Sign differs` = ifelse(x$sensitivity$sign_disagree, "yes", ""),
                    check.names = FALSE)
  print(out, row.names = FALSE)
  if (!is.null(x$contrasts)) {
    cat("\nEstimator sensitivity: group contrasts\n\n")
    cn <- x$contrasts
    out <- data.frame(Contrast = cn$contrast,
                      round(cn[, x$methods, drop = FALSE], digits),
                      `Means 95% CI` = sprintf("[%.*f, %.*f]", digits, cn$means_ci_low,
                                               digits, cn$means_ci_high),
                      Range = round(cn$range, digits),
                      `Sign differs` = ifelse(cn$sign_disagree, "yes", ""),
                      check.names = FALSE)
    names(out)[seq_along(x$methods) + 1] <- method_labels[x$methods]
    names(out)[names(out) == "Means 95% CI"] <- sprintf("Means %g%% CI", 100 * s$level)
    print(out, row.names = FALSE)
  }
  if (length(x$excluded)) cat("\nExcluded groups:", paste(x$excluded, collapse = ", "), "\n")
  if (isTRUE(x$n_dropped > 0)) {
    cat(sprintf("\n%s of %s rows were dropped for missing or non-finite values in the outcome,\ntreatment, group or covariates; all estimates and diagnostics use the remaining rows.\n",
                format(x$n_dropped, big.mark = ","), format(x$n_original, big.mark = ",")))
  }
  if (NROW(x$failures)) {
    cat("\nEstimator failures (estimates set to NA for the affected groups):\n")
    for (i in seq_len(nrow(x$failures))) {
      f <- x$failures[i, ]
      cat(sprintf("  %s, fold %d (%s): %s\n", method_labels[[f$method]], f$fold, f$groups, f$reason))
    }
  }
  cat("\nRange = largest minus smallest estimate across estimators. It measures sensitivity",
      "\nto estimator choice, not which estimate is closer to the true effect.\n")
  invisible(x)
}

#' @export
plot.hte_compare <- function(x, type = c("gates", "contrasts"), main = NULL, ...) {
  type <- match.arg(type)
  if (type == "gates") {
    est <- x$gates[x$sensitivity$group, , drop = FALSE]
    i <- match(rownames(est), x$means$group)
    lo <- x$means$ci_low[i]
    hi <- x$means$ci_high[i]
    lab <- paste0(rownames(est), " (n=", x$sensitivity$n, ")")
    xlab <- "Group treatment effect"
  } else {
    if (is.null(x$contrasts)) stop("No contrasts were computed.", call. = FALSE)
    est <- as.matrix(x$contrasts[, x$methods, drop = FALSE])
    lo <- x$contrasts$means_ci_low
    hi <- x$contrasts$means_ci_high
    lab <- x$contrasts$contrast
    xlab <- "Difference in treatment effects"
  }
  if (x$settings$standardize) xlab <- paste(xlab, "(control-group SD)")

  cols <- c(means = "#1F4E79", ols = "#8C1D40", cf = "#222222", bart = "#7F7F7F")
  pchs <- c(means = 16, ols = 17, cf = 18, bart = 15)
  methods <- colnames(est)
  y <- rev(seq_len(nrow(est)))
  xlim <- range(c(0, est, lo, hi), na.rm = TRUE)
  pad <- max(diff(xlim) * 0.05, 0.01)
  xlim <- xlim + c(-pad, pad)

  old <- graphics::par(mar = c(7, max(6, max(nchar(lab)) * 0.6), 3, 1), las = 1)
  on.exit(graphics::par(old))
  graphics::plot(NA, xlim = xlim, ylim = c(0.5, length(y) + 0.5),
                 yaxt = "n", ylab = "", xlab = xlab, main = main, ...)
  graphics::axis(2, at = y, labels = lab)
  graphics::abline(v = 0, col = "gray60", lty = 2)
  graphics::segments(lo, y, hi, y, col = "gray70", lwd = 2)

  offsets <- if (length(methods) == 1) 0 else seq(0.15, -0.15, length.out = length(methods))
  for (j in seq_along(methods)) {
    graphics::points(est[, j], y + offsets[j], col = cols[methods[j]],
                     pch = pchs[methods[j]], cex = 1.3)
  }
  # Legend below the plot so it never covers the estimates.
  graphics::legend(mean(graphics::par("usr")[1:2]), graphics::grconvertY(0.01, "ndc", "user"),
                   legend = c(method_labels[methods], "Means CI"),
                   col = c(cols[methods], "gray70"),
                   pch = c(pchs[methods], NA), lty = c(rep(NA, length(methods)), 1),
                   xjust = 0.5, yjust = 0, horiz = TRUE, xpd = NA, bty = "n", cex = 0.8)
  invisible(x)
}
