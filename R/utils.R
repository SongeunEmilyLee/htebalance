# Small helper functions

method_labels <- c(means = "Subgroup means", ols = "Interacted OLS",
                   cf = "Causal forest", bart = "BART")

check_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop("`", arg, "` must be a single column name.", call. = FALSE)
  }
}

check_number <- function(x, arg, lower = -Inf, upper = Inf, whole = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) && x > lower && x < upper
  if (ok && whole) ok <- x == round(x)
  if (!ok) stop("`", arg, "` is out of range.", call. = FALSE)
}

check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

check_names <- function(x, arg) {
  if (is.null(x)) return(invisible())
  if (!is.character(x) || !length(x) || anyNA(x) || any(x == "") || anyDuplicated(x)) {
    stop("`", arg, "` must be a character vector of distinct column names.", call. = FALSE)
  }
}

# Check estimator names.
check_methods <- function(methods) {
  if (!is.character(methods) || length(methods) == 0L || anyNA(methods) ||
      any(methods == "")) {
    stop("`methods` must be a non-empty character vector of method names.", call. = FALSE)
  }
  dup <- unique(methods[duplicated(methods)])
  if (length(dup)) {
    stop("`methods` lists the same method more than once: ", paste(dup, collapse = ", "),
         ".", call. = FALSE)
  }
  bad <- setdiff(methods, names(method_labels))
  if (length(bad)) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "), ". Choose from: ",
         paste(names(method_labels), collapse = ", "), call. = FALSE)
  }
  c("means", setdiff(methods, "means"))
}

# Welch degrees of freedom.
welch_df <- function(v, n) {
  num <- sum(v)^2
  den <- sum(v^2 / pmax(n - 1, 1))
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(max(sum(n) - length(n), 1))
  max(num / den, 1)
}

# Degrees of freedom for a difference in means.
df_dim <- function(n, p) {
  as.numeric(mapply(function(ni, pi) {
    n1 <- ni * pi; n0 <- ni * (1 - pi)
    welch_df(c(1 / n1, 1 / n0), c(n1, n0))
  }, n, p))
}

# Keep complete cases and count excluded rows.
clean_data <- function(data, group, treatment, outcome = NULL, covariates = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  check_name(group, "group")
  check_name(treatment, "treatment")
  if (!is.null(outcome)) check_name(outcome, "outcome")
  check_names(covariates, "covariates")
  cols <- c(group, treatment, outcome, covariates)
  if (anyDuplicated(cols)) {
    stop("`group`, `treatment`, `outcome` and `covariates` must name different columns.",
         call. = FALSE)
  }
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("Column(s) not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # Inf is treated as missing.
  bad_val <- function(v) is.na(v) | (is.numeric(v) & !is.finite(v))
  na_by_col <- vapply(cols, function(v) sum(bad_val(data[[v]])), numeric(1))
  ok <- !Reduce(`|`, lapply(cols, function(v) bad_val(data[[v]])))
  n_dropped <- sum(!ok)
  if (n_dropped > 0) {
    why <- na_by_col[na_by_col > 0]
    message(n_dropped, " of ", nrow(data), " row(s) dropped because of missing or ",
            "non-finite values in: ",
            paste(sprintf("%s (%d)", names(why), why), collapse = ", "),
            ". All results use the remaining ", sum(ok), " rows.")
  }
  g_all <- as.character(data[[group]])
  g_lev <- sort(unique(g_all[!is.na(g_all)]))
  sample_info <- data.frame(
    group = g_lev,
    n_original = as.integer(table(factor(g_all, levels = g_lev))),
    n_analysis = as.integer(table(factor(g_all[ok], levels = g_lev))),
    stringsAsFactors = FALSE)
  sample_info$n_dropped <- sample_info$n_original - sample_info$n_analysis
  d <- data[ok, , drop = FALSE]

  Tr <- d[[treatment]]
  if (is.logical(Tr)) Tr <- as.integer(Tr)
  if (is.factor(Tr) || is.character(Tr) || !all(Tr %in% c(0, 1))) {
    vals <- sort(unique(as.character(d[[treatment]])))
    stop("`", treatment, "` must be coded 0/1 or TRUE/FALSE. Found: ",
         paste(vals[seq_len(min(5, length(vals)))], collapse = ", "), call. = FALSE)
  }
  if (length(unique(Tr)) < 2) {
    stop("`", treatment, "` has only one value; both arms are needed.", call. = FALSE)
  }

  G <- as.factor(d[[group]])
  empty <- levels(G)[table(G) == 0]
  if (length(empty)) message("Dropping empty group level(s): ", paste(empty, collapse = ", "), ".")
  G <- droplevels(G)
  if (nlevels(G) < 2) stop("`", group, "` has fewer than two groups.", call. = FALSE)

  Y <- NULL
  if (!is.null(outcome)) {
    Y <- d[[outcome]]
    if (is.logical(Y)) Y <- as.numeric(Y)
    if (!is.numeric(Y)) stop("`", outcome, "` must be numeric or logical.", call. = FALSE)
  }

  list(d = d, G = G, T = Tr, Y = Y, n = nrow(d), n_dropped = n_dropped,
       n_original = nrow(data), na_by_col = na_by_col, sample = sample_info)
}

# Scale effects by the control-group SD.
control_sd <- function(Y, Tr) {
  s <- stats::sd(Y[Tr == 0])
  if (!is.finite(s) || s == 0) {
    stop("The outcome has no variation in the control group; cannot standardize. ",
         "Use `standardize = FALSE`.", call. = FALSE)
  }
  s
}

# Make folds within each group and treatment arm.
make_folds <- function(G, Tr, J) {
  strata <- interaction(G, Tr, drop = TRUE)
  f <- integer(length(G))
  for (s in levels(strata)) {
    idx <- which(strata == s)
    f[idx] <- if (length(idx) == 1) sample.int(J, 1) else sample(rep_len(seq_len(J), length(idx)))
  }
  f
}

# Estimate a treatment effect for each group.
dim_by_group <- function(Y, Tr, G, level = 0.95) {
  out <- lapply(levels(G), function(g) {
    y1 <- Y[G == g & Tr == 1]; y0 <- Y[G == g & Tr == 0]
    n1 <- length(y1); n0 <- length(y0)
    est <- if (n1 > 0 && n0 > 0) mean(y1) - mean(y0) else NA_real_
    if (n1 > 1 && n0 > 1) {
      v1 <- stats::var(y1) / n1; v0 <- stats::var(y0) / n0
      se <- sqrt(v1 + v0)
      df <- if (se > 0) (v1 + v0)^2 / (v1^2 / (n1 - 1) + v0^2 / (n0 - 1)) else n1 + n0 - 2
      q <- stats::qt(1 - (1 - level) / 2, df)
    } else {
      se <- NA_real_; q <- NA_real_; df <- NA_real_
    }
    data.frame(group = g, estimate = est, se = se, df = df,
               ci_low = est - q * se, ci_high = est + q * se)
  })
  do.call(rbind, out)
}
