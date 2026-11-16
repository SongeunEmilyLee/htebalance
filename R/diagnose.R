#' Diagnose subgroup imbalance and precision
#'
#' Check subgroup sizes, allocation, and precision.
#' @export
diagnose_imbalance <- function(data, group, treatment, outcome = NULL,
                               small_n = 50, small_share = 0.05, min_cell = 10,
                               J = 5, alpha = 0.05, power = 0.80,
                               standardize = TRUE, target_effect = NULL,
                               contrasts = "reference", reference = NULL,
                               outcome_sd = NULL, covariates = NULL) {
  check_number(small_n, "small_n", lower = 0)
  check_number(small_share, "small_share", lower = 0, upper = 1)
  check_number(min_cell, "min_cell", lower = 0)
  check_number(J, "J", lower = 1, whole = TRUE)
  check_number(alpha, "alpha", lower = 0, upper = 1)
  check_number(power, "power", lower = 0, upper = 1)
  if (!is.null(target_effect)) check_number(target_effect, "target_effect", lower = 0)
  if (!is.null(outcome_sd)) check_number(outcome_sd, "outcome_sd", lower = 0)
  check_flag(standardize, "standardize")

  cd <- clean_data(data, group, treatment, outcome, covariates)
  G <- cd$G
  Tr <- cd$T
  N <- cd$n

  n <- as.integer(table(G))
  n1 <- as.integer(tapply(Tr == 1, G, sum))
  n0 <- n - n1
  n_eff <- ifelse(n1 > 0 & n0 > 0, 4 * n1 * n0 / (n1 + n0), 0)

  out <- data.frame(group = levels(G), n = n, share = n / N,
                    n_treated = n1, n_control = n0, n_eff = n_eff,
                    n_per_fold = n / J, stringsAsFactors = FALSE)

  # MDE = (t_{1-alpha/2, df} + t_{power, df}) x SE
  mde_mult <- function(df) stats::qt(1 - alpha / 2, df) + stats::qt(power, df)
  scale_y <- NA_real_
  if (is.null(outcome)) {
    # No outcome: SD = 1 in both arms.
    ok_cells <- n1 > 0 & n0 > 0 & n1 + n0 > 2
    out$se <- ifelse(n1 > 0 & n0 > 0, sqrt(1 / pmax(n1, 1) + 1 / pmax(n0, 1)), NA_real_)
    out$df <- ifelse(ok_cells, df_dim(pmax(n, 1), ifelse(n > 0, n1 / pmax(n, 1), 0.5)),
                     NA_real_)
    out$mde <- mde_mult(out$df) * out$se
    se_basis <- "assumed"
  } else {
    Y <- cd$Y
    if (standardize) {
      scale_y <- control_sd(Y, Tr)
      Y <- Y / scale_y
    }
    dm <- dim_by_group(Y, Tr, G, level = 1 - alpha)
    out$estimate <- dm$estimate
    out$se <- dm$se
    out$ci_low <- dm$ci_low
    out$ci_high <- dm$ci_high
    out$df <- dm$df
    out$mde <- mde_mult(out$df) * out$se
    se_basis <- "estimated"
    # Estimate SD for power calculations.
    out$sd <- sapply(levels(G), function(g) {
      y1 <- Y[G == g & Tr == 1]; y0 <- Y[G == g & Tr == 0]
      if (length(y1) < 2 || length(y0) < 2) return(NA_real_)
      sqrt(((length(y1) - 1) * stats::var(y1) + (length(y0) - 1) * stats::var(y0)) /
             (length(y1) + length(y0) - 2))
    })
    sd_all <- sqrt(sum(tapply(Y, interaction(G, Tr), function(v) sum((v - mean(v))^2)), na.rm = TRUE) /
                     (length(Y) - nlevels(interaction(G, Tr, drop = TRUE))))
  }

  # Sample size needed for a target effect
  plan_g <- NULL
  plan_c <- NULL
  sd_source <- NULL
  if (!is.null(target_effect)) {
    if (!is.null(outcome_sd)) {
      sigma <- rep(outcome_sd, nrow(out)); sd_source <- "user"
    } else if (se_basis == "estimated") {
      sigma <- out$sd
      sigma[is.na(sigma)] <- sd_all
      sd_source <- "estimated"
    } else {
      sigma <- rep(1, nrow(out)); sd_source <- "assumed"
    }
    plan_g <- plan_groups(out, sigma, target_effect, alpha, power)
    pairs <- make_pairs(contrasts, reference, out$group, stats::setNames(out$n, out$group))
    if (length(pairs)) {
      plan_c <- do.call(rbind, lapply(pairs, function(pr)
        plan_contrast(pr[1], pr[2], out, plan_g, target_effect, alpha, power)))
    }
  }
  out$sd <- NULL

  out$small_n <- out$n < small_n
  out$small_share <- out$share < small_share
  out$thin_cell <- pmin(n1, n0) < min_cell
  out <- out[order(out$n, out$group), ]
  rownames(out) <- NULL

  attr(out, "N") <- N
  attr(out, "n_dropped") <- cd$n_dropped
  attr(out, "n_original") <- cd$n_original
  attr(out, "sample") <- cd$sample
  attr(out, "settings") <- list(small_n = small_n, small_share = small_share,
                                min_cell = min_cell, J = J, alpha = alpha,
                                power = power, se_basis = se_basis,
                                outcome = outcome, standardize = standardize,
                                scale_y = scale_y, target_effect = target_effect,
                                sd_source = sd_source, covariates = covariates)
  if (!is.null(plan_g)) {
    attr(out, "plan_groups") <- plan_g[match(out$group, plan_g$group), , drop = FALSE]
    rownames(attr(out, "plan_groups")) <- NULL
  }
  attr(out, "plan_contrasts") <- plan_c
  class(out) <- c("hte_imbalance", "data.frame")
  out
}

#' @export
print.hte_imbalance <- function(x, digits = 2, ...) {
  s <- attr(x, "settings")
  cat("Subgroup diagnostics (N =", attr(x, "N"), ")\n\n")
  if (attr(x, "n_dropped") > 0) {
    cat(attr(x, "n_dropped"), "rows removed for missing values.\n\n")
  }

  tab <- data.frame(Group = x$group, N = x$n,
                    Share = sprintf("%.1f%%", 100 * x$share),
                    Treated = x$n_treated, Control = x$n_control,
                    `Eff. n` = round(x$n_eff),
                    `Per fold` = round(x$n_per_fold, 1), check.names = FALSE)
  print(tab, row.names = FALSE)

  for (i in seq_len(nrow(x))) {
    flags <- c(if (x$small_n[i]) "small n",
               if (x$small_share[i]) "small share",
               if (x$thin_cell[i]) "small treatment cell")
    if (length(flags)) cat(" ", x$group[i], ": ", paste(flags, collapse = ", "), "\n", sep = "")
  }

  cat("\nDifference in means and precision\n\n")
  effects <- data.frame(Group = x$group, check.names = FALSE)
  if (s$se_basis == "estimated") {
    effects$Estimate <- round(x$estimate, digits)
    effects$SE <- round(x$se, digits)
    effects[[sprintf("%g%% CI", 100 * (1 - s$alpha))]] <-
      sprintf("[%.*f, %.*f]", digits, x$ci_low, digits, x$ci_high)
  } else {
    effects$`SE (assumed)` <- round(x$se, digits)
  }
  effects$MDE <- round(x$mde, digits)
  print(effects, row.names = FALSE)

  pg <- attr(x, "plan_groups")
  pc <- attr(x, "plan_contrasts")
  if (!is.null(pg)) {
    fmt <- function(z) ifelse(is.infinite(z), "not reachable", format(z, trim = TRUE))
    cat("\nSample size planning (target effect =", s$target_effect, ")\n\n")
    print(data.frame(Group = pg$group,
                     Power = sprintf("%.0f%%", 100 * pg$power_now),
                     Needed = fmt(pg$n_needed), Additional = fmt(pg$n_additional)),
          row.names = FALSE)
    if (!is.null(pc)) {
      cat("\nGroup contrasts\n\n")
      print(data.frame(Contrast = pc$contrast,
                       Power = sprintf("%.0f%%", 100 * pc$power_now),
                       `A needed` = fmt(pc$a_needed_if_b_fixed),
                       `A additional` = fmt(pc$a_additional_if_b_fixed),
                       `Both needed` = paste(fmt(pc$a_needed_if_both_grow),
                                             fmt(pc$b_needed_if_both_grow), sep = ", "),
                       check.names = FALSE), row.names = FALSE)
    }
  }
  cat("\nEff. n = 4 * treated * control / total.\n")
  cat("MDE uses", 100 * s$power, "% power and", 100 * s$alpha, "% alpha.\n")
  if (s$se_basis == "assumed") cat("SE and MDE assume an outcome SD of 1.\n")
  if (s$se_basis == "estimated" && s$standardize)
    cat("Effects are in control-group SD units.\n")
  invisible(x)
}
