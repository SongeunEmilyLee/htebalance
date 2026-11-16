#' Summarize subgroup treatment effects
#' @export
hte_report <- function(x, groups = NULL, digits = 2, file = NULL) {
  if (!inherits(x, "hte_compare")) stop("Use the result of compare_hte().", call. = FALSE)

  s <- x$settings
  dg <- x$diagnostics
  f <- function(z) formatC(z, digits = digits, format = "f")
  pct <- function(z) sprintf("%.1f%%", 100 * z)

  if (is.null(groups)) groups <- c(x$excluded, x$sensitivity$group)
  bad <- setdiff(groups, c(x$excluded, x$sensitivity$group))
  if (length(bad)) stop("Unknown groups: ", paste(bad, collapse = ", "), call. = FALSE)

  paras <- character(0)
  for (g in groups) {
    d <- dg[dg$group == g, ]
    p <- sprintf("%s: %s respondents (%s), with %s treated and %s control.",
                 g, format(d$n, big.mark = ","), pct(d$share),
                 d$n_treated, d$n_control)

    if (g %in% x$excluded) {
      paras <- c(paras, paste(p, "Not estimated: fewer than two observations in one arm."))
      next
    }

    m <- x$means[x$means$group == g, ]
    p <- c(p, sprintf("Difference in means: %s (SE %s, %g%% CI %s to %s).",
                      f(m$estimate), f(m$se), 100 * s$level, f(m$ci_low), f(m$ci_high)))
    p <- c(p, sprintf("MDE: %s (power = %g, alpha = %g).",
                      f(d$mde), attr(dg, "settings")$power, attr(dg, "settings")$alpha))

    est <- x$gates[g, ]
    est <- est[is.finite(est)]
    if (length(est) > 1) {
      p <- c(p, sprintf("Across estimators: %s to %s (range %s).",
                        f(min(est)), f(max(est)), f(diff(range(est)))))
      if (any(est > 0) && any(est < 0)) p <- c(p, "The estimates have different signs.")
    }
    paras <- c(paras, paste(p, collapse = " "))
  }

  tab <- NULL
  if (!is.null(x$contrasts)) {
    cn <- x$contrasts[x$contrasts$group_a %in% groups |
                      x$contrasts$group_b %in% groups, , drop = FALSE]
    if (nrow(cn)) {
      tab <- data.frame(Contrast = cn$contrast, check.names = FALSE)
      for (m in x$methods) tab[[method_labels[[m]]]] <- round(cn[[m]], digits)
      tab[[sprintf("Means %g%% CI", 100 * s$level)]] <-
        sprintf("[%s, %s]", f(cn$means_ci_low), f(cn$means_ci_high))
      tab$Range <- round(cn$range, digits)
      tab$`Sign differs` <- ifelse(cn$sign_disagree, "yes", "")
    }
  }

  note <- c(
    if (x$n_dropped > 0) sprintf("%d rows were removed for missing values.", x$n_dropped),
    if (nrow(x$failures)) "Some estimates are missing due to model failures; see $failures.",
    if (s$standardize) "Effects are in control-group SD units." else "Effects are in outcome units.",
    "Confidence intervals apply only to difference-in-means estimates.",
    "Differences across estimators show sensitivity, not estimation error."
  )
  out <- structure(list(text = paras, contrasts = tab, note = paste(note, collapse = " ")),
                   class = "hte_report")
  if (!is.null(file)) writeLines(format(out), file)
  out
}

#' @export
print.hte_report <- function(x, width = 0.9 * getOption("width"), ...) {
  for (p in x$text) cat(strwrap(p, width = width), "", sep = "\n")
  if (!is.null(x$contrasts)) {
    cat("Group contrasts\n\n")
    print(x$contrasts, row.names = FALSE)
    cat("\n")
  }
  cat(strwrap(x$note, width = width), sep = "\n")
  invisible(x)
}

#' @export
format.hte_report <- function(x, ...) {
  tab <- if (is.null(x$contrasts)) character(0) else
    c("Group contrasts", "", utils::capture.output(print(x$contrasts, row.names = FALSE)), "")
  c(paste0(x$text, "\n"), tab, x$note)
}
