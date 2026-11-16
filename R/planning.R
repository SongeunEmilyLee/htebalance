# Power and sample size calculations.

# Standard error of a difference in means.
se_dim <- function(n, p, sigma) sigma / sqrt(n * p * (1 - p))

power_at <- function(effect, se, df, alpha) {
  crit <- stats::qt(1 - alpha / 2, df)
  ncp <- abs(effect) / se
  out <- stats::pt(ncp - crit, df) + stats::pt(-ncp - crit, df)
  out[is.na(se) | is.na(df) | !(se > 0)] <- NA_real_
  out
}

# Find the smallest sample size with enough power.
solve_n <- function(se_fun, df_fun, effect, alpha, power, lower, upper = 1e8) {
  f <- function(x) power_at(effect, se_fun(x), df_fun(x), alpha) - power
  if (!is.finite(f(upper)) || f(upper) < 0) return(Inf)
  if (f(lower) >= 0) return(lower)
  stats::uniroot(f, c(lower, upper), tol = 1e-6)$root
}

# Power and sample size for each group.
plan_groups <- function(out, sigma, target, alpha, power) {
  p <- ifelse(out$n > 0, out$n_treated / out$n, NA_real_)
  p[!is.finite(p) | p <= 0 | p >= 1] <- 0.5   # empty arm: plan for 1:1 assignment
  se_now <- se_dim(out$n, p, sigma)
  df_now <- df_dim(pmax(out$n, 1), p)
  needed <- mapply(function(pi, si) {
    if (is.na(si)) return(NA_real_)
    ceiling(solve_n(function(m) se_dim(m, pi, si), function(m) df_dim(m, pi),
                    target, alpha, power, lower = 4))
  }, p, sigma)
  data.frame(group = out$group, treated_share = p, sd = sigma,
             power_now = power_at(target, se_now, df_now, alpha),
             n_needed = needed, n_additional = pmax(needed - out$n, 0),
             stringsAsFactors = FALSE)
}

# For a contrast, grow group A or grow both groups.
plan_contrast <- function(a, b, out, gp, target, alpha, power) {
  ra <- match(a, out$group)
  rb <- match(b, out$group)
  na <- out$n[ra]
  nb <- out$n[rb]
  pa <- gp$treated_share[ra]
  pb <- gp$treated_share[rb]
  sa <- gp$sd[ra]
  sb <- gp$sd[rb]
  se_c <- function(ma, mb) sqrt(se_dim(ma, pa, sa)^2 + se_dim(mb, pb, sb)^2)
  df_c <- function(ma, mb) {
    cells <- c(ma * pa, ma * (1 - pa), mb * pb, mb * (1 - pb))
    welch_df(c(sa, sa, sb, sb)^2 / cells, cells)
  }

  pw_now <- power_at(target, se_c(na, nb), df_c(na, nb), alpha)

  xa <- solve_n(function(m) se_c(m, nb), function(m) df_c(m, nb),
                target, alpha, power, lower = 4)
  na_needed_a <- if (is.finite(xa)) ceiling(xa) else Inf

  k <- solve_n(function(s) se_c(s * na, s * nb), function(s) df_c(s * na, s * nb),
               target, alpha, power, lower = 1)
  na_both <- if (k == 1) na else ceiling(k * na)
  nb_both <- if (k == 1) nb else ceiling(k * nb)

  data.frame(contrast = paste(a, "-", b), group_a = a, group_b = b,
             power_now = pw_now,
             a_needed_if_b_fixed = na_needed_a,
             a_additional_if_b_fixed = if (is.finite(na_needed_a)) max(na_needed_a - na, 0) else Inf,
             a_needed_if_both_grow = na_both, b_needed_if_both_grow = nb_both,
             stringsAsFactors = FALSE)
}
