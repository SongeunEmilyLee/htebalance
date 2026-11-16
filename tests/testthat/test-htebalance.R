toy <- function(seed = 1) {
  set.seed(seed)
  n <- c(White = 600, Black = 200, Latino = 100, Asian = 40)
  race <- rep(names(n), n)
  treat <- unlist(lapply(n, function(m) sample(rep(0:1, length.out = m))))
  age <- rnorm(sum(n))
  eff <- c(White = 0.2, Black = 0.5, Latino = 0.3, Asian = -0.2)[race]
  data.frame(race, treat, age, y = 0.5 * age + treat * eff + rnorm(sum(n)))
}

# diagnose_imbalance -----------------------------------------------------------

test_that("design-only diagnostics are correct", {
  d <- diagnose_imbalance(toy(), "race", "treat")
  expect_s3_class(d, "hte_imbalance")
  expect_equal(d$group, c("Asian", "Latino", "Black", "White"))
  expect_equal(d$n_treated[1] + d$n_control[1], 40)
  expect_equal(d$n_eff[1], 40)                        # balanced arms
  expect_equal(d$se[1], sqrt(1 / 20 + 1 / 20))
  expect_equal(d$mde[1], (qt(0.975, 38) + qt(0.8, 38)) * d$se[1])
  expect_equal(d$small_n, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(d$small_share, c(TRUE, FALSE, FALSE, FALSE))
})

test_that("effective n falls when arms are unbalanced", {
  x <- toy(); a <- which(x$race == "Asian")
  x$treat[a] <- c(rep(1, 30), rep(0, 10))
  d <- diagnose_imbalance(x, "race", "treat")
  expect_equal(d$n_eff[d$group == "Asian"], 4 * 30 * 10 / 40)
})

test_that("outcome-based estimates match hand calculation", {
  x <- toy()
  d <- diagnose_imbalance(x, "race", "treat", outcome = "y", standardize = FALSE)
  a <- x[x$race == "Asian", ]
  expect_equal(d$estimate[1], mean(a$y[a$treat == 1]) - mean(a$y[a$treat == 0]))
  expect_equal(d$se[1], sqrt(var(a$y[a$treat == 1]) / 20 + var(a$y[a$treat == 0]) / 20))
  expect_true(d$ci_low[1] < d$estimate[1] && d$estimate[1] < d$ci_high[1])
  ds <- diagnose_imbalance(x, "race", "treat", outcome = "y")
  expect_equal(ds$estimate[1], d$estimate[1] / sd(x$y[x$treat == 0]))
})

test_that("inputs are validated", {
  x <- toy()
  x2 <- x; x2$race[1:3] <- NA
  expect_message(diagnose_imbalance(x2, "race", "treat"), "3 of 940 row")
  x3 <- x; x3$treat <- ifelse(x3$treat == 1, "T", "C")
  expect_error(diagnose_imbalance(x3, "race", "treat"), "0/1")
  x4 <- x; x4$treat <- 1
  expect_error(diagnose_imbalance(x4, "race", "treat"), "only one value")
  x5 <- x; x5$race <- factor(x5$race, levels = c(unique(x$race), "Other"))
  expect_message(diagnose_imbalance(x5, "race", "treat"), "Other")
  expect_error(diagnose_imbalance(x, "ethnicity", "treat"), "not found")
  expect_error(diagnose_imbalance(x, "race", "treat", small_share = 5), "out of range")
})

# compare_hte -------------------------------------------------------------------

test_that("subgroup means in compare_hte match diagnose_imbalance", {
  x <- toy()
  r <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"))
  d <- diagnose_imbalance(x, "race", "treat", outcome = "y")
  expect_equal(unname(r$gates["Asian", "means"]), d$estimate[d$group == "Asian"])
  expect_false(anyNA(r$gates))
})

test_that("interacted OLS without covariates stays close to subgroup means", {
  x <- toy()
  r <- compare_hte(x, "y", "treat", "race", methods = c("means", "ols"), J = 2)
  # With only group dummies, each fold's OLS fit is the training-fold difference in
  # means; the cross-fitted GATE is the average over folds, so it is close, not equal.
  expect_true(all(abs(r$gates[, "ols"] - r$gates[, "means"]) < 0.5))
})

test_that("contrasts are differences of GATEs", {
  x <- toy()
  r <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"),
                   contrasts = list(c("Asian", "White"), c("Black", "Latino")))
  expect_equal(r$contrasts$contrast, c("Asian - White", "Black - Latino"))
  expect_equal(r$contrasts$ols[1], unname(r$gates["Asian", "ols"] - r$gates["White", "ols"]))
  r2 <- compare_hte(x, "y", "treat", "race", methods = "means")
  expect_true(all(r2$contrasts$group_b == "White"))       # default reference = largest group
})

test_that("a group with an empty arm is excluded, not filled with pooled values", {
  x <- toy(); x$treat[x$race == "Asian"] <- 0
  expect_warning(r <- compare_hte(x, "y", "treat", "race", methods = c("means", "ols")),
                 "Asian")
  expect_false("Asian" %in% rownames(r$gates))
  expect_equal(r$excluded, "Asian")
})

test_that("report states the facts it is given", {
  x <- toy()
  r <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"))
  rep <- hte_report(r, groups = "Asian")
  expect_s3_class(rep, "hte_report")
  expect_true(grepl("40 respondents", rep$text[1]))
  expect_true(grepl("20 treated and 20 control", rep$text[1]))
  expect_true(grepl("MDE:", rep$text[1]))
  expect_equal(rep$contrasts$Contrast, "Asian - White")
  expect_true(grepl("sensitivity, not estimation error", rep$note))
})

# Validation against known truth and invariances ----------------------------------

test_that("subgroup means recover the true group effects in large samples", {
  set.seed(42)
  n <- c(White = 20000, Black = 20000, Latino = 20000, Asian = 20000)
  truth <- c(White = 0.20, Black = 0.45, Latino = 0.30, Asian = -0.10)
  race <- rep(names(n), n)
  treat <- unlist(lapply(n, function(m) sample(rep(0:1, length.out = m))))
  y <- treat * truth[race] + rnorm(sum(n))
  r <- compare_hte(data.frame(y, treat, race), "y", "treat", "race",
                   methods = c("means", "ols"), standardize = FALSE)
  expect_true(all(abs(r$gates[names(truth), "means"] - truth) < 0.05))
  expect_true(all(abs(r$gates[names(truth), "ols"] - truth) < 0.05))
})

test_that("standardizing divides every linear estimate by the control-group SD", {
  x <- toy()
  raw <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"),
                     standardize = FALSE)
  std <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"))
  s0 <- sd(x$y[x$treat == 0])
  expect_equal(std$gates, raw$gates / s0)
  expect_equal(std$contrasts$means, raw$contrasts$means / s0)
  expect_equal(std$means$se, raw$means$se / s0)
})

test_that("the same seed reproduces the results", {
  x <- toy()
  a <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"), seed = 7)
  b <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"), seed = 7)
  expect_equal(a$gates, b$gates)
  expect_equal(a$tau, b$tau)
  skip_if_not_installed("grf")
  a <- compare_hte(x, "y", "treat", "race", "age", methods = "cf", seed = 7,
                   cf_args = list(num.trees = 200))
  b <- compare_hte(x, "y", "treat", "race", "age", methods = "cf", seed = 7,
                   cf_args = list(num.trees = 200))
  expect_equal(a$gates, b$gates)
})

test_that("every contrast equals the difference of the corresponding GATEs", {
  x <- toy()
  r <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"),
                   contrasts = "all")
  expect_equal(nrow(r$contrasts), choose(4, 2))
  for (i in seq_len(nrow(r$contrasts))) {
    a <- r$contrasts$group_a[i]; b <- r$contrasts$group_b[i]
    for (m in c("means", "ols")) {
      expect_equal(r$contrasts[[m]][i], unname(r$gates[a, m] - r$gates[b, m]))
    }
  }
})

test_that("an empty arm triggers a warning in compare_hte and a flag in diagnose_imbalance", {
  x <- toy(); x$treat[x$race == "Latino"] <- 1
  expect_warning(compare_hte(x, "y", "treat", "race", methods = "means"), "Latino")
  d <- diagnose_imbalance(x, "race", "treat", outcome = "y")
  expect_true(d$thin_cell[d$group == "Latino"])
  expect_true(is.na(d$estimate[d$group == "Latino"]))
})

test_that("the MDE matches the effect an exact two-sample t-test detects with the stated power", {
  # Equal arms, common SD: the Welch df equals n1 + n0 - 2 when the sample SDs match,
  # so the t-approximation MDE should be within 0.5% of power.t.test (noncentral t).
  set.seed(3)
  sigma <- 1.3; m <- 20
  y1 <- as.numeric(scale(rnorm(m))) * sigma; y0 <- as.numeric(scale(rnorm(m))) * sigma
  x <- data.frame(g = rep(c("A", "B"), each = 2 * m),
                  treat = rep(rep(1:0, each = m), 2), y = c(y1, y0, y1, y0))
  d <- diagnose_imbalance(x, "g", "treat", outcome = "y", standardize = FALSE)
  exact <- power.t.test(n = m, sd = sigma, power = 0.8)$delta
  expect_equal(d$mde[1], exact, tolerance = 0.005)
})

# Precision planning ---------------------------------------------------------------

test_that("power at the MDE equals the target power", {
  d <- diagnose_imbalance(toy(), "race", "treat", target_effect = 0.2)
  r <- d[d$group == "Asian", ]
  # The MDE ignores the far rejection tail; power_at includes it, so allow a tiny gap.
  expect_equal(htebalance:::power_at(r$mde, r$se, r$df, 0.05), 0.8, tolerance = 1e-4)
})

test_that("required n is the smallest group size reaching the target power", {
  d <- diagnose_imbalance(toy(), "race", "treat", target_effect = 0.2)
  pg <- attr(d, "plan_groups")
  n_req <- pg$n_needed[pg$group == "Asian"]
  pw <- function(n) htebalance:::power_at(0.2, htebalance:::se_dim(n, 0.5, 1), n - 2, 0.05)
  expect_true(pw(n_req) >= 0.8)
  expect_true(pw(n_req - 1) < 0.8)
  # Matches the textbook approximation for 1:1 assignment, SD = 1, d = 0.2: about 788.
  expect_true(abs(n_req - 788) <= 3)
  expect_equal(pg$n_additional[pg$group == "Asian"], n_req - 40)
})

test_that("planning uses the estimated SD or a user-supplied SD", {
  x <- toy()
  d1 <- diagnose_imbalance(x, "race", "treat", outcome = "y", target_effect = 0.2)
  d2 <- diagnose_imbalance(x, "race", "treat", outcome = "y", target_effect = 0.2, outcome_sd = 1)
  expect_equal(attr(d1, "settings")$sd_source, "estimated")
  expect_equal(attr(d2, "plan_groups")$sd, rep(1, 4))
})

test_that("contrast planning: fixed reference can be unreachable, growing both cannot", {
  # Large target relative to the reference group's own precision: adding to A alone
  # cannot push the contrast SE below what the target needs.
  n <- c(White = 60, Asian = 40)
  x <- data.frame(race = rep(names(n), n),
                  treat = unlist(lapply(n, function(m) rep(0:1, length.out = m))))
  d <- diagnose_imbalance(x, "race", "treat", target_effect = 0.3)
  pc <- attr(d, "plan_contrasts")
  expect_equal(pc$contrast, "Asian - White")
  expect_true(is.infinite(pc$a_needed_if_b_fixed))
  expect_true(is.finite(pc$a_needed_if_both_grow))
  # With a big reference group, adding to A alone works and needs more A than "both grow".
  n2 <- c(White = 3000, Asian = 40)
  x2 <- data.frame(race = rep(names(n2), n2),
                   treat = unlist(lapply(n2, function(m) rep(0:1, length.out = m))))
  pc2 <- attr(diagnose_imbalance(x2, "race", "treat", target_effect = 0.3), "plan_contrasts")
  expect_true(is.finite(pc2$a_needed_if_b_fixed))
  se_c <- sqrt(htebalance:::se_dim(pc2$a_needed_if_b_fixed, 0.5, 1)^2 +
                 htebalance:::se_dim(3000, 0.5, 1)^2)
  expect_true(htebalance:::power_at(0.3, se_c, pc2$a_needed_if_b_fixed + 3000 - 4, 0.05) >= 0.8)
})

# Edge cases added before the public release ----------------------------------------

test_that("diagnostics and estimates use the same rows when covariates are missing", {
  x <- toy()
  a <- which(x$race == "Asian")
  x$age[a[1:10]] <- NA                       # 10 of 40 Asian respondents lack age
  r <- suppressMessages(compare_hte(x, "y", "treat", "race", "age",
                                    methods = c("means", "ols")))
  dg <- r$diagnostics
  expect_equal(dg$n[dg$group == "Asian"], 30)
  expect_equal(r$sensitivity$n[r$sensitivity$group == "Asian"], 30)
  expect_equal(dg$estimate[dg$group == "Asian"], unname(r$gates["Asian", "means"]))
  expect_equal(dg$se[dg$group == "Asian"], r$means$se[r$means$group == "Asian"])
  expect_equal(r$n_dropped, 10)
  expect_equal(r$sample$n_original[r$sample$group == "Asian"], 40)
  expect_equal(r$sample$n_analysis[r$sample$group == "Asian"], 30)
  # Same standardizing SD in both places.
  expect_equal(attr(dg, "settings")$scale_y, r$settings$scale_y)
  rep <- hte_report(r, groups = "Asian")
  expect_true(grepl("30 respondents", rep$text[1]))
  expect_true(grepl("10 rows were removed", rep$note))
  # diagnose_imbalance(covariates = ) gives the same sample on its own.
  d2 <- suppressMessages(diagnose_imbalance(x, "race", "treat", "y", covariates = "age"))
  expect_equal(d2$n, dg$n)
})

test_that("the standardizing SD does not change when a group is excluded", {
  x <- toy(); x$treat[x$race == "Asian"] <- 0
  r <- suppressWarnings(compare_hte(x, "y", "treat", "race", methods = "means"))
  dg <- r$diagnostics
  expect_equal(r$settings$scale_y, sd(x$y[x$treat == 0]))
  expect_equal(dg$estimate[dg$group == "Black"], unname(r$gates["Black", "means"]))
})

test_that("infinite values are dropped like missing values", {
  x <- toy(); x$y[1] <- Inf; x$age[2] <- -Inf
  expect_message(r <- compare_hte(x, "y", "treat", "race", "age", methods = "means"),
                 "2 of 940")
  expect_equal(sum(r$sensitivity$n), 938)
})

test_that("`methods` is validated", {
  x <- toy()
  expect_error(compare_hte(x, "y", "treat", "race", methods = c("means", "means", "ols")),
               "more than once")
  expect_error(compare_hte(x, "y", "treat", "race", methods = c("ols", NA)), "non-empty")
  expect_error(compare_hte(x, "y", "treat", "race", methods = character(0)), "non-empty")
  expect_error(compare_hte(x, "y", "treat", "race", methods = 1), "non-empty")
  expect_error(compare_hte(x, "y", "treat", "race", methods = "lasso"), "Unknown")
  r <- compare_hte(x, "y", "treat", "race", methods = "ols")
  expect_equal(r$methods, c("means", "ols"))
})

test_that("other arguments are validated", {
  x <- toy()
  expect_error(compare_hte(x, "y", "treat", "race", methods = "means", standardize = NA), "TRUE or FALSE")
  expect_error(compare_hte(x, "y", "treat", "race", methods = "means", J = 1), "out of range")
  expect_error(compare_hte(x, "y", "treat", "race", methods = "means", J = 2.5), "out of range")
  expect_error(compare_hte(x, "y", "treat", "race", methods = "means", covariates = c("age", "age")), "distinct")
  expect_error(compare_hte(x, "y", "treat", "race", methods = "means", covariates = "race"), "different columns")
  expect_error(compare_hte(x, "y", "treat", "race", methods = "ols", propensity = 1.2),
               "out of range")
  expect_error(compare_hte(x, "y", "treat", "race", methods = "ols", propensity = "p"),
               "not found")
  x$p <- 0.5; x$p[1] <- NA
  expect_error(compare_hte(x, "y", "treat", "race", methods = "ols", propensity = "p"),
               "strictly between")
  x$p <- "half"
  expect_error(compare_hte(x, "y", "treat", "race", methods = "ols", propensity = "p"),
               "numeric")
})

test_that("thin training folds are flagged and recorded", {
  x <- toy(); a <- which(x$race == "Asian")
  x$treat[a] <- c(1, 1, rep(0, 38))         # 2 treated Asian respondents, J = 5
  expect_warning(r <- compare_hte(x, "y", "treat", "race", methods = c("means", "ols")),
                 "fewer than 2 treated or control")
  fc <- r$fold_cells
  expect_equal(nrow(fc), 5 * 4)
  expect_true(min(fc$train_treated[fc$group == "Asian"]) < 2)
  expect_true(all(fc$train_treated[fc$group == "Asian"] >= 1))   # still identified
  expect_false(is.na(r$gates["Asian", "ols"]))
})

test_that("a failing estimator is reported with its reason, not silently dropped", {
  x <- toy()
  # Force a fold failure by making OLS unidentifiable: a covariate that equals the
  # treatment indicator within the training data.
  x$tcopy <- x$treat
  w <- character(0)
  r <- withCallingHandlers(
    compare_hte(x, "y", "treat", "race", "tcopy", methods = c("means", "ols")),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("failed in fold", w)))
  expect_true(nrow(r$failures) > 0)
  expect_true(all(r$failures$method == "ols"))
  expect_true(all(is.na(r$gates[, "ols"])))
  expect_false(anyNA(r$gates[, "means"]))
  out <- capture.output(print(r))
  expect_true(any(grepl("Estimator failures", out)))
  expect_true(grepl("model failures", hte_report(r)$note))
})

test_that("assignment probabilities that differ by group are handled", {
  set.seed(11)
  n <- c(White = 800, Asian = 200)
  race <- rep(names(n), n)
  p <- c(White = 0.5, Asian = 0.8)[race]
  treat <- rbinom(sum(n), 1, p)
  y <- treat * c(White = 0.2, Asian = 0.6)[race] + rnorm(sum(n))
  x <- data.frame(y, treat, race, p)
  r1 <- compare_hte(x, "y", "treat", "race", methods = c("means", "ols"))
  r2 <- compare_hte(x, "y", "treat", "race", methods = c("means", "ols"), propensity = "p")
  expect_equal(r1$gates[, "means"], r2$gates[, "means"])
  d <- diagnose_imbalance(x, "race", "treat")
  expect_true(d$n_eff[d$group == "Asian"] < 200)
  skip_if_not_installed("grf")
  r3 <- compare_hte(x, "y", "treat", "race", methods = c("means", "cf"), propensity = "p",
                    cf_args = list(num.trees = 200))
  expect_false(anyNA(r3$gates))
})

test_that("contrast intervals use Welch-Satterthwaite t, not the normal", {
  x <- toy()
  r <- compare_hte(x, "y", "treat", "race", methods = "means")
  cn <- r$contrasts[r$contrasts$contrast == "Asian - White", ]
  ma <- r$means[r$means$group == "Asian", ]; mw <- r$means[r$means$group == "White", ]
  df <- (ma$se^2 + mw$se^2)^2 / (ma$se^4 / ma$df + mw$se^4 / mw$df)
  expect_equal(cn$means_df, df)
  expect_equal(cn$means_ci_high - cn$means, qt(0.975, df) * cn$means_se)
})

test_that("plot() runs for both types", {
  x <- toy()
  r <- compare_hte(x, "y", "treat", "race", "age", methods = c("means", "ols"))
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  expect_invisible(plot(r))
  expect_invisible(plot(r, type = "contrasts"))
  r0 <- compare_hte(x, "y", "treat", "race", methods = "means", contrasts = "none")
  expect_error(plot(r0, type = "contrasts"), "No contrasts")
})

# Planning edge cases --------------------------------------------------------------

test_that("'both grow' never proposes fewer respondents than the current sample", {
  # Already well powered: both groups large relative to the target.
  n <- c(White = 4000, Asian = 2000)
  x <- data.frame(race = rep(names(n), n),
                  treat = unlist(lapply(n, function(m) rep(0:1, length.out = m))))
  pc <- attr(diagnose_imbalance(x, "race", "treat", target_effect = 0.3), "plan_contrasts")
  expect_true(pc$power_now > 0.8)
  expect_equal(pc$a_needed_if_both_grow, 2000)
  expect_equal(pc$b_needed_if_both_grow, 4000)
  expect_equal(pc$a_additional_if_b_fixed, 0)
  # Under-powered: both grow and both increase.
  pc2 <- attr(diagnose_imbalance(toy(), "race", "treat", target_effect = 0.3,
                                 contrasts = list(c("Asian", "White"))), "plan_contrasts")
  expect_true(pc2$a_needed_if_both_grow > 40 && pc2$b_needed_if_both_grow > 600)
})

test_that("contrast planning uses Welch df and reaches the target power with unbalanced arms", {
  n <- c(White = 1000, Asian = 60)
  x <- data.frame(race = rep(names(n), n))
  x$treat <- c(rep(0:1, length.out = 1000), rep(c(1, 0, 0), length.out = 60))  # Asian 1:2
  d <- diagnose_imbalance(x, "race", "treat", target_effect = 0.3)
  pc <- attr(d, "plan_contrasts"); pg <- attr(d, "plan_groups")
  pa <- pg$treated_share[pg$group == "Asian"]
  m <- pc$a_needed_if_b_fixed
  se <- sqrt(htebalance:::se_dim(m, pa, 1)^2 + htebalance:::se_dim(1000, 0.5, 1)^2)
  cells <- c(m * pa, m * (1 - pa), 500, 500)
  df <- htebalance:::welch_df(1 / cells, cells)
  expect_true(df < m + 1000 - 4)                  # Welch is more conservative
  expect_true(htebalance:::power_at(0.3, se, df, 0.05) >= 0.8)
  se1 <- sqrt(htebalance:::se_dim(m - 1, pa, 1)^2 + htebalance:::se_dim(1000, 0.5, 1)^2)
  cells1 <- c((m - 1) * pa, (m - 1) * (1 - pa), 500, 500)
  expect_true(htebalance:::power_at(0.3, se1, htebalance:::welch_df(1 / cells1, cells1), 0.05) < 0.8)
})

test_that("group planning: n needed is a minimum; additional is never negative", {
  n <- c(White = 3000, Asian = 2000)
  x <- data.frame(race = rep(names(n), n),
                  treat = unlist(lapply(n, function(m) rep(0:1, length.out = m))))
  pg <- attr(diagnose_imbalance(x, "race", "treat", target_effect = 0.2), "plan_groups")
  expect_true(all(pg$power_now > 0.8))
  expect_true(all(pg$n_needed < n[pg$group]))
  expect_equal(pg$n_additional, c(0, 0))
})

test_that("t-approximation power stays close to exact noncentral-t power", {
  # Documented accuracy of the approximation used for MDE and planning.
  exact <- function(effect, se, df, alpha) {
    crit <- qt(1 - alpha / 2, df); ncp <- effect / se
    1 - pt(crit, df, ncp) + pt(-crit, df, ncp)
  }
  for (df in c(4, 8, 18, 38, 98, 498)) {
    se <- 1
    mde <- (qt(0.975, df) + qt(0.8, df)) * se
    approx <- htebalance:::power_at(mde, se, df, 0.05)
    expect_lt(abs(approx - exact(mde, se, df, 0.05)), if (df < 5) 0.015 else 0.002)
  }
})

# Machine-learning estimators on small subgroups -------------------------------------

small_toy <- function(seed = 2) {
  set.seed(seed)
  n <- c(White = 400, Black = 100, Asian = 12)   # Asian: 6 treated, 6 control
  race <- rep(names(n), n)
  treat <- unlist(lapply(n, function(m) sample(rep(0:1, length.out = m))))
  age <- rnorm(sum(n))
  y <- 0.5 * age + treat * c(White = 0.2, Black = 0.5, Asian = -0.2)[race] + rnorm(sum(n))
  data.frame(y, treat, race, age)
}

test_that("causal forest runs on a small subgroup or reports why it cannot", {
  skip_if_not_installed("grf")
  r <- suppressWarnings(compare_hte(small_toy(), "y", "treat", "race", "age",
                                    methods = c("means", "cf"), J = 3,
                                    cf_args = list(num.trees = 200)))
  ok <- !is.na(r$gates["Asian", "cf"])
  expect_true(ok || nrow(r$failures) > 0)
  if (ok) expect_true(is.finite(r$gates["Asian", "cf"]))
  expect_equal(r$sensitivity$n[r$sensitivity$group == "Asian"], 12)
})

test_that("BART runs on a small subgroup or reports why it cannot", {
  skip_if_not_installed("dbarts")
  r <- suppressWarnings(compare_hte(small_toy(), "y", "treat", "race", "age",
                                    methods = c("means", "bart"), J = 3,
                                    bart_args = list(ndpost = 100, nskip = 50)))
  ok <- !is.na(r$gates["Asian", "bart"])
  expect_true(ok || nrow(r$failures) > 0)
  if (ok) expect_true(is.finite(r$gates["Asian", "bart"]))
})

test_that("a failing ML fit is caught per fold and recorded", {
  skip_if_not_installed("grf")
  # An invalid grf argument makes every fold fail; nothing else should break.
  w <- character(0)
  r <- withCallingHandlers(
    compare_hte(small_toy(), "y", "treat", "race", "age", methods = c("means", "cf"), J = 3,
                cf_args = list(num.trees = 200, sample.fraction = 2)),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_equal(nrow(r$failures), 3)
  expect_true(all(is.na(r$gates[, "cf"])))
  expect_false(anyNA(r$gates[, "means"]))
})
