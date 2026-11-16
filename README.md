# htebalance

An R package for checking subgroup imbalance and comparing heterogeneous
treatment effect (HTE) estimates across estimators in randomized experiments.

**Paper:** Songeun Emily Lee, "Estimating Heterogeneous Treatment Effects under
Subgroup Imbalance." Working paper.
[Link](https://songeunemilylee.com)

## Why this package

Political scientists often ask whether a treatment affects White, Black,
Latino, and Asian respondents differently. Experimental samples rarely give
these groups equal information. A study can have thousands of respondents and
only a few dozen from one racial group, and random assignment then splits that
small group again into treatment and control.

This creates two separate problems:

- **Small group size.** A group with few respondents gives a noisy estimate of
  its own treatment effect.
- **Small group share.** A group that is a small part of the sample has little
  influence when a model is fit or selected on the whole sample.

In the paper, I reanalyzed U.S. experiments on race and ethnicity published in
the *APSR*, *AJPS*, and *JOP* since 2015. Small racial subgroups are common: in
42 percent of the 89 experimental comparisons with at least two racial groups,
the smallest group had fewer than 50 respondents.

![Racial subgroup imbalance in political science experiments](man/figures/paper-imbalance.png)

I then re-estimated each group's treatment effect with four estimators:
subgroup means, interacted OLS, a causal forest, and BART. The respondents,
treatment, outcome, and estimand stay the same; only the estimator changes. The
smaller the group, the more the estimates disagree.

![Estimator sensitivity is larger in smaller racial subgroups](man/figures/paper-sensitivity.png)

| Respondents in the racial group | Median range across the four estimators |
|---|---|
| Fewer than 50 | 0.41 SD |
| 400 or more | 0.05 SD |

*287 group effects from 85 treatment-versus-control comparisons. Effects are in
control-group standard deviations.*

The groups with the least information are the groups for which the choice of
estimator matters most. A confidence interval from one estimator does not show
this, because it describes uncertainty under that estimator only.

`htebalance` puts these checks into three functions so researchers can run them
on their own experiments before reporting subgroup effects.

## Installation

```r
# install.packages("remotes")
remotes::install_github("SongeunEmilyLee/htebalance")
```

Causal forests require `grf`, and BART requires `dbarts`. Both are optional.

```r
install.packages(c("grf", "dbarts"))
```

## Step-by-step example

The example uses simulated data so it runs as is.

### Step 1. Prepare the data

You need one row per respondent, with an outcome, a 0/1 treatment indicator,
and a group variable.

```r
library(htebalance)

set.seed(10)
n <- c(White = 1500, Black = 300, Latino = 150, Asian = 40)
race <- rep(names(n), n)
age <- rnorm(sum(n), 45, 15)
treat <- unlist(lapply(n, function(m) sample(rep(0:1, length.out = m))))
effect <- c(White = 0.2, Black = 0.45, Latino = 0.3, Asian = -0.1)[race]
y <- 0.02 * (age - 45) + treat * effect + rnorm(sum(n))
data <- data.frame(y, treat, race, age)
```

### Step 2. Check how much information each group has

```r
diagnose_imbalance(data, "race", "treat", "y")
```

```
  Group    N Share Treated Control Eff. n Per fold
  Asian   40  2.0%      20      20     40        8
 Latino  150  7.5%      75      75    150       30
  Black  300 15.1%     150     150    300       60
  White 1500 75.4%     750     750   1500      300
 Asian: small n, small share

  Group Estimate   SE        95% CI  MDE
  Asian     0.17 0.36 [-0.57, 0.91] 1.05
 Latino     0.38 0.18  [0.03, 0.73] 0.50
  Black     0.28 0.11  [0.06, 0.50] 0.31
  White     0.26 0.05  [0.16, 0.36] 0.14
```

With 40 respondents, the Asian group can only detect an effect of about one
control-group SD with 80 percent power. Add `target_effect = 0.2` to see how
many respondents each group would need.

### Step 3. Compare estimators

```r
res <- compare_hte(data, "y", "treat", "race",
                   covariates = "age",
                   methods = c("means", "ols", "cf", "bart"))
res
```

```
  Group    n Subgroup means Interacted OLS Causal forest BART Range
  Asian   40           0.17           0.13          0.30 0.26  0.17
 Latino  150           0.38           0.39          0.31 0.37  0.08
  Black  300           0.28           0.33          0.30 0.31  0.05
  White 1500           0.26           0.25          0.26 0.25  0.02
```

`Range` is the largest minus the smallest estimate. It falls from 0.17 in the
Asian group to 0.02 in the White group. Between-group contrasts (here, each
group minus White) are reported in the same way.

If `grf` or `dbarts` is not installed, use `methods = c("means", "ols")`.

### Step 4. Plot the results

```r
plot(res)                      # group effects
plot(res, type = "contrasts")  # differences between groups
```

![Group effects by estimator](man/figures/example-gates.png)

Points are estimates from each estimator. Gray lines are 95 percent confidence
intervals for subgroup means only.

### Step 5. Summarize the results

```r
hte_report(res, groups = "Asian")
```

```
Asian: 40 respondents (2.0%), with 20 treated and 20 control.
Difference in means: 0.17 (SE 0.36, 95% CI -0.57 to 0.91). MDE: 1.05
(power = 0.8, alpha = 0.05). Across estimators: 0.13 to 0.30 (range 0.17).
```

## Functions

| Function | What it reports |
|---|---|
| `diagnose_imbalance()` | Group sizes and shares, treated and control counts, effective sample size, difference in means, SE, minimum detectable effect, and (with `target_effect`) the sample size needed |
| `compare_hte()` | Group treatment effects and between-group contrasts from subgroup means, interacted OLS, causal forests, and BART, with the range across estimators |
| `hte_report()` | A short written summary of each group |

## Interpretation

- The package is designed for randomized experiments.
- Small group size and small group share are reported separately.
- Differences across estimators show how sensitive a conclusion is to the
  method used. They do not show which estimate is more accurate, because the
  true effects are unknown in real data.
- The paper also proposes STAGE, an ensemble rule that chooses estimator
  weights using held-out error on prespecified group effects. STAGE is not part
  of this release.

## Citation

Lee, Songeun Emily. "Estimating Heterogeneous Treatment Effects under Subgroup
Imbalance." Working paper.

## License

MIT
