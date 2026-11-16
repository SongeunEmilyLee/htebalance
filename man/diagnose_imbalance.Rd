\name{diagnose_imbalance}
\alias{diagnose_imbalance}
\alias{print.hte_imbalance}
\title{Subgroup composition and precision}
\description{
Reports how much information an experiment holds about each prespecified
subgroup: its size, its share of the sample, treated and control cell counts,
the effective sample size, and (when an outcome is supplied) the within-group
difference in means with its standard error, confidence interval, and minimum
detectable effect. Small size and small share are flagged separately.
}
\usage{
diagnose_imbalance(data, group, treatment, outcome = NULL,
                   small_n = 50, small_share = 0.05, min_cell = 10, J = 5,
                   alpha = 0.05, power = 0.80, standardize = TRUE,
                   target_effect = NULL, contrasts = "reference",
                   reference = NULL, outcome_sd = NULL, covariates = NULL)

\method{print}{hte_imbalance}(x, digits = 2, ...)
}
\arguments{
\item{data}{A data frame.}
\item{group}{Name of the subgroup column, e.g. \code{"race"}.}
\item{treatment}{Name of the treatment column, coded 0/1 or \code{TRUE}/\code{FALSE}.}
\item{outcome}{Optional name of a numeric outcome column. If \code{NULL}, the
standard error and minimum detectable effect assume an outcome SD of 1 in both arms.}
\item{small_n}{Flag groups with fewer respondents than this.}
\item{small_share}{Flag groups below this share of the sample.}
\item{min_cell}{Flag groups with fewer respondents than this in either arm.}
\item{J}{Number of cross-fitting folds you plan to use; used only to report
respondents per held-out fold.}
\item{alpha}{Two-sided significance level for the confidence interval and MDE.}
\item{power}{Target power for the MDE.}
\item{standardize}{If \code{TRUE} and an outcome is given, divide the outcome by
the control-group standard deviation so estimates are in control-group SD units.}
\item{target_effect}{Optional smallest effect of interest (in the same units as
the estimates: control-group SDs when \code{standardize = TRUE}). When supplied,
the function reports the power to detect it now, the total group size needed for
the target \code{power}, and the same for between-group contrasts.}
\item{contrasts}{Contrasts to plan for when \code{target_effect} is given:
\code{"reference"} (each group minus the reference group), \code{"all"},
\code{"none"}, or a list of group pairs such as \code{list(c("Asian", "White"))}.}
\item{reference}{Reference group for \code{contrasts = "reference"}; defaults to
the largest group.}
\item{outcome_sd}{Optional outcome SD to assume for planning, overriding the SD
estimated from the data (or the default of 1 when no outcome is supplied).}
\item{covariates}{Optional covariate names. They are not used in any calculation;
rows with missing values in them are dropped so that the diagnostics describe
the same sample that \code{\link{compare_hte}} estimates on with these
covariates. \code{compare_hte} does this automatically.}
\item{x}{An \code{hte_imbalance} object.}
\item{digits}{Digits to print.}
\item{...}{Ignored.}
}
\details{
Rows with a missing or non-finite group, treatment, outcome, or covariate value
are dropped with a message giving the count per column, and empty factor levels
are removed. All counts refer to the remaining rows; the original and analysed
group sizes are in the \code{sample} attribute.

The effective sample size is \eqn{4 n_T n_C / (n_T + n_C)}. It equals the group
size when the arms are balanced and falls as they become unbalanced.

With an outcome, the standard error is the Neyman variance estimator
\eqn{\sqrt{s_1^2/n_1 + s_0^2/n_0}} and the interval uses Welch degrees of
freedom. The minimum detectable effect is
\eqn{(t_{1-\alpha/2,\,df} + t_{power,\,df}) \times SE}, using the Welch degrees of
freedom (or, when no outcome is supplied, Welch degrees of freedom under equal
SDs in both arms, which equal \eqn{n_1 + n_0 - 2} for balanced arms). This is the
usual approximation to the power of a two-sided t-test; it is not the exact
noncentral-t calculation. At the MDE the approximate power differs from the
exact noncentral-t power by about 0.01 at df = 4 and by less than 0.002 for
df of 8 or more (checked in the package tests). It is computed from
the standard error estimated in this sample, so it reflects this study's outcome
variance and design and is not a general threshold for groups of a given size.
It is a planning quantity, not a test of the observed effect.

With \code{target_effect}, planning uses the same t approximation as the MDE. It
assumes each group's pooled within-arm outcome SD estimated here (or
\code{outcome_sd}, or 1 without an outcome), each group's current treated
share, and the same SD in future samples. Degrees of freedom are
Welch-Satterthwaite, computed from the planning SDs and cell sizes. \code{n_needed}
is the minimum total group size that reaches the target power and can be below
the current size; \code{n_additional} is never negative. For a contrast A - B,
the standard error is \eqn{\sqrt{SE_A^2 + SE_B^2}}, the degrees of freedom are
Welch-Satterthwaite over the four group-by-arm cells, and two scenarios are
reported: adding respondents to A only while B stays fixed, and growing A and B
together in their current proportion. The second scenario never proposes fewer
respondents than the current sample; if the target power is already reached it
returns the current sizes. Adding to A alone can be unable to reach the target,
because B's own sampling error sets a floor on the contrast's standard error;
the table then reports "not reachable".

The flags are screening rules chosen by the user. They do not establish whether
any estimate is valid.
}
\value{
A data frame of class \code{hte_imbalance}, one row per group, sorted from
smallest to largest. With \code{target_effect}, the attributes
\code{plan_groups} (power now, \code{n_needed}, \code{n_additional}) and
\code{plan_contrasts} (power now; group A size needed with B fixed; A and B
sizes needed when both grow) hold the planning results.
}
\examples{
set.seed(1)
n <- c(White = 1500, Black = 300, Latino = 150, Asian = 40)
experiment <- data.frame(race = rep(names(n), n),
                         treat = unlist(lapply(n, function(m) rep(0:1, length.out = m))))
experiment$y <- rnorm(nrow(experiment)) + 0.2 * experiment$treat

diagnose_imbalance(experiment, group = "race", treatment = "treat")
diagnose_imbalance(experiment, group = "race", treatment = "treat", outcome = "y")

# How many respondents would detect a 0.2 SD effect, alone and against White?
diagnose_imbalance(experiment, "race", "treat", outcome = "y", target_effect = 0.2,
                   contrasts = list(c("Asian", "White")))
}
