\name{compare_hte}
\alias{compare_hte}
\alias{print.hte_compare}
\alias{plot.hte_compare}
\title{Compare group treatment effects across estimators}
\description{
Estimates the average treatment effect for each prespecified group with several
estimators that borrow different amounts of information across groups, along
with between-group contrasts, and reports how much the estimates differ. In real
data the true effects are unknown, so the differences measure sensitivity to
estimator choice, not estimation error.

The function is intended for randomized experiments, where the group-level
difference in means is an unbiased estimate of each group's average treatment
effect. With observational data, supplying \code{propensity} does not make the
subgroup-means column a causal estimate, and none of the output should be read
as one without a separate identification argument.
}
\usage{
compare_hte(data, outcome, treatment, group, covariates = NULL,
            methods = c("means", "ols", "cf"), contrasts = "reference",
            reference = NULL, propensity = NULL, standardize = TRUE,
            J = 5, seed = 1, level = 0.95, cf_args = list(), bart_args = list(),
            min_train_cell = 2)

\method{print}{hte_compare}(x, digits = 2, ...)

\method{plot}{hte_compare}(x, type = c("gates", "contrasts"), main = NULL, ...)
}
\arguments{
\item{data}{A data frame.}
\item{outcome, treatment, group}{Column names. Treatment must be 0/1 or logical.}
\item{covariates}{Optional character vector of pretreatment covariates used by
the adjusted estimators.}
\item{methods}{Any of \code{"means"} (within-group difference in means),
\code{"ols"} (interacted OLS), \code{"cf"} (causal forest, needs \pkg{grf}),
and \code{"bart"} (BART, needs \pkg{dbarts}). Each method may appear only once.
Subgroup means are always included (added first if not listed).}
\item{contrasts}{\code{"reference"} (each group minus a reference group),
\code{"all"} (every pair), \code{"none"}, or a list of group pairs such as
\code{list(c("Asian", "White"))}.}
\item{reference}{Reference group for \code{contrasts = "reference"}; defaults
to the largest group.}
\item{propensity}{Known treatment probability: \code{NULL} uses the treated share
within each group, a single number applies to everyone, or a column name gives
individual probabilities (numeric, strictly between 0 and 1, no missing values
in the analysis sample). Used by the causal forest.}
\item{standardize}{Divide the outcome by the control-group SD.}
\item{J}{Number of cross-fitting folds for the adjusted estimators.}
\item{seed}{Random seed for the fold split.}
\item{level}{Confidence level for the subgroup-means intervals.}
\item{cf_args, bart_args}{Lists of extra arguments passed to
\code{grf::causal_forest} and \code{dbarts::bart}.}
\item{min_train_cell}{Warn when, in any cross-fitting fold, a group has fewer
than this many treated or control respondents in the training data. The
per-fold counts are returned in \code{fold_cells}.}
\item{x}{An \code{hte_compare} object.}
\item{digits}{Digits to print.}
\item{type}{Plot group effects (\code{"gates"}) or between-group contrasts (\code{"contrasts"}).}
\item{main}{Optional plot title. By default no title is drawn, so the figure can
be captioned in a manuscript.}
\item{...}{Passed to \code{plot}.}
}
\details{
All results use one analysis sample: rows with a missing or non-finite value in
the outcome, treatment, group, or any covariate are dropped first (with a
message giving the count per column), and the diagnostics in
\code{diagnostics} are computed on the same rows. The original and analysed
group sizes are returned in \code{sample}. With \code{standardize = TRUE} the
outcome is divided by the control-group SD of this analysis sample, before any
group is excluded, so the scale matches \code{\link{diagnose_imbalance}}.

Subgroup means use only treated and control respondents within each group and
are computed on the full analysis sample. Interacted OLS regresses the outcome on
treatment, group indicators, and covariates, with treatment interacted with
both. The causal forest and BART are estimated with the same group indicators
and covariates. The adjusted estimators are cross-fitted: folds are stratified
by group and treatment arm, each respondent's effect is predicted from a model
fit without that respondent's fold, and the predictions are averaged within
group.

A group with fewer than two respondents in either arm is excluded from
estimation with a warning; its effect is not replaced with a pooled value.
Before fitting, the training data of every fold are checked: a group with fewer
than \code{min_train_cell} treated or control respondents in some training fold
triggers a warning, because its model-based estimates then depend mostly on
other groups. If an estimator fails in a fold (an error, or missing or
non-finite predictions), the groups in that held-out fold receive \code{NA} for
that estimator, and the method, fold, groups, and error message are recorded in
\code{failures} and printed.

The model-based estimates are averages of cross-fitted individual predictions
within each group. Causal forests target treatment-effect heterogeneity
directly, while BART here models the outcome with treatment as a predictor;
both regularize, and small groups can be shrunk toward the others. The
comparison therefore describes how sensitive a group's estimate is to the
choice of estimator. It is not a test of which estimator is more accurate or of
bias for a particular group.

Intervals are reported for subgroup means only. Contrast intervals treat the
two groups' estimates as independent and use Welch-Satterthwaite degrees of
freedom combining the two groups' Welch degrees of freedom.

In the plot, points are the estimates from each estimator (estimator
sensitivity) and gray intervals are 95\% confidence intervals for subgroup
means only (sampling uncertainty). The gray intervals are not intervals for the
other estimators.

\code{sign_disagree} is \code{TRUE} whenever at least one estimate is positive
and another negative, including when an estimate is very close to zero. Read it
together with the estimates themselves.
}
\references{
For formal tests of whether machine-learning group estimates differ from the
experimental difference in means, and for calibration or bias correction, see
Leng, Y. and Dimmery, D. (2024), "Calibration of Heterogeneous Treatment Effects
in Randomized Experiments," \emph{Information Systems Research} 35(4):1721--1742,
and Persson, J. et al. (2026), "Detecting and Mitigating Group Bias in
Heterogeneous Treatment Effects," arXiv:2602.20383. \code{compare_hte} reports
differences across estimators descriptively and does not test or correct them.
}
\value{
An object of class \code{hte_compare}: a list with \code{gates} (groups by
estimators), \code{means} (difference in means with SE and CI),
\code{sensitivity} (range and sign disagreement by group), \code{contrasts},
\code{diagnostics} (the \code{\link{diagnose_imbalance}} table), cross-fitted
respondent-level predictions in \code{tau}, \code{sample} (original and
analysed group sizes), \code{n_original} and \code{n_dropped},
\code{fold_cells} (training-fold cell counts by group), \code{failures}
(estimator failures by fold), and the settings used.
}
\examples{
set.seed(1)
n <- c(White = 600, Black = 200, Latino = 100, Asian = 40)
race <- rep(names(n), n)
treat <- unlist(lapply(n, function(m) sample(rep(0:1, length.out = m))))
age <- rnorm(sum(n))
y <- 0.5 * age + treat * c(White = 0.2, Black = 0.5, Latino = 0.3, Asian = -0.2)[race] +
  rnorm(sum(n))
experiment <- data.frame(y, treat, race, age)

res <- compare_hte(experiment, "y", "treat", "race", covariates = "age",
                   methods = c("means", "ols"))
res
plot(res)
plot(res, type = "contrasts")
}
